#include "pe_options.h"

int main(int argc, char **argv)
{
  GtError *error;
  GtLogger *logger;
  GtQueue *streams;
  GtNodeStream *current_stream, *last_stream, *refrgff3, *predgff3, *tempstream;
  GtNodeVisitor *rpt;
  char *start_time;

  gt_lib_init();
  start_time = pe_get_start_time();

  // Parse command-line options
  ParsEvalOptions options;
  pe_set_option_defaults(&options);
  error = gt_error_new();
  pe_parse_options(argc, argv, &options, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[ParsEval] error: %s", gt_error_get(error));
    return 1;
  }
  int numfiles = argc - optind;
  if(numfiles != 2)
  {
    fprintf(stderr, "[ParsEval] error: must provide two GFF3 files as input");
    pe_print_usage(stderr);
    return 1;
  }

  logger = gt_logger_new(true, "", stderr);
  streams = gt_queue_new();


  //----- Set up the node processing stream -----//
  //---------------------------------------------//

  refrgff3 = gt_gff3_in_stream_new_unsorted(1, &options.refrfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)refrgff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)refrgff3);
  gt_queue_add(streams, refrgff3);
  if(options.refrlabel != NULL)
  {
    GtStr *label = gt_str_new_cstr(options.refrlabel);
    GtNodeVisitor *nv = gt_set_source_visitor_new(label);
    tempstream = gt_visitor_stream_new(refrgff3, nv);
    gt_queue_add(streams, tempstream);
    refrgff3 = tempstream;
    gt_str_delete(label);
  }
  tempstream = agn_gene_stream_new(refrgff3, logger);
  gt_queue_add(streams, tempstream);
  refrgff3 = tempstream;

  predgff3 = gt_gff3_in_stream_new_unsorted(1, &options.predfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)predgff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)predgff3);
  gt_queue_add(streams, predgff3);
  if(options.predlabel != NULL)
  {
    GtStr *label = gt_str_new_cstr(options.predlabel);
    GtNodeVisitor *nv = gt_set_source_visitor_new(label);
    tempstream = gt_visitor_stream_new(predgff3, nv);
    gt_queue_add(streams, tempstream);
    predgff3 = tempstream;
    gt_str_delete(label);
  }
  tempstream = agn_gene_stream_new(predgff3, logger);
  gt_queue_add(streams, tempstream);
  predgff3 = tempstream;

  current_stream = agn_locus_stream_new_pairwise(refrgff3, predgff3, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  /* FIXME I don't understand why this is needed, but without it memory is
   * leaked; if it's not in this precise location, no memory is leaked but
   * segfaults result */
  current_stream = agn_node_delete_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(options.filters != NULL)
  {
    current_stream = agn_locus_filter_stream_new(last_stream, options.filters);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  switch(options.outfmt)
  {
    case TEXTMODE:
      if(options.summary_only)
        rpt = agn_compare_report_text_new(NULL, logger);
      else
        rpt = agn_compare_report_text_new(options.outfile, logger);
      break;
    case HTMLMODE:
      if(options.graphics)
      {
        rpt = agn_compare_report_html_new(options.outfilename, &options.pngdata,
                                          logger);
      }
      else
        rpt = agn_compare_report_html_new(options.outfilename, NULL, logger);
      break;
    case CSVMODE:
      fprintf(stderr, "error: CSV output mode support temporarily "
              "unavailable\n");
      return 1;
      break;
    default:
      fprintf(stderr, "error: unknown output format\n");
      return 1;
      break;
  }
  current_stream = gt_visitor_stream_new(last_stream, rpt);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;


  //----- Execute the node processing stream -----//
  //----------------------------------------------//

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
    fprintf(stderr, "[ParsEval] error: %s", gt_error_get(error));

  if(options.outfmt == TEXTMODE)
  {
    pe_summary_header(&options, options.outfile, start_time, argc, argv);
    agn_compare_report_text_create_summary((AgnCompareReportText *)rpt,
                                           options.outfile);
  }
  else if(options.outfmt == HTMLMODE)
  {
    // FIXME Place the next 4-6 lines elsewhere?
    const char *refrlabel = options.refrfile;
    if(options.refrlabel)
      refrlabel = options.refrlabel;
    const char *predlabel = options.predfile;
    if(options.predlabel)
      predlabel = options.predlabel;
    agn_compare_report_html_create_summary((AgnCompareReportHTML *)rpt, argc,
                                           argv, refrlabel, predlabel,
                                           start_time);
  }

  // Free memory and terminate
  gt_free(start_time);
  pe_free_option_memory(&options);
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_queue_delete(streams);
  gt_logger_delete(logger);
  gt_error_delete(error);
  gt_lib_clean();
  return 0;
}
