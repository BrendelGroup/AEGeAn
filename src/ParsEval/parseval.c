/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include "pe_options.h"
#include "pe_utils.h"

int main(int argc, char **argv)
{
  GtError *error;
  GtLogger *logger;
  GtQueue *streams;
  GtNodeStream *current_stream, *last_stream;
  GtNodeVisitor *rpt;
  PeHtmlOverviewData odata;
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

  const char * infiles[] = { options.refrfile, options.predfile };
  current_stream = gt_gff3_in_stream_new_unsorted(2, infiles);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = gt_sort_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;
  
  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, 0);
  agn_locus_stream_skip_empty_loci((AgnLocusStream *)current_stream);
  agn_locus_stream_label_pairwise((AgnLocusStream *)current_stream,
                                  options.refrfile, options.predfile);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(gt_array_size(options.filters) > 0)
  {
    current_stream = agn_locus_filter_stream_new(last_stream, options.filters);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  switch(options.outfmt)
  {
    case TEXTMODE:
      if(options.summary_only)
        rpt = agn_compare_report_text_new(NULL, false, logger);
      else
        rpt = agn_compare_report_text_new(options.outfile,options.gff3,logger);
      break;
    case HTMLMODE:
      if(options.graphics)
      {
        rpt = agn_compare_report_html_new(options.outfilename, options.gff3,
                                          &options.pngdata, logger);
      }
      else
      {
        rpt = agn_compare_report_html_new(options.outfilename, options.gff3,
                                          NULL, logger);
      }
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
    odata.refrlabel = options.refrfile;
    if(options.refrlabel)
      odata.refrlabel = options.refrlabel;
    odata.predlabel = options.predfile;
    if(options.predlabel)
      odata.predlabel = options.predlabel;
    odata.start_time = start_time;
    odata.argc = argc;
    odata.argv = argv;

    agn_compare_report_html_set_overview_func((AgnCompareReportHTML *)rpt,
                                              pe_summary_html_overview,
                                              &odata);
    agn_compare_report_html_create_summary((AgnCompareReportHTML *)rpt);
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
