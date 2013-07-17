#include <ctype.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include "AgnCanonNodeVisitor.h"
#include "AgnLocusIndex.h"
#include "AgnGeneLocus.h"
#include "AgnUtils.h"
#include "PeReports.h"

typedef struct
{
  GtHashmap *comp_evals;
  GtHashmap *locus_summaries;
  PeOptions *options;
  FILE *seqfile;
} PeAnalysisData;

char *pe_get_start_time()
{
  time_t start_time;
  struct tm *start_time_info;
  time(&start_time);
  start_time_info = localtime(&start_time);

  char timestr[128];
  strftime(timestr, 128, "%d %b %Y, %I:%M%p", start_time_info);
  return gt_cstr_dup(timestr);
}

void pe_seqid_check(const char *seqid, AgnLogger *logger)
{
  size_t n = strlen(seqid);
  int i;
  for(i = 0; i < n; i++)
  {
    char c = seqid[i];
    if(!isalnum(c) && c != '.' && c != '-' && c != '_')
    {
      agn_logger_log_error(logger, "seqid '%s' contains illegal characters; "
                           "only alphanumeric characters and . and _ and - are "
                           "allowed.", seqid);
      return;
    }
  }
}

unsigned long pe_load_and_parse_loci(AgnLocusIndex **locusindexp, GtArray **locip,
                                     GtStrArray **seqidsp, PeOptions *options,
                                     AgnLogger *logger)
{
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin loading data and parsing loci\n", stderr);

  // Load loci into memory
  AgnLocusIndex *locusindex = agn_locus_index_new(false);
  unsigned long total = agn_locus_index_parse_pairwise_disk(locusindex,
                            options->refrfile, options->predfile,
                            options->numprocs, &options->filters, logger);

  // Collect IDs of all sequences annotated by input files
  GtStrArray *seqids = agn_locus_index_seqids(locusindex);
  GtArray *loci = gt_array_new( sizeof(GtArray *) );
  int i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    pe_seqid_check(seqid, logger);

    GtArray *seq_loci = agn_locus_index_get(locusindex, seqid);
    gt_array_sort(seq_loci,(GtCompare)agn_gene_locus_array_compare);
    gt_array_add(loci, seq_loci);
  }

  *locusindexp = locusindex;
  *locip = loci;
  *seqidsp = seqids;

  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished loading data and "
                          "parsing loci (%ld.%06ld seconds)\n", stderr);
  gt_timer_delete(timer);
  return total;
}

GtArray *pe_prep_output(GtStrArray *seqids, PeOptions *options)
{
  GtArray *seqfiles = gt_array_new( sizeof(FILE *) );

  if(strcmp(options->outfmt, "csv") == 0)
    pe_print_csv_header(options->outfile);

  unsigned long i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    FILE *seqfile = NULL;
    if(!options->summary_only)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      char filename[512];
      if(options->html)
      {
        char dircmd[512];
        sprintf(dircmd, "mkdir %s/%s", options->outfilename, seqid);
        if(system(dircmd) != 0)
        {
          fprintf(stderr, "error: could not open directory '%s/%s'\n",
                  options->outfilename, seqid);
          exit(1);
        }
        if(options->debug)
          fprintf(stderr, "debug: opening directory '%s'\n", dircmd);

        sprintf(dircmd, "ln -s ../LICENSE %s/%s/LICENSE", options->outfilename,
                seqid);
        if(system(dircmd) != 0)
        {
          fputs("warning: could not create symbolic link to LICENSE\n", stderr);
        }

        // Create summary page for this sequence
        sprintf(filename, "%s/%s/index.html", options->outfilename, seqid);
        if(options->debug)
          fprintf(stderr, "debug: opening outfile '%s'\n", filename);
        seqfile = agn_fopen(filename, "w", stderr);
        pe_print_seqfile_header(seqfile, seqid);
        gt_array_add(seqfiles, seqfile);
      }
      else
      {
        sprintf(filename, "%s.%s", options->outfilename, seqid);
        if(options->debug)
          fprintf(stderr, "debug: opening temp outfile '%s'\n", filename);
        seqfile = agn_fopen(filename, "w", stderr);
        gt_array_add(seqfiles, seqfile);
      }
    }
  }

  return seqfiles;
}

void pe_pre_analysis(AgnGeneLocus *locus, PeAnalysisData *data)
{
  unsigned long npairs = agn_gene_locus_enumerate_clique_pairs(locus);
  if(data->options->complimit != 0 && npairs > data->options->complimit)
  {
    if(data->options->debug)
    {
      fprintf(stderr, "debug: locus %s[%lu, %lu] contains %lu transcript (or "
              "transcript clique) pairs, exceeds the limit of %d, moving on\n",
              agn_gene_locus_get_seqid(locus), agn_gene_locus_get_start(locus),
              agn_gene_locus_get_end(locus), npairs, data->options->complimit);
    }
  }
  else
  {
    PeCompEvaluation *compeval = gt_hashmap_get(data->comp_evals, locus);
    compeval->counts.num_loci++;
    compeval->counts.refr_genes += agn_gene_locus_num_refr_genes(locus);
    compeval->counts.pred_genes += agn_gene_locus_num_pred_genes(locus);
    unsigned long rt = agn_gene_locus_num_refr_transcripts(locus);
    unsigned long pt = agn_gene_locus_num_pred_transcripts(locus);
    compeval->counts.refr_transcripts += rt;
    compeval->counts.pred_transcripts += pt;
    if(agn_gene_locus_num_refr_genes(locus) == 0)
      compeval->counts.unique_pred++;
    else if(agn_gene_locus_num_pred_genes(locus) == 0)
      compeval->counts.unique_refr++;
  }
}

void pe_post_analysis(AgnGeneLocus *locus, PeAnalysisData *data)
{
  PeCompEvaluation *compeval = gt_hashmap_get(data->comp_evals, locus);
  GtArray *pairs = agn_gene_locus_pairs_to_report(locus);
  if(gt_array_size(pairs) > 0)
    agn_gene_locus_aggregate_results(locus, compeval);

  if(!data->options->summary_only)
  {
    pe_gene_locus_print_results(locus, data->seqfile, data->options);

    AgnGeneLocusSummary *locsum = gt_hashmap_get(data->locus_summaries, locus);
    locsum->start = agn_gene_locus_get_start(locus);
    locsum->end = agn_gene_locus_get_end(locus);
    locsum->length = agn_gene_locus_get_length(locus);
    locsum->refrtrans = agn_gene_locus_num_refr_transcripts(locus);
    locsum->predtrans = agn_gene_locus_num_pred_transcripts(locus);
    locsum->reported = gt_array_size(pairs);
    locsum->counts = compeval->counts;

#ifndef WITHOUT_CAIRO
    if(data->options->locus_graphics)
    {
      AgnGeneLocusPngMetadata metadata;
      pe_gene_locus_get_png_filename(locus, metadata.filename,
                                     data->options->outfilename);
      metadata.graphic_width = pe_gene_locus_get_graphic_width(locus);
      sprintf(metadata.stylefile, "%s/pe.style", data->options->data_path);
      metadata.refrfile = data->options->refrfile;
      metadata.predfile = data->options->predfile;
      metadata.refrlabel = data->options->refrlabel;
      metadata.predlabel = data->options->predlabel;
      metadata.track_order_func = pe_track_order;

      agn_gene_locus_print_png(locus, &metadata);
    }
#endif
  }

  agn_gene_locus_delete(locus);
}

void pe_comp_eval(AgnLocusIndex *locusindex, GtHashmap **comp_evalsp,
                  GtHashmap **locus_summariesp, GtStrArray *seqids,
                  GtArray *seqfiles, GtArray *loci, PeOptions *options)
{
  GtHashmap *comp_evals;
  GtHashmap *locus_summaries;

  AgnLogger *logger = agn_logger_new();
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin comparative analysis\n", stderr);
  locus_summaries = gt_hashmap_new(GT_HASH_DIRECT, NULL, (GtFree)gt_free_mem);
  comp_evals = gt_hashmap_new(GT_HASH_DIRECT, NULL, (GtFree)gt_free_mem);

  int i, j;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *seqloci = *(GtArray **)gt_array_get(loci, i);
    PeAnalysisData analysis_data;
    analysis_data.seqfile = *(FILE **)gt_array_get(seqfiles, i);
    analysis_data.options = options;

    for(j = 0; j < gt_array_size(seqloci); j++)
    {
      AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(seqloci, j);
      PeCompEvaluation *compeval = gt_malloc( sizeof(PeCompEvaluation) );
      pe_comp_evaluation_init(compeval);
      gt_hashmap_add(comp_evals, locus, compeval);
      AgnGeneLocusSummary *locsum = gt_malloc( sizeof(AgnGeneLocusSummary) );
      agn_gene_locus_summary_init(locsum);
      gt_hashmap_add(locus_summaries, locus, locsum);
    }
    analysis_data.comp_evals = comp_evals;
    analysis_data.locus_summaries = locus_summaries;

    agn_locus_index_comparative_analysis(locusindex, seqid, options->numprocs,
                                     (AgnLocusIndexVisitFunc)pe_pre_analysis,
                                     (AgnLocusIndexVisitFunc)pe_post_analysis,
                                     &analysis_data, logger);
    if(options->debug)
      agn_logger_print_status(logger, stderr, "comparative analysis");
    if(agn_logger_has_error(logger))
    {
      agn_logger_print_error(logger, stderr, "issues with comparative analysis "
                             "of sequence '%s'", seqid);
      exit(1);
    }
  }

  *comp_evalsp = comp_evals;
  *locus_summariesp = locus_summaries;

  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished comparative "
                          "analysis (%ld.%06ld seconds)\n", stderr);
  agn_logger_delete(logger);
}

void pe_agg_results(PeCompEvaluation *overall_eval, GtArray **seqlevel_evalsp,
                    GtArray *loci, GtArray *seqfiles, GtHashmap *comp_evals,
                    GtHashmap *locus_summaries, PeOptions *options)
{
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin aggregating locus-level results\n", stderr);

  pe_comp_evaluation_init(overall_eval);
  GtArray *seqlevel_evals = gt_array_new( sizeof(PeCompEvaluation) );
  unsigned long i;
  for(i = 0; i < gt_array_size(seqfiles); i++)
  {
    FILE *seqfile = *(FILE **)gt_array_get(seqfiles, i);
    PeCompEvaluation seqeval;
    pe_comp_evaluation_init(&seqeval);

    int j;
    GtArray *seqloci = *(GtArray **)gt_array_get(loci, i);
    for(j = 0; j < gt_array_size(seqloci); j++)
    {
      AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(seqloci, j);
      PeCompEvaluation *eval = gt_hashmap_get(comp_evals, locus);
      pe_comp_evaluation_combine(&seqeval, eval);
      pe_comp_evaluation_combine(overall_eval, eval);
      AgnGeneLocusSummary *locsum = gt_hashmap_get(locus_summaries, locus);
      if(options->html && !options->summary_only)
      {
        pe_print_locus_to_seqfile(seqfile, locsum->start, locsum->end,
                                  locsum->length, locsum->refrtrans,
                                  locsum->predtrans, &locsum->counts);
      }
    }
    gt_array_add(seqlevel_evals, seqeval);
    gt_array_delete(seqloci);
  }

  *seqlevel_evalsp = seqlevel_evals;

  gt_hashmap_delete(locus_summaries);
  gt_array_delete(loci);
  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished aggregating locus-"
                          "level results (%ld.%06ld seconds)\n", stderr);
}

// Main method
int main(int argc, char * const argv[])
{
  // Initialize ParsEval
  gt_lib_init();
  char *start_time_str = pe_get_start_time();
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin ParsEval\n", stderr);

  // Parse command-line arguments
  PeOptions options;
  pe_set_option_defaults(&options);
  int optind = pe_parse_options(argc, argv, &options);
  options.refrfile = argv[optind];
  options.predfile = argv[optind + 1];
  if(options.refrfile == NULL || options.predfile == NULL)
  {
    fprintf(stderr, "[ParsEval] error: could not parse input filenames\n");
    return EXIT_FAILURE;
  }

  // Load data into memory and parse loci
  AgnLocusIndex *locusindex;
  GtArray *loci;
  GtStrArray *seqids;
  AgnLogger *logger = agn_logger_new();
  unsigned long totalloci = pe_load_and_parse_loci(&locusindex, &loci, &seqids,
                                                   &options, logger);
  bool haderror = agn_logger_print_all(logger, stderr, NULL);
  if(haderror)
    return EXIT_FAILURE;

  // Stop now if there are no loci
  if(totalloci == 0)
  {
    fprintf(stderr, "[ParsEval] Warning: found no loci to analyze");
    gt_timer_stop(timer);
    gt_timer_show_formatted(timer,"[ParsEval] ParsEval complete! (total "
                            "runtime: %ld.%06ld seconds)\n\n", stderr );
    fclose(options.outfile);
    gt_timer_delete(timer);
    if(gt_lib_clean() != 0)
    {
      fputs("error: issue cleaning GenomeTools library\n", stderr);
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  // Prep report output
  GtArray *seqfiles = pe_prep_output(seqids, &options);

  // Comparative analysis of loci
  GtHashmap *comp_evals, *locus_summaries;
  pe_comp_eval(locusindex, &comp_evals, &locus_summaries, seqids, seqfiles,
               loci, &options);

  // Aggregate locus-level results at the sequence level and overall
  PeCompEvaluation overall_eval;
  GtArray *seqlevel_evals;
  pe_agg_results(&overall_eval, &seqlevel_evals, loci, seqfiles, comp_evals,
                 locus_summaries, &options);

  // Print summary statistics, combine output
  GtTimer *timer_short = gt_timer_new();
  gt_timer_start(timer_short);
  fputs("[ParsEval] Begin printing summary, combining output\n", stderr);

  pe_print_summary(start_time_str, argc, argv, seqids, &overall_eval,
                   seqlevel_evals, options.outfile, &options );
  if(!options.summary_only)
  {
    unsigned long i;
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      FILE *seqfile = *(FILE **)gt_array_get(seqfiles, i);
      if(options.html)
        pe_print_seqfile_footer(seqfile);
      fclose(seqfile);
    }
  }
  gt_array_delete(seqfiles);

  if(options.outfile == stdout)
    fflush(stdout);
  else
    fclose(options.outfile);
  if(!options.html && !options.summary_only)
  {
    int result;
    unsigned long i;
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      char command[1024];
      sprintf( command, "cat %s.%s >> %s && rm %s.%s", options.outfilename, seqid,
               options.outfilename, options.outfilename, seqid );
      if(options.outfile == stdout)
      {
        sprintf( command, "cat %s.%s && rm %s.%s", options.outfilename, seqid, options.outfilename,
                 seqid );
      }
      if(options.debug)
        fprintf(stderr, "debug: merging output files: %s\n", command);
      result = system(command);
      if(result)
      {
        fprintf(stderr, "[ParsEval] error: issue merging GFF3 files: %s\n", command);
        exit(1);
      }
    }
  }
  gt_timer_stop(timer_short);
  gt_timer_show_formatted( timer_short, "[ParsEval] Finished printing summary, combining output "
                           "(%ld.%06ld seconds)\n", stderr );
  if(options.outfile == stdout)
    fclose(stdout);


  // All done!
  gt_timer_stop(timer);
  gt_timer_show_formatted( timer,"[ParsEval] ParsEval complete! (total runtime: %ld.%06ld "
                           "seconds)\n\n", stderr );

  // Free up remaining memory
  gt_hashmap_delete(comp_evals);
  agn_logger_delete(logger);
  agn_locus_index_delete(locusindex);
  gt_array_delete(seqlevel_evals);
  gt_timer_delete(timer);
  gt_timer_delete(timer_short); // FIXME remove
  if(gt_lib_clean() != 0)
  {
    fputs("error: issue cleaning GenomeTools library\n", stderr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
