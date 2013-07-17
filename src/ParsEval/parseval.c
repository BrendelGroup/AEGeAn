#include <string.h>
#include "AgnCanonNodeVisitor.h"
#include "AgnLocusIndex.h"
#include "AgnGeneLocus.h"
#include "AgnUtils.h"
#include "PeReports.h"

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
  pe_print_summary(start_time_str, argc, argv, seqids, &overall_eval,
                   seqlevel_evals, options.outfile, &options);
  pe_print_combine_output(seqids, seqfiles, &options);


  // All done!
  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] ParsEval complete! (total runtime:"
                          " %ld.%06ld seconds)\n\n", stderr );

  // Free up remaining memory
  gt_free(start_time_str);
  gt_array_delete(seqfiles);
  gt_hashmap_delete(comp_evals);
  agn_logger_delete(logger);
  agn_locus_index_delete(locusindex);
  gt_array_delete(seqlevel_evals);
  gt_timer_delete(timer);
  if(gt_lib_clean() != 0)
  {
    fputs("error: issue cleaning GenomeTools library\n", stderr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
