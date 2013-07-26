#ifndef PARSEVAL_COMPARATIVE_EVALUATION
#define PARSEVAL_COMPARATIVE_EVALUATION

#include "AgnCliquePair.h"
#include "AgnComparEval.h"
#include "AgnLocusIndex.h"
#include "PeOptions.h"

/**
 * Supplemental data structure used during comparative analysis of gene loci.
 */
typedef struct
{
  GtHashmap *comp_evals;
  GtHashmap *locus_summaries;
  PeOptions *options;
  FILE *seqfile;
} PeAnalysisData;

/**
 * Perform comparison of annotations (reference vs prediction) for each locus in
 * the given locus index.
 *
 * @param[in]  locusindex          data structure in which all loci are stored
 * @param[out] comp_evalsp         pointer to a hashmap for storing comparison
 *                                 results; key=memory address for locus,
 *                                 value=pointer to AgnCompEvaluation struct
 * @param[out] locus_summariesp    pointer to a hashmap for storing locus
 *                                 summaries; key=memory address for locus,
 *                                 value=pointer to AgnGeneLocusSummary struct
 * @param[in]  seqids              list of sequence IDs
 * @param[in]  seqfiles            list of sequence-specific output files (for
 *                                 HTML mode only)
 * @param[in]  loci                array of arrays, each sub-array containing
 *                                 loci for a given sequence, sorted by position
 * @param[in]  options             ParsEval options
 */
void pe_comparative_analysis(AgnLocusIndex *locusindex, GtHashmap **comp_evalsp,
                             GtHashmap **locus_summariesp, GtStrArray *seqids,
                             GtArray *seqfiles, GtArray *loci,
                             PeOptions *options);

/**
 * Collect information from the given locus prior to comparative analysis.
 *
 * @param[in]  locus    locus object
 * @param[out] data     data structures to store the needed information
 */
void pe_pre_analysis(AgnGeneLocus *locus, PeAnalysisData *data);

/**
 * Collect information from the given locus following comparative analysis.
 *
 * @param[in]  locus    locus object
 * @param[out] data     data structures to store the needed information
 */
void pe_post_analysis(AgnGeneLocus *locus, PeAnalysisData *data);

#endif
