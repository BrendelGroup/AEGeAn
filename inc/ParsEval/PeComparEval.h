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
 * Each transcript clique pair that is compared is classified as one of the
 * following.
 *   - perfect match
 *   - perfect match with mislabeled UTRs
 *   - CDS match
 *   - exon structure match
 *   - UTR structure match
 *   - non-match
 *
 * When reporting the results of a comparative analysis, it may be useful to (as
 * is done by ParsEval) show some basic information about clique pairs that fall
 * under each classification category. The counts in this struct are necessary
 * to calculate those summary characteristics.
 */
typedef struct
{
  unsigned long total_length;
  unsigned long transcript_count;
  unsigned long refr_cds_length;
  unsigned long pred_cds_length;
  unsigned long refr_exon_count;
  unsigned long pred_exon_count;
} PeCompResultDesc;

/**
 * This struct is used to aggregate characteristics for all of the
 * classification categories.
 */
typedef struct
{
  PeCompResultDesc perfect_matches;
  PeCompResultDesc perfect_mislabeled;
  PeCompResultDesc cds_matches;
  PeCompResultDesc exon_matches;
  PeCompResultDesc utr_matches;
  PeCompResultDesc non_matches;
} PeCompResultSummary;

/**
 * This struct provides a convenient way to manage the counts, stats, and
 * results corresponding to one or more comparisons.
 */
typedef struct
{
  AgnCompSummary counts;
  AgnComparison stats;
  PeCompResultSummary results;
} PeCompEvaluation;

/**
 * Add information about a clique pair to a set of aggregate characteristics.
 *
 * @param[in]  pair               the clique pair
 * @param[out] characteristics    a set of aggregate characteristics
 */
void pe_clique_pair_record_characteristics(AgnCliquePair *pair,
                                           PeCompResultDesc *desc);

/**
 * Perform comparison of annotations (reference vs prediction) for each locus in
 * the given locus index.
 *
 * @param[in]  locusindex          data structure in which all loci are stored
 * @param[out] comp_evalsp         pointer to a hashmap for storing comparison
 *                                 results; key=memory address for locus,
 *                                 value=pointer to PeCompEvaluation struct
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
 * Take values from one data set and add them to the other.
 *
 * @param[out] data           stats, counts, and descriptions relevant to
 *                            comparative analysis
 * @param[in]  data_to_add    a data set which will be added to the first
 */
void pe_comp_evaluation_combine(PeCompEvaluation *data,
                                 PeCompEvaluation *data_to_add);

/**
 * Initialize to default values.
 *
 * @param[in] data    the data
 */
void pe_comp_evaluation_init(PeCompEvaluation *data);


/**
 * Take values from one description and add them to the other.
 *
 * @param[out] desc           a set of descriptive values relevant to
 *                            comparative analysis
 * @param[in]  desc_to_add    a set of descriptive values which will be added to
 *                            the first
 */
void pe_comp_result_summary_combine(PeCompResultSummary *desc,
                                    PeCompResultSummary *desc_to_add);

/**
 * Initialize to default values.
 *
 * @param[in] results    the descriptive values
 */
void pe_comp_result_summary_init(PeCompResultSummary *desc);

/**
 * Take the counts from one description and add them to a larger aggregate set
 * of counts.
 *
 * @param[out] desc           a set of aggregate descriptive counts
 * @param[in]  desc_to_add    a smaller set of counts which will be added to the
 *                            first set
 */
//void pe_classification_characteristics_combine
void pe_comp_result_desc_combine(PeCompResultDesc *desc,
                                 PeCompResultDesc *desc_to_add);

/**
 * Initialize characteristics to default values.
 *
 * @param[in] characteristics    the counts/stats
 */
//void pe_classification_characteristics_init
void pe_comp_result_desc_init(PeCompResultDesc *desc);

/**
 *
 */
void pe_pre_analysis(AgnGeneLocus *locus, PeAnalysisData *data);

/**
 *
 */
void pe_post_analysis(AgnGeneLocus *locus, PeAnalysisData *data);

#endif
