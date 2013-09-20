#ifndef AEGEAN_COMPARATIVE_EVALUATION
#define AEGEAN_COMPARATIVE_EVALUATION

#include <stdio.h>
#include "genometools.h"
#include "AgnLogger.h"

/**
 * @module AgnComparEval
 *
 * A collection of data structures and functions for comparative evaluation of
 * annotations.
 */ //;

/**
 * @type This struct is used to aggregate counts and statistics regarding the
 * nucleotide-level comparison and analysis of gene structure.
 */
struct AgnCompStatsScaled
{
  GtUword tp;
  GtUword fn;
  GtUword fp;
  GtUword tn;
  double  mc;
  double  cc;
  double  sn;
  double  sp;
  double  f1;
  double  ed;
  char    mcs[7];
  char    ccs[7];
  char    sns[7];
  char    sps[7];
  char    f1s[7];
  char    eds[16];
};
typedef struct AgnCompStatsScaled AgnCompStatsScaled;

/**
 * @type This struct is used to aggregate counts and statistics regarding the
 * structural-level comparison (i.e., at the level of whole CDS segments, whole
 * exons, and whole UTRs) and analysis of gene structure.
 */
struct AgnCompStatsBinary
{
  GtUword correct;
  GtUword missing;
  GtUword wrong;
  double  sn;
  double  sp;
  double  f1;
  double  ed;
  char    sns[7];
  char    sps[7];
  char    f1s[7];
  char    eds[16];
};
typedef struct AgnCompStatsBinary AgnCompStatsBinary;

/**
 * @type This struct contains various counts to be reported in the summary
 * report.
 */
struct AgnCompSummary
{
  unsigned int unique_refr;
  unsigned int unique_pred;
  GtUword      refr_genes;
  GtUword      pred_genes;
  GtUword      refr_transcripts;
  GtUword      pred_transcripts;
  GtUword      num_loci;
  unsigned int num_comparisons;
  unsigned int num_perfect;
  unsigned int num_mislabeled;
  unsigned int num_cds_match;
  unsigned int num_exon_match;
  unsigned int num_utr_match;
  unsigned int non_match;
};
typedef struct AgnCompSummary AgnCompSummary;

/**
 * @type This struct aggregates all the counts and stats that go into a
 * comparison, including structural-level and nucleotide-level counts and stats.
 */
struct AgnComparison
{
  AgnCompStatsScaled cds_nuc_stats;
  AgnCompStatsScaled utr_nuc_stats;
  AgnCompStatsBinary cds_struc_stats;
  AgnCompStatsBinary exon_struc_stats;
  AgnCompStatsBinary utr_struc_stats;
  GtUword overall_matches;
  GtUword overall_length;
  double overall_identity;
  double tolerance;
};
typedef struct AgnComparison AgnComparison;

/**
 * @type Each transcript clique pair that is compared is classified as one of
 * the following: perfect match; perfect match with mislabeled UTRs; CDS match;
 * exon structure match; UTR structure match; non-match. When reporting the
 * results of a comparative analysis, it may be useful to (as is done by
 * ParsEval) show some basic information about clique pairs that fall under each
 * classification category. The counts in this struct are necessary to calculate
 * those summary characteristics.
 */
struct AgnCompResultDesc
{
  GtUword total_length;
  GtUword transcript_count;
  GtUword refr_cds_length;
  GtUword pred_cds_length;
  GtUword refr_exon_count;
  GtUword pred_exon_count;
};
typedef struct AgnCompResultDesc AgnCompResultDesc;

/**
 * @type This struct is used to aggregate characteristics for all of the
 * classification categories.
 */
struct AgnCompResultSummary
{
  AgnCompResultDesc perfect_matches;
  AgnCompResultDesc perfect_mislabeled;
  AgnCompResultDesc cds_matches;
  AgnCompResultDesc exon_matches;
  AgnCompResultDesc utr_matches;
  AgnCompResultDesc non_matches;
};
typedef struct AgnCompResultSummary AgnCompResultSummary;

/**
 * @type This struct provides a convenient way to manage the counts, stats, and
 * results corresponding to one or more comparisons.
 */
struct AgnCompEvaluation
{
  AgnCompSummary counts;
  AgnComparison stats;
  AgnCompResultSummary results;
};
typedef struct AgnCompEvaluation AgnCompEvaluation;

/**
 * @type This struct contains a list of filters to be used in determining which
 * loci should be included/excluded in a comparative analysis.
 */
struct AgnCompareFilters
{
  GtUword LocusLengthUpperLimit;
  GtUword LocusLengthLowerLimit;
  GtUword MinReferenceGeneModels;
  GtUword MaxReferenceGeneModels;
  GtUword MinPredictionGeneModels;
  GtUword MaxPredictionGeneModels;
  GtUword MinReferenceTranscriptModels;
  GtUword MaxReferenceTranscriptModels;
  GtUword MinPredictionTranscriptModels;
  GtUword MaxPredictionTranscriptModels;
  GtUword MinTranscriptsPerReferenceGeneModel;
  GtUword MaxTranscriptsPerReferenceGeneModel;
  GtUword MinTranscriptsPerPredictionGeneModel;
  GtUword MaxTranscriptsPerPredictionGeneModel;
  GtUword MinReferenceExons;
  GtUword MaxReferenceExons;
  GtUword MinPredictionExons;
  GtUword MaxPredictionExons;
  GtUword MinReferenceCDSLength;
  GtUword MaxReferenceCDSLength;
  GtUword MinPredictionCDSLength;
  GtUword MaxPredictionCDSLength;
};
typedef struct AgnCompareFilters AgnCompareFilters;

/**
 * @function Take values from one data set and add them to the other.
 */
void agn_comp_evaluation_combine(AgnCompEvaluation *data,
                                 AgnCompEvaluation *data_to_add);

/**
 * @function Initialize to default values.
 */
void agn_comp_evaluation_init(AgnCompEvaluation *data);


/**
 * @function Take values from one description and add them to the other.
 */
void agn_comp_result_summary_combine(AgnCompResultSummary *desc,
                                     AgnCompResultSummary *desc_to_add);

/**
 * @function Initialize to default values.
 */
void agn_comp_result_summary_init(AgnCompResultSummary *desc);

/**
 * @function Take the counts from one description and add them to a larger
 * aggregate set of counts.
 */
void agn_comp_result_desc_combine(AgnCompResultDesc *desc,
                                  AgnCompResultDesc *desc_to_add);

/**
 * @function Initialize characteristics to default values.
 */
void agn_comp_result_desc_init(AgnCompResultDesc *desc);

/**
 * @function Take one set of values and add them to the other.
 */
void agn_comp_summary_combine(AgnCompSummary *s1, AgnCompSummary *s2);

/**
 * @function Initialize to default values.
 */
void agn_comp_summary_init(AgnCompSummary *summary);

/**
 * @function Take stats from one comparison and add them to the other.
 */
void agn_comparison_combine(AgnComparison *c1, AgnComparison *c2);

/**
 * @function Initialize comparison stats to default values.
 */
void agn_comparison_init(AgnComparison *comparison);

/**
 * @function Initialize filters to default values.
 */
void agn_compare_filters_init(AgnCompareFilters *filters);

/**
 * @function Parse the filter configuration file (from ``instream``) to set the
 * filters appropriately.
 */
void agn_compare_filters_parse(AgnCompareFilters *filters, FILE *instream,
                               AgnLogger *logger);

/**
 * @function Initialize comparison counts/stats to default values.
 */
void agn_comp_stats_binary_init(AgnCompStatsBinary *stats);

/**
 * @function Calculate stats from the given counts.
 */
void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats);

/**
 * @function Initialize comparison counts/stats to default values.
 */
void agn_comp_stats_scaled_init(AgnCompStatsScaled *stats);

/**
 * @function Calculate stats from the given counts.
 */
void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats);

#endif
