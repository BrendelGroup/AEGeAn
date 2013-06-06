#ifndef AEGEAN_COMPARATIVE_EVALUATION
#define AEGEAN_COMPARATIVE_EVALUATION

#include <stdio.h>

/**
 * This struct is used to aggregate counts and statistics regarding the
 * nucleotide-level comparison and analysis of gene structure.
 */
typedef struct
{
  unsigned long tp;
  unsigned long fn;
  unsigned long fp;
  unsigned long tn;
  double        mc;
  double        cc;
  double        sn;
  double        sp;
  double        f1;
  double        ed;
  char          mcs[7];
  char          ccs[7];
  char          sns[7];
  char          sps[7];
  char          f1s[7];
  char          eds[16];
} AgnCompStatsScaled;

/**
 * This struct is used to aggregate counts and statistics regarding the
 * structural-level comparison (i.e., at the level of whole CDS segments, whole
 * exons, and whole UTRs) and analysis of gene structure.
 */
typedef struct
{
  unsigned long correct;
  unsigned long missing;
  unsigned long wrong;
  double        sn;
  double        sp;
  double        f1;
  double        ed;
  char          sns[7];
  char          sps[7];
  char          f1s[7];
  char          eds[16];
} AgnCompStatsBinary;

/**
 * This struct contains various counts to be reported in the summary report.
 */
typedef struct
{
  unsigned int unique_refr;
  unsigned int unique_pred;
  unsigned long refr_genes;
  unsigned long pred_genes;
  unsigned long refr_transcripts;
  unsigned long pred_transcripts;
  unsigned long num_loci;
  unsigned int num_comparisons;
  unsigned int num_perfect;
  unsigned int num_mislabeled;
  unsigned int num_cds_match;
  unsigned int num_exon_match;
  unsigned int num_utr_match;
  unsigned int non_match;
} AgnComparisonCounts;

/**
 * This struct aggregates all the counts and stats that go into a comparison,
 * including structural-level and nucleotide-level counts and stats.
 */
typedef struct
{
  AgnCompStatsScaled cds_nuc_stats;
  AgnCompStatsScaled utr_nuc_stats;
  AgnCompStatsBinary cds_struc_stats;
  AgnCompStatsBinary exon_struc_stats;
  AgnCompStatsBinary utr_struc_stats;
  unsigned long overall_matches;
  unsigned long overall_length;
  double overall_identity;
  double tolerance;
} AgnComparisonStats;

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
} AgnCompareClassDescription;

/**
 * This struct is used to aggregate characteristics for all of the
 * classification categories.
 */
typedef struct
{
  AgnCompareClassDescription perfect_matches;
  AgnCompareClassDescription perfect_mislabeled;
  AgnCompareClassDescription cds_matches;
  AgnCompareClassDescription exon_matches;
  AgnCompareClassDescription utr_matches;
  AgnCompareClassDescription non_matches;
} AgnCompareClassAggregateDescription;

/**
 * This struct provides a convenient way to manage the counts, stats, and
 * results corresponding to one or more comparisons.
 */
typedef struct
{
  AgnComparisonCounts counts;
  AgnComparisonStats stats;
  AgnCompareClassAggregateDescription results;
} AgnSummaryData;

/**
 * This struct contains a list of filters to be used in determining which loci
 * should be included/excluded in a comparative analysis.
 */
typedef struct
{
  unsigned long LocusLengthUpperLimit;
  unsigned long LocusLengthLowerLimit;
  unsigned long MinReferenceGeneModels;
  unsigned long MaxReferenceGeneModels;
  unsigned long MinPredictionGeneModels;
  unsigned long MaxPredictionGeneModels;
  unsigned long MinReferenceTranscriptModels;
  unsigned long MaxReferenceTranscriptModels;
  unsigned long MinPredictionTranscriptModels;
  unsigned long MaxPredictionTranscriptModels;
  unsigned long MinTranscriptsPerReferenceGeneModel;
  unsigned long MaxTranscriptsPerReferenceGeneModel;
  unsigned long MinTranscriptsPerPredictionGeneModel;
  unsigned long MaxTranscriptsPerPredictionGeneModel;
  unsigned long MinReferenceExons;
  unsigned long MaxReferenceExons;
  unsigned long MinPredictionExons;
  unsigned long MaxPredictionExons;
  unsigned long MinReferenceCDSLength;
  unsigned long MaxReferenceCDSLength;
  unsigned long MinPredictionCDSLength;
  unsigned long MaxPredictionCDSLength;
} AgnCompareFilters;

/**
 * Take values from one description and add them to the other.
 *
 * @param[out] desc           a set of descriptive values relevant to
 *                            comparative analysis
 * @param[in]  desc_to_add    a set of descriptive values which will be added to
 *                            the first
 */
void agn_compare_class_agg_desc_combine
(
  AgnCompareClassAggregateDescription *desc,
  AgnCompareClassAggregateDescription *desc_to_add
);

/**
 * Initialize to default values.
 *
 * @param[in] results    the descriptive values
 */
void agn_compare_class_agg_desc_init(AgnCompareClassAggregateDescription *desc);

/**
 * Take the counts from one description and add them to a larger aggregate set
 * of counts.
 *
 * @param[out] desc           a set of aggregate descriptive counts
 * @param[in]  desc_to_add    a smaller set of counts which will be added to the
 *                            first set
 */
//void pe_classification_characteristics_combine
void agn_compare_class_description_combine
(
  AgnCompareClassDescription *desc,
  AgnCompareClassDescription *desc_to_add
);

/**
 * Initialize characteristics to default values.
 *
 * @param[in] characteristics    the counts/stats
 */
//void pe_classification_characteristics_init
void agn_compare_class_description_init(AgnCompareClassDescription *desc);

/**
 * Take one set of counts and add them to the other.
 *
 * @param[out] counts           a set of counts relevant to comparative analysis
 * @param[in]  counts_to_add    a smaller set of counts which will be added to
 *                              the first set
 */
void agn_comparison_counts_combine( AgnComparisonCounts *counts,
                                    AgnComparisonCounts *counts_to_add );

/**
 * Initialize counts to default values.
 *
 * @param[in] counts    the counts
 */
void agn_comparison_counts_init(AgnComparisonCounts *counts);

/**
 * Take one set of stats and add them to the other.
 *
 * @param[out] stats           a set of stats relevant to comparative analysis
 * @param[in]  stats_to_add    a smaller set of stats which will be added to
 *                             the first set
 */
void agn_comparison_stats_combine( AgnComparisonStats *stats,
                                   AgnComparisonStats *stats_to_add );

/**
 * Initialize comparison stats to default values.
 *
 * @param[in] stats    the comparison stats
 */
void agn_comparison_stats_init(AgnComparisonStats *stats);

/**
 * Initialize filters to default values.
 *
 * @param[in] filters    the filters
 */
void agn_compare_filters_init(AgnCompareFilters *filters);

/**
 * Parse the filter configuration file to set the filters appropriately.
 *
 * @param[out] filters     set of filters to be set
 * @param[in]  instream    pointer to the filter config file
 */
void agn_compare_filters_parse(AgnCompareFilters *filters, FILE *instream);

/**
 * Calculate stats from the given counts.
 *
 * @param[in] stats    pointer to the counts and the stats
 */
void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats);

/**
 * Calculate stats from the given counts.
 *
 * @param[in] stats    pointer to the counts and the stats
 */
void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats);

/**
 * Take values from one data set and add them to the other.
 *
 * @param[out] data           stats, counts, and descriptions relevant to
 *                            comparative analysis
 * @param[in]  data_to_add    a data set which will be added to the first
 */
void agn_summary_data_combine(AgnSummaryData *data, AgnSummaryData *data_to_add);

/**
 * Initialize to default values.
 *
 * @param[in] data    the data
 */
void agn_summary_data_init(AgnSummaryData *data);

#endif
