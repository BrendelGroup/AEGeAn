#ifndef AEGEAN_COMPARISON
#define AEGEAN_COMPARISON

#include <stdbool.h>
#include <stdio.h>
#include "core/types_api.h"

/**
 * @module AgnComparison
 *
 * Data structures and functions related to comparative assessment of
 * gene/transcript annotations.
 */ //;

/**
 * @type This struct is used to aggregate counts and statistics regarding the
 * structural-level comparison (i.e., at the level of whole CDS segments, whole
 * exons, and whole UTRs) and analysis of gene structure. See header file for
 * details.
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
 * @type This struct is used to aggregate counts and statistics regarding the
 * nucleotide-level comparison and analysis of gene structure. See header file
 * for details.
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
 * @type This struct aggregates all the counts and stats that go into a
 * comparison, including structural-level and nucleotide-level counts and stats.
 * See header file for details.
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
};
typedef struct AgnComparison AgnComparison;

/**
 * @type This enumerated type refers to all the possible outcomes when
 * annotations from two different sources are compared:
 * ``AGN_COMP_CLASS_UNCLASSIFIED``, ``AGN_COMP_CLASS_PERFECT_MATCH``,
 * ``AGN_COMP_CLASS_MISLABELED``,   ``AGN_COMP_CLASS_CDS_MATCH``,
 * ``AGN_COMP_CLASS_EXON_MATCH``,   ``AGN_COMP_CLASS_UTR_MATCH``, and
 * ``AGN_COMP_CLASS_NON_MATCH``.
 */
enum AgnCompClassification
{
  AGN_COMP_CLASS_UNCLASSIFIED,
  AGN_COMP_CLASS_PERFECT_MATCH,
  AGN_COMP_CLASS_MISLABELED,
  AGN_COMP_CLASS_CDS_MATCH,
  AGN_COMP_CLASS_EXON_MATCH,
  AGN_COMP_CLASS_UTR_MATCH,
  AGN_COMP_CLASS_NON_MATCH
};
typedef enum AgnCompClassification AgnCompClassification;

/**
 * @type This struct contains various counts to be reported in the summary
 * report.
 */
struct AgnCompInfo
{
  GtUword num_loci;
  GtUword unique_refr_loci;
  GtUword unique_pred_loci;
  GtUword refr_genes;
  GtUword pred_genes;
  GtUword refr_transcripts;
  GtUword pred_transcripts;
  GtUword num_comparisons;
};
typedef struct AgnCompInfo AgnCompInfo;

/**
 * @type When reporting the results of a comparative analysis, it may be useful
 * to (as is done by ParsEval) show some basic information about clique pairs
 * that fall under each classification category. The counts in this struct are
 * necessary to calculate those summary characteristics.
 */
struct AgnCompClassDesc
{
  GtUword comparison_count;
  GtUword total_length;
  GtUword refr_cds_length;
  GtUword pred_cds_length;
  GtUword refr_exon_count;
  GtUword pred_exon_count;
};
typedef struct AgnCompClassDesc AgnCompClassDesc;

/**
 * @type This struct is used to aggregate descriptions for all of the
 * classification categories.
 */
struct AgnCompClassSummary
{
  AgnCompClassDesc perfect_matches;
  AgnCompClassDesc perfect_mislabeled;
  AgnCompClassDesc cds_matches;
  AgnCompClassDesc exon_matches;
  AgnCompClassDesc utr_matches;
  AgnCompClassDesc non_matches;
};
typedef struct AgnCompClassSummary AgnCompClassSummary;

/**
 * @type Aggregate various data related to comparison of annotations.
 */
struct AgnComparisonData
{
  AgnCompClassSummary summary;
  AgnCompInfo info;
  AgnComparison stats;
};
typedef struct AgnComparisonData AgnComparisonData;

/**
 * @function Function used to combine similarity stats from many different
 * comparisons into a single aggregate summary.
 */
void agn_comparison_aggregate(AgnComparison *agg_cmp, AgnComparison *cmp);

/**
 * @function Add counts and stats from ``data`` to ``agg_data``.
 */
void agn_comparison_data_aggregate(AgnComparisonData *agg_data,
                                   AgnComparisonData *data);

/**
 * @function Initialize counts and stats to default values.
 */
void agn_comparison_data_init(AgnComparisonData *data);

/**
 * @function Initialize comparison stats to default values.
 */
void agn_comparison_init(AgnComparison *comparison);

/**
 * @function Print the comparison stats to the given file.
 */
void agn_comparison_print(AgnComparison *stats, FILE *outstream);

/**
 * @function Calculate stats from the given counts.
 */
void agn_comparison_resolve(AgnComparison *comparison);

/**
 * @function Returns true if c1 and c2 contain identical values, false
 * otherwise.
 */
bool agn_comparison_test(AgnComparison *c1, AgnComparison *c2);

/**
 * @function Add values from ``desc`` to ``agg_desc``.
 */
void agn_comp_class_desc_aggregate(AgnCompClassDesc *agg_desc,
                                   AgnCompClassDesc *desc);

/**
 * @function Initialize to default values.
 */
void agn_comp_class_desc_init(AgnCompClassDesc *desc);

/**
 * @function Add values from ``summ`` to ``agg_summ``.
 */
void agn_comp_class_summary_aggregate(AgnCompClassSummary *agg_summ,
                                      AgnCompClassSummary *summ);

/**
 * @function Initialize to default values.
 */
void agn_comp_class_summary_init(AgnCompClassSummary *summ);

/**
 * @function Add values from ``info`` to ``agg_info``.
 */
void agn_comp_info_aggregate(AgnCompInfo *agg_info, AgnCompInfo *info);

/**
 * @function Initialize to default values.
 */
void agn_comp_info_init(AgnCompInfo *info);

/**
 * @function Function used to combine similarity stats from many different
 * comparisons into a single aggregate summary.
 */
void agn_comp_stats_binary_aggregate(AgnCompStatsBinary *agg_stats,
                                     AgnCompStatsBinary *stats);

/**
 * @function Initialize comparison counts/stats to default values.
 */
void agn_comp_stats_binary_init(AgnCompStatsBinary *stats);

/**
 * @function Print the comparison stats to the given file.
 */
void agn_comp_stats_binary_print(AgnCompStatsBinary *stats, FILE *outstream);

/**
 * @function Calculate stats from the given counts.
 */
void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats);

/**
 * @function Returns true if s1 and s2 contain identical values, false
 * otherwise.
 */
bool agn_comp_stats_binary_test(AgnCompStatsBinary *s1, AgnCompStatsBinary *s2);

/**
 * @function Function used to combine similarity stats from many different
 * comparisons into a single aggregate summary.
 */
void agn_comp_stats_scaled_aggregate(AgnCompStatsScaled *agg_stats,
                                     AgnCompStatsScaled *stats);

/**
 * @function Initialize comparison counts/stats to default values.
 */
void agn_comp_stats_scaled_init(AgnCompStatsScaled *stats);

/**
 * @function Print the comparison stats to the given file.
 */
void agn_comp_stats_scaled_print(AgnCompStatsScaled *stats, FILE *outstream);

/**
 * @function Calculate stats from the given counts.
 */
void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats);

/**
 * @function Returns true if s1 and s2 contain identical values, false
 * otherwise.
 */
bool agn_comp_stats_scaled_test(AgnCompStatsScaled *s1, AgnCompStatsScaled *s2);

#endif
