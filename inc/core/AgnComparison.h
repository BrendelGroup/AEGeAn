#ifndef AEGEAN_COMPARISON
#define AEGEAN_COMPARISON

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
 * @function Function used to combine similarity stats from many different
 * comparisons into a single aggregate summary. 
 */
void agn_comparison_aggregate(AgnComparison *agg_cmp, AgnComparison *cmp);

/**
 * @function Initialize comparison stats to default values.
 */
void agn_comparison_init(AgnComparison *comparison);

/**
 * @function Calculate stats from the given counts.
 */
void agn_comparison_resolve(AgnComparison *comparison);

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
 * @function Calculate stats from the given counts.
 */
void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats);

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
 * @function Calculate stats from the given counts.
 */
void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats);

#endif
