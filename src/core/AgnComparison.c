#include <math.h>
#include <stdio.h>
#include "AgnComparison.h"

void agn_comparison_aggregate(AgnComparison *a, AgnComparison *b)
{
  agn_comp_stats_scaled_aggregate(&a->cds_nuc_stats,    &b->cds_nuc_stats);
  agn_comp_stats_scaled_aggregate(&a->utr_nuc_stats,    &b->utr_nuc_stats);
  agn_comp_stats_binary_aggregate(&a->cds_struc_stats,  &b->cds_struc_stats);
  agn_comp_stats_binary_aggregate(&a->exon_struc_stats, &b->exon_struc_stats);
  agn_comp_stats_binary_aggregate(&a->utr_struc_stats,  &b->utr_struc_stats);
  a->overall_matches += b->overall_matches;
  a->overall_length  += b->overall_length;
}

void agn_comparison_init(AgnComparison *comparison)
{
  agn_comp_stats_scaled_init(&comparison->cds_nuc_stats);
  agn_comp_stats_scaled_init(&comparison->utr_nuc_stats);
  agn_comp_stats_binary_init(&comparison->cds_struc_stats);
  agn_comp_stats_binary_init(&comparison->exon_struc_stats);
  agn_comp_stats_binary_init(&comparison->utr_struc_stats);

  comparison->overall_matches  = 0;
  comparison->overall_length   = 0;
}

void agn_comparison_resolve(AgnComparison *comparison)
{
  agn_comp_stats_scaled_resolve(&comparison->cds_nuc_stats);
  agn_comp_stats_scaled_resolve(&comparison->utr_nuc_stats);
  agn_comp_stats_binary_resolve(&comparison->cds_struc_stats);
  agn_comp_stats_binary_resolve(&comparison->exon_struc_stats);
  agn_comp_stats_binary_resolve(&comparison->utr_struc_stats);
}

void agn_comp_stats_binary_aggregate(AgnCompStatsBinary *a,
                                     AgnCompStatsBinary *b)
{
  a->correct += b->correct;
  a->missing += b->missing;
  a->wrong   += b->wrong;
}

void agn_comp_stats_binary_init(AgnCompStatsBinary *stats)
{
  stats->correct = 0;
  stats->missing = 0;
  stats->wrong   = 0;
  stats->sn      = 0.0;
  stats->sp      = 0.0;
  stats->f1      = 0.0;
  stats->ed      = 0.0;
}

void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats)
{
  double correct = (double)stats->correct;
  double missing = (double)stats->missing;
  double wrong = (double)stats->wrong;

  // Sensitivity
  stats->sn = correct / (correct + missing);
  if(isnan(stats->sn))
    sprintf(stats->sns, "--");
  else
    sprintf(stats->sns, "%.3lf", stats->sn);

  // Specificity
  stats->sp = correct / (correct + wrong);
  if(isnan(stats->sp))
    sprintf(stats->sps, "--");
  else
    sprintf(stats->sps, "%.3lf", stats->sp);

  // F1 score
  double precision = correct/(correct + wrong);
  double recall = correct/(correct + missing);
  stats->f1 = (2.0 * precision * recall) / (precision + recall);
  if(isnan(stats->f1))
    sprintf(stats->f1s, "--");
  else
    sprintf(stats->f1s, "%.3lf", stats->f1);

  // Annotation edit distance
  double congruency = (stats->sn+stats->sp)*0.5;
  stats->ed = 1 - congruency;
  if(isnan(stats->ed))
    sprintf(stats->eds, "--");
  else
    sprintf(stats->eds, "%.3lf", stats->ed);
}

void agn_comp_stats_scaled_aggregate(AgnCompStatsScaled *a,
                                     AgnCompStatsScaled *b)
{
  a->tp += b->tp;
  a->fn += b->fn;
  a->fp += b->fp;
  a->tn += b->tn;
}

void agn_comp_stats_scaled_init(AgnCompStatsScaled *stats)
{
  stats->tp = 0;
  stats->fn = 0;
  stats->fp = 0;
  stats->tn = 0;
  stats->mc = 0.0;
  stats->cc = 0.0;
  stats->sn = 0.0;
  stats->sp = 0.0;
  stats->f1 = 0.0;
  stats->ed = 0.0;
}

void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats)
{
  double tp = (double)stats->tp;
  double fn = (double)stats->fn;
  double fp = (double)stats->fp;
  double tn = (double)stats->tn;

  // Simple matching coefficient
  stats->mc = (tp + tn) / (tp + fn + fp + tn);
  if(isnan(stats->mc))
    sprintf(stats->mcs, "--");
  else
    sprintf(stats->mcs, "%.3lf", stats->mc);

  // Correlation coefficient
  stats->cc = ((tp*tn)-(fn*fp)) / pow(((tp+fn)*(tn+fp)*(tp+fp)*(tn+fn)), 0.5);
  if(isnan(stats->cc))
    sprintf(stats->ccs, "--");
  else
    sprintf(stats->ccs, "%.3lf", stats->cc);

  // Sensitivity
  stats->sn = tp/(tp + fn);
  if(isnan(stats->sn))
    sprintf(stats->sns, "--");
  else
    sprintf(stats->sns, "%.3lf", stats->sn);

  // Specificity
  stats->sp = tn/(tn + fp);
  if(isnan(stats->sp))
    sprintf(stats->sps, "--");
  else
    sprintf(stats->sps, "%.3lf", stats->sp);

  // F1 score
  double precision = tp/(tp + fp);
  double recall = tp/(tp + fn);
  stats->f1 = (2.0 * precision * recall) / (precision + recall);
  if(isnan(stats->f1))
    sprintf(stats->f1s, "--");
  else
    sprintf(stats->f1s, "%.3lf", stats->f1);

  // Annotation edit distance
  double congruency = (stats->sn+stats->sp)*0.5;
  stats->ed = 1 - congruency;
  if(isnan(stats->ed))
    sprintf(stats->eds, "--");
  else
    sprintf(stats->eds, "%.3lf", stats->ed);
}