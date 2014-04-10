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

void agn_comparison_data_aggregate(AgnComparisonData *agg_data,
                                   AgnComparisonData *data)
{
  agn_comp_class_summary_aggregate(&agg_data->summary, &data->summary);
  agn_comp_info_aggregate(&agg_data->info, &data->info);
  agn_comparison_aggregate(&agg_data->stats, &data->stats);
}

void agn_comparison_data_init(AgnComparisonData *data)
{
  agn_comp_class_summary_init(&data->summary);
  agn_comp_info_init(&data->info);
  agn_comparison_init(&data->stats);
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

void agn_comparison_print(AgnComparison *stats, FILE *outstream)
{
  fprintf(outstream, "CDS nucleotide stats: ");
  agn_comp_stats_scaled_print(&stats->cds_nuc_stats, outstream);
  fprintf(outstream, "UTR nucleotide stats: ");
  agn_comp_stats_scaled_print(&stats->utr_nuc_stats, outstream);
  fprintf(outstream, "CDS structure stats: ");
  agn_comp_stats_binary_print(&stats->cds_struc_stats, outstream);
  fprintf(outstream, "exon structure stats: ");
  agn_comp_stats_binary_print(&stats->exon_struc_stats, outstream);
  fprintf(outstream, "UTR structure stats: ");
  agn_comp_stats_binary_print(&stats->utr_struc_stats, outstream);
  fprintf(outstream, "Overall matches: %lu, overall length: %lu\n",
          stats->overall_matches, stats->overall_length);
}

void agn_comparison_resolve(AgnComparison *comparison)
{
  agn_comp_stats_scaled_resolve(&comparison->cds_nuc_stats);
  agn_comp_stats_scaled_resolve(&comparison->utr_nuc_stats);
  agn_comp_stats_binary_resolve(&comparison->cds_struc_stats);
  agn_comp_stats_binary_resolve(&comparison->exon_struc_stats);
  agn_comp_stats_binary_resolve(&comparison->utr_struc_stats);
}

bool agn_comparison_test(AgnComparison *c1, AgnComparison *c2)
{
  return agn_comp_stats_scaled_test(&c1->cds_nuc_stats, &c2->cds_nuc_stats) &&
         agn_comp_stats_scaled_test(&c1->utr_nuc_stats, &c2->utr_nuc_stats) &&
     agn_comp_stats_binary_test(&c1->cds_struc_stats, &c2->cds_struc_stats) &&
     agn_comp_stats_binary_test(&c1->exon_struc_stats, &c2->exon_struc_stats) &&
     agn_comp_stats_binary_test(&c1->utr_struc_stats, &c2->utr_struc_stats);
}

void agn_comp_class_desc_aggregate(AgnCompClassDesc *agg_desc,
                                   AgnCompClassDesc *desc)
{
  agg_desc->comparison_count = desc->comparison_count;
  agg_desc->total_length = desc->total_length;
  agg_desc->refr_cds_length = desc->refr_cds_length;
  agg_desc->pred_cds_length = desc->pred_cds_length;
  agg_desc->refr_exon_count = desc->refr_exon_count;
  agg_desc->pred_exon_count = desc->pred_exon_count;
}

void agn_comp_class_desc_init(AgnCompClassDesc *desc)
{
  desc->comparison_count = 0;
  desc->total_length = 0;
  desc->refr_cds_length = 0;
  desc->pred_cds_length = 0;
  desc->refr_exon_count = 0;
  desc->pred_exon_count = 0;
}

void agn_comp_class_summary_aggregate(AgnCompClassSummary *agg_summ,
                                      AgnCompClassSummary *summ)
{
  agn_comp_class_desc_aggregate(&agg_summ->perfect_matches,
                                &summ->perfect_matches);
  agn_comp_class_desc_aggregate(&agg_summ->perfect_mislabeled,
                                &summ->perfect_mislabeled);
  agn_comp_class_desc_aggregate(&agg_summ->cds_matches, &summ->cds_matches);
  agn_comp_class_desc_aggregate(&agg_summ->exon_matches, &summ->exon_matches);
  agn_comp_class_desc_aggregate(&agg_summ->utr_matches, &summ->utr_matches);
  agn_comp_class_desc_aggregate(&agg_summ->non_matches, &summ->non_matches);
}

void agn_comp_class_summary_init(AgnCompClassSummary *summ)
{
  agn_comp_class_desc_init(&summ->perfect_matches);
  agn_comp_class_desc_init(&summ->perfect_mislabeled);
  agn_comp_class_desc_init(&summ->cds_matches);
  agn_comp_class_desc_init(&summ->exon_matches);
  agn_comp_class_desc_init(&summ->utr_matches);
  agn_comp_class_desc_init(&summ->non_matches);
}

void agn_comp_info_aggregate(AgnCompInfo *agg_info, AgnCompInfo *info)
{
  agg_info->num_loci = info->num_loci;
  agg_info->unique_refr_loci = info->unique_refr_loci;
  agg_info->unique_pred_loci = info->unique_pred_loci;
  agg_info->refr_genes = info->refr_genes;
  agg_info->pred_genes = info->pred_genes;
  agg_info->refr_transcripts = info->refr_transcripts;
  agg_info->pred_transcripts = info->pred_transcripts;
  agg_info->num_comparisons = info->num_comparisons;
}

void agn_comp_info_init(AgnCompInfo *info)
{
  info->num_loci = 0;
  info->unique_refr_loci = 0;
  info->unique_pred_loci = 0;
  info->refr_genes = 0;
  info->pred_genes = 0;
  info->refr_transcripts = 0;
  info->pred_transcripts = 0;
  info->num_comparisons = 0;
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

void agn_comp_stats_binary_print(AgnCompStatsBinary *stats, FILE *outstream)
{
  fprintf(outstream, "%lu, %lu, %lu, %.6lf, %.6lf, %.6lf, %.6lf\n",
          stats->correct, stats->missing, stats->wrong, stats->sn, stats->sp,
          stats->f1, stats->ed);
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

bool agn_comp_stats_binary_test(AgnCompStatsBinary *s1, AgnCompStatsBinary *s2)
{
  bool counts = s1->correct == s2->correct && s1->missing == s2->missing &&
                s1->wrong == s2->wrong;
  bool sn = (isnan(s1->sn) && isnan(s2->sn)) || fabs(s1->sn - s2->sn) < 0.0001;
  bool sp = (isnan(s1->sp) && isnan(s2->sp)) || fabs(s1->sp - s2->sp) < 0.0001;
  bool f1 = (isnan(s1->f1) && isnan(s2->f1)) || fabs(s1->f1 - s2->f1) < 0.0001;
  bool ed = (isnan(s1->ed) && isnan(s2->ed)) || fabs(s1->ed - s2->ed) < 0.0001;
  return counts && sn && sp && f1 && ed;
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

void agn_comp_stats_scaled_print(AgnCompStatsScaled *stats, FILE *outstream)
{
  fprintf(outstream, "%lu, %lu, %lu, %lu, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, "
          "%.6lf\n", stats->tp, stats->fn, stats->fp, stats->tn, stats->mc,
          stats->cc, stats->sn, stats->sp, stats->f1, stats->ed);
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
  stats->sp = tp/(tp + fp);
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

bool agn_comp_stats_scaled_test(AgnCompStatsScaled *s1, AgnCompStatsScaled *s2)
{
  bool counts = s1->tp == s2->tp && s1->fn == s2->fn && s1->fp == s2->fp &&
                s1->tn == s2->tn;
  bool mc = (isnan(s1->mc) && isnan(s2->mc)) || fabs(s1->mc - s2->mc) < 0.0001;
  bool cc = (isnan(s1->cc) && isnan(s2->cc)) || fabs(s1->cc - s2->cc) < 0.0001;
  bool sn = (isnan(s1->sn) && isnan(s2->sn)) || fabs(s1->sn - s2->sn) < 0.0001;
  bool sp = (isnan(s1->sp) && isnan(s2->sp)) || fabs(s1->sp - s2->sp) < 0.0001;
  bool f1 = (isnan(s1->f1) && isnan(s2->f1)) || fabs(s1->f1 - s2->f1) < 0.0001;
  bool ed = (isnan(s1->ed) && isnan(s2->ed)) || fabs(s1->ed - s2->ed) < 0.0001;
  return counts && mc && cc && sn && sp && f1 && ed;
}
