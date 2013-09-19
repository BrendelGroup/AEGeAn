#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "AgnComparEval.h"

void agn_comp_evaluation_combine(AgnCompEvaluation *data,
                                 AgnCompEvaluation *data_to_add)
{
  agn_comp_summary_combine(&data->counts, &data_to_add->counts);
  agn_comparison_combine(&data->stats, &data_to_add->stats);
  agn_comp_result_summary_combine(&data->results, &data_to_add->results);
}

void agn_comp_evaluation_init(AgnCompEvaluation *data)
{
  agn_comp_summary_init(&data->counts);
  agn_comparison_init(&data->stats);
  agn_comp_result_summary_init(&data->results);
}

void agn_comp_result_summary_combine(AgnCompResultSummary *desc,
                                     AgnCompResultSummary *desc_to_add)
{
  agn_comp_result_desc_combine(&desc->perfect_matches,
                               &desc_to_add->perfect_matches);
  agn_comp_result_desc_combine(&desc->perfect_mislabeled,
                               &desc_to_add->perfect_mislabeled);
  agn_comp_result_desc_combine(&desc->cds_matches,
                               &desc_to_add->cds_matches);
  agn_comp_result_desc_combine(&desc->exon_matches,
                               &desc_to_add->exon_matches);
  agn_comp_result_desc_combine(&desc->utr_matches,
                               &desc_to_add->utr_matches);
  agn_comp_result_desc_combine(&desc->non_matches,
                               &desc_to_add->non_matches);
}

void agn_comp_result_summary_init(AgnCompResultSummary *desc)
{
  agn_comp_result_desc_init(&desc->perfect_matches);
  agn_comp_result_desc_init(&desc->perfect_mislabeled);
  agn_comp_result_desc_init(&desc->cds_matches);
  agn_comp_result_desc_init(&desc->exon_matches);
  agn_comp_result_desc_init(&desc->utr_matches);
  agn_comp_result_desc_init(&desc->non_matches);
}

void agn_comp_result_desc_combine(AgnCompResultDesc *desc,
                                  AgnCompResultDesc *desc_to_add)
{
  desc->total_length     += desc_to_add->total_length;
  desc->transcript_count += desc_to_add->transcript_count;
  desc->refr_cds_length  += desc_to_add->refr_cds_length;
  desc->pred_cds_length  += desc_to_add->pred_cds_length;
  desc->refr_exon_count  += desc_to_add->refr_exon_count;
  desc->pred_exon_count  += desc_to_add->pred_exon_count;
}

void agn_comp_result_desc_init(AgnCompResultDesc *desc)
{
  desc->total_length = 0;
  desc->transcript_count = 0;
  desc->refr_cds_length = 0;
  desc->pred_cds_length = 0;
  desc->refr_exon_count = 0;
  desc->pred_exon_count = 0;
}

void agn_comp_summary_combine(AgnCompSummary *s1, AgnCompSummary *s2)
{
  s1->unique_refr      += s2->unique_refr;
  s1->unique_pred      += s2->unique_pred;
  s1->refr_genes       += s2->refr_genes;
  s1->pred_genes       += s2->pred_genes;
  s1->refr_transcripts += s2->refr_transcripts;
  s1->pred_transcripts += s2->pred_transcripts;
  s1->num_loci         += s2->num_loci;
  s1->num_comparisons  += s2->num_comparisons;
  s1->num_perfect      += s2->num_perfect;
  s1->num_mislabeled   += s2->num_mislabeled;
  s1->num_cds_match    += s2->num_cds_match;
  s1->num_exon_match   += s2->num_exon_match;
  s1->num_utr_match    += s2->num_utr_match;
  s1->non_match        += s2->non_match;
}

void agn_comp_summary_init(AgnCompSummary *summary)
{
  summary->unique_refr = 0;
  summary->unique_pred = 0;
  summary->refr_genes = 0;
  summary->pred_genes = 0;
  summary->refr_transcripts = 0;
  summary->pred_transcripts = 0;
  summary->num_loci = 0;
  summary->num_comparisons = 0;
  summary->num_perfect = 0;
  summary->num_mislabeled = 0;
  summary->num_cds_match = 0;
  summary->num_exon_match = 0;
  summary->num_utr_match = 0;
  summary->non_match = 0;
}

void agn_comparison_combine(AgnComparison *c1, AgnComparison *c2)
{
  c1->cds_nuc_stats.tp += c2->cds_nuc_stats.tp;
  c1->cds_nuc_stats.fn += c2->cds_nuc_stats.fn;
  c1->cds_nuc_stats.fp += c2->cds_nuc_stats.fp;
  c1->cds_nuc_stats.tn += c2->cds_nuc_stats.tn;

  c1->utr_nuc_stats.tp += c2->utr_nuc_stats.tp;
  c1->utr_nuc_stats.fn += c2->utr_nuc_stats.fn;
  c1->utr_nuc_stats.fp += c2->utr_nuc_stats.fp;
  c1->utr_nuc_stats.tn += c2->utr_nuc_stats.tn;

  c1->cds_struc_stats.correct += c2->cds_struc_stats.correct;
  c1->cds_struc_stats.missing += c2->cds_struc_stats.missing;
  c1->cds_struc_stats.wrong   += c2->cds_struc_stats.wrong  ;

  c1->exon_struc_stats.correct += c2->exon_struc_stats.correct;
  c1->exon_struc_stats.missing += c2->exon_struc_stats.missing;
  c1->exon_struc_stats.wrong   += c2->exon_struc_stats.wrong  ;

  c1->utr_struc_stats.correct += c2->utr_struc_stats.correct;
  c1->utr_struc_stats.missing += c2->utr_struc_stats.missing;
  c1->utr_struc_stats.wrong   += c2->utr_struc_stats.wrong  ;

  c1->overall_matches += c2->overall_matches;
  c1->overall_length += c2->overall_length;
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
  comparison->overall_identity = 0.0;
  comparison->tolerance        = 0.0;
}

void agn_compare_filters_init(AgnCompareFilters *filters)
{
  filters->LocusLengthUpperLimit = 0;
  filters->LocusLengthLowerLimit = 0;
  filters->MinReferenceGeneModels = 0;
  filters->MaxReferenceGeneModels = 0;
  filters->MinPredictionGeneModels = 0;
  filters->MaxPredictionGeneModels = 0;
  filters->MinReferenceTranscriptModels = 0;
  filters->MaxReferenceTranscriptModels = 0;
  filters->MinPredictionTranscriptModels = 0;
  filters->MaxPredictionTranscriptModels = 0;
  filters->MinTranscriptsPerReferenceGeneModel = 0;
  filters->MaxTranscriptsPerReferenceGeneModel = 0;
  filters->MinTranscriptsPerPredictionGeneModel = 0;
  filters->MaxTranscriptsPerPredictionGeneModel = 0;
  filters->MinReferenceExons = 0;
  filters->MaxReferenceExons = 0;
  filters->MinPredictionExons = 0;
  filters->MaxPredictionExons = 0;
  filters->MinReferenceCDSLength = 0;
  filters->MaxReferenceCDSLength = 0;
  filters->MinPredictionCDSLength = 0;
  filters->MaxPredictionCDSLength = 0;
}

void agn_compare_filters_parse(AgnCompareFilters *filters, FILE *instream,
                               AgnLogger *logger)
{
  char buffer[256];
  while(fgets(buffer, 255, instream))
  {
    if(strlen(buffer) > 0 && buffer[0] != '\n' && buffer[0] != '#')
    {
      char *key = strtok(buffer, " :\n");
      char *value = strtok(NULL, " :\n");

      if(strcmp(key, "LocusLengthUpperLimit") == 0)
      {
        GtUword val = atol(value);
        filters->LocusLengthUpperLimit = val;
      }
      else if(strcmp(key, "LocusLengthLowerLimit") == 0)
      {
        GtUword val = atol(value);
        filters->LocusLengthLowerLimit = val;
      }
      else if(strcmp(key, "MinReferenceGeneModels") == 0)
      {
        GtUword val = atol(value);
        filters->MinReferenceGeneModels = val;
      }
      else if(strcmp(key, "MaxReferenceGeneModels") == 0)
      {
        GtUword val = atol(value);
        filters->MaxReferenceGeneModels = val;
      }
      else if(strcmp(key, "MinPredictionGeneModels") == 0)
      {
        GtUword val = atol(value);
        filters->MinPredictionGeneModels = val;
      }
      else if(strcmp(key, "MaxPredictionGeneModels") == 0)
      {
        GtUword val = atol(value);
        filters->MaxPredictionGeneModels = val;
      }
      else if(strcmp(key, "MinReferenceTranscriptModels") == 0)
      {
        GtUword val = atol(value);
        filters->MinReferenceTranscriptModels = val;
      }
      else if(strcmp(key, "MaxReferenceTranscriptModels") == 0)
      {
        GtUword val = atol(value);
        filters->MaxReferenceTranscriptModels = val;
      }
      else if(strcmp(key, "MinPredictionTranscriptModels") == 0)
      {
        GtUword val = atol(value);
        filters->MinPredictionTranscriptModels = val;
      }
      else if(strcmp(key, "MaxPredictionTranscriptModels") == 0)
      {
        GtUword val = atol(value);
        filters->MaxPredictionTranscriptModels = val;
      }
      else if(strcmp(key, "MinTranscriptsPerReferenceGeneModel") == 0)
      {
        GtUword val = atol(value);
        filters->MinTranscriptsPerReferenceGeneModel = val;
      }
      else if(strcmp(key, "MaxTranscriptsPerReferenceGeneModel") == 0)
      {
        GtUword val = atol(value);
        filters->MaxTranscriptsPerReferenceGeneModel = val;
      }
      else if(strcmp(key, "MinTranscriptsPerPredictionGeneModel") == 0)
      {
        GtUword val = atol(value);
        filters->MinTranscriptsPerPredictionGeneModel = val;
      }
      else if(strcmp(key, "MaxTranscriptsPerPredictionGeneModel") == 0)
      {
        GtUword val = atol(value);
        filters->MaxTranscriptsPerPredictionGeneModel = val;
      }
      else if(strcmp(key, "MinReferenceExons") == 0)
      {
        GtUword val = atol(value);
        filters->MinReferenceExons = val;
      }
      else if(strcmp(key, "MaxReferenceExons") == 0)
      {
        GtUword val = atol(value);
        filters->MaxReferenceExons = val;
      }
      else if(strcmp(key, "MinPredictionExons") == 0)
      {
        GtUword val = atol(value);
        filters->MinPredictionExons = val;
      }
      else if(strcmp(key, "MaxPredictionExons") == 0)
      {
        GtUword val = atol(value);
        filters->MaxPredictionExons = val;
      }
      else if(strcmp(key, "MinReferenceCDSLength") == 0)
      {
        GtUword val = atol(value);
        filters->MinReferenceCDSLength = val;
      }
      else if(strcmp(key, "MaxReferenceCDSLength") == 0)
      {
        GtUword val = atol(value);
        filters->MaxReferenceCDSLength = val;
      }
      else if(strcmp(key, "MinPredictionCDSLength") == 0)
      {
        GtUword val = atol(value);
        filters->MinPredictionCDSLength = val;
      }
      else if(strcmp(key, "MaxPredictionCDSLength") == 0)
      {
        GtUword val = atol(value);
        filters->MaxPredictionCDSLength = val;
      }
      else
      {
        agn_logger_log_error(logger, "unrecognized filter option '%s'\n", key);
        return;
      }
    }
  }
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
