#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "AgnComparEval.h"

void agn_comp_summary_combine( AgnCompSummary *counts,
                                    AgnCompSummary *counts_to_add )
{
  counts->unique_refr += counts_to_add->unique_refr;
  counts->unique_pred += counts_to_add->unique_pred;
  counts->refr_genes += counts_to_add->refr_genes;
  counts->pred_genes += counts_to_add->pred_genes;
  counts->refr_transcripts += counts_to_add->refr_transcripts;
  counts->pred_transcripts += counts_to_add->pred_transcripts;
  counts->num_loci += counts_to_add->num_loci;
  counts->num_comparisons += counts_to_add->num_comparisons;
  counts->num_perfect += counts_to_add->num_perfect;
  counts->num_mislabeled += counts_to_add->num_mislabeled;
  counts->num_cds_match += counts_to_add->num_cds_match;
  counts->num_exon_match += counts_to_add->num_exon_match;
  counts->num_utr_match += counts_to_add->num_utr_match;
  counts->non_match += counts_to_add->non_match;
}

void agn_comp_summary_init(AgnCompSummary *counts)
{
  counts->unique_refr = 0;
  counts->unique_pred = 0;
  counts->refr_genes = 0;
  counts->pred_genes = 0;
  counts->refr_transcripts = 0;
  counts->pred_transcripts = 0;
  counts->num_loci = 0;
  counts->num_comparisons = 0;
  counts->num_perfect = 0;
  counts->num_mislabeled = 0;
  counts->num_cds_match = 0;
  counts->num_exon_match = 0;
  counts->num_utr_match = 0;
  counts->non_match = 0;
}

void agn_comparison_combine( AgnComparison *stats,
                                   AgnComparison *stats_to_add )
{
  stats->cds_nuc_stats.tp += stats_to_add->cds_nuc_stats.tp;
  stats->cds_nuc_stats.fn += stats_to_add->cds_nuc_stats.fn;
  stats->cds_nuc_stats.fp += stats_to_add->cds_nuc_stats.fp;
  stats->cds_nuc_stats.tn += stats_to_add->cds_nuc_stats.tn;

  stats->utr_nuc_stats.tp += stats_to_add->utr_nuc_stats.tp;
  stats->utr_nuc_stats.fn += stats_to_add->utr_nuc_stats.fn;
  stats->utr_nuc_stats.fp += stats_to_add->utr_nuc_stats.fp;
  stats->utr_nuc_stats.tn += stats_to_add->utr_nuc_stats.tn;

  stats->cds_struc_stats.correct += stats_to_add->cds_struc_stats.correct;
  stats->cds_struc_stats.missing += stats_to_add->cds_struc_stats.missing;
  stats->cds_struc_stats.wrong   += stats_to_add->cds_struc_stats.wrong  ;

  stats->exon_struc_stats.correct += stats_to_add->exon_struc_stats.correct;
  stats->exon_struc_stats.missing += stats_to_add->exon_struc_stats.missing;
  stats->exon_struc_stats.wrong   += stats_to_add->exon_struc_stats.wrong  ;

  stats->utr_struc_stats.correct += stats_to_add->utr_struc_stats.correct;
  stats->utr_struc_stats.missing += stats_to_add->utr_struc_stats.missing;
  stats->utr_struc_stats.wrong   += stats_to_add->utr_struc_stats.wrong  ;

  stats->overall_matches += stats_to_add->overall_matches;
  stats->overall_length += stats_to_add->overall_length;
}

void agn_comparison_init(AgnComparison *stats)
{
  stats->cds_nuc_stats.tp = 0;
  stats->cds_nuc_stats.fn = 0;
  stats->cds_nuc_stats.fp = 0;
  stats->cds_nuc_stats.tn = 0;
  stats->cds_nuc_stats.mc = 0.0;
  stats->cds_nuc_stats.cc = 0.0;
  stats->cds_nuc_stats.sn = 0.0;
  stats->cds_nuc_stats.sp = 0.0;
  stats->cds_nuc_stats.f1 = 0.0;
  stats->cds_nuc_stats.ed = 0.0;

  stats->utr_nuc_stats.tp = 0;
  stats->utr_nuc_stats.fn = 0;
  stats->utr_nuc_stats.fp = 0;
  stats->utr_nuc_stats.tn = 0;
  stats->utr_nuc_stats.mc = 0.0;
  stats->utr_nuc_stats.cc = 0.0;
  stats->utr_nuc_stats.sn = 0.0;
  stats->utr_nuc_stats.sp = 0.0;
  stats->utr_nuc_stats.f1 = 0.0;
  stats->utr_nuc_stats.ed = 0.0;

  stats->cds_struc_stats.correct = 0;
  stats->cds_struc_stats.missing = 0;
  stats->cds_struc_stats.wrong   = 0;
  stats->cds_struc_stats.sn      = 0.0;
  stats->cds_struc_stats.sp      = 0.0;
  stats->cds_struc_stats.f1      = 0.0;
  stats->cds_struc_stats.ed      = 0.0;

  stats->exon_struc_stats.correct = 0;
  stats->exon_struc_stats.missing = 0;
  stats->exon_struc_stats.wrong   = 0;
  stats->exon_struc_stats.sn      = 0.0;
  stats->exon_struc_stats.sp      = 0.0;
  stats->exon_struc_stats.f1      = 0.0;
  stats->exon_struc_stats.ed      = 0.0;

  stats->utr_struc_stats.correct = 0;
  stats->utr_struc_stats.missing = 0;
  stats->utr_struc_stats.wrong   = 0;
  stats->utr_struc_stats.sn      = 0.0;
  stats->utr_struc_stats.sp      = 0.0;
  stats->utr_struc_stats.f1      = 0.0;
  stats->utr_struc_stats.ed      = 0.0;

  stats->overall_matches = 0;
  stats->overall_length = 0;
  stats->overall_identity = 0.0;
  stats->tolerance = 0.0;
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

void agn_compare_filters_parse(AgnCompareFilters *filters, FILE *instream)
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
        unsigned long val = atol(value);
        filters->LocusLengthUpperLimit = val;
      }
      else if(strcmp(key, "LocusLengthLowerLimit") == 0)
      {
        unsigned long val = atol(value);
        filters->LocusLengthLowerLimit = val;
      }
      else if(strcmp(key, "MinReferenceGeneModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MinReferenceGeneModels = val;
      }
      else if(strcmp(key, "MaxReferenceGeneModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxReferenceGeneModels = val;
      }
      else if(strcmp(key, "MinPredictionGeneModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MinPredictionGeneModels = val;
      }
      else if(strcmp(key, "MaxPredictionGeneModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxPredictionGeneModels = val;
      }
      else if(strcmp(key, "MinReferenceTranscriptModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MinReferenceTranscriptModels = val;
      }
      else if(strcmp(key, "MaxReferenceTranscriptModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxReferenceTranscriptModels = val;
      }
      else if(strcmp(key, "MinPredictionTranscriptModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MinPredictionTranscriptModels = val;
      }
      else if(strcmp(key, "MaxPredictionTranscriptModels") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxPredictionTranscriptModels = val;
      }
      else if(strcmp(key, "MinTranscriptsPerReferenceGeneModel") == 0)
      {
        unsigned long val = atol(value);
        filters->MinTranscriptsPerReferenceGeneModel = val;
      }
      else if(strcmp(key, "MaxTranscriptsPerReferenceGeneModel") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxTranscriptsPerReferenceGeneModel = val;
      }
      else if(strcmp(key, "MinTranscriptsPerPredictionGeneModel") == 0)
      {
        unsigned long val = atol(value);
        filters->MinTranscriptsPerPredictionGeneModel = val;
      }
      else if(strcmp(key, "MaxTranscriptsPerPredictionGeneModel") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxTranscriptsPerPredictionGeneModel = val;
      }
      else if(strcmp(key, "MinReferenceExons") == 0)
      {
        unsigned long val = atol(value);
        filters->MinReferenceExons = val;
      }
      else if(strcmp(key, "MaxReferenceExons") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxReferenceExons = val;
      }
      else if(strcmp(key, "MinPredictionExons") == 0)
      {
        unsigned long val = atol(value);
        filters->MinPredictionExons = val;
      }
      else if(strcmp(key, "MaxPredictionExons") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxPredictionExons = val;
      }
      else if(strcmp(key, "MinReferenceCDSLength") == 0)
      {
        unsigned long val = atol(value);
        filters->MinReferenceCDSLength = val;
      }
      else if(strcmp(key, "MaxReferenceCDSLength") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxReferenceCDSLength = val;
      }
      else if(strcmp(key, "MinPredictionCDSLength") == 0)
      {
        unsigned long val = atol(value);
        filters->MinPredictionCDSLength = val;
      }
      else if(strcmp(key, "MaxPredictionCDSLength") == 0)
      {
        unsigned long val = atol(value);
        filters->MaxPredictionCDSLength = val;
      }
      else
      {
        fprintf(stderr, "error: unrecognized filter option '%s'\n", key);
        exit(1);
      }
    }
  }
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
