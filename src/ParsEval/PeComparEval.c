#include "PeComparEval.h"

void pe_comp_result_summary_combine(PeCompResultSummary *desc,
                                    PeCompResultSummary *desc_to_add)
{
  pe_comp_result_desc_combine(&desc->perfect_matches,
                              &desc_to_add->perfect_matches);
  pe_comp_result_desc_combine(&desc->perfect_mislabeled,
                              &desc_to_add->perfect_mislabeled);
  pe_comp_result_desc_combine(&desc->cds_matches,
                              &desc_to_add->cds_matches);
  pe_comp_result_desc_combine(&desc->exon_matches,
                              &desc_to_add->exon_matches);
  pe_comp_result_desc_combine(&desc->utr_matches,
                              &desc_to_add->utr_matches);
  pe_comp_result_desc_combine(&desc->non_matches,
                              &desc_to_add->non_matches);
}

void pe_comp_result_summary_init(PeCompResultSummary *desc)
{
  pe_comp_result_desc_init(&desc->perfect_matches);
  pe_comp_result_desc_init(&desc->perfect_mislabeled);
  pe_comp_result_desc_init(&desc->cds_matches);
  pe_comp_result_desc_init(&desc->exon_matches);
  pe_comp_result_desc_init(&desc->utr_matches);
  pe_comp_result_desc_init(&desc->non_matches);
}

void pe_comp_result_desc_combine(PeCompResultDesc *desc,
                                 PeCompResultDesc *desc_to_add)
{
  desc->total_length     += desc_to_add->total_length;
  desc->transcript_count += desc_to_add->transcript_count;
  desc->refr_cds_length  += desc_to_add->refr_cds_length;
  desc->pred_cds_length  += desc_to_add->pred_cds_length;
  desc->refr_exon_count  += desc_to_add->refr_exon_count;
  desc->pred_exon_count  += desc_to_add->pred_exon_count;
}

void pe_comp_result_desc_init(PeCompResultDesc *desc)
{
  desc->total_length = 0;
  desc->transcript_count = 0;
  desc->refr_cds_length = 0;
  desc->pred_cds_length = 0;
  desc->refr_exon_count = 0;
  desc->pred_exon_count = 0;
}

void agn_summary_data_combine(AgnSummaryData *data, AgnSummaryData *data_to_add)
{
  agn_comparison_counts_combine(&data->counts, &data_to_add->counts);
  agn_comparison_stats_combine(&data->stats, &data_to_add->stats);
  pe_comp_result_summary_combine(&data->results, &data_to_add->results);
}

void agn_summary_data_init(AgnSummaryData *data)
{
  agn_comparison_counts_init(&data->counts);
  agn_comparison_stats_init(&data->stats);
  pe_comp_result_summary_init(&data->results);
}
