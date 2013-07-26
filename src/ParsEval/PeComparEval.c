#include "PeComparEval.h"
#include "PeReports.h"

void pe_clique_pair_record_characteristics(AgnCliquePair *pair,
                                           PeCompResultDesc *desc)
{
  AgnTranscriptClique *refr = agn_clique_pair_get_refr_clique(pair);
  AgnTranscriptClique *pred = agn_clique_pair_get_refr_clique(pair);

  desc->transcript_count += 1;
  desc->total_length += agn_clique_pair_length(pair);
  desc->refr_cds_length += agn_transcript_clique_cds_length(refr);
  desc->pred_cds_length += agn_transcript_clique_cds_length(pred);
  desc->refr_exon_count += agn_transcript_clique_num_exons(refr);
  desc->pred_exon_count += agn_transcript_clique_num_exons(pred);
}

void pe_comparative_analysis(AgnLocusIndex *locusindex, GtHashmap **comp_evalsp,
                             GtHashmap **locus_summariesp, GtStrArray *seqids,
                             GtArray *seqfiles, GtArray *loci,
                             PeOptions *options)
{
  GtHashmap *comp_evals;
  GtHashmap *locus_summaries;

  AgnLogger *logger = agn_logger_new();
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin comparative analysis\n", stderr);
  locus_summaries = gt_hashmap_new(GT_HASH_DIRECT, NULL, (GtFree)gt_free_mem);
  comp_evals = gt_hashmap_new(GT_HASH_DIRECT, NULL, (GtFree)gt_free_mem);

  int i, j;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *seqloci = *(GtArray **)gt_array_get(loci, i);
    PeAnalysisData analysis_data;
    analysis_data.seqfile = *(FILE **)gt_array_get(seqfiles, i);
    analysis_data.options = options;

    for(j = 0; j < gt_array_size(seqloci); j++)
    {
      AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(seqloci, j);
      PeCompEvaluation *compeval = gt_malloc( sizeof(PeCompEvaluation) );
      pe_comp_evaluation_init(compeval);
      gt_hashmap_add(comp_evals, locus, compeval);
      AgnGeneLocusSummary *locsum = gt_malloc( sizeof(AgnGeneLocusSummary) );
      agn_gene_locus_summary_init(locsum);
      gt_hashmap_add(locus_summaries, locus, locsum);
    }
    analysis_data.comp_evals = comp_evals;
    analysis_data.locus_summaries = locus_summaries;

    agn_locus_index_comparative_analysis(locusindex, seqid, options->numprocs,
                                     (AgnLocusIndexVisitFunc)pe_pre_analysis,
                                     (AgnLocusIndexVisitFunc)pe_post_analysis,
                                     &analysis_data, logger);
    if(options->debug)
      agn_logger_print_status(logger, stderr, "comparative analysis");
    if(agn_logger_has_error(logger))
    {
      agn_logger_print_error(logger, stderr, "issues with comparative analysis "
                             "of sequence '%s'", seqid);
      exit(1);
    }
  }

  *comp_evalsp = comp_evals;
  *locus_summariesp = locus_summaries;

  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished comparative "
                          "analysis (%ld.%06ld seconds)\n", stderr);
  agn_logger_delete(logger);
  gt_timer_delete(timer);
}

void pe_comp_evaluation_combine(PeCompEvaluation *data,
                                 PeCompEvaluation *data_to_add)
{
  agn_comp_summary_combine(&data->counts, &data_to_add->counts);
  agn_comparison_combine(&data->stats, &data_to_add->stats);
  pe_comp_result_summary_combine(&data->results, &data_to_add->results);
}

void pe_comp_evaluation_init(PeCompEvaluation *data)
{
  agn_comp_summary_init(&data->counts);
  agn_comparison_init(&data->stats);
  pe_comp_result_summary_init(&data->results);
}

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

void pe_pre_analysis(AgnGeneLocus *locus, PeAnalysisData *data)
{
  unsigned long npairs = agn_gene_locus_num_clique_pairs(locus);
  if(data->options->complimit != 0 && npairs > data->options->complimit)
  {
    if(data->options->debug)
    {
      fprintf(stderr, "debug: locus %s[%lu, %lu] contains %lu transcript (or "
              "transcript clique) pairs, exceeds the limit of %d, moving on\n",
              agn_gene_locus_get_seqid(locus), agn_gene_locus_get_start(locus),
              agn_gene_locus_get_end(locus), npairs, data->options->complimit);
    }
  }
  else
  {
    PeCompEvaluation *compeval = gt_hashmap_get(data->comp_evals, locus);
    compeval->counts.num_loci++;
    compeval->counts.refr_genes += agn_gene_locus_num_refr_genes(locus);
    compeval->counts.pred_genes += agn_gene_locus_num_pred_genes(locus);
    unsigned long rt = agn_gene_locus_num_refr_transcripts(locus);
    unsigned long pt = agn_gene_locus_num_pred_transcripts(locus);
    compeval->counts.refr_transcripts += rt;
    compeval->counts.pred_transcripts += pt;
    if(agn_gene_locus_num_refr_genes(locus) == 0)
      compeval->counts.unique_pred++;
    else if(agn_gene_locus_num_pred_genes(locus) == 0)
      compeval->counts.unique_refr++;
  }
}

void pe_post_analysis(AgnGeneLocus *locus, PeAnalysisData *data)
{
  PeCompEvaluation *compeval = gt_hashmap_get(data->comp_evals, locus);
  GtArray *pairs = agn_gene_locus_pairs_to_report(locus);
  if(gt_array_size(pairs) > 0)
    pe_gene_locus_aggregate_results(locus, compeval);

  if(!data->options->summary_only)
  {
    pe_gene_locus_print_results(locus, data->seqfile, data->options);

    AgnGeneLocusSummary *locsum = gt_hashmap_get(data->locus_summaries, locus);
    locsum->start = agn_gene_locus_get_start(locus);
    locsum->end = agn_gene_locus_get_end(locus);
    locsum->length = agn_gene_locus_get_length(locus);
    locsum->refrtrans = agn_gene_locus_num_refr_transcripts(locus);
    locsum->predtrans = agn_gene_locus_num_pred_transcripts(locus);
    locsum->reported = gt_array_size(pairs);
    locsum->counts = compeval->counts;

#ifndef WITHOUT_CAIRO
    if(data->options->locus_graphics)
    {
      AgnGeneLocusPngMetadata metadata;
      pe_gene_locus_get_png_filename(locus, metadata.filename,
                                     data->options->outfilename);
      metadata.graphic_width = pe_gene_locus_get_graphic_width(locus);
      sprintf(metadata.stylefile, "%s/pe.style", data->options->data_path);
      metadata.refrfile = data->options->refrfile;
      metadata.predfile = data->options->predfile;
      metadata.refrlabel = data->options->refrlabel;
      metadata.predlabel = data->options->predlabel;
      metadata.track_order_func = pe_track_order;

      agn_gene_locus_print_png(locus, &metadata);
    }
#endif
  }

  agn_gene_locus_delete(locus);
}
