#include <math.h>
#include <string.h>
#include "AgnCompareReportText.h"


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Print an overview of reference and prediction annotations for the
 * summary report.
 */
static void compare_report_text_annot_summary(AgnCompInfo *info,
                                              FILE *outstream);

/**
 * @function Print overall comparison statistics for a particular class of
 * feature comparisons in the summary report.
 */
static void compare_report_text_comp_class_summary(AgnCompClassDesc *summ,
                                                   GtUword num_comparisons,
                                                   const char *label,
                                                   FILE *outstream);

/**
 * @function Create a report for each locus.
 */
static void compare_report_text_locus_handler(AgnLocus *locus, void *data);

/**
 * @function Print locus report header.
 */
static void compare_report_text_locus_header(AgnLocus *locus, FILE *outstream);

/**
 * @function Print gene IDs for locus report header.
 */
static void compare_report_text_locus_gene_ids(AgnLocus *locus,FILE *outstream);

/**
 * @function Print a report of nucleotide-level structure comparison for the
 * given clique pair.
 */
static void compare_report_text_pair_nucleotide(FILE *outstream,
                                                AgnCliquePair *pair);

/**
 * @function Print a report of feature-level structure comparison for the
 * given clique pair.
 */
static void compare_report_text_pair_structure(FILE *outstream,
                                               AgnCompStatsBinary *stats,
                                               const char *label,
                                               const char *units);

/**
 * @function Print a comparison report for the given clique pair.
 */
static void compare_report_text_print_pair(AgnCliquePair *pair,FILE *outstream);

/**
 * @function Print a breakdown of characteristics of loci that fall into a
 * particular comparison class.
 */
static void compare_report_text_summary_struc(FILE *outstream,
                                              AgnCompStatsBinary *stats,
                                              const char *label,
                                              const char *units);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_compare_report_text_create_summary(AgnCompareReportText *rpt,
                                            FILE *outstream)
{
  AgnComparisonData *data = agn_compare_report_data(rpt);
  GtStrArray *seqids = agn_compare_report_seqids(rpt);
  fprintf(outstream, "  Sequences compared\n");
  GtUword i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    fprintf(outstream, "    %s\n", gt_str_array_get(seqids, i));
  }
  fputs("\n", outstream);

  compare_report_text_annot_summary(&data->info, outstream);

  fprintf(outstream, "  Total comparisons........................%lu\n",
          data->info.num_comparisons);

  compare_report_text_comp_class_summary(&data->summary.perfect_matches,
                                         data->info.num_comparisons,
                                         "perfect matches", outstream);
  compare_report_text_comp_class_summary(&data->summary.perfect_mislabeled,
                                         data->info.num_comparisons,
                                        "perfect matches with mislabeled UTRs",
                                         outstream);
  compare_report_text_comp_class_summary(&data->summary.cds_matches,
                                         data->info.num_comparisons,
                                         "CDS structure matches", outstream);
  compare_report_text_comp_class_summary(&data->summary.exon_matches,
                                         data->info.num_comparisons,
                                         "exon structure matches", outstream);
  compare_report_text_comp_class_summary(&data->summary.utr_matches,
                                         data->info.num_comparisons,
                                         "UTR structure matches", outstream);
  compare_report_text_comp_class_summary(&data->summary.non_matches,
                                         data->info.num_comparisons,
                                         "non-matches", outstream);
  fputs("\n", outstream);

  compare_report_text_summary_struc(outstream, &data->stats.cds_struc_stats,
                                    "CDS", "CDS segments");
  compare_report_text_summary_struc(outstream, &data->stats.exon_struc_stats,
                                    "Exon", "exons");
  compare_report_text_summary_struc(outstream, &data->stats.utr_struc_stats,
                                    "UTR", "UTR segments");

  double identity = (double)data->stats.overall_matches /
                    (double)data->stats.overall_length;
  fprintf( outstream, "  %-30s   %-10s   %-10s   %-10s\n",
           "Nucleotide-level comparison", "CDS", "UTRs", "Overall" );
  fprintf( outstream, "    %-30s %-10s   %-10s   %-.3lf\n",
           "Matching coefficient:", data->stats.cds_nuc_stats.mcs,
           data->stats.utr_nuc_stats.mcs, identity);
  fprintf( outstream, "    %-30s %-10s   %-10s   %-10s\n",
           "Correlation coefficient:", data->stats.cds_nuc_stats.ccs,
           data->stats.utr_nuc_stats.ccs, "--");
  fprintf( outstream, "    %-30s %-10s   %-10s   %-10s\n", "Sensitivity:",
           data->stats.cds_nuc_stats.sns, data->stats.utr_nuc_stats.sns, "--");
  fprintf( outstream, "    %-30s %-10s   %-10s   %-10s\n", "Specificity:",
           data->stats.cds_nuc_stats.sps, data->stats.utr_nuc_stats.sps, "--");
  fprintf( outstream, "    %-30s %-10s   %-10s   %-10s\n", "F1 Score:",
           data->stats.cds_nuc_stats.f1s, data->stats.utr_nuc_stats.f1s, "--");
  fprintf( outstream, "    %-30s %-10s   %-10s   %-10s\n", "Annotation edit distance:",
           data->stats.cds_nuc_stats.eds, data->stats.utr_nuc_stats.eds, "--");
}

GtNodeVisitor *agn_compare_report_text_new(FILE *outstream, GtLogger *logger)
{
  GtNodeVisitor *rpt;
  GtArray *filters;

  filters = gt_array_new( sizeof(AgnLocusFilter) );
  rpt = agn_compare_report_new(filters, logger);
  agn_compare_report_set_locus_callback((AgnCompareReport *)rpt,
                                        compare_report_text_locus_handler,
                                        outstream);
  gt_array_delete(filters);
  return rpt;
}

static void compare_report_text_annot_summary(AgnCompInfo *info,
                                              FILE *outstream)
{
  GtUword numnotshared = info->unique_refr_loci +
                         info->unique_pred_loci;
  fprintf(outstream, "  Gene loci................................%lu\n"
                     "    shared.................................%lu\n"
                     "    unique to reference....................%lu\n"
                     "    unique to prediction...................%lu\n\n",
           info->num_loci, info->num_loci - numnotshared,
           info->unique_refr_loci, info->unique_pred_loci);

  fprintf(outstream, "  Reference annotations\n"
                     "    genes..................................%lu\n"
                     "      average per locus....................%.3f\n"
                     "    transcripts............................%lu\n"
                     "      average per locus....................%.3f\n"
                     "      average per gene.....................%.3f\n\n",
           info->refr_genes, (float)info->refr_genes / (float)info->num_loci,
           info->refr_transcripts,
           (float)info->refr_transcripts / (float)info->num_loci,
           (float)info->refr_transcripts / (float)info->refr_genes );

  fprintf(outstream, "  Prediction annotations\n"
                     "    genes..................................%lu\n"
                     "      average per locus....................%.3f\n"
                     "    transcripts............................%lu\n"
                     "      average per locus....................%.3f\n"
                     "      average per gene.....................%.3f\n\n",
           info->pred_genes, (float)info->pred_genes / (float)info->num_loci,
           info->pred_transcripts,
           (float)info->pred_transcripts / (float)info->num_loci,
           (float)info->pred_transcripts / (float)info->pred_genes );
}

static void compare_report_text_comp_class_summary(AgnCompClassDesc *summ,
                                                   GtUword num_comparisons,
                                                   const char *label,
                                                   FILE *outstream)
{
  GtUword numclass     = summ->comparison_count;
  float perc_class     = 100.0 * (float)numclass /
                         (float)num_comparisons;
  float mean_len       = (float)summ->total_length /
                         (float)summ->comparison_count;
  float mean_refr_exon = (float)summ->refr_exon_count /
                         (float)summ->comparison_count;
  float mean_pred_exon = (float)summ->pred_exon_count /
                         (float)summ->comparison_count;
  float mean_refr_cds  = (float)summ->refr_cds_length / 3 /
                         (float)summ->comparison_count;
  float mean_pred_cds  = (float)summ->pred_cds_length / 3 /
                         (float)summ->comparison_count;

  char header[128];
  sprintf(header, "    .......................................%lu (%.1f%%)\n",
          numclass, perc_class);
  strncpy(header + 4, label, strlen(label));
  fputs(header, outstream);

  if(numclass == 0)
    return;

  fprintf(outstream, "      avg. length..........................%.2f bp\n"
                     "      avg. # refr exons....................%.2f\n"
                     "      avg. # pred exons....................%.2f\n"
                     "      avg. refr CDS length.................%.2f aa\n"
                     "      avg. pred CDS length.................%.2f aa\n",
          mean_len,mean_refr_exon,mean_pred_exon,mean_refr_cds,mean_pred_cds);
}

static void compare_report_text_locus_handler(AgnLocus *locus, void *data)
{
  GtArray *pairs2report, *unique;
  GtUword i;
  FILE *outstream = data;

  compare_report_text_locus_header(locus, outstream);
  pairs2report = agn_locus_pairs_to_report(locus);
  if(pairs2report == NULL || gt_array_size(pairs2report) == 0)
  {
    fprintf(outstream,
            "     |\n"
            "     | No comparisons were performed for this locus.\n"
            "     |\n");
    return;
  }

  for(i = 0; i < gt_array_size(pairs2report); i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(pairs2report, i);
    compare_report_text_print_pair(pair, outstream);
  }

  unique = agn_locus_get_unique_refr_cliques(locus);
  if(gt_array_size(unique) > 0)
  {
    fprintf(outstream,
            "     |\n"
            "     |--- Unmatched reference transcripts ---\n"
            "     |---------------------------------------\n"
            "     |\n");
    for(i = 0; i < gt_array_size(unique); i++)
    {
      GtUword j;
      AgnTranscriptClique **clique = gt_array_get(unique, i);
      GtArray *ids = agn_transcript_clique_ids(*clique);
      for(j = 0; j < gt_array_size(ids); j++)
      {
        const char *id = *(const char **)gt_array_get(ids, j);
        fprintf(outstream, "     |    %s\n", id);
      }
      gt_array_delete(ids);
    }
    fputs("     |\n", outstream);
  }

  unique = agn_locus_get_unique_pred_cliques(locus);
  if(gt_array_size(unique) > 0)
  {
    fprintf(outstream,
            "     |\n"
            "     |--- Unmatched prediction transcripts ---\n"
            "     |----------------------------------------\n"
            "     |\n");
    for(i = 0; i < gt_array_size(unique); i++)
    {
      GtUword j;
      AgnTranscriptClique **clique = gt_array_get(unique, i);
      GtArray *ids = agn_transcript_clique_ids(*clique);
      for(j = 0; j < gt_array_size(ids); j++)
      {
        const char *id = *(const char **)gt_array_get(ids, j);
        fprintf(outstream, "     |    %s\n", id);
      }
      gt_array_delete(ids);
    }
    fputs("     |\n", outstream);
  }
  fputs("\n", outstream);
}

static void compare_report_text_locus_header(AgnLocus *locus, FILE *outstream)
{
  GtRange range = gt_genome_node_get_range(locus);
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  fprintf(outstream,
          "|-------------------------------------------------\n"
          "|---- Locus: seqid=%s range=%lu-%lu\n"
          "|-------------------------------------------------\n"
          "|\n",
          gt_str_get(seqid), range.start, range.end);
  compare_report_text_locus_gene_ids(locus, outstream);
  fprintf(outstream,
          "|\n"
          "|----------\n");
}

static void compare_report_text_locus_gene_ids(AgnLocus *locus, FILE *outstream)
{
  GtArray *genes;

  genes = agn_locus_refr_gene_ids(locus);
  fprintf(outstream, "|  reference genes:\n");
  if(gt_array_size(genes) == 0)
    fprintf(outstream, "|    None!\n");
  while(gt_array_size(genes) > 0)
  {
    const char **geneid = gt_array_pop(genes);
    fprintf(outstream, "|    %s\n", *geneid);
  }
  fprintf(outstream, "|\n");
  gt_array_delete(genes);

  genes = agn_locus_pred_gene_ids(locus);
  fprintf(outstream, "|  prediction genes:\n");
  if(gt_array_size(genes) == 0)
    fprintf(outstream, "|    None!\n");
  while(gt_array_size(genes) > 0)
  {
    const char **geneid = gt_array_pop(genes);
    fprintf(outstream, "|    %s\n", *geneid);
  }
  gt_array_delete(genes);
}

static void compare_report_text_pair_nucleotide(FILE *outstream,
                                                AgnCliquePair *pair)
{
  AgnComparison *pairstats = agn_clique_pair_get_stats(pair);
  double identity = (double)pairstats->overall_matches /
                    (double)pairstats->overall_length;
  if(pairstats->overall_matches == pairstats->overall_length)
    fprintf(outstream, "     |    Gene structures match perfectly!\n");
  else
  {
    fprintf(outstream, "     |  %-30s   %-10s %-10s %-10s\n",
            "Nucleotide-level comparison", "CDS", "UTRs", "Overall" );
    fprintf(outstream, "     |    %-30s %-10s %-10s %.3lf\n",
            "Matching coefficient:", pairstats->cds_nuc_stats.mcs,
            pairstats->utr_nuc_stats.mcs, identity);
    fprintf(outstream, "     |    %-30s %-10s %-10s %-10s\n",
            "Correlation coefficient:", pairstats->cds_nuc_stats.ccs,
            pairstats->utr_nuc_stats.ccs, "--");
    fprintf(outstream, "     |    %-30s %-10s %-10s %-10s\n", "Sensitivity:",
            pairstats->cds_nuc_stats.sns, pairstats->utr_nuc_stats.sns, "--");
    fprintf(outstream, "     |    %-30s %-10s %-10s %-10s\n", "Specificity:",
            pairstats->cds_nuc_stats.sps, pairstats->utr_nuc_stats.sps, "--");
    fprintf(outstream, "     |    %-30s %-10s %-10s %-10s\n", "F1 Score:",
            pairstats->cds_nuc_stats.f1s, pairstats->utr_nuc_stats.f1s, "--");
    fprintf(outstream, "     |    %-30s %-10s %-10s %-10s\n",
            "Annotation edit distance:", pairstats->cds_nuc_stats.eds,
            pairstats->utr_nuc_stats.eds, "--");
  }

  fprintf(outstream, "     |\n");
}

static void compare_report_text_pair_structure(FILE *outstream,
                                               AgnCompStatsBinary *stats,
                                               const char *label,
                                               const char *units)
{
  fprintf(outstream, "     |  %s structure comparison\n", label);
  if(stats->missing == 0 && stats->wrong == 0)
  {
    fprintf(outstream,
            "     |    %lu reference  %s\n"
            "     |    %lu prediction %s\n"
            "     |    %s structures match perfectly!\n",
            stats->correct, units, stats->correct, units, label);
  }
  else
  {
    fprintf(outstream,
            "     |    %lu reference %s\n"
            "     |        %lu match prediction\n"
            "     |        %lu don't match prediction\n"
            "     |    %lu prediction %s\n"
            "     |        %lu match reference\n"
            "     |        %lu don't match reference\n",
             stats->correct + stats->missing, units,
             stats->correct, stats->missing,
             stats->correct + stats->wrong, units,
             stats->correct, stats->wrong);
    fprintf(outstream,
            "     |    %-30s %-10s\n" 
            "     |    %-30s %-10s\n" 
            "     |    %-30s %-10s\n" 
            "     |    %-30s %-10s\n",
            "Sensitivity:", stats->sns, "Specificity:", stats->sps,
            "F1 Score:", stats->f1s, "Annotation edit distance:",stats->eds);
  }
  fprintf(outstream, "     |\n");
}


static void compare_report_text_print_pair(AgnCliquePair *pair, FILE *outstream)
{
  GtArray *tids;
  AgnTranscriptClique *refrclique, *predclique;

  fprintf(outstream,
          "     |\n"
          "     |--------------------------\n"
          "     |---- Begin comparison ----\n"
          "     |--------------------------\n"
          "     |\n");

  refrclique = agn_clique_pair_get_refr_clique(pair);
  tids = agn_transcript_clique_ids(refrclique);
  fprintf(outstream, "     |  reference transcripts:\n");
  while(gt_array_size(tids) > 0)
  {
    const char **tid = gt_array_pop(tids);
    fprintf(outstream, "     |    %s\n", *tid);
  }
  gt_array_delete(tids);
  predclique = agn_clique_pair_get_pred_clique(pair);
  tids = agn_transcript_clique_ids(predclique);
  fprintf(outstream, "     |  prediction transcripts:\n");
  while(gt_array_size(tids) > 0)
  {
    const char **tid = gt_array_pop(tids);
    fprintf(outstream, "     |    %s\n", *tid);
  }
  gt_array_delete(tids);
  fprintf(outstream, "     |\n");

  fprintf(outstream, "     | reference GFF3:\n");
  agn_transcript_clique_to_gff3(refrclique, outstream, " | ");
  fprintf(outstream, "     | prediction GFF3:\n");
  agn_transcript_clique_to_gff3(predclique, outstream, " | ");
  fprintf(outstream, " |\n");

  AgnComparison *pairstats = agn_clique_pair_get_stats(pair);
  compare_report_text_pair_structure(outstream, &pairstats->cds_struc_stats,
                                     "CDS", "CDS segments");
  compare_report_text_pair_structure(outstream, &pairstats->exon_struc_stats,
                                     "Exon", "exons");
  compare_report_text_pair_structure(outstream, &pairstats->utr_struc_stats,
                                     "UTR", "UTR segments");
  compare_report_text_pair_nucleotide(outstream, pair);  

  fprintf(outstream,
          "     |\n"
          "     |--------------------------\n"
          "     |----- End comparison -----\n"
          "     |--------------------------\n");
}

static void compare_report_text_summary_struc(FILE *outstream,
                                              AgnCompStatsBinary *stats,
                                              const char *label,
                                              const char *units)
{
  char buffer[128];
  fprintf(outstream, "  %s structure comparison\n", label);

  GtUword refrcnt = stats->correct + stats->missing;
  sprintf(buffer, "    reference .............................%lu\n", refrcnt);
  strncpy(buffer + 14, units, strlen(units));
  fputs(buffer, outstream);
  fprintf(outstream,
          "      match prediction.....................%lu (%.1f%%)\n",
          stats->correct, (float)stats->correct / (float)refrcnt * 100 );
  fprintf(outstream,
          "      don't match prediction...............%lu (%.1f%%)\n",
          stats->missing, (float)stats->missing / (float)refrcnt * 100 );

  GtUword predcnt = stats->correct + stats->wrong;
  sprintf(buffer, "    prediction ............................%lu\n", predcnt);
  strncpy(buffer + 15, units, strlen(units));
  fputs(buffer, outstream);
  fprintf(outstream,
          "      match reference......................%lu (%.1f%%)\n",
          stats->correct, (float)stats->correct / (float)predcnt * 100 );
  fprintf(outstream,
          "      don't match reference................%lu (%.1f%%)\n",
          stats->wrong, (float)stats->wrong / (float)predcnt * 100 );

  fprintf(outstream, "    Sensitivity............................%.3lf\n"
                     "    Specificity............................%.3lf\n"
                     "    F1 Score...............................%.3lf\n"
                     "    Annotation edit distance...............%.3lf\n\n",
          stats->sn, stats->sp, stats->f1, stats->ed);
}
