#include <math.h>
#include "AgnCompareReportText.h"


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

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
 * @function FIXME
 */
static void compare_report_text_pair_nucleotide(FILE *outstream,
                                                AgnCliquePair *pair);

/**
 * @function FIXME
 */
static void compare_report_text_pair_structure(FILE *outstream,
                                               AgnCompStatsBinary *stats,
                                               const char *label,
                                               const char *units);

/**
 * @function FIXME
 */
static void compare_report_text_print_pair(AgnCliquePair *pair,FILE *outstream);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_compare_report_text_create_summary(AgnCompareReportText *rpt,
                                            FILE *outstream)
{
  GT_UNUSED AgnComparisonData *data = agn_compare_report_data(rpt);
  GtStrArray *seqids = agn_compare_report_seqids(rpt);
  fprintf(outstream, "  Sequences compared\n");
  GtUword i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    fprintf(outstream, "    %s\n", gt_str_array_get(seqids, i));
  }

  GtUword numnotshared = data->info.unique_refr_loci +
                         data->info.unique_pred_loci;
  fprintf(outstream, "\n  Gene loci................................%lu\n"
                     "    shared.................................%lu\n"
                     "    unique to reference....................%lu\n"
                     "    unique to prediction...................%lu\n\n",
           data->info.num_loci, data->info.num_loci - numnotshared,
           data->info.unique_refr_loci, data->info.unique_pred_loci);

  fprintf(outstream, "  Reference annotations\n"
                     "    genes..................................%lu\n"
                     "      average per locus....................%.3f\n"
                     "    transcripts............................%lu\n"
                     "      average per locus....................%.3f\n"
                     "      average per gene.....................%.3f\n\n",
           data->info.refr_genes,
           (float)data->info.refr_genes / (float)data->info.num_loci,
           data->info.refr_transcripts,
           (float)data->info.refr_transcripts / (float)data->info.num_loci,
           (float)data->info.refr_transcripts / (float)data->info.refr_genes );

  fprintf(outstream, "  Prediction annotations\n"
                     "    genes..................................%lu\n"
                     "      average per locus....................%.3f\n"
                     "    transcripts............................%lu\n"
                     "      average per locus....................%.3f\n"
                     "      average per gene.....................%.3f\n\n",
           data->info.pred_genes,
           (float)data->info.pred_genes / (float)data->info.num_loci,
           data->info.pred_transcripts,
           (float)data->info.pred_transcripts / (float)data->info.num_loci,
           (float)data->info.pred_transcripts / (float)data->info.pred_genes );

  fprintf(outstream, "  Total comparisons........................%lu\n",
          data->info.num_comparisons);

  GtUword numperfect = data->summary.perfect_matches.comparison_count;
  fprintf( outstream, "    perfect matches........................%lu (%.1f%%)\n",
           numperfect,
           ((float)numperfect / (float)data->info.num_comparisons)*100.0);
  if(numperfect > 0)
  {
    fprintf(outstream, "      avg. length..........................%.2lf bp\n",
            (double)data->summary.perfect_matches.total_length /
            (double)data->summary.perfect_matches.comparison_count );
    fprintf(outstream, "      avg. # refr exons....................%.2lf\n",
            (double)data->summary.perfect_matches.refr_exon_count /
            (double)data->summary.perfect_matches.comparison_count );
    fprintf(outstream, "      avg. # pred exons....................%.2lf\n",
            (double)data->summary.perfect_matches.pred_exon_count /
            (double)data->summary.perfect_matches.comparison_count );
    fprintf(outstream, "      avg. refr CDS length.................%.2lf aa\n",
            (double)data->summary.perfect_matches.refr_cds_length / 3 /
            (double)data->summary.perfect_matches.comparison_count );
    fprintf(outstream, "      avg. pred CDS length.................%.2lf aa\n",
            (double)data->summary.perfect_matches.pred_cds_length / 3 /
            (double)data->summary.perfect_matches.comparison_count );
  }

  GtUword nummislabeled = data->summary.perfect_mislabeled.comparison_count;
  fprintf( outstream, "    perfect matches with mislabeled UTRs...%lu (%.1f%%)\n",
           nummislabeled,
           ((float)nummislabeled / (float)data->info.num_comparisons)*100.0);
  if(nummislabeled > 0)
  {
    fprintf(outstream, "      avg. length..........................%.2lf bp\n",
            (double)data->summary.perfect_mislabeled.total_length /
            (double)data->summary.perfect_mislabeled.comparison_count );
    fprintf(outstream, "      avg. # refr exons....................%.2lf\n",
            (double)data->summary.perfect_mislabeled.refr_exon_count /
            (double)data->summary.perfect_mislabeled.comparison_count );
    fprintf(outstream, "      avg. # pred exons....................%.2lf\n",
            (double)data->summary.perfect_mislabeled.pred_exon_count /
            (double)data->summary.perfect_mislabeled.comparison_count );
    fprintf(outstream, "      avg. refr CDS length.................%.2lf aa\n",
            (double)data->summary.perfect_mislabeled.refr_cds_length / 3 /
            (double)data->summary.perfect_mislabeled.comparison_count );
    fprintf(outstream, "      avg. pred CDS length.................%.2lf aa\n",
            (double)data->summary.perfect_mislabeled.pred_cds_length / 3 /
            (double)data->summary.perfect_mislabeled.comparison_count );
  }

  GtUword numcdsmatch = data->summary.cds_matches.comparison_count;
  fprintf( outstream, "    CDS structure matches..................%lu (%.1f%%)\n",
           numcdsmatch,
           ((float)numcdsmatch / (float)data->info.num_comparisons)*100.0);
  if(numcdsmatch > 0)
  {
    fprintf(outstream, "      avg. length..........................%.2lf bp\n",
            (double)data->summary.cds_matches.total_length /
            (double)data->summary.cds_matches.comparison_count );
    fprintf(outstream, "      avg. # refr exons....................%.2lf\n",
            (double)data->summary.cds_matches.refr_exon_count /
            (double)data->summary.cds_matches.comparison_count );
    fprintf(outstream, "      avg. # pred exons....................%.2lf\n",
            (double)data->summary.cds_matches.pred_exon_count /
            (double)data->summary.cds_matches.comparison_count );
    fprintf(outstream, "      avg. refr CDS length.................%.2lf aa\n",
            (double)data->summary.cds_matches.refr_cds_length / 3 /
            (double)data->summary.cds_matches.comparison_count );
    fprintf(outstream, "      avg. pred CDS length.................%.2lf aa\n",
            (double)data->summary.cds_matches.pred_cds_length / 3 /
            (double)data->summary.cds_matches.comparison_count );
  }

  GtUword numexonmatch = data->summary.exon_matches.comparison_count;
  fprintf( outstream, "    Exon structure matches.................%lu (%.1f%%)\n",
           numexonmatch,
           ((float)numexonmatch / (float)data->info.num_comparisons)*100.0);
  if(numexonmatch > 0)
  {
    fprintf(outstream, "      avg. length..........................%.2lf bp\n",
            (double)data->summary.exon_matches.total_length /
            (double)data->summary.exon_matches.comparison_count );
    fprintf(outstream, "      avg. # refr exons....................%.2lf\n",
            (double)data->summary.exon_matches.refr_exon_count /
            (double)data->summary.exon_matches.comparison_count );
    fprintf(outstream, "      avg. # pred exons....................%.2lf\n",
            (double)data->summary.exon_matches.pred_exon_count /
            (double)data->summary.exon_matches.comparison_count );
    fprintf(outstream, "      avg. refr CDS length.................%.2lf aa\n",
            (double)data->summary.exon_matches.refr_cds_length / 3 /
            (double)data->summary.exon_matches.comparison_count );
    fprintf(outstream, "      avg. pred CDS length.................%.2lf aa\n",
            (double)data->summary.exon_matches.pred_cds_length / 3 /
            (double)data->summary.exon_matches.comparison_count );
  }

  GtUword numutrmatch = data->summary.utr_matches.comparison_count;
  fprintf( outstream, "    UTR structure matches..................%lu (%.1f%%)\n",
           numutrmatch,
           ((float)numutrmatch / (float)data->info.num_comparisons)*100.0);
  if(numutrmatch > 0)
  {
    fprintf(outstream, "      avg. length..........................%.2lf bp\n",
            (double)data->summary.utr_matches.total_length /
            (double)data->summary.utr_matches.comparison_count );
    fprintf(outstream, "      avg. # refr exons....................%.2lf\n",
            (double)data->summary.utr_matches.refr_exon_count /
            (double)data->summary.utr_matches.comparison_count );
    fprintf(outstream, "      avg. # pred exons....................%.2lf\n",
            (double)data->summary.utr_matches.pred_exon_count /
            (double)data->summary.utr_matches.comparison_count );
    fprintf(outstream, "      avg. refr CDS length.................%.2lf aa\n",
            (double)data->summary.utr_matches.refr_cds_length / 3 /
            (double)data->summary.utr_matches.comparison_count );
    fprintf(outstream, "      avg. pred CDS length.................%.2lf aa\n",
            (double)data->summary.utr_matches.pred_cds_length / 3 /
            (double)data->summary.utr_matches.comparison_count );
  }

  GtUword numnonmatch = data->summary.non_matches.comparison_count;
  fprintf( outstream, "    non-matches............................%lu (%.1f%%)\n",
           numnonmatch,
           ((float)numnonmatch / (float)data->info.num_comparisons)*100.0);
  if(numnonmatch > 0)
  {
    fprintf(outstream, "      avg. length..........................%.2lf bp\n",
            (double)data->summary.non_matches.total_length /
            (double)data->summary.non_matches.comparison_count );
    fprintf(outstream, "      avg. # refr exons....................%.2lf\n",
            (double)data->summary.non_matches.refr_exon_count /
            (double)data->summary.non_matches.comparison_count );
    fprintf(outstream, "      avg. # pred exons....................%.2lf\n",
            (double)data->summary.non_matches.pred_exon_count /
            (double)data->summary.non_matches.comparison_count );
    fprintf(outstream, "      avg. refr CDS length.................%.2lf aa\n",
            (double)data->summary.non_matches.refr_cds_length / 3 /
            (double)data->summary.non_matches.comparison_count );
    fprintf(outstream, "      avg. pred CDS length.................%.2lf aa\n",
            (double)data->summary.non_matches.pred_cds_length / 3 /
            (double)data->summary.non_matches.comparison_count );
  }
  fputs("\n", outstream);

  fprintf( outstream, "  CDS structure comparison\n" );
  fprintf( outstream, "    reference CDS segments.................%lu\n",
           data->stats.cds_struc_stats.correct + data->stats.cds_struc_stats.missing );
  fprintf( outstream, "      match prediction.....................%lu (%.1f%%)\n",
           data->stats.cds_struc_stats.correct,
           ((float)data->stats.cds_struc_stats.correct /(float)
           (data->stats.cds_struc_stats.correct+data->stats.cds_struc_stats.missing))*100 );
  fprintf( outstream, "      don't match prediction...............%lu (%.1f%%)\n",
           data->stats.cds_struc_stats.missing,
           ((float)data->stats.cds_struc_stats.missing/(float)
           (data->stats.cds_struc_stats.correct+data->stats.cds_struc_stats.missing))*100 );
  fprintf( outstream, "    prediction CDS segments................%lu\n",
           data->stats.cds_struc_stats.correct + data->stats.cds_struc_stats.wrong );
  fprintf( outstream, "      match reference......................%lu (%.1f%%)\n",
           data->stats.cds_struc_stats.correct,
           ((float)data->stats.cds_struc_stats.correct/(float)
           (data->stats.cds_struc_stats.correct+data->stats.cds_struc_stats.wrong))*100 );
  fprintf( outstream, "      don't match reference................%lu (%.1f%%)\n",
           data->stats.cds_struc_stats.wrong,
           ((float)data->stats.cds_struc_stats.wrong/(float)
           (data->stats.cds_struc_stats.correct+data->stats.cds_struc_stats.wrong))*100 );
  fprintf( outstream, "    Sensitivity............................%.3lf\n",
           data->stats.cds_struc_stats.sn );
  fprintf( outstream, "    Specificity............................%.3lf\n",
           data->stats.cds_struc_stats.sp );
  fprintf( outstream, "    F1 Score...............................%.3lf\n",
           data->stats.cds_struc_stats.f1 );
  fprintf( outstream, "    Annotation edit distance...............%.3lf\n\n",
           data->stats.cds_struc_stats.ed );

  fprintf( outstream, "  Exon structure comparison\n");
  fprintf( outstream, "    reference exons........................%lu\n",
           data->stats.exon_struc_stats.correct + data->stats.exon_struc_stats.missing );
  fprintf( outstream, "      match prediction.....................%lu (%.1f%%)\n",
           data->stats.exon_struc_stats.correct,
           ((float)data->stats.exon_struc_stats.correct/(float)
           (data->stats.exon_struc_stats.correct+data->stats.exon_struc_stats.missing))*100 );
  fprintf( outstream, "      don't match prediction...............%lu (%.1f%%)\n",
           data->stats.exon_struc_stats.missing,
           ((float)data->stats.exon_struc_stats.missing/(float)
           (data->stats.exon_struc_stats.correct+data->stats.exon_struc_stats.missing))*100 );
  fprintf( outstream, "    prediction exons.......................%lu\n",
           data->stats.exon_struc_stats.correct + data->stats.exon_struc_stats.wrong );
  fprintf( outstream, "      match reference......................%lu (%.1f%%)\n",
           data->stats.exon_struc_stats.correct,
           ((float)data->stats.exon_struc_stats.correct/(float)
           (data->stats.exon_struc_stats.correct+data->stats.exon_struc_stats.wrong))*100 );
  fprintf( outstream, "      don't match reference................%lu (%.1f%%)\n",
           data->stats.exon_struc_stats.wrong,
           ((float)data->stats.exon_struc_stats.wrong/(float)
           (data->stats.exon_struc_stats.correct+data->stats.exon_struc_stats.wrong))*100 );
  fprintf( outstream, "    Sensitivity............................%.3lf\n",
           data->stats.exon_struc_stats.sn );
  fprintf( outstream, "    Specificity............................%.3lf\n",
           data->stats.exon_struc_stats.sp );
  fprintf( outstream, "    F1 Score...............................%.3lf\n",
           data->stats.exon_struc_stats.f1 );
  fprintf( outstream, "    Annotation edit distance...............%.3lf\n\n",
           data->stats.exon_struc_stats.ed );

  fprintf( outstream, "  UTR structure comparison\n");
  fprintf( outstream, "    reference UTR segments.................%lu\n",
           data->stats.utr_struc_stats.correct + data->stats.utr_struc_stats.missing );
  if(data->stats.utr_struc_stats.correct + data->stats.utr_struc_stats.missing > 0)
  {
    fprintf( outstream, "      match prediction.....................%lu (%.1f%%)\n",
             data->stats.utr_struc_stats.correct,
             ((float)data->stats.utr_struc_stats.correct/(float)
             (data->stats.utr_struc_stats.correct+data->stats.utr_struc_stats.missing))*100 );
    fprintf( outstream, "      don't match prediction...............%lu (%.1f%%)\n",
             data->stats.utr_struc_stats.missing,
             ((float)data->stats.utr_struc_stats.missing/(float)
             (data->stats.utr_struc_stats.correct+data->stats.utr_struc_stats.missing))*100 );
  }
  fprintf( outstream, "    prediction UTR segments................%lu\n",
           data->stats.utr_struc_stats.correct + data->stats.utr_struc_stats.wrong );
  if(data->stats.utr_struc_stats.correct + data->stats.utr_struc_stats.wrong > 0)
  {
    fprintf( outstream, "      match reference......................%lu (%.1f%%)\n",
             data->stats.utr_struc_stats.correct,
             ((float)data->stats.utr_struc_stats.correct/(float)
             (data->stats.utr_struc_stats.correct+data->stats.utr_struc_stats.wrong))*100 );
    fprintf( outstream, "      don't match reference................%lu (%.1f%%)\n",
             data->stats.utr_struc_stats.wrong,
             ((float)data->stats.utr_struc_stats.wrong/(float)
             (data->stats.utr_struc_stats.correct+data->stats.utr_struc_stats.wrong))*100 );
  }
  fprintf( outstream, "    Sensitivity............................%s\n",
           data->stats.utr_struc_stats.sns );
  fprintf( outstream, "    Specificity............................%s\n",
           data->stats.utr_struc_stats.sps );
  fprintf( outstream, "    F1 Score...............................%s\n",
           data->stats.utr_struc_stats.f1s );
  fprintf( outstream, "    Annotation edit distance...............%s\n\n",
           data->stats.utr_struc_stats.eds );

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

static void compare_report_text_locus_handler(AgnLocus *locus, void *data)
{
  FILE *outstream = data;
  compare_report_text_locus_header(locus, outstream);
  /*
       implement method here
   */
}

static void compare_report_text_locus_header(AgnLocus *locus, FILE *outstream)
{
  GtArray *pairs2report;
  GtUword i;
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

  pairs2report = agn_locus_pairs_to_report(locus);
  if(pairs2report == NULL || gt_array_size(pairs2report) == 0)
  {
    fprintf(outstream,
            "     |\n"
            "     | No comparisons were performed for this locus.\n"
            "     |\n\n");
    return;
  }

  for(i = 0; i < gt_array_size(pairs2report); i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(pairs2report, i);
    compare_report_text_print_pair(pair, outstream);
  }
  fputs("\n", outstream);
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
  double tolerance = agn_clique_pair_tolerance(pair);
  if(fabs(identity - 1.0) < tolerance)
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

  fprintf(outstream, " | reference GFF3:\n");
  agn_transcript_clique_to_gff3(refrclique, outstream, " | ");
  fprintf(outstream, " | prediction GFF3:\n");
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
