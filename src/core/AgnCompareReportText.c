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
  /*
       implement method here
   */
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
