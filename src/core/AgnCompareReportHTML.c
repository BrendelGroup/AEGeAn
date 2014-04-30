#include <math.h>
#include <string.h>
#include "AgnCompareReportHTML.h"
#include "AgnVersion.h"

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------
typedef struct
{
  GtRange range;
  GtUword numrefrmrnas;
  GtUword numpredmrnas;
  GtUword numperfect;
  GtUword nummislabeled;
  GtUword numcdsmatch;
  GtUword numexonmatch;
  GtUword numutrmatch;
  GtUword numnonmatch;
} LocusSummary;

typedef struct
{
  GtArray *locus_summaries;
  const char *outdir;
} CallbackData;


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Print an overview of reference and prediction annotations for the
 * summary report.
 */
static void compare_report_html_annot_summary(AgnCompInfo *info,
                                              FILE *outstream);

/**
 * @function Print overall comparison statistics for a particular class of
 * feature comparisons in the summary report.
 */
static void compare_report_html_comp_class_summary(AgnCompClassDesc *summ,
                                                   GtUword num_comparisons,
                                                   const char *label,
                                                   FILE *outstream);

/**
 * @function Write the HTML footer to the given output stream.
 */
static void compare_report_html_footer(FILE *outstream);

/**
 * @function Create a report for each locus.
 */
static void compare_report_html_locus_handler(AgnLocus *locus, void *data);

/**
 * @function Print locus report header.
 */
static void compare_report_html_locus_header(AgnLocus *locus, FILE *outstream);

/**
 * @function Print gene IDs for locus report header.
 */
static void compare_report_html_locus_gene_ids(AgnLocus *locus,FILE *outstream);

/**
 * @function Print a report of nucleotide-level structure comparison for the
 * given clique pair.
 */
static void compare_report_html_pair_nucleotide(FILE *outstream,
                                                AgnCliquePair *pair);

/**
 * @function Print a report of feature-level structure comparison for the
 * given clique pair.
 */
static void compare_report_html_pair_structure(FILE *outstream,
                                               AgnCompStatsBinary *stats,
                                               const char *label,
                                               const char *units);

/**
 * @function Print a comparison report for the given clique pair.
 */
static void compare_report_html_print_pair(AgnCliquePair *pair, FILE *outstream,
                                           GtUword k);

/**
 * @function FIXME
 */
static void compare_report_html_sequence_handler(const AgnComparisonData *cd,
                                                 const char *seqid, void *data);

/**
 * @function Print a breakdown of characteristics of loci that fall into a
 * particular comparison class.
 */
static void compare_report_html_summary_struc(FILE *outstream,
                                              AgnCompStatsBinary *stats,
                                              const char *label,
                                              const char *units);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_compare_report_html_create_summary(AgnCompareReportHTML *rpt,
                                            const char *outdir)
{
  AgnComparisonData *data = agn_compare_report_data(rpt);
  GtStrArray *seqids = agn_compare_report_seqids(rpt);
  
  char filename[1024];
  sprintf(filename, "%s/index.html", outdir);
  FILE *outstream = fopen(filename, "w");
  if(outstream == NULL)
  {
    fprintf(stderr, "error: unable to open output file '%s'\n", filename);
    exit(1);
  }
  
  fprintf(outstream, "  Sequences compared\n");
  GtUword i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    fprintf(outstream, "    %s\n", gt_str_array_get(seqids, i));
  }
  fputs("\n", outstream);

  compare_report_html_annot_summary(&data->info, outstream);

  fprintf(outstream, "  Total comparisons........................%lu\n",
          data->info.num_comparisons);

  compare_report_html_comp_class_summary(&data->summary.perfect_matches,
                                         data->info.num_comparisons,
                                         "perfect matches", outstream);
  compare_report_html_comp_class_summary(&data->summary.perfect_mislabeled,
                                         data->info.num_comparisons,
                                        "perfect matches with mislabeled UTRs",
                                         outstream);
  compare_report_html_comp_class_summary(&data->summary.cds_matches,
                                         data->info.num_comparisons,
                                         "CDS structure matches", outstream);
  compare_report_html_comp_class_summary(&data->summary.exon_matches,
                                         data->info.num_comparisons,
                                         "exon structure matches", outstream);
  compare_report_html_comp_class_summary(&data->summary.utr_matches,
                                         data->info.num_comparisons,
                                         "UTR structure matches", outstream);
  compare_report_html_comp_class_summary(&data->summary.non_matches,
                                         data->info.num_comparisons,
                                         "non-matches", outstream);
  fputs("\n", outstream);

  compare_report_html_summary_struc(outstream, &data->stats.cds_struc_stats,
                                    "CDS", "CDS segments");
  compare_report_html_summary_struc(outstream, &data->stats.exon_struc_stats,
                                    "Exon", "exons");
  compare_report_html_summary_struc(outstream, &data->stats.utr_struc_stats,
                                    "UTR", "UTR segments");

  double identity = (double)data->stats.overall_matches /
                    (double)data->stats.overall_length;
  fprintf(outstream, "  %-30s   %-10s   %-10s   %-10s\n",
          "Nucleotide-level comparison", "CDS", "UTRs", "Overall" );
  fprintf(outstream, "    %-30s %-10s   %-10s   %-.3lf\n",
          "Matching coefficient:", data->stats.cds_nuc_stats.mcs,
          data->stats.utr_nuc_stats.mcs, identity);
  fprintf(outstream, "    %-30s %-10s   %-10s   %-10s\n",
          "Correlation coefficient:", data->stats.cds_nuc_stats.ccs,
          data->stats.utr_nuc_stats.ccs, "--");
  fprintf(outstream, "    %-30s %-10s   %-10s   %-10s\n", "Sensitivity:",
          data->stats.cds_nuc_stats.sns, data->stats.utr_nuc_stats.sns, "--");
  fprintf(outstream, "    %-30s %-10s   %-10s   %-10s\n", "Specificity:",
          data->stats.cds_nuc_stats.sps, data->stats.utr_nuc_stats.sps, "--");
  fprintf(outstream, "    %-30s %-10s   %-10s   %-10s\n", "F1 Score:",
          data->stats.cds_nuc_stats.f1s, data->stats.utr_nuc_stats.f1s, "--");
  fprintf(outstream, "    %-30s %-10s   %-10s   %-10s\n", "Annotation edit distance:",
          data->stats.cds_nuc_stats.eds, data->stats.utr_nuc_stats.eds, "--");
}

GtNodeVisitor *agn_compare_report_html_new(const char *outdir, GtLogger *logger)
{
  GtNodeVisitor *rpt = agn_compare_report_new(logger);
  CallbackData cd = { NULL, outdir };
  agn_compare_report_set_locus_callback((AgnCompareReport *)rpt,
                                        compare_report_html_locus_handler, &cd);
  agn_compare_report_set_sequence_callback((AgnCompareReport *)rpt,
                                           compare_report_html_sequence_handler,
                                           &cd);
  return rpt;
}

static void compare_report_html_annot_summary(AgnCompInfo *info,
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

static void compare_report_html_comp_class_summary(AgnCompClassDesc *summ,
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

static void compare_report_html_footer(FILE *outstream)
{
  char shortversion[11];
  strncpy(shortversion, AEGEAN_VERSION, 10);
  shortversion[10] = '\0';

  fprintf(outstream,
          "      <p class=\"footer\">\n"
          "        Generated by <a href=\"%s\">AEGeAn version %s</a>.<br />\n"
          "        Copyright © %s <a href=\"http://aegean.readthedocs.org/en/"
          "latest/contrib.html\">AEGeAn authors</a>.<br />\n"
          "        See <a href=\"LICENSE\">LICENSE</a> for details."
          "      </p>\n", AEGEAN_LINK, shortversion, AEGEAN_COPY_DATE);
}

static void compare_report_html_locus_handler(AgnLocus *locus, void *data)
{
  GtArray *pairs2report, *unique;
  GtUword i;
  CallbackData *cdata = data;
  LocusSummary s;

  if(cdata->locus_summaries == NULL)
    cdata->locus_summaries = gt_array_new( sizeof(LocusSummary) );
  s.range = gt_genome_node_get_range(locus);
  s.numrefrmrnas = agn_locus_num_refr_mrnas(locus);
  s.numpredmrnas = agn_locus_num_pred_mrnas(locus);
  s.numperfect = 0;
  s.nummislabeled = 0;
  s.numcdsmatch = 0;
  s.numexonmatch = 0;
  s.numutrmatch = 0;
  s.numnonmatch = 0;
  
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  char filename[1024];
  sprintf(filename, "%s/%s/%lu-%lu.html", cdata->outdir, gt_str_get(seqid),
          s.range.start, s.range.end);
  FILE *outstream = fopen(filename, "w");
  if(outstream == NULL)
  {
    fprintf(stderr, "error: unable to open output file '%s'\n", filename);
    exit(1);
  }
  
  compare_report_html_locus_header(locus, outstream);
  pairs2report = agn_locus_pairs_to_report(locus);
  if(pairs2report == NULL || gt_array_size(pairs2report) == 0)
  {
    fputs("      <p>No comparisons were performed for this locus.</p>\n\n",
          outstream);
    return;
  }

  fputs("      <h2 class=\"bottomspace\">Comparisons</h2>\n", outstream);
  for(i = 0; i < gt_array_size(pairs2report); i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(pairs2report, i);
    AgnCompClassification cls = agn_clique_pair_classify(pair);
    switch(cls)
    {
      case AGN_COMP_CLASS_PERFECT_MATCH: s.numperfect++;    break;
      case AGN_COMP_CLASS_MISLABELED:    s.nummislabeled++; break;
      case AGN_COMP_CLASS_CDS_MATCH:     s.numcdsmatch++;   break;
      case AGN_COMP_CLASS_EXON_MATCH:    s.numexonmatch++;  break;
      case AGN_COMP_CLASS_UTR_MATCH:     s.numutrmatch++;   break;
      case AGN_COMP_CLASS_NON_MATCH:     s.numnonmatch++;   break;
      case AGN_COMP_CLASS_UNCLASSIFIED:  /* nothing */      break;
      default:                                              break;
    }
    compare_report_html_print_pair(pair, outstream, i);
  }
  gt_array_add(cdata->locus_summaries, s);

  unique = agn_locus_get_unique_refr_cliques(locus);
  if(gt_array_size(unique) > 0)
  {
    fputs("      <h2>Unmatched reference transcripts</h2>\n"
          "      <ul>\n",
          outstream);
    for(i = 0; i < gt_array_size(unique); i++)
    {
      GtUword j;
      AgnTranscriptClique **clique = gt_array_get(unique, i);
      GtArray *ids = agn_transcript_clique_ids(*clique);
      for(j = 0; j < gt_array_size(ids); j++)
      {
        const char *id = *(const char **)gt_array_get(ids, j);
        fprintf(outstream, "        <li>%s</li>\n", id);
      }
      gt_array_delete(ids);
    }
    fputs("      </ul>\n\n", outstream);
  }

  unique = agn_locus_get_unique_pred_cliques(locus);
  if(gt_array_size(unique) > 0)
  {
    fputs("      <h2>Unmatched prediction transcripts</h2>\n"
          "      <ul>\n",
          outstream);
    for(i = 0; i < gt_array_size(unique); i++)
    {
      GtUword j;
      AgnTranscriptClique **clique = gt_array_get(unique, i);
      GtArray *ids = agn_transcript_clique_ids(*clique);
      for(j = 0; j < gt_array_size(ids); j++)
      {
        const char *id = *(const char **)gt_array_get(ids, j);
        fprintf(outstream, "        <li>%s</li>\n", id);
      }
      gt_array_delete(ids);
    }
    fputs("      </ul>\n\n", outstream);
  }
  fputs("\n", outstream);

  compare_report_html_footer(outstream);
  fputs("    </div>\n"
        "  </body>\n"
        "</html>",
        outstream);
  fclose(outstream);
}

static void compare_report_html_locus_header(AgnLocus *locus, FILE *outstream)
{
  GtRange range = gt_genome_node_get_range(locus);
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtArray *pairs2report = agn_locus_pairs_to_report(locus);
  GtUword numpairs = gt_array_size(pairs2report);

  fprintf( outstream,
           "<!doctype html>\n"
           "<html lang=\"en\">\n"
           "  <head>\n"
           "    <meta charset=\"utf-8\" />\n"
           "    <title>ParsEval: Locus at %s[%lu, %lu]</title>\n"
           "    <link rel=\"stylesheet\" type=\"text/css\" "
           "href=\"../parseval.css\" />\n",
           gt_str_get(seqid), range.start, range.end);

  if(numpairs > 0)
  {
    GtUword i;
    fputs("    <script type=\"text/javascript\""
          " src=\"../vendor/mootools-core-1.3.2-full-nocompat-yc.js\"></script>\n"
          "    <script type=\"text/javascript\" src=\"../vendor/mootools-more-1.3.2.1.js\"></script>\n"
          "    <script type=\"text/javascript\">\n"
          "window.addEvent('domready', function() {\n"
          "  var status =\n"
          "  {\n"
          "    'true': \"(hide details)\",\n"
          "    'false': \"(show details)\",\n"
          "  }\n",
          outstream);
    for(i = 0; i < numpairs; i++)
    {
      fprintf(outstream,
              "  var compareWrapper%lu = new Fx.Slide('compare_wrapper_%lu');\n"
              "  compareWrapper%lu.hide();\n"
              "  $('toggle_compare_%lu').addEvent('click', function(event){\n"
              "    event.stop();\n"
              "    compareWrapper%lu.toggle();\n"
              "  });\n"
              "  compareWrapper%lu.addEvent('complete', function() {\n"
              "    $('toggle_compare_%lu').set('text', status[compareWrapper%lu.open]);\n"
              "  });\n",
              i, i, i, i, i, i, i, i);
    }
    fputs("});\n"
          "    </script>\n",
          outstream);
  }

  fprintf(outstream,
          "  </head>\n"
          "  <body>\n"
          "    <div id=\"content\">\n"
          "      <h1>Locus at %s[%lu, %lu]</h1>\n"
          "      <p><a href=\"index.html\">⇐ Back to %s loci</a></p>\n\n",
          gt_str_get(seqid), range.start, range.end, gt_str_get(seqid) );

  compare_report_html_locus_gene_ids(locus, outstream);
}

static void compare_report_html_locus_gene_ids(AgnLocus *locus, FILE *outstream)
{
  GtUword i;
  GtArray *refr_genes = agn_locus_refr_gene_ids(locus);
  GtArray *pred_genes = agn_locus_pred_gene_ids(locus);
  
  fputs("      <h2>Gene annotations</h2>\n"
        "      <table>\n"
        "        <tr><th>Reference</th><th>Prediction</th></tr>\n",
        outstream);
  for(i = 0; i < gt_array_size(refr_genes) || i < gt_array_size(pred_genes); i++)
  {
    fputs("        <tr>", outstream);
    if(i < gt_array_size(refr_genes))
    {
      const char *gid = *(const char **)gt_array_get(refr_genes, i);
      fprintf(outstream, "<td>%s</td>", gid);
    }
    else
    {
      if(i == 0) fputs("<td>None</td>", outstream);
      else       fputs("<td>&nbsp;</td>", outstream);
    }

    if(i < gt_array_size(pred_genes))
    {
      const char *gid = *(const char **)gt_array_get(pred_genes, i);
      fprintf(outstream, "<td>%s</td>", gid);
    }
    else
    {
      if(i == 0) fputs("<td>None</td>", outstream);
      else       fputs("<td>&nbsp;</td>", outstream);
    }
    fputs("</tr>\n", outstream);
  }
  fputs("      </table>\n\n", outstream);
  gt_array_delete(refr_genes);
  gt_array_delete(pred_genes);
}

static void compare_report_html_pair_nucleotide(FILE *outstream,
                                                AgnCliquePair *pair)
{
  AgnComparison *pairstats = agn_clique_pair_get_stats(pair);
  double identity = (double)pairstats->overall_matches /
                    (double)pairstats->overall_length;
  if(pairstats->overall_matches == pairstats->overall_length)
    fputs("        <h3>Gene structures match perfectly!</h3>\n", outstream);
  else
  {
    fprintf(outstream,
            "        <h3>Nucleotide-level comparison</h3>\n"
            "        <table class=\"table_wide table_extra_indent\">\n"
            "          <tr><td>&nbsp;</td><th>CDS</th><th>UTRs</th><th>Overall</th></tr>\n"
            "          <tr><th class=\"left-align\">matching coefficient</th><td>%-10s</td><td>%-10s</td><td>%.3f</td></tr>\n"
            "          <tr><th class=\"left-align\">correlation coefficient</th><td>%-10s</td><td>%-10s</td><td>--</td></tr>\n"
            "          <tr><th class=\"left-align\">sensitivity</th><td>%-10s</td><td>%-10s</td><td>--</td></tr>\n"
            "          <tr><th class=\"left-align\">specificity</th><td>%-10s</td><td>%-10s</td><td>--</td></tr>\n"
            "          <tr><th class=\"left-align\">F1 Score</th><td>%-10s</td><td>%-10s</td><td>--</td></tr>\n"
            "          <tr><th class=\"left-align\">Annotation edit distance</th><td>%-10s</td><td>%-10s</td><td>--</td></tr>\n"
            "        </table>\n",
            pairstats->cds_nuc_stats.mcs, pairstats->utr_nuc_stats.mcs,
            identity, pairstats->cds_nuc_stats.ccs,
            pairstats->utr_nuc_stats.ccs, pairstats->cds_nuc_stats.sns,
            pairstats->utr_nuc_stats.sns, pairstats->cds_nuc_stats.sps,
            pairstats->utr_nuc_stats.sps, pairstats->cds_nuc_stats.f1s,
            pairstats->utr_nuc_stats.f1s, pairstats->cds_nuc_stats.eds,
            pairstats->utr_nuc_stats.eds);
  }
}

static void compare_report_html_pair_structure(FILE *outstream,
                                               AgnCompStatsBinary *stats,
                                               const char *label,
                                               const char *units)
{
  fprintf(outstream,
          "        <h3>%s structure comparison</h3>\n"
          "        <table class=\"table_normal table_extra_indent\">\n",
          label);
  if(stats->missing == 0 && stats->wrong == 0)
  {
    fprintf(outstream,
            "          <tr><td>reference %s</td><td>%lu</td></tr>\n"
            "          <tr><td>prediction %s</td><td>%lu</td></tr>\n"
            "          <tr><th class=\"left-align\" colspan=\"2\">%s structures"
            " match perfectly!</th></tr>\n",
            units, stats->correct, units, stats->correct, units);
  }
  else
  {
    fprintf(outstream,
            "          <tr><td>reference %s</td><td>%lu</td></tr>\n"
            "          <tr class=\"cell_small\"><td class=\"cell_indent\">match"
            " prediction</td><td>%lu</td></tr>\n"
            "          <tr class=\"cell_small\"><td class=\"cell_indent\">don't"
            " match prediction</td><td>%lu</td></tr>\n"
            "          <tr><td>prediction %s</td><td>%lu</td></tr>\n"
            "          <tr class=\"cell_small\"><td class=\"cell_indent\">match"
            " reference</td><td>%lu</td></tr>\n"
            "          <tr class=\"cell_small\"><td class=\"cell_indent\">don't"
            " match reference</td><td>%lu</td></tr>\n",
             units, stats->correct + stats->missing,
             stats->correct, stats->missing,
             units, stats->correct + stats->wrong,
             stats->correct, stats->wrong);
    fprintf(outstream,
            "          <tr><td>Sensitivity</td><td>%-10s</td></tr>\n"
            "          <tr><td>Specificity</td><td>%-10s</td></tr>\n"
            "          <tr><td>F1 score</td><td>%-10s</td></tr>\n"
            "          <tr><td>Annotation edit distance</td><td>%-10s</td></tr>\n",
            stats->sns, stats->sps, stats->f1s, stats->eds);
  }
  fputs("        </table>\n\n", outstream);
}


static void compare_report_html_print_pair(AgnCliquePair *pair, FILE *outstream,
                                           GtUword k)
{
  AgnTranscriptClique *refrclique = agn_clique_pair_get_refr_clique(pair);
  AgnTranscriptClique *predclique = agn_clique_pair_get_pred_clique(pair);

  fprintf(outstream,
          "      <h3 class=\"compare-header\">Comparison"
          "<a id=\"toggle_compare_%lu\" href=\"#\">(show details)</a></h3>\n",
          k);
  fprintf(outstream, "      <div id=\"compare_wrapper_%lu\""
                     "class=\"compare-wrapper\">\n", k);

  fputs("        <h3>Reference GFF3</h3>\n"
        "        <pre class=\"gff3 refr\">\n", outstream);
  agn_transcript_clique_to_gff3(refrclique, outstream, NULL);
  fputs("</pre>\n", outstream);
  fputs("        <h3>Prediction GFF3</h3>\n"
        "        <pre class=\"gff3 pred\">\n", outstream);
  agn_transcript_clique_to_gff3(predclique, outstream, NULL);
  fputs("</pre>\n", outstream);

  AgnComparison *pairstats = agn_clique_pair_get_stats(pair);
  compare_report_html_pair_structure(outstream, &pairstats->cds_struc_stats,
                                     "CDS", "CDS segments");
  compare_report_html_pair_structure(outstream, &pairstats->exon_struc_stats,
                                     "Exon", "exons");
  compare_report_html_pair_structure(outstream, &pairstats->utr_struc_stats,
                                     "UTR", "UTR segments");
  compare_report_html_pair_nucleotide(outstream, pair);  

  fputs("      </div>\n\n", outstream);
}

static void compare_report_html_sequence_handler(const AgnComparisonData *cd,
                                                 const char *seqid, void *data)
{
  CallbackData *cdata = data;
  // blah blah
  gt_array_delete(cdata->locus_summaries);
}

static void compare_report_html_summary_struc(FILE *outstream,
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
