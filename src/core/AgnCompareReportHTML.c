#include <string.h>
#include "core/hashmap_api.h"
#include "AgnComparison.h"
#include "AgnCompareReportHTML.h"
#include "AgnLocus.h"
#include "AgnVersion.h"

#define compare_report_html_cast(GV)\
        gt_node_visitor_cast(compare_report_html_class(), GV)

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

struct AgnCompareReportHTML
{
  const GtNodeVisitor parent_instance;
  AgnLocusPngMetadata *pngdata;
  AgnComparisonData data;
  GtStrArray *seqids;
  const char *outdir;
  GtHashmap *seqdata;
  FILE *seqfilestream;
  GtLogger *logger;
};

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
 * @function Implement the GtNodeVisitor interface.
 */
static const GtNodeVisitorClass *compare_report_html_class();

/**
 * @function Print overall comparison statistics for a particular class of
 * feature comparisons in the summary report.
 */
static void compare_report_html_comp_class_summary(AgnCompClassDesc *summ,
                                                   GtUword num_comparisons,
                                                   const char *label,
                                                   const char *desc,
                                                   FILE *outstream);

/**
 * @function HTML footer for comparison reports.
 */
static void compare_report_html_footer(FILE *outstream);

/**
 * @function Free memory used by this node visitor.
 */
static void compare_report_html_free(GtNodeVisitor *nv);

/**
 * @function Create a report for each locus.
 */
static void compare_report_html_locus_handler(AgnCompareReportHTML *rpt,
                                              AgnLocus *locus);

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
 * @function Add the locus' information to the sequence summary page
 */
static void
compare_report_html_print_locus_to_seqfile(AgnCompareReportHTML *rpt,
                                           AgnLocus *locus);

/**
 * @function Print a comparison report for the given clique pair.
 */
static void compare_report_html_print_pair(AgnCliquePair *pair, FILE *outstream,
                                           GtUword k);

/**
 * @function FIXME
 */
static void compare_report_html_seqfile_header(AgnCompareReportHTML *rpt,
                                               const char *seqid);

/**
 * @function FIXME
 */
static void compare_report_html_seqfile_footer(AgnCompareReportHTML *rpt);

/**
 * @function Print a breakdown of characteristics of loci that fall into a
 * particular comparison class.
 */
static void compare_report_html_summary_struc(FILE *outstream,
                                              AgnCompStatsBinary *stats,
                                              const char *label,
                                              const char *units);

/**
 * @function Process feature nodes.
 */
static int compare_report_html_visit_feature_node(GtNodeVisitor *nv,
                                                  GtFeatureNode *fn,
                                                  GtError *error);

/**
 * @function Process region nodes.
 */
static int compare_report_html_visit_region_node(GtNodeVisitor *nv,
                                                 GtRegionNode *rn,
                                                 GtError *error);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_compare_report_html_create_summary(AgnCompareReportHTML *rpt,
                                            int argc, char **argv,
                                            const char *refrlabel,
                                            const char *predlabel,
                                            const char *start_time)
{
  AgnComparisonData *data = &rpt->data;
  
  char filename[1024];
  sprintf(filename, "%s/index.html", rpt->outdir);
  FILE *outstream = fopen(filename, "w");
  if(outstream == NULL)
  {
    fprintf(stderr, "error: unable to open output file '%s'\n", filename);
    exit(1);
  }

  // Print header
  fputs( "<!doctype html>\n"
         "<html lang=\"en\">\n"
         "  <head>\n"
         "    <meta charset=\"utf-8\" />\n"
         "    <title>ParsEval Summary</title>\n"
         "    <link rel=\"stylesheet\" type=\"text/css\" href=\"parseval.css\" />\n"
         "    <script type=\"text/javascript\" language=\"javascript\" src=\"vendor/jquery.js\"></script>\n"
         "    <script type=\"text/javascript\" language=\"javascript\" src=\"vendor/jquery.dataTables.js\"></script>\n"
         "    <script type=\"text/javascript\">\n"
         "      $(document).ready(function() {\n"
         "        $('#seqlist').dataTable( {\n"
         "          \"sScrollY\": \"400px\",\n"
         "          \"bPaginate\": false,\n"
         "          \"bScrollCollapse\": true,\n"
         "          \"bSort\": false,\n"
         "          \"bFilter\": false,\n"
         "          \"bInfo\": false\n"
         "        });\n"
         "      } );\n"
         "    </script>\n"
         "  </head>\n", outstream );

  // Print runtime summary
  fprintf(outstream,
          "  <body>\n"
          "    <div id=\"content\">\n"
          "      <h1>ParsEval Summary</h1>\n"
          "      <pre class=\"command\">\n"
          "Started:                %s\n"
          "Reference annotations:  %s\n"
          "Prediction annotations: %s\n"
          "Executing command:      ",
          start_time, refrlabel, predlabel);
  int x;
  for(x = 0; x < argc; x++)
  {
    fprintf(outstream, "%s ", argv[x]);
  }
  fprintf(outstream, "</pre>\n\n");

  fputs("      <h2>Sequences compared</h2>\n"
        "      <p class=\"indent\">Click on a sequence ID below to see "
        "comparison results for individual loci.</p>\n", outstream);
  fputs("      <table id=\"seqlist\" class=\"indent\">\n"
        "        <thead>\n"
        "          <tr>\n"
        "            <th>Sequence</th>\n"
        "            <th>Refr genes</th>\n"
        "            <th>Pred genes</th>\n"
        "            <th>Loci</th>\n"
        "          </tr>\n"
        "        </thead>\n"
        "        <tbody>\n", outstream);
  GtUword i;
  for(i = 0; i < gt_str_array_size(rpt->seqids); i++)
  {
    const char *seqid = gt_str_array_get(rpt->seqids, i);
    AgnComparisonData *seqdat = gt_hashmap_get(rpt->seqdata, seqid);
    agn_assert(data != NULL);
    fprintf(outstream,
            "        <tr><td><a href=\"%s/index.html\">%s</a></td>"
            "<td>%lu</td><td>%lu</td><td>%lu</td></tr>\n",
            seqid, seqid, seqdat->info.refr_genes, seqdat->info.pred_genes,
            seqdat->info.num_loci);
  }
  fputs("        </tbody>\n\n"
        "      </table>\n\n", outstream);

  compare_report_html_annot_summary(&data->info, outstream);

  fprintf(outstream, "      <h2>Comparisons</h2>\n"
                     "      <table class=\"comparisons\">\n"
                     "        <tr><th>Total comparisons</th><th>%lu</th></tr>\n",
          data->info.num_comparisons);
  compare_report_html_comp_class_summary(
      &data->summary.perfect_matches, data->info.num_comparisons,
      "perfect matches", "Prediction transcripts (exons, coding sequences, and"
      "UTRs) line up perfectly with reference transcripts.", outstream);
  compare_report_html_comp_class_summary(
      &data->summary.perfect_mislabeled, data->info.num_comparisons,
      "perfect matches with mislabeled UTRs", "5'/3' orientation of UTRs is "
      "reversed between reference and prediction, but a perfect match in all "
      "other aspects.", outstream);
  compare_report_html_comp_class_summary(
      &data->summary.cds_matches, data->info.num_comparisons,
      "CDS structure matches", "Not a perfect match, but prediction coding "
      "sequence(s) line up perfectly with reference coding sequence(s).",
      outstream);
  compare_report_html_comp_class_summary(
      &data->summary.exon_matches, data->info.num_comparisons,
      "exon structure matches", "Not a perfect match or CDS match, but "
      "prediction exon structure is identical to reference exon structure.",
      outstream);
  compare_report_html_comp_class_summary(
      &data->summary.utr_matches, data->info.num_comparisons,
      "UTR structure matches", "Not a perfect match, CDS match, or exon "
      "structure match, but prediction UTRs line up perfectly with reference "
      "UTRs.", outstream);
  compare_report_html_comp_class_summary(
      &data->summary.non_matches, data->info.num_comparisons, "non-matches",
      "Differences in CDS, exon, and UTR structure.", outstream);
  fputs("      </table>\n\n", outstream);

  compare_report_html_summary_struc(outstream, &data->stats.cds_struc_stats,
                                    "CDS", "CDS segments");
  compare_report_html_summary_struc(outstream, &data->stats.exon_struc_stats,
                                    "Exon", "exons");
  compare_report_html_summary_struc(outstream, &data->stats.utr_struc_stats,
                                    "UTR", "UTR segments");
  
  double identity = (double)data->stats.overall_matches /
                    (double)data->stats.overall_length;
  fprintf(outstream,
          "      <h3>Nucleotide-level comparison</h3>\n"
          "      <table class=\"table_wide table_extra_indent\">\n"
          "        <tr><th>&nbsp;</th><th>CDS</th><th>UTRs</th><th>Overall</th></tr>\n"
          "        <tr><th class=\"left-align\">matching coefficient</th><td>%s</td><td>%s</td><td>%.3lf</td></tr>\n"
          "        <tr><th class=\"left-align\">correlation coefficient</th><td>%s</td><td>%s</td><td>--</td></tr>\n"
          "        <tr><th class=\"left-align\">sensitivity</th><td>%s</td><td>%s</td><td>--</td></tr>\n"
          "        <tr><th class=\"left-align\">specificity</th><td>%s</td><td>%s</td><td>--</td></tr>\n"
          "        <tr><th class=\"left-align\">F1 score</th><td>%s</td><td>%s</td><td>--</td></tr>\n"
          "        <tr><th class=\"left-align\">annotation edit distance</th><td>%s</td><td>%s</td><td>--</td></tr>\n"
          "      </table>\n\n",
          data->stats.cds_nuc_stats.mcs, data->stats.utr_nuc_stats.mcs, identity,
          data->stats.cds_nuc_stats.ccs, data->stats.utr_nuc_stats.ccs,
          data->stats.cds_nuc_stats.sns, data->stats.utr_nuc_stats.sns,
          data->stats.cds_nuc_stats.sps, data->stats.utr_nuc_stats.sps,
          data->stats.cds_nuc_stats.f1s, data->stats.utr_nuc_stats.f1s,
          data->stats.cds_nuc_stats.eds, data->stats.utr_nuc_stats.eds);
  
  compare_report_html_footer(outstream);
}

GtNodeVisitor *agn_compare_report_html_new(const char *outdir,
                                           AgnLocusPngMetadata *pngdata,
                                           GtLogger *logger)
{
  GtNodeVisitor *nv = gt_node_visitor_create(compare_report_html_class());
  AgnCompareReportHTML *rpt = compare_report_html_cast(nv);
  rpt->pngdata = pngdata;
  agn_comparison_data_init(&rpt->data);
  rpt->seqids = gt_str_array_new();
  rpt->outdir = outdir;
  rpt->seqdata = gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
  rpt->seqfilestream = NULL;
  rpt->logger = logger;

  return nv;
}

static void compare_report_html_annot_summary(AgnCompInfo *info,
                                              FILE *outstream)
{
  GtUword numnotshared = info->unique_refr_loci +
                         info->unique_pred_loci;
  fprintf(outstream,
          "      <h2>Gene loci <span class=\"tooltip\">[?]<span class=\"tooltip_text\">If a gene "
          "annotation overlaps with another gene annotation, those annotations are associated "
          "with the same gene locus. See <a target=\"_blank\" "
          "href=\"http://aegean.readthedocs.org/en/refactor/loci.html\">"
          "this page</a> for a formal definition of a locus annotation.</span></span></h2>\n"
          "      <table class=\"table_normal\">\n"
          "        <tr><td>shared</td><td>%lu</td></tr>\n"
          "        <tr><td>unique to reference</td><td>%lu</td></tr>\n"
          "        <tr><td>unique to prediction</td><td>%lu</td></tr>\n"
          "        <tr><th class=\"right-align\">Total</th><td>%lu</td></tr>\n"
          "      </table>\n\n",
          info->num_loci - numnotshared, info->unique_refr_loci,
          info->unique_pred_loci, info->num_loci);

  fprintf(outstream,
          "      <h2>Reference annotations</h2>\n"
          "      <table class=\"table_normal\">\n"
          "        <tr><td>genes</td><td>%lu</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average per locus</td><td>%.3f</td></tr>\n"
          "        <tr><td>transcripts</td><td>%lu</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average per locus</td><td>%.3f</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average per gene</td><td>%.3f</td></tr>\n"
          "      </table>\n\n",
          info->refr_genes, (float)info->refr_genes / (float)info->num_loci,
          info->refr_transcripts,
          (float)info->refr_transcripts / (float)info->num_loci,
          (float)info->refr_transcripts / (float)info->refr_genes);

  fprintf(outstream,
          "      <h2>Prediction annotations</h2>\n"
          "      <table class=\"table_normal\">\n"
          "        <tr><td>genes</td><td>%lu</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average per locus</td><td>%.3f</td></tr>\n"
          "        <tr><td>transcripts</td><td>%lu</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average per locus</td><td>%.3f</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average per gene</td><td>%.3f</td></tr>\n"
          "      </table>\n\n",
          info->pred_genes, (float)info->pred_genes / (float)info->num_loci,
          info->pred_transcripts,
          (float)info->pred_transcripts / (float)info->num_loci,
          (float)info->pred_transcripts / (float)info->pred_genes);
}

static const GtNodeVisitorClass *compare_report_html_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnCompareReportHTML),
                                    compare_report_html_free, NULL,
                                    compare_report_html_visit_feature_node,
                                    compare_report_html_visit_region_node,
                                    NULL, NULL);
  }
  return nvc;
}

static void compare_report_html_comp_class_summary(AgnCompClassDesc *summ,
                                                   GtUword num_comparisons,
                                                   const char *label,
                                                   const char *desc,
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

  fprintf(outstream,
          "        <tr><td>%s <span class=\"tooltip\">"
          "<span class=\"small_tooltip\">[?]</span>"
          "<span class=\"tooltip_text\">%s</span></span></td>"
          "<td>%lu (%.1f%%)</tr>\n", label, desc, numclass, perc_class);

  if(numclass == 0)
    return;

  fprintf(outstream,
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average length</td><td>%.2lf bp</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average # refr exons</td><td>%.2lf</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average # pred exons</td><td>%.2lf</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average refr CDS length</td><td>%.2lf aa</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">average pred CDS length</td><td>%.2lf aa</td></tr>\n",
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

static void compare_report_html_free(GtNodeVisitor *nv)
{
  AgnCompareReportHTML *rpt;
  agn_assert(nv);

  rpt = compare_report_html_cast(nv);
  gt_str_array_delete(rpt->seqids);
  {
    compare_report_html_seqfile_footer(rpt);
    fclose(rpt->seqfilestream);
  }
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


static void compare_report_html_locus_handler(AgnCompareReportHTML *rpt,
                                              AgnLocus *locus)
{
  GtArray *pairs2report, *unique;
  GtUword i;
  
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtRange rng = gt_genome_node_get_range(locus);
  compare_report_html_print_locus_to_seqfile(rpt, locus);
  AgnComparisonData *seqdat = gt_hashmap_get(rpt->seqdata, gt_str_get(seqid));
  agn_locus_data_aggregate(locus, seqdat);

  char filename[1024];
  sprintf(filename, "%s/%s/%lu-%lu.html", rpt->outdir, gt_str_get(seqid),
          rng.start, rng.end);
  FILE *outstream = fopen(filename, "w");
  if(outstream == NULL)
  {
    fprintf(stderr, "error: unable to open output file '%s'\n", filename);
    exit(1);
  }
  
  compare_report_html_locus_header(locus, outstream);
  if(rpt->pngdata != NULL)
  {
    agn_locus_print_png(locus, rpt->pngdata);
    fprintf(outstream,
            "      <div class='graphic'>\n"
            "        <a href='%s_%lu-%lu.png'><img src='%s_%lu-%lu.png' /></a>\n"
            "      </div>\n\n",
            gt_str_get(seqid), rng.start, rng.end,
            gt_str_get(seqid), rng.start, rng.end);
  }
  
  pairs2report = agn_locus_pairs_to_report(locus);
  if(pairs2report == NULL || gt_array_size(pairs2report) == 0)
  {
    fputs("      <p>No comparisons were performed for this locus.</p>\n\n",
          outstream);
    compare_report_html_footer(outstream);
    return;
  }

  fputs("      <h2 class=\"bottomspace\">Comparisons</h2>\n", outstream);
  for(i = 0; i < gt_array_size(pairs2report); i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(pairs2report, i);
    compare_report_html_print_pair(pair, outstream, i);
  }

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

static void
compare_report_html_print_locus_to_seqfile(AgnCompareReportHTML *rpt,
                                           AgnLocus *locus)
{
  GtArray *pairs2report;
  GtUword i;
  GtRange rng;
  unsigned numperfect = 0, nummislabeled = 0, numcdsmatch = 0, numexonmatch = 0,
           numutrmatch = 0, numnonmatch = 0;
  char sstart[64], send[64], slength[64];

  rng = gt_genome_node_get_range(locus);
  pairs2report = agn_locus_pairs_to_report(locus);
  for(i = 0; i < gt_array_size(pairs2report); i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(pairs2report, i);
    AgnCompClassification cls = agn_clique_pair_classify(pair);
    switch(cls)
    {
      case AGN_COMP_CLASS_PERFECT_MATCH: numperfect++;    break;
      case AGN_COMP_CLASS_MISLABELED:    nummislabeled++; break;
      case AGN_COMP_CLASS_CDS_MATCH:     numcdsmatch++;   break;
      case AGN_COMP_CLASS_EXON_MATCH:    numexonmatch++;  break;
      case AGN_COMP_CLASS_UTR_MATCH:     numutrmatch++;   break;
      case AGN_COMP_CLASS_NON_MATCH:     numnonmatch++;   break;
      case AGN_COMP_CLASS_UNCLASSIFIED:  /* nothing */    break;
      default:                                            break;
    }
  }
  
  agn_sprintf_comma(rng.start, sstart);
  agn_sprintf_comma(rng.end, send);
  agn_sprintf_comma(gt_range_length(&rng), slength);
  fprintf(rpt->seqfilestream,
          "        <tr>\n"
          "          <td><a href=\"%lu-%lu.html\">(+)</a></td>\n"
          "          <td>%s</td>\n"
          "          <td>%s</td>\n"
          "          <td>%s</td>\n"
          "          <td>%lu / %lu</td>\n"
          "          <td>\n",
          rng.start, rng.end, sstart, send, slength,
          agn_locus_num_refr_mrnas(locus),
          agn_locus_num_pred_mrnas(locus));
  if(numperfect > 0)
  {
    fprintf(rpt->seqfilestream,
            "            <a class=\"pointer left20\" title=\"Perfect "
            "matches at this locus\">[P]</a> %u\n", numperfect);
  }
  if(nummislabeled > 0)
  {
    fprintf(rpt->seqfilestream,
            "            <a class=\"pointer left20\" title=\"Perfect "
            "matches at this locus with mislabeled UTRs\">[M]</a> %u\n",
            nummislabeled);
  }
  if(numcdsmatch > 0)
  {
    fprintf(rpt->seqfilestream,
            "            <a class=\"pointer left20\" title=\"CDS "
            "matches at this locus\">[C]</a> %u\n", numcdsmatch);
  }
  if(numexonmatch > 0)
  {
    fprintf(rpt->seqfilestream,
            "            <a class=\"pointer left20\" title=\"Exon "
            "structure matches at this locus\">[E]</a> %u\n",
            numexonmatch);
  }
  if(numutrmatch > 0)
  {
    fprintf(rpt->seqfilestream,
            "            <a class=\"pointer left20\" title=\"UTR "
            "matches at this locus\">[U]</a> %u\n", numutrmatch);
  }
  if(numnonmatch > 0)
  {
     fprintf(rpt->seqfilestream, "            <a class=\"pointer left20\" "
             "title=\"Non-matches at this locus\">[N]</a> %u\n",
             numnonmatch);
  }
  fprintf(rpt->seqfilestream, "          </td>\n"
          "        </tr>\n");
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

static void compare_report_html_seqfile_header(AgnCompareReportHTML *rpt,
                                               const char *seqid)
{
  fprintf(rpt->seqfilestream,
          "<!doctype html>\n"
          "<html lang=\"en\">\n"
          "  <head>\n"
          "    <meta charset=\"utf-8\" />\n"
          "    <title>ParsEval: Loci for %s</title>\n"
          "    <link rel=\"stylesheet\" type=\"text/css\" href=\"../parseval.css\" />\n"
          "    <script type=\"text/javascript\" language=\"javascript\" src=\"../vendor/jquery.js\"></script>\n"
          "    <script type=\"text/javascript\" language=\"javascript\" src=\"../vendor/jquery.dataTables.js\"></script>\n"
          "    <script type=\"text/javascript\">\n"
          "      $(document).ready(function() {\n"
          "        $('#locus_table').dataTable( {\n"
          "          \"sScrollY\": \"400px\",\n"
          "          \"bPaginate\": false,\n"
          "          \"bScrollCollapse\": true,\n"
          "          \"bSort\": false,\n"
          "          \"bFilter\": false,\n"
          "          \"bInfo\": false\n"
          "        });\n"
          "      } );\n"
          "    </script>\n"
          "  </head>\n"
          "  <body>\n"
          "    <div id=\"content\">\n"
          "      <h1>Loci for %s</h1>\n"
          "      <p><a href=\"../index.html\">⇐ Back to summary</a></p>\n\n"
          "      <p class=\"indent\">\n"
          "        Below is a list of all loci identified for sequence <strong>%s</strong>.\n"
          "        Click on the <a>(+)</a> symbol for a report of the complete comparative analysis corresponding to each locus.\n"
          "      </p>\n\n"
          "      <table class=\"loci\" id=\"locus_table\">\n"
          "        <thead>\n"
          "          <tr>\n"
          "            <th>&nbsp;</th>\n"
          "            <th>Start</th>\n"
          "            <th>End</th>\n"
          "            <th>Length</th>\n"
          "            <th>#Trans</th>\n"
          "            <th>Comparisons</th>\n"
          "          </tr>\n"
          "        </thead>\n"
          "        <tbody>\n",
          seqid, seqid, seqid);
}

static void compare_report_html_seqfile_footer(AgnCompareReportHTML *rpt)
{
  fputs("        </tbody>\n", rpt->seqfilestream);
  fputs("      </table>\n\n", rpt->seqfilestream);
  compare_report_html_footer(rpt->seqfilestream);
  fputs("    </div>\n", rpt->seqfilestream);
  fputs("  </body>\n", rpt->seqfilestream);
  fputs("</html>\n", rpt->seqfilestream);
}

static void compare_report_html_summary_struc(FILE *outstream,
                                              AgnCompStatsBinary *stats,
                                              const char *label,
                                              const char *units)
{
  fprintf(outstream, "      <h3>%s structure comparison</h3>\n", label);
  
  GtUword refrcnt = stats->correct + stats->missing;
  GtUword predcnt = stats->correct + stats->wrong;
  fprintf(outstream,
          "      <table class=\"table_normal table_extra_indent\">\n"
          "        <tr><td>reference %s</td><td>%lu</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">match prediction</td><td>%lu (%.1f%%)</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">don't match prediction</td><td>%lu (%.1f%%)</td></tr>\n"
          "        <tr><td>prediction %s</td><td>%lu</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">match prediction</td><td>%lu (%.1f%%)</td></tr>\n"
          "        <tr class=\"cell_small\"><td class=\"cell_indent\">don't match prediction</td><td>%lu (%.1f%%)</td></tr>\n"
          "        <tr><td>Sensitivity</td><td>%.3f</td></tr>\n"
          "        <tr><td>Specificity</td><td>%.3f</td></tr>\n"
          "        <tr><td>F1 score</td><td>%.3f</td></tr>\n"
          "        <tr><td>Annotation edit distance</td><td>%.3f</td></tr>\n"
          "      </table>\n\n",
          units, stats->correct + stats->missing,
          stats->correct, (float)stats->correct / (float)refrcnt * 100,
          stats->missing, (float)stats->missing / (float)refrcnt * 100,
          units, stats->correct + stats->wrong,
          stats->correct, (float)stats->correct / (float)predcnt * 100,
          stats->wrong, (float)stats->wrong / (float)predcnt * 100,
          stats->sn, stats->sp, stats->f1, stats->ed);
}

static int compare_report_html_visit_feature_node(GtNodeVisitor *nv,
                                                  GtFeatureNode *fn,
                                                  GtError *error)
{
  AgnCompareReportHTML *rpt;
  AgnLocus *locus;

  gt_error_check(error);
  agn_assert(nv && fn && gt_feature_node_has_type(fn, "locus"));

  rpt = compare_report_html_cast(nv);
  locus = (AgnLocus *)fn;
  agn_locus_comparative_analysis(locus, rpt->logger);
  agn_locus_data_aggregate(locus, &rpt->data);
  compare_report_html_locus_handler(rpt, locus);

  return 0;
}

static int compare_report_html_visit_region_node(GtNodeVisitor *nv,
                                                 GtRegionNode *rn,
                                                 GtError *error)
{
  AgnCompareReportHTML *rpt;
  AgnComparisonData *data;
  GtStr *seqidstr;
  const char *seqid;
  char seqfilename[512];

  gt_error_check(error);
  agn_assert(nv && rn);

  rpt = compare_report_html_cast(nv);
  seqidstr = gt_genome_node_get_seqid((GtGenomeNode *)rn);
  seqid = gt_cstr_dup(gt_str_get(seqidstr));
  gt_str_array_add(rpt->seqids, seqidstr);
  data = gt_malloc( sizeof(AgnComparisonData) );
  agn_comparison_data_init(data);
  gt_hashmap_add(rpt->seqdata, (char *)seqid, data);

  char seqdircmd[AGN_MAX_FILENAME_SIZE];
  sprintf(seqdircmd, "mkdir %s/%s", rpt->outdir, seqid);
  if(system(seqdircmd))
  {
    fprintf(stderr, "error: could not create directory %s/%s\n", rpt->outdir,
            seqid);
    exit(1);
  }

  if(rpt->seqfilestream)
  {
    compare_report_html_seqfile_footer(rpt);
    fclose(rpt->seqfilestream);
  }
  sprintf(seqfilename, "%s/%s/index.html", rpt->outdir, seqid);
  rpt->seqfilestream = fopen(seqfilename, "w");
  if(!rpt->seqfilestream)
  {
    fprintf(stderr, "error: unable to open %s\n", seqfilename);
    exit(1);
  }
  compare_report_html_seqfile_header(rpt, seqid);

  return 0;
}
