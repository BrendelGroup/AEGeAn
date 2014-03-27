#include "AgnCliquePair.h"
#include "AgnCompareTextReportVisitor.h"
#include "AgnLocus.h"
#include "AgnTranscriptClique.h"

#define compare_text_report_visitor_cast(GV)\
        gt_node_visitor_cast(compare_text_report_visitor_class(), GV)

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

struct AgnCompareTextReportVisitor
{
  const GtNodeVisitor parent_instance;
  GtUword max_locus_transcripts;
  FILE *fp_reports;
  FILE *fp_summary;
  GtLogger *logger;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *compare_text_report_visitor_class();

/**
 * @function For each locus, perform comparative analysis and report the
 * comparison statistics.
 */
static int ctrv_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                   GtError *error);

/**
 * @function FIXME
 */
static void locus_report_print_geneids(AgnCompareTextReportVisitor *v,
                                       AgnLocus *locus);

/**
 * @function FIXME
 */
static void locus_report_print_pair(AgnCompareTextReportVisitor *v,
                                    AgnCliquePair *pair);

/**
 * @function FIXME
 */
static void print_locus_report(AgnCompareTextReportVisitor *v, AgnLocus *locus);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream *agn_compare_text_report_stream_new(GtNodeStream *instream,
                                                 FILE *reports,
                                                 FILE *summary,
                                                 GtLogger *logger)
{
  GtNodeVisitor*
    nv = agn_compare_text_report_visitor_new(reports, summary, logger);
  GtNodeStream *ns = gt_visitor_stream_new(instream, nv);
  return ns;
}

GtNodeVisitor *agn_compare_text_report_visitor_new(FILE *reports,
                                                   FILE *summary,
                                                   GtLogger *logger)
{
  GtNodeVisitor*
    nv = gt_node_visitor_create(compare_text_report_visitor_class());
  AgnCompareTextReportVisitor *v = compare_text_report_visitor_cast(nv);
  v->max_locus_transcripts = 0;
  v->fp_reports = reports;
  v->fp_summary = summary;
  v->logger = logger;
  return nv;
}

void agn_compare_text_report_visitor_trans_max(AgnCompareTextReportVisitor *v,
                                               GtUword max_locus_transcripts)
{
  v->max_locus_transcripts = max_locus_transcripts;
}

static const GtNodeVisitorClass *compare_text_report_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnCompareTextReportVisitor), NULL,
                                    NULL, ctrv_visit_feature_node, NULL, NULL,
                                    NULL);
  }
  return nvc;
}

static int ctrv_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                   GtError *error)
{
  AgnCompareTextReportVisitor *v = compare_text_report_visitor_cast(nv);
  gt_error_check(error);

  gt_assert(gt_feature_node_has_type(fn, "locus"));
  AgnLocus *locus = (AgnLocus *)fn;
  agn_locus_comparative_analysis(locus, v->max_locus_transcripts, 0, v->logger);
  print_locus_report(v, locus);

  return 0;
}

static void locus_report_print_geneids(AgnCompareTextReportVisitor *v,
                                       AgnLocus *locus)
{
  GtArray *genes;

  genes = agn_locus_refr_gene_ids(locus);
  fprintf(v->fp_reports, "  |  reference genes:\n");
  if(gt_array_size(genes) == 0)
    fprintf(v->fp_reports, "|    None!\n");
  while(gt_array_size(genes) > 0)
  {
    const char **geneid = gt_array_pop(genes);
    fprintf(v->fp_reports, "|    %s\n", *geneid);
  }
  fprintf(v->fp_reports, "|\n");
  gt_array_delete(genes);

  genes = agn_locus_pred_gene_ids(locus);
  fprintf(v->fp_reports, "  |  prediction genes:\n");
  if(gt_array_size(genes) == 0)
    fprintf(v->fp_reports, "|    None!\n");
  while(gt_array_size(genes) > 0)
  {
    const char **geneid = gt_array_pop(genes);
    fprintf(v->fp_reports, "|    %s\n", *geneid);
  }
  fprintf(v->fp_reports, "|\n");
  gt_array_delete(genes);
}

static void locus_report_print_pair(AgnCompareTextReportVisitor *v,
                                    AgnCliquePair *pair)
{
  
}

static void print_locus_report(AgnCompareTextReportVisitor *v, AgnLocus *locus)
{
  GtRange range = gt_genome_node_get_range(locus);
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  fprintf(v->fp_reports,
          "|-------------------------------------------------\n"
          "|---- Locus: %s_%lu-%lu\n"
          "|-------------------------------------------------\n"
          "|\n",
          gt_str_get(seqid), range.start, range.end);

  locus_report_print_geneids(v, locus);
  fprintf(v->fp_reports,
          "|\n"
          "|----------\n");

  GtArray *pairs2report = agn_locus_pairs_to_report(locus);
  if(pairs2report == NULL || gt_array_size(pairs2report) == 0)
  {
    fprintf(v->fp_reports,
            "     |\n"
            "     | No comparisons were performed for this locus.\n"
            "     |\n");
    return;
  }

  GtUword i;
  for(i = 0; i < gt_array_size(pairs2report); i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(pairs2report, i);
    fprintf(v->fp_reports,
            "     |\n"
            "     |--------------------------\n"
            "     |---- Begin comparison ----\n"
            "     |--------------------------\n"
            "     |\n");
    locus_report_print_pair(v, pair);
  }
}
