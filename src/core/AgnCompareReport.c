#include "AgnCompareReport.h"
#include "AgnTypecheck.h"

#define compare_report_visitor_cast(GV)\
        gt_node_visitor_cast(compare_report_visitor_class(), GV)

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

struct AgnCompareReport
{
  const GtNodeVisitor parent_instance;
  AgnComparisonData data;
  GtStrArray *seqids;
  GtHashmap *seqdata;
  const char *last_seqid;
  GtLogger *logger;
  AgnCompareReportLocusFunc locusfunc;
  void *locusfuncdata;
  AgnCompareReportSequenceFunc sequencefunc;
  void *sequencefuncdata;
};

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Free class memory when parent destructor is called.
 */
static void compare_report_free(GtNodeVisitor *nv);

/**
 * @function Aggregate comparison statistics for ``locus`` either at the
 * sequence level or over all sequences.
 */
void compare_report_record_locus_stats(AgnComparisonData *data,AgnLocus *locus);

/**
 * @function Implement the ``GtNodeVisitor`` interface.
 */
static const GtNodeVisitorClass *compare_report_visitor_class();

/**
 * @function Callback function to apply to each feature node (AgnLocus object)
 * in the node stream. Performs a comparative analysis of the two alternative
 * sources of annotation at that locus.
 */
static int compare_report_visit_feature_node(GtNodeVisitor *nv,
                                             GtFeatureNode *fn, GtError *error);

/**
 * @function Callback function to apply to each region node in the node stream.
 * Stores sequence IDs and sequence-level comparison statistics.
 */
static int compare_report_visit_region_node(GtNodeVisitor *nv,
                                            GtRegionNode *rn, GtError *error);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

AgnComparisonData *agn_compare_report_data(AgnCompareReport *rpt)
{
  gt_assert(rpt);
  return &rpt->data;
}

GtNodeVisitor *agn_compare_report_new(GtLogger *logger)
{
  GtNodeVisitor *nv = gt_node_visitor_create(compare_report_visitor_class());
  AgnCompareReport *rpt = compare_report_visitor_cast(nv);
  agn_comparison_data_init(&rpt->data);
  rpt->seqids = gt_str_array_new();
  rpt->last_seqid = NULL;
  rpt->seqdata = gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
  rpt->logger = logger;
  rpt->locusfunc = NULL;
  rpt->locusfuncdata = NULL;
  rpt->sequencefunc = NULL;
  rpt->sequencefuncdata = NULL;

  return nv;
}

GtStrArray *agn_compare_report_seqids(AgnCompareReport *rpt)
{
  return rpt->seqids;
}

void agn_compare_report_set_locus_callback(AgnCompareReport *rpt,
                                           AgnCompareReportLocusFunc func,
                                           void *data)
{
  gt_assert(rpt);
  rpt->locusfunc = func;
  rpt->locusfuncdata = data;
}

void agn_compare_report_set_sequence_callback(AgnCompareReport *rpt,
                                              AgnCompareReportSequenceFunc func,
                                              void *data)
{
  gt_assert(rpt);
  rpt->sequencefunc = func;
  rpt->sequencefuncdata = data;
}

static void compare_report_free(GtNodeVisitor *nv)
{
  AgnCompareReport *rpt;
  gt_assert(nv);

  rpt = compare_report_visitor_cast(nv);
  if(rpt->last_seqid && rpt->sequencefunc)
  {
    AgnComparisonData *seqdata = gt_hashmap_get(rpt->seqdata, rpt->last_seqid);
    rpt->sequencefunc(seqdata, rpt->last_seqid, rpt->sequencefuncdata);
  }

  gt_str_array_delete(rpt->seqids);
  gt_hashmap_delete(rpt->seqdata);
}

void compare_report_record_locus_stats(AgnComparisonData *data, AgnLocus *locus)
{
  GtUword numrefrgenes, numpredgenes;
  GtUword i;
  gt_assert(data && locus);

  data->info.num_loci++;
  numrefrgenes = agn_locus_num_refr_genes(locus);
  numpredgenes = agn_locus_num_pred_genes(locus);
  data->info.refr_genes += numrefrgenes;
  data->info.pred_genes += numpredgenes;
  data->info.refr_transcripts += agn_locus_num_refr_mrnas(locus);
  data->info.pred_transcripts += agn_locus_num_pred_mrnas(locus);
  if(numrefrgenes > 0 && numpredgenes == 0)
    data->info.unique_refr_loci++;
  if(numpredgenes > 0 && numrefrgenes == 0)
    data->info.unique_pred_loci++;

  agn_locus_comparison_aggregate(locus, &data->stats);
  agn_comparison_resolve(&data->stats);

  GtRange locusrange = gt_genome_node_get_range(locus);
  GtArray *pairs2report = agn_locus_pairs_to_report(locus);
  data->info.num_comparisons += gt_array_size(pairs2report);
  for(i = 0; i < gt_array_size(pairs2report); i++)
  {
    AgnCliquePair **pair = gt_array_get(pairs2report, i);
    AgnCompClassification c = agn_clique_pair_classify(*pair);
    AgnCompClassDesc *desc;
    switch(c)
    {
      case AGN_COMP_CLASS_PERFECT_MATCH:
        desc = &data->summary.perfect_matches;
        break;
      case AGN_COMP_CLASS_MISLABELED:
        desc = &data->summary.perfect_mislabeled;
        break;
      case AGN_COMP_CLASS_CDS_MATCH:
        desc = &data->summary.cds_matches;
        break;
      case AGN_COMP_CLASS_EXON_MATCH:
        desc = &data->summary.exon_matches;
        break;
      case AGN_COMP_CLASS_UTR_MATCH:
        desc = &data->summary.utr_matches;
        break;
      case AGN_COMP_CLASS_NON_MATCH:
        desc = &data->summary.non_matches;
        break;
      default:
        desc = NULL;
        fprintf(stderr, "error: unknow comp classification %d\n", c);
        break;
    }
    desc->comparison_count++;
    desc->total_length += gt_range_length(&locusrange);

    AgnTranscriptClique *rclique = agn_clique_pair_get_refr_clique(*pair);
    AgnTranscriptClique *pclique = agn_clique_pair_get_pred_clique(*pair);
    desc->refr_cds_length += agn_transcript_clique_cds_length(rclique);
    desc->pred_cds_length += agn_transcript_clique_cds_length(pclique);
    GtFeatureNode *rcliquefn = gt_feature_node_cast(rclique);
    GtFeatureNode *pcliquefn = gt_feature_node_cast(pclique);
    desc->refr_exon_count += agn_typecheck_count(rcliquefn, agn_typecheck_exon);
    desc->pred_exon_count += agn_typecheck_count(pcliquefn, agn_typecheck_exon);
  }
}

static const GtNodeVisitorClass *compare_report_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnCompareReport),
                                    compare_report_free, NULL,
                                    compare_report_visit_feature_node,
                                    compare_report_visit_region_node,
                                    NULL, NULL);
  }
  return nvc;
}

static int compare_report_visit_feature_node(GtNodeVisitor *nv,
                                             GtFeatureNode *fn, GtError *error)
{
  AgnCompareReport *rpt;
  AgnComparisonData *seqdata;
  AgnLocus *locus;
  GtStr *seqid;

  gt_error_check(error);
  gt_assert(nv && fn && gt_feature_node_has_type(fn, "locus"));

  rpt = compare_report_visitor_cast(nv);
  locus = (AgnLocus *)fn;

  agn_locus_comparative_analysis(locus, rpt->logger);
  compare_report_record_locus_stats(&rpt->data, locus);

  seqid = gt_genome_node_get_seqid((GtGenomeNode *)fn);
  seqdata = gt_hashmap_get(rpt->seqdata, gt_str_get(seqid));
  compare_report_record_locus_stats(seqdata, locus);

  if(rpt->locusfunc != NULL)
    rpt->locusfunc(locus, rpt->locusfuncdata);

  return 0;
}

static int compare_report_visit_region_node(GtNodeVisitor *nv,
                                            GtRegionNode *rn, GtError *error)
{
  AgnCompareReport *rpt;
  AgnComparisonData *data;
  GtStr *seqidstr;
  const char *seqid;

  gt_error_check(error);
  gt_assert(nv && rn);

  rpt = compare_report_visitor_cast(nv);
  seqidstr = gt_genome_node_get_seqid((GtGenomeNode *)rn);
  seqid = gt_cstr_dup(gt_str_get(seqidstr));
  gt_str_array_add(rpt->seqids, seqidstr);
  data = gt_malloc( sizeof(AgnComparisonData) );
  agn_comparison_data_init(data);
  gt_hashmap_add(rpt->seqdata, (char *)seqid, data);

  if(rpt->last_seqid && rpt->sequencefunc)
  {
    AgnComparisonData *seqdata = gt_hashmap_get(rpt->seqdata, seqid);
    rpt->sequencefunc(seqdata, seqid, rpt->sequencefuncdata);
  }
  rpt->last_seqid = seqid;

  return 0;
}
