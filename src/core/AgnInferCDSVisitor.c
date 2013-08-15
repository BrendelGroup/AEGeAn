#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "AgnInferCDSVisitor.h"
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

#define agn_infer_cds_visitor_cast(GV)\
        gt_node_visitor_cast(agn_infer_cds_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnInferCDSVisitor
{
  const GtNodeVisitor parent_instance;
  AgnLogger *logger;
};


//----------------------------------------------------------------------------//
// Prototype of private function(s)
//----------------------------------------------------------------------------//

/**
 * Cast a node visitor object as a AgnInferCDSVisitor
 *
 * @returns    a node visitor object cast as a AgnInferCDSVisitor
 */
const GtNodeVisitorClass* agn_infer_cds_visitor_class();

/**
 * Destructor for the AgnInferCDSVisitor class
 *
 * @param[in] nv    the node visitor object
 */
static void infer_cds_visitor_free(GtNodeVisitor *nv);

/**
 * Identify any mRNA subfeatures associated with this top-level feature and
 * apply the CDS inference procedure.
 *
 * @param[in]  nv       a node visitor
 * @param[in]  fn       node representing a top-level GFF3 feature entry
 * @param[out] error    stores error messages
 * @returns             0 for success
 */
static int infer_cds_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                        GtError *error);

/**
 * Infer CDS for any mRNAs that have no CDS specified but do have exons and
 * start/stop codons explicitly specified.
 * 
 * @param[in]  nv      a node visitor
 * @param[in]  mrna    an mRNA feature
 */
void infer_cds_visit_mrna(GtNodeVisitor *nv, GtFeatureNode *mrna);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
const GtNodeVisitorClass* agn_infer_cds_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInferCDSVisitor),
                                    infer_cds_visitor_free,
                                    NULL,
                                    infer_cds_visit_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

static void infer_cds_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED AgnInferCDSVisitor *aiv = agn_infer_cds_visitor_cast(nv);
}

GtNodeVisitor* agn_infer_cds_visitor_new(AgnLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(agn_infer_cds_visitor_class());
  AgnInferCDSVisitor *v = agn_infer_cds_visitor_cast(nv);
  v->logger = logger;
  return nv;
}

static int infer_cds_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                        GtError *error)
{
  gt_error_check(error);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(agn_gt_feature_node_is_mrna_feature(current))
      infer_cds_visit_mrna(nv, current);
  }
  gt_feature_node_iterator_delete(iter);
  return 0;
}

void infer_cds_visit_mrna(GtNodeVisitor *nv, GtFeatureNode *mrna)
{
  AgnInferCDSVisitor *v = agn_infer_cds_visitor_cast(nv);

  GtArray *cds   = agn_gt_feature_node_children_of_type(mrna, "CDS");
  GtArray *exons = agn_gt_feature_node_children_of_type(mrna, "exon");
  GtArray *starts = agn_gt_feature_node_children_of_type(mrna, "start_codon");
  GtArray *stops  = agn_gt_feature_node_children_of_type(mrna, "stop_codon");
  GtFeatureNode *start_codon, *stop_codon;

  bool cdsexplicit      = gt_array_size(cds)   > 0;
  bool exonsexplicit    = gt_array_size(exons) > 0;
  bool start_codon_sane = gt_array_size(starts) == 1 &&
                          (start_codon = gt_array_get(starts, 0)) != NULL;
  bool stop_codon_sane  = gt_array_size(stops)  == 1 &&
                          (stop_codon  = gt_array_get(stops,  0)) != NULL;

  if(cdsexplicit || !exonsexplicit || !start_codon_sane || !stop_codon_sane)
    return;

  GtRange left_codon_range, right_codon_range;
  left_codon_range  = gt_genome_node_get_range((GtGenomeNode *)start_codon);
  right_codon_range = gt_genome_node_get_range((GtGenomeNode *)stop_codon);
  if(gt_feature_node_get_strand(mrna) == GT_STRAND_REVERSE)
  {
    left_codon_range  = gt_genome_node_get_range((GtGenomeNode *)stop_codon);
    right_codon_range = gt_genome_node_get_range((GtGenomeNode *)start_codon);
  }
  unsigned long i;
  for(i = 0; i < gt_array_size(exons); i++)
  {
    GtFeatureNode *exon = *(GtFeatureNode **)gt_array_get(exons, i);
    GtGenomeNode *exon_gn = (GtGenomeNode *)exon;
    GtRange exon_range = gt_genome_node_get_range(exon_gn);
    GtStrand exon_strand = gt_feature_node_get_strand(exon);

    GtRange cdsrange;
    bool exon_includes_cds = agn_infer_cds_range_from_exon_and_codons(
                                 &exon_range, &left_codon_range,
                                 &right_codon_range, &cdsrange);
    if(exon_includes_cds)
    {
      GtGenomeNode *cdsfeat = gt_feature_node_new
      (
        gt_genome_node_get_seqid(exon_gn), "CDS", cdsrange.start,
        cdsrange.end, exon_strand
      );
      gt_feature_node_add_child(mrna, (GtFeatureNode *)cdsfeat);
      gt_feature_node_add_attribute((GtFeatureNode *)cdsfeat, "Parent",
                                    gt_feature_node_get_attribute(mrna, "ID"));
      gt_array_add(cds, cdsfeat);
    }
  }

  if(gt_array_size(cds) == 0)
  {
    agn_logger_log_error(v->logger, "error inferring CDS from codons and exons "
                         "for mRNA '%s' (line %u)",
                         gt_feature_node_get_attribute(mrna, "ID"),
                         gt_genome_node_get_line_number((GtGenomeNode *)mrna));
  }
}
