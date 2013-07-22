#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "AgnSimpleNodeVisitor.h"
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

#define agn_simple_node_visitor_cast(GV)\
        gt_node_visitor_cast(agn_simple_node_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnSimpleNodeVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureIndex *index;
  const char *type;
  AgnLogger *logger;
};


//----------------------------------------------------------------------------//
// Prototype of private function(s)
//----------------------------------------------------------------------------//

/**
 * Cast a node visitor object as a AgnSimpleNodeVisitor
 *
 * @returns    a node visitor object cast as a AgnSimpleNodeVisitor
 */
const GtNodeVisitorClass* agn_simple_node_visitor_class();

/**
 * Destructor for the AgnSimpleNodeVisitor class
 *
 * @param[in] nv    the node visitor object
 */
static void agn_simple_node_visitor_free(GtNodeVisitor *nv);

/**
 * Delete any comment nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  cn       node representing a GFF3 comment entry
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success
 */
static int agn_simple_node_visitor_visit_comment_node(GtNodeVisitor *nv,
                                                     GtCommentNode *cn,
                                                     GT_UNUSED GtError *error);

/**
 * Delete any EOF nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  en       node representing the end of a GFF3 file
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success
 */
static int agn_simple_node_visitor_visit_eof_node(GtNodeVisitor *nv,
                                                 GtEOFNode *en,
                                                 GT_UNUSED GtError *error);

/**
 * Validate and store any feature nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  fn       node representing a GFF3 feature entry
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success
 */
static int agn_simple_node_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                     GtFeatureNode *fn,
                                                     GT_UNUSED GtError *error);

/**
 * Store any region nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  rn       node representing a GFF3 sequence-region entry
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success
 */
static int agn_simple_node_visitor_visit_region_node(GtNodeVisitor *nv,
                                                    GtRegionNode *rn,
                                                    GT_UNUSED GtError *error);

/**
 * Delete any sequence nodes encountered while loading the data
 * At some point these will be handled differently, but for now this is it
 *
 * @param[in]  nv       a node visitor
 * @param[in]  sn       node representing a Fasta sequence
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success
 */
static int agn_simple_node_visitor_visit_sequence_node(GtNodeVisitor *nv,
                                                      GtSequenceNode *sn,
                                                      GT_UNUSED GtError *error);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
const GtNodeVisitorClass* agn_simple_node_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnSimpleNodeVisitor),
                                    agn_simple_node_visitor_free,
                                    agn_simple_node_visitor_visit_comment_node,
                                    agn_simple_node_visitor_visit_feature_node,
                                    agn_simple_node_visitor_visit_region_node,
                                    agn_simple_node_visitor_visit_sequence_node,
                                    agn_simple_node_visitor_visit_eof_node);
  }
  return nvc;
}

static void agn_simple_node_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED AgnSimpleNodeVisitor *aiv = agn_simple_node_visitor_cast(nv);
}

GtNodeVisitor* agn_simple_node_visitor_new(GtFeatureIndex *index,
                                           const char *type,
                                           AgnLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(agn_simple_node_visitor_class());
  AgnSimpleNodeVisitor *v = agn_simple_node_visitor_cast(nv);

  v->index = index;
  v->type = type;
  v->logger = logger;

  return nv;
}

static int agn_simple_node_visitor_visit_comment_node(GtNodeVisitor *nv,
                                                      GtCommentNode *cn,
                                                      GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)cn);
  return 0;
}

static int agn_simple_node_visitor_visit_eof_node(GtNodeVisitor *nv,
                                                 GtEOFNode *en,
                                                 GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)en);
  return 0;
}

static int agn_simple_node_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                      GtFeatureNode *fn,
                                                      GT_UNUSED GtError *error)
{
  AgnSimpleNodeVisitor *v;
  gt_error_check(error);
  v = agn_simple_node_visitor_cast(nv);

  if(gt_feature_node_has_type(fn, v->type))
  {
    gt_feature_index_add_feature_node(v->index, fn, error);
    return 0;
  }

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    if(gt_feature_node_has_type(feature, v->type))
    {
      gt_genome_node_ref((GtGenomeNode *)feature);
      if(feature != fn)
        agn_gt_feature_node_remove_tree(fn, feature);
      gt_feature_index_add_feature_node(v->index, feature, error);
    }
  }
  gt_feature_node_iterator_delete(iter);
  gt_genome_node_delete((GtGenomeNode *)fn);
  return 0;
}

static int agn_simple_node_visitor_visit_region_node(GtNodeVisitor *nv,
                                                     GtRegionNode *rn,
                                                     GT_UNUSED GtError *error)
{
  AgnSimpleNodeVisitor *v;
  gt_error_check(error);
  v = agn_simple_node_visitor_cast(nv);

  if(gt_feature_index_add_region_node(v->index, rn, error))
  {
    GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)rn);
    agn_logger_log_error(v->logger, "unable to add region node for seqid "
                         "'%s': %s", seqid, gt_error_get(error));
    gt_error_unset(error);
    return 1;
  }

  return 0;
}

static int agn_simple_node_visitor_visit_sequence_node(GtNodeVisitor *nv,
                                                       GtSequenceNode *sn,
                                                       GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)sn);
  return 0;
}
