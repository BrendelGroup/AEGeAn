#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "extended/feature_node_rep.h"
#include "extended/feature_node.h"
#include "extended/node_visitor.h"
#include "AgnCanonNodeVisitor.h"
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

#define agn_canon_node_visitor_cast(GV)\
        gt_node_visitor_cast(agn_canon_node_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnCanonNodeVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureIndex *index;
  AgnGeneValidator *validator;
  AgnLogger *logger;
};


//----------------------------------------------------------------------------//
// Prototype of private function(s)
//----------------------------------------------------------------------------//

/**
 * Cast a node visitor object as a AgnCanonNodeVisitor
 *
 * @returns    a node visitor object cast as a AgnCanonNodeVisitor
 */
const GtNodeVisitorClass* agn_canon_node_visitor_class();

/**
 * Destructor for the AgnCanonNodeVisitor class
 *
 * @param[in] nv    the node visitor object
 */
static void agn_canon_node_visitor_free(GtNodeVisitor *nv);

/**
 * Delete any comment nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  cn       node representing a GFF3 comment entry
 * @param[out] error    error object to which error messages, if any, will be
 *                      written
 * @returns             0 for success
 */
static int agn_canon_node_visitor_visit_comment_node(GtNodeVisitor *nv,
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
static int agn_canon_node_visitor_visit_eof_node(GtNodeVisitor *nv,
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
static int agn_canon_node_visitor_visit_feature_node(GtNodeVisitor *nv,
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
static int agn_canon_node_visitor_visit_region_node(GtNodeVisitor *nv,
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
static int agn_canon_node_visitor_visit_sequence_node(GtNodeVisitor *nv,
                                                      GtSequenceNode *sn,
                                                      GT_UNUSED GtError *error);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
const GtNodeVisitorClass* agn_canon_node_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnCanonNodeVisitor),
                                    agn_canon_node_visitor_free,
                                    agn_canon_node_visitor_visit_comment_node,
                                    agn_canon_node_visitor_visit_feature_node,
                                    agn_canon_node_visitor_visit_region_node,
                                    agn_canon_node_visitor_visit_sequence_node,
                                    agn_canon_node_visitor_visit_eof_node);
  }
  return nvc;
}

static void agn_canon_node_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED AgnCanonNodeVisitor *aiv = agn_canon_node_visitor_cast(nv);
}

GtNodeVisitor* agn_canon_node_visitor_new(GtFeatureIndex *index,
                                          AgnGeneValidator *validator,
                                          AgnLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(agn_canon_node_visitor_class());
  AgnCanonNodeVisitor *v = agn_canon_node_visitor_cast(nv);

  v->index = index;
  v->validator = validator;
  v->logger = logger;

  return nv;
}

static int agn_canon_node_visitor_visit_comment_node(GtNodeVisitor *nv,
                                                     GtCommentNode *cn,
                                                     GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)cn);
  return 0;
}

static int agn_canon_node_visitor_visit_eof_node(GtNodeVisitor *nv,
                                                 GtEOFNode *en,
                                                 GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)en);
  return 0;
}

static int agn_canon_node_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                     GtFeatureNode *fn,
                                                     GT_UNUSED GtError *error)
{
  AgnCanonNodeVisitor *v;
  gt_error_check(error);
  v = agn_canon_node_visitor_cast(nv);

  GtArray *features = gt_array_new( sizeof(GtFeatureNode *) );
  agn_gt_feature_node_resolve_pseudo_node(fn, features);
  while(gt_array_size(features) > 0)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(features);
    bool isvalid = agn_gene_validator_validate_gene(v->validator,fn,v->logger);
    if(isvalid)
    {
      bool adderror = gt_feature_index_add_feature_node(v->index, fn, error);
      if(adderror)
      {
        agn_logger_log_error(v->logger, "%s", gt_error_get(error));
        gt_error_unset(error);
        return 1;
      }
    }
    else
    {
      gt_genome_node_delete((GtGenomeNode *)fn);
      if(agn_logger_has_error(v->logger))
        return 1;
    }
  }
  gt_array_delete(features);

  return 0;
}

static int agn_canon_node_visitor_visit_region_node(GtNodeVisitor *nv,
                                                    GtRegionNode *rn,
                                                    GT_UNUSED GtError *error)
{
  AgnCanonNodeVisitor *v;
  gt_error_check(error);
  v = agn_canon_node_visitor_cast(nv);

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

static int agn_canon_node_visitor_visit_sequence_node(GtNodeVisitor *nv,
                                                      GtSequenceNode *sn,
                                                      GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)sn);
  return 0;
}
