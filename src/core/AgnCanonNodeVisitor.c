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
// Perhaps add an AgnError/AgnLogger object?
struct AgnCanonNodeVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureIndex *index;
  AgnGeneValidator *validator;
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
                                          AgnGeneValidator *validator)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(agn_canon_node_visitor_class());
  AgnCanonNodeVisitor *v = agn_canon_node_visitor_cast(nv);

  v->index = index;
  v->validator = validator;

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
  AgnLogger *logger = agn_logger_new();

  if(gt_feature_node_is_pseudo(fn))
  {
    GtArray *features = gt_array_new( sizeof(GtFeatureNode *) );
    agn_gt_feature_node_resolve_pseudo_node(fn, features);
    while(gt_array_size(features) > 0)
    {
      GtFeatureNode *fn2add = *(GtFeatureNode **)gt_array_pop(features);
      if(agn_gene_validator_validate_gene(v->validator, fn2add, logger))
      {
        if(gt_feature_index_add_feature_node(v->index, fn2add, error))
        {
          return 1;
        }
      }
      else
      {
        gt_genome_node_delete((GtGenomeNode *)fn2add);
        // FIXME this should not be handled here
        bool haderror = agn_logger_print_all(logger, stderr,
                            "[ParsEval] validating gene '%s'",
                            gt_feature_node_get_attribute(fn2add, "ID"));
        if(haderror)
        {
          gt_error_set(error, "    fatal error(s)");
          return -1;
        }
      }
    }
  }
  else
  {
    if(agn_gene_validator_validate_gene(v->validator, fn, logger))
    {
      if(gt_feature_index_add_feature_node(v->index, fn, error))
      {
        return 1;
      }
    }
    else
    {
      gt_genome_node_delete((GtGenomeNode *)fn);
      // FIXME this should not be handled here
      bool haderror = agn_logger_print_all(logger, stderr,
                          "[ParsEval] validating gene '%s'",
                          gt_feature_node_get_attribute(fn, "ID"));
      if(haderror)
      {
        gt_error_set(error, "    fatal error(s)");
        return -1;
      }
    }
  }
  agn_logger_delete(logger);

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
    // FIXME this should not be handled here
    fprintf(stderr, "[ParsEval] error: issue loading GFF3 data into memory: %s",
            gt_error_get(error));
    exit(1);
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
