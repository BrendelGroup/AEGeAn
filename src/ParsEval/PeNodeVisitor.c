#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "extended/feature_node_rep.h"
#include "extended/feature_node.h"
#include "extended/node_visitor.h"
#include "AgnGtExtensions.h"
#include "PeNodeVisitor.h"
#include "AgnUtils.h"

#define pe_node_visitor_cast(GV)\
        gt_node_visitor_cast(pe_node_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct PeNodeVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureIndex *index;
  AgnGeneValidator *validator;
};


//----------------------------------------------------------------------------//
// Prototype of private function(s)
//----------------------------------------------------------------------------//

/**
 * Cast a node visitor object as a PeNodeVisitor
 *
 * @returns    a node visitor object cast as a PeNodeVisitor
 */
const GtNodeVisitorClass* pe_node_visitor_class();

/**
 * Destructor for the PeNodeVisitor class
 *
 * @param[in] nv    the node visitor object
 */
static void pe_node_visitor_free(GtNodeVisitor *nv);

/**
 * Delete any comment nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  cn       node representing a GFF3 comment entry
 * @param[out] error    error object to which error messages, if any, will be written
 * @returns             0 for success
 */
static int pe_node_visitor_visit_comment_node(GtNodeVisitor *nv, GtCommentNode *cn, GT_UNUSED GtError *error);

/**
 * Delete any EOF nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  en       node representing the end of a GFF3 file
 * @param[out] error    error object to which error messages, if any, will be written
 * @returns             0 for success
 */
static int pe_node_visitor_visit_eof_node(GtNodeVisitor *nv, GtEOFNode *en, GT_UNUSED GtError *error);

/**
 * Validate and store any feature nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  fn       node representing a GFF3 feature entry
 * @param[out] error    error object to which error messages, if any, will be written
 * @returns             0 for success
 */
static int pe_node_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GT_UNUSED GtError *error);

/**
 * Store any region nodes encountered while loading the data
 *
 * @param[in]  nv       a node visitor
 * @param[in]  rn       node representing a GFF3 sequence-region entry
 * @param[out] error    error object to which error messages, if any, will be written
 * @returns             0 for success
 */
static int pe_node_visitor_visit_region_node(GtNodeVisitor *nv, GtRegionNode *rn, GT_UNUSED GtError *error);

/**
 * Delete any sequence nodes encountered while loading the data
 * At some point these will be handled differently, but for now this is it
 *
 * @param[in]  nv       a node visitor
 * @param[in]  sn       node representing a Fasta sequence
 * @param[out] error    error object to which error messages, if any, will be written
 * @returns             0 for success
 */
static int pe_node_visitor_visit_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn, GT_UNUSED GtError *error);


//------------------------------------------------------------------------------------------------//
// Method implementations
//------------------------------------------------------------------------------------------------//
const GtNodeVisitorClass* pe_node_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new( sizeof (PeNodeVisitor),
                                     pe_node_visitor_free,
                                     pe_node_visitor_visit_comment_node,
                                     pe_node_visitor_visit_feature_node,
                                     pe_node_visitor_visit_region_node,
                                     pe_node_visitor_visit_sequence_node,
                                     pe_node_visitor_visit_eof_node );
  }
  return nvc;
}

static void pe_node_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED PeNodeVisitor *aiv = pe_node_visitor_cast(nv);
}

GtNodeVisitor* pe_node_visitor_new(GtFeatureIndex *index, AgnGeneValidator *validator)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(pe_node_visitor_class());
  PeNodeVisitor *v = pe_node_visitor_cast(nv);

  v->index = index;
  v->validator = validator;

  return nv;
}

static int pe_node_visitor_visit_comment_node(GtNodeVisitor *nv, GtCommentNode *cn, GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)cn);
  return 0;
}

static int pe_node_visitor_visit_eof_node(GtNodeVisitor *nv, GtEOFNode *en, GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)en);
  return 0;
}

static int pe_node_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GT_UNUSED GtError *error)
{
  PeNodeVisitor *v;
  gt_error_check(error);
  v = pe_node_visitor_cast(nv);
  AgnError *myerror = agn_error_new();

  if(gt_feature_node_is_pseudo(fn))
  {
    GtArray *features = gt_array_new( sizeof(GtFeatureNode *) );
    agn_gt_feature_node_resolve_pseudo_node(fn, features);
    while(gt_array_size(features) > 0)
    {
      GtFeatureNode *fn2add = *(GtFeatureNode **)gt_array_pop(features);
      if( agn_gene_validator_validate_gene(v->validator, fn2add, myerror) )
      {
        if(gt_feature_index_add_feature_node(v->index, fn2add, error))
        {
          return 1;
        }
      }
      else
      {
        gt_genome_node_delete((GtGenomeNode *)fn2add);
        if(agn_error_is_set(myerror))
        {
          agn_error_print( myerror, stderr,
                           "[ParsEval] issue validating gene '%s'",
                           gt_feature_node_get_attribute(fn2add, "ID") );
          if(agn_error_is_fatal(myerror))
          {
            gt_error_set(error, "    fatal error");
            return -1;
          }
          agn_error_unset(myerror);
        }
      }
    }
  }
  else
  {
    if( agn_gene_validator_validate_gene(v->validator, fn, myerror) )
    {
      if(gt_feature_index_add_feature_node(v->index, fn, error))
      {
        return 1;
      }
    }
    else
    {
      gt_genome_node_delete((GtGenomeNode *)fn);
      if(agn_error_is_set(myerror))
      {
        agn_error_print( myerror, stderr,
                         "[ParsEval] issue validating gene '%s'",
                         gt_feature_node_get_attribute(fn, "ID") );
        if(agn_error_is_fatal(myerror))
        {
          gt_error_set(error, "    fatal error");
          return -1;
        }
        agn_error_unset(myerror);
      }
    }
  }
  agn_error_delete(myerror);

  return 0;
}

static int pe_node_visitor_visit_region_node(GtNodeVisitor *nv, GtRegionNode *rn, GT_UNUSED GtError *error)
{
  PeNodeVisitor *v;
  gt_error_check(error);
  v = pe_node_visitor_cast(nv);

  if(gt_feature_index_add_region_node(v->index, rn, error))
  {
    fprintf(stderr, "[ParsEval] error: issue loading GFF3 data into memory: %s", gt_error_get(error));
    exit(1);
  }

  return 0;
}

static int pe_node_visitor_visit_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn, GT_UNUSED GtError *error)
{
  gt_error_check(error);
  gt_genome_node_delete((GtGenomeNode *)sn);
  return 0;
}
