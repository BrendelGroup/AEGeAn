/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include "extended/visitor_stream_api.h"
#include "AgnInspectVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define inspect_visitor_cast(GV)\
        gt_node_visitor_cast(inspect_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//

struct AgnInspectVisitor
{
  const GtNodeVisitor parent_instance;
  FILE *outstream;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *inspect_visitor_class();

/**
 * @function Print node info
 */
static int
visit_comment_node(GtNodeVisitor *nv, GtCommentNode  *cn, GtError *error);
static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode  *fn, GtError *error);
static int
visit_region_node(GtNodeVisitor *nv, GtRegionNode   *rn, GtError *error);
static int
visit_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn, GtError *error);
static int
visit_eof_node(GtNodeVisitor *nv, GtEOFNode      *en, GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream*
agn_inspect_stream_new(GtNodeStream *in, FILE *outstream)
{
  GtNodeVisitor *nv = agn_inspect_visitor_new(outstream);
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_inspect_visitor_new(FILE *outstream)
{
  GtNodeVisitor *nv = gt_node_visitor_create(inspect_visitor_class());
  AgnInspectVisitor *v = inspect_visitor_cast(nv);
  v->outstream = outstream;
  return nv;
}

static const GtNodeVisitorClass *inspect_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInspectVisitor), NULL,
                                    visit_comment_node, visit_feature_node,
                                    visit_region_node,  visit_sequence_node,
                                    visit_eof_node);
  }
  return nvc;
}

static int
visit_comment_node(GtNodeVisitor *nv, GtCommentNode  *cn, GtError *error)
{
  AgnInspectVisitor *iv = inspect_visitor_cast(nv);
  fprintf(iv->outstream, "[AgnInspectVisitor] comment node: '%s'\n", 
          gt_comment_node_get_comment(cn));
  return 0;
}

static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode  *fn, GtError *error)
{
  AgnInspectVisitor *iv = inspect_visitor_cast(nv);
  const char *feat_type = gt_feature_node_get_type(fn);
  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)fn);
  GtRange range = gt_genome_node_get_range((GtGenomeNode *)fn);
  fprintf(iv->outstream, "[AgnInspectVisitor] feature node: %s at %s[%lu, %lu] %p\n",
          feat_type, gt_str_get(seqid), range.start, range.end, fn);
  return 0;
}

static int
visit_region_node(GtNodeVisitor *nv, GtRegionNode   *rn, GtError *error)
{
  AgnInspectVisitor *iv = inspect_visitor_cast(nv);
  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)rn);
  GtRange range = gt_genome_node_get_range((GtGenomeNode *)rn);
  fprintf(iv->outstream, "[AgnInspectVisitor] region node: %s[%lu, %lu] %p\n",
          gt_str_get(seqid), range.start, range.end, rn);
  return 0;
}

static int
visit_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn, GtError *error)
{
  AgnInspectVisitor *iv = inspect_visitor_cast(nv);
  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)sn);
  fprintf(iv->outstream, "[AgnInspectVisitor] sequence node: %s\n",
          gt_str_get(seqid));
  return 0;
}

static int
visit_eof_node(GtNodeVisitor *nv, GtEOFNode      *en, GtError *error)
{
  AgnInspectVisitor *iv = inspect_visitor_cast(nv);
  fprintf(iv->outstream, "[AgnInspectVisitor] EOF node\n");
  return 0;
}

