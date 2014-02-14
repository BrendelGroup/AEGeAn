#include "extended/visitor_stream_api.h"
#include "AgnNodeDeleteVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *node_delete_visitor_class();

/**
 * @function Call ``gt_genome_node_delete`` on all feature nodes.
 */
static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_node_delete_stream_new(GtNodeStream *in)
{
  GtNodeVisitor *nv = agn_node_delete_visitor_new();
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_node_delete_visitor_new()
{
  return gt_node_visitor_create(node_delete_visitor_class());
}

static const GtNodeVisitorClass *node_delete_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (GtNodeVisitor), NULL, NULL,
                                    visit_feature_node, NULL, NULL, NULL);
  }
  return nvc;
}

static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error)
{
  gt_genome_node_delete((GtGenomeNode *)fn);
  return 0;
}
