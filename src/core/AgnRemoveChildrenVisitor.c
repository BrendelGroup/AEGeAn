#include "AgnRemoveChildrenVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *remove_children_visitor_class();

/**
 * @function Generate data for unit testing.
 */
static void remove_children_visitor_test_data(GtQueue *queue);

/**
 * @function If ``fn`` is not a pseudofeature, remove all of its children.
 */
static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_remove_children_stream_new(GtNodeStream *in)
{
  GtNodeVisitor *nv = agn_remove_children_visitor_new();
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_remove_children_visitor_new(GtLogger *logger)
{
  return gt_node_visitor_create(remove_children_visitor_class());
}

bool agn_remove_children_visitor_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  remove_children_visitor_test_data(queue);

  gt_queue_delete(queue);
  return false;
}

static const GtNodeVisitorClass *remove_children_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (GtNodeVisitor), NULL, NULL,
                                    visit_feature_node, NULL, NULL, NULL);
  }
  return nvc;
}

static void remove_children_visitor_test_data(GtQueue *queue)
{
}

static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error)
{
  gt_error_check(error);
  if(gt_feature_node_is_pseudo(fn))
    return 0;

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *child;
  for(child  = gt_feature_node_iterator_next(iter);
      child != NULL;
      child  = gt_feature_node_iterator_next(iter))
  {
    agn_feature_node_remove_tree(fn, child);
  }
  gt_feature_node_iterator_delete(iter);
  return 0;
}
