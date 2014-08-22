/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
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
  bool test1 = true;
  while(gt_queue_size(queue) > 0)
  {
    GtFeatureNode *fn = gt_queue_get(queue);
    GtUword nc = gt_feature_node_number_of_children(fn);
    test1 = test1 && nc == 0 && gt_feature_node_has_type(fn, "gene");
    gt_genome_node_delete((GtGenomeNode *)fn);
  }
  agn_unit_test_result(test, "Test 1: grape (Gaze)", test);

  gt_queue_delete(queue);
  return agn_unit_test_success(test);
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
  GtError *error = gt_error_new();
  const char *file = "data/gff3/grape-refr.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtNodeStream *rcs = agn_remove_children_stream_new(gff3in);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(rcs, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnRemoveChildrenVisitor::remove_children_visitor_test_"
            "data] error processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(rcs);
  gt_node_stream_delete(arraystream);

  while(gt_array_size(feats) > 0)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(feats);
    gt_queue_add(queue, fn);
  }
  gt_array_delete(feats);
  gt_error_delete(error);
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
