#include <string.h>
#include "core/array_api.h"
#include "AgnPseudogeneFixVisitor.h"
#include "AgnTypecheck.h"


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *pseudogene_fix_visitor_class();

/**
 * @function Generate data for unit testing.
 */
static void pseudogene_fix_visitor_test_data(GtQueue *queue);

/**
 * @function For any gene feature with attribute 'pseudo=true', set type to
 * 'pseudogene'.
 */
static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_pseudogene_fix_stream_new(GtNodeStream *in)
{
  GtNodeVisitor *nv = agn_pseudogene_fix_visitor_new();
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_pseudogene_fix_visitor_new(GtLogger *logger)
{
  GtNodeVisitor *nv = gt_node_visitor_create(pseudogene_fix_visitor_class());
  //AgnPseudogeneFixVisitor *v = pseudogene_fix_visitor_cast(nv);
  return nv;
}

bool agn_pseudogene_fix_visitor_unit_test(AgnUnitTest *test)
{
  const char *type;
  GtFeatureNode *fn;
  GtGenomeNode *gn;
  GtQueue *queue = gt_queue_new();
  GtRange frange;
  pseudogene_fix_visitor_test_data(queue);

  gn = gt_queue_get(queue);
  fn = gt_feature_node_cast(gn);
  frange = gt_genome_node_get_range(gn);
  type = gt_feature_node_get_type(fn);
//fprintf(stderr, "%s[%lu, %lu]\n", type, frange.start, frange.end);
  bool t1 = frange.start == 1 && frange.end == 1000 && strcmp(type,"gene") == 0;
  agn_unit_test_result(test, "Gene, no pseudo attribute", t1);
  gt_genome_node_delete(gn);

  gn = gt_queue_get(queue);
  fn = gt_feature_node_cast(gn);
  frange = gt_genome_node_get_range(gn);
  type = gt_feature_node_get_type(fn);
  bool t2 = (frange.start == 2000 && frange.end == 3000 &&
             strcmp(type, "gene") == 0);
  agn_unit_test_result(test, "Gene, bogus pseudo attribute", t2);
  gt_genome_node_delete(gn);

  gn = gt_queue_get(queue);
  fn = gt_feature_node_cast(gn);
  frange = gt_genome_node_get_range(gn);
  type = gt_feature_node_get_type(fn);
  bool t3 = (frange.start == 4000 && frange.end == 5000 &&
             strcmp(type, "pseudogene") == 0);
  agn_unit_test_result(test, "Corrected pseudogene", t3);
  gt_genome_node_delete(gn);

  gn = gt_queue_get(queue);
  fn = gt_feature_node_cast(gn);
  frange = gt_genome_node_get_range(gn);
  type = gt_feature_node_get_type(fn);
  bool t4 = (frange.start == 6000 && frange.end == 7000 &&
             strcmp(type, "pseudogene") == 0);
  agn_unit_test_result(test, "Uncorrected pseudogene", t4);
  gt_genome_node_delete(gn);

  GtGenomeNode *locus = gt_queue_get(queue);
  GtFeatureNode *locusfn = gt_feature_node_cast(locus);
  GtArray *genes = agn_typecheck_select(locusfn, agn_typecheck_gene);
  bool t5a = gt_array_size(genes) == 1;
  if(t5a)
  {
    gn = *(GtGenomeNode **)gt_array_pop(genes);
    frange = gt_genome_node_get_range(gn);
    t5a = frange.start == 10000 && frange.end == 11500;
  }
  gt_array_delete(genes);
  genes = agn_typecheck_select(locusfn, agn_typecheck_pseudogene);
  bool t5b = gt_array_size(genes) == 1;
  if(t5b)
  {
    gn = *(GtGenomeNode **)gt_array_pop(genes);
    frange = gt_genome_node_get_range(gn);
    t5b = frange.start == 11000 && frange.end == 12000;
  }
  gt_array_delete(genes);
  agn_unit_test_result(test, "iLocus with gene & pseudogene", t5a && t5b);
  gt_genome_node_delete(locus);

  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static const GtNodeVisitorClass *pseudogene_fix_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (GtNodeStream), NULL, NULL,
                                    visit_feature_node, NULL, NULL, NULL);
  }
  return nvc;
}

static void pseudogene_fix_visitor_test_data(GtQueue *queue)
{
  GtArray *source, *dest;
  GtFeatureNode *fn1, *fn2;
  GtGenomeNode *gn1, *gn2;
  GtNodeStream *in, *fix, *out;
  GtStr *seqid = gt_str_new_cstr("chr");
  source = gt_array_new( sizeof(GtGenomeNode *) );
  dest = gt_array_new( sizeof(GtGenomeNode *) );

  gn1 = gt_feature_node_new(seqid, "gene", 1, 1000, GT_STRAND_REVERSE);
  gt_array_add(source, gn1);

  gn1 = gt_feature_node_new(seqid, "gene", 2000, 3000, GT_STRAND_FORWARD);
  fn1 = gt_feature_node_cast(gn1);
  gt_feature_node_add_attribute(fn1, "pseudo", "false");
  gt_array_add(source, gn1);

  gn1 = gt_feature_node_new(seqid, "gene", 4000, 5000, GT_STRAND_FORWARD);
  fn1 = gt_feature_node_cast(gn1);
  gt_feature_node_add_attribute(fn1, "pseudo", "true");
  gt_array_add(source, gn1);

  gn1 = gt_feature_node_new(seqid, "pseudogene", 6000, 7000, GT_STRAND_FORWARD);
  gt_array_add(source, gn1);

  gn1 = gt_feature_node_new(seqid, "locus", 10000, 12000, GT_STRAND_BOTH);
  fn1 = gt_feature_node_cast(gn1);
  gn2 = gt_feature_node_new(seqid, "gene", 10000, 11500, GT_STRAND_REVERSE);
  fn2 = gt_feature_node_cast(gn2);
  gt_feature_node_add_child(fn1, fn2);
  gn2 = gt_feature_node_new(seqid, "gene", 11000, 12000, GT_STRAND_REVERSE);
  fn2 = gt_feature_node_cast(gn2);
  gt_feature_node_add_attribute(fn2, "pseudo", "true");
  gt_feature_node_add_child(fn1, fn2);
  gt_array_add(source, gn1);

  gt_str_delete(seqid);

  GtUword n;
  GtError *error = gt_error_new();
  in = gt_array_in_stream_new(source, &n, error);
  fix = agn_pseudogene_fix_stream_new(in);
  out = gt_array_out_stream_new(fix, dest, error);
  int pullresult = gt_node_stream_pull(out, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnPsuedogeneFixVisitor::pseudogene_fix_visitor_test_"
            "data] error processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(in);
  gt_node_stream_delete(fix);
  gt_node_stream_delete(out);
  gt_error_delete(error);

  gt_array_reverse(dest);
  while(gt_array_size(dest) > 0)
  {
    gn1 = *(GtGenomeNode **)gt_array_pop(dest);
    gt_queue_add(queue, gn1);
  }
  gt_array_delete(dest);
  gt_array_delete(source);
}

static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error)
{
  gt_error_check(error);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_typecheck_gene(current))
      continue;

    const char *attrvalue = gt_feature_node_get_attribute(current, "pseudo");
    if(attrvalue && strcmp(attrvalue, "true") == 0)
      gt_feature_node_set_type(current, "pseudogene");
  }
  gt_feature_node_iterator_delete(iter);
  return 0;
}
