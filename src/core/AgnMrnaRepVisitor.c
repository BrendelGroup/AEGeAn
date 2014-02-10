#include <string.h>
#include "core/array_api.h"
#include "AgnFilterStream.h"
#include "AgnMrnaRepVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Determine the length of an mRNA's coding sequence.
 */
static GtUword mrna_rep_cds_length(GtFeatureNode *mrna);

/**
 * @function Remove feature ``fn`` and all its subfeatures from ``root``.
 */
static void mrna_rep_remove_tree(GtFeatureNode *root, GtFeatureNode *fn);

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *mrna_rep_visitor_class();

/**
 * @function Generate data for unit testing.
 */
static void mrna_rep_visitor_test_data(GtQueue *queue);

/**
 * @function Identify any mRNA subfeatures associated with this top-level
 * feature and apply the CDS inference procedure.
 */
static int
mrna_rep_visit_feature_node(GtNodeVisitor *nv,GtFeatureNode *fn,GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_mrna_rep_stream_new(GtNodeStream *in)
{
  GtNodeVisitor *nv = agn_mrna_rep_visitor_new();
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_mrna_rep_visitor_new(GtLogger *logger)
{
  return gt_node_visitor_create(mrna_rep_visitor_class());
}

bool agn_mrna_rep_visitor_unit_test(AgnUnitTest *test)
{
  mrna_rep_visitor_test_data(NULL);
  return false;
}

static GtUword mrna_rep_cds_length(GtFeatureNode *mrna)
{
  GtUword totallength = 0;
  const char *cdsid = NULL;
  GtArray *cds_segments = agn_typecheck_select(mrna, agn_typecheck_cds);
  while(gt_array_size(cds_segments) > 0)
  {
    GtGenomeNode **segment = gt_array_pop(cds_segments);
    GtFeatureNode *segmentfn = gt_feature_node_cast(*segment);
    const char *fid = gt_feature_node_get_attribute(segmentfn, "ID");
    if(fid)
    {
      if(cdsid == NULL)
        cdsid = fid;
      else
        gt_assert(strcmp(cdsid, fid) == 0);
    }
    totallength += gt_genome_node_get_length(*segment);
  }
  gt_array_delete(cds_segments);
  return totallength;
}

static void mrna_rep_remove_tree(GtFeatureNode *root, GtFeatureNode *fn)
{
  gt_assert(root && fn);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *child;
  for(child = gt_feature_node_iterator_next(iter);
      child != NULL;
      child = gt_feature_node_iterator_next(iter))
  {
    mrna_rep_remove_tree(root, child);
    gt_feature_node_remove_leaf(fn, child);
  }
  gt_feature_node_iterator_delete(iter);
  gt_feature_node_remove_leaf(root, fn);
}

static const GtNodeVisitorClass *mrna_rep_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (GtNodeVisitor), NULL, NULL,
                                    mrna_rep_visit_feature_node, NULL, NULL,
                                    NULL);
  }
  return nvc;
}

static void mrna_rep_visitor_test_data(GtQueue *queue)
{
  /*GtError *error = gt_error_new();
  const char *file = "data/gff3/grape-codons.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtNodeStream *rep_stream = agn_mrna_rep_stream_new(gff3in);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(rep_stream, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnMrnaRepVisitor::mrna_rep_visitor_test_data] error "
            "processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(rep_stream);
  gt_node_stream_delete(arraystream);
  gt_array_sort(feats, (GtCompare)agn_genome_node_compare);
  gt_array_reverse(feats);
  while(gt_array_size(feats) > 0)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(feats);
    gt_queue_add(queue, fn);
  }
  gt_array_delete(feats);
  gt_error_delete(error);*/
}

static int
mrna_rep_visit_feature_node(GtNodeVisitor *nv,GtFeatureNode *fn,GtError *error)
{
  gt_error_check(error);

  GtUword i,j;
  GtArray *genes = agn_typecheck_select(fn, agn_typecheck_gene);
  for(i = 0; i < gt_array_size(genes); i++)
  {
    GtGenomeNode **gene = gt_array_get(genes, i);
    GtFeatureNode *genefn = gt_feature_node_cast(*gene);
    GtArray *mrnas = agn_typecheck_select(genefn, agn_typecheck_mrna);
    if(gt_array_size(mrnas) <= 1)
    {
      gt_array_delete(mrnas);
      continue;
    }

    // Find the longest mRNA
    GtFeatureNode *longest_mrna = NULL;
    GtUword longest_length = 0;
    for(j = 0; j < gt_array_size(mrnas); j++)
    {
      GtGenomeNode **mrna = gt_array_get(mrnas, j);
      GtFeatureNode *mrnafn = gt_feature_node_cast(*mrna);
      GtUword length = mrna_rep_cds_length(mrnafn);
      if(longest_length == 0 || length > longest_length)
      {
        longest_mrna = mrnafn;
        longest_length = length;
      }
    }

    // Now, remove all other mRNAs
    for(j = 0; j < gt_array_size(mrnas); j++)
    {
      GtFeatureNode *mrna = *(GtFeatureNode **)gt_array_get(mrnas, j);
      if(mrna != longest_mrna)
        mrna_rep_remove_tree(fn, mrna);
    }
    gt_array_delete(mrnas);
  }
  gt_array_delete(genes);
  
  return 0;
}
