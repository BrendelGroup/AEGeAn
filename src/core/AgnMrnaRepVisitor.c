/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include <string.h>
#include "core/array_api.h"
#include "AgnFilterStream.h"
#include "AgnMrnaRepVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define mrna_rep_visitor_cast(GV)\
        gt_node_visitor_cast(mrna_rep_visitor_class(), GV)


//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnMrnaRepVisitor
{
  const GtNodeVisitor parent_instance;
  char *parenttype;
  FILE *mapstream;
  bool useacc;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *mrna_rep_visitor_class();

/**
 * @function Release memory.
 */
static void mrna_rep_visitor_free(GtNodeVisitor *nv);

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

GtNodeStream* agn_mrna_rep_stream_new(GtNodeStream *in, FILE *mapstream,
                                      bool useacc)
{
  GtNodeVisitor *nv = agn_mrna_rep_visitor_new(mapstream);
  if(useacc)
  {
    AgnMrnaRepVisitor *v = mrna_rep_visitor_cast(nv);
    agn_mrna_rep_visitor_use_accession(v);
  }
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_mrna_rep_visitor_new(FILE *mapstream)
{
  GtNodeVisitor *nv = gt_node_visitor_create(mrna_rep_visitor_class());
  AgnMrnaRepVisitor *v = mrna_rep_visitor_cast(nv);
  v->parenttype = gt_cstr_dup("gene");
  v->mapstream = mapstream;
  v->useacc = false;
  return nv;
}

void agn_mrna_rep_visitor_set_parent_type(AgnMrnaRepVisitor *v,
                                          const char *type)
{
  gt_free(v->parenttype);
  v->parenttype = gt_cstr_dup(type);
}

void agn_mrna_rep_visitor_use_accession(AgnMrnaRepVisitor *v)
{
  agn_assert(v);
  v->useacc = true;
}

bool agn_mrna_rep_visitor_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  mrna_rep_visitor_test_data(queue);
  agn_assert(gt_queue_size(queue) == 1);

  GtGenomeNode *gene = gt_queue_get(queue);
  GtFeatureNode *genefn = gt_feature_node_cast(gene);
  GtArray *mrnas = agn_typecheck_select(genefn, agn_typecheck_mrna);
  bool test1 = gt_array_size(mrnas) == 1;
  if(test1)
  {
    GtGenomeNode **mrna = gt_array_get(mrnas, 0);
    GtFeatureNode *mrnafn = gt_feature_node_cast(*mrna);
    GtUword cdslength = agn_mrna_cds_length(mrnafn);
    GtRange range = gt_genome_node_get_range(*mrna);
    test1 = cdslength == 738 && range.start == 5928 && range.end == 8737;
  }
  agn_unit_test_result(test, "TAIR10: test 1", test1);
  gt_genome_node_delete(gene);
  gt_array_delete(mrnas);

  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static const GtNodeVisitorClass *mrna_rep_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnMrnaRepVisitor),
                                    mrna_rep_visitor_free, NULL,
                                    mrna_rep_visit_feature_node, NULL, NULL,
                                    NULL);
  }
  return nvc;
}

static void mrna_rep_visitor_free(GtNodeVisitor *nv)
{
  AgnMrnaRepVisitor *v = mrna_rep_visitor_cast(nv);
  gt_free(v->parenttype);
}

static void mrna_rep_visitor_test_data(GtQueue *queue)
{
  GtError *error = gt_error_new();
  const char *file = "data/gff3/tair-altsplice.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtNodeStream *rep_stream = agn_mrna_rep_stream_new(gff3in, NULL, false);
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
  gt_error_delete(error);
}

static int
mrna_rep_visit_feature_node(GtNodeVisitor *nv,GtFeatureNode *fn,GtError *error)
{
  gt_error_check(error);
  AgnMrnaRepVisitor *v = mrna_rep_visitor_cast(nv);

  GtUword i,j;
  GtArray *genes = agn_typecheck_select_str(fn, v->parenttype);
  for(i = 0; i < gt_array_size(genes); i++)
  {
    GtGenomeNode **gene = gt_array_get(genes, i);
    GtFeatureNode *genefn = gt_feature_node_cast(*gene);
    const char *gid = gt_feature_node_get_attribute(genefn, "accession");
    if(gid == NULL || v->useacc == false)
      gid = gt_feature_node_get_attribute(genefn, "Name");
    if(gid == NULL || v->useacc == false)
      gid = gt_feature_node_get_attribute(genefn, "ID");
    GtArray *mrnas = agn_typecheck_select(genefn, agn_typecheck_mrna);
    if(gt_array_size(mrnas) <= 1)
    {
      if(v->mapstream != NULL && gt_array_size(mrnas) == 1)
      {
        GtFeatureNode **mrna = gt_array_pop(mrnas);
        const char *mid = gt_feature_node_get_attribute(*mrna, "accession");
        if(mid == NULL || v->useacc == false)
        {
          mid = gt_feature_node_get_attribute(*mrna, "ID");
          if(mid == NULL)
            mid = gt_feature_node_get_attribute(*mrna, "Parent");
        }
        fprintf(v->mapstream, "%s\t%s\n", gid, mid);
      }
      gt_array_delete(mrnas);
      continue;
    }

    // Find the longest mRNA
    const char *longest_id = NULL;
    GtFeatureNode *longest_mrna = NULL;
    GtUword longest_length = 0;
    for(j = 0; j < gt_array_size(mrnas); j++)
    {
      GtGenomeNode **mrna = gt_array_get(mrnas, j);
      GtFeatureNode *mrnafn = gt_feature_node_cast(*mrna);
      GtUword length = agn_mrna_cds_length(mrnafn);
      const char *mid = gt_feature_node_get_attribute(mrnafn, "accession");
      if(mid == NULL || v->useacc == false)
      {
        mid = gt_feature_node_get_attribute(mrnafn, "ID");
        if(mid == NULL)
          mid = gt_feature_node_get_attribute(mrnafn, "Parent");
      }
      bool preempt = (length == longest_length && strcmp(mid, longest_id) < 0);
      if(longest_length == 0 || length > longest_length || preempt)
      {
        longest_mrna = mrnafn;
        longest_id = mid;
        longest_length = length;
      }
    }
    if(v->mapstream != NULL)
      fprintf(v->mapstream, "%s\t%s\n", gid, longest_id);

    // Now, remove all other mRNAs
    for(j = 0; j < gt_array_size(mrnas); j++)
    {
      GtFeatureNode *mrna = *(GtFeatureNode **)gt_array_get(mrnas, j);
      if(mrna != longest_mrna)
        agn_feature_node_remove_tree(fn, mrna);
    }
    gt_array_delete(mrnas);
  }
  gt_array_delete(genes);

  return 0;
}
