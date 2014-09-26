/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include "AgnASInspectVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define as_inspect_visitor_cast(GV)\
        gt_node_visitor_cast(as_inspect_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definitions
//----------------------------------------------------------------------------//

struct AgnASInspectVisitor
{
  const GtNodeVisitor parent_instance;
  GtArray *as_genes;
  GtArray *se_events;
  FILE *report;
};

typedef struct
{
  GtFeatureNode *gene;
  GtRange left;
  GtRange skipped;
  GtRange right;
} SkippedExonEvent;


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function FIXME
 */
static void as_inspect_free(GtNodeVisitor *nv);

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *as_inspect_visitor_class();

/**
 * @function FIXME
 */
static GtUword check_skipped_exons(AgnASInspectVisitor *v, GtFeatureNode *gene,
                                   GtArray *mrnas);

/**
 * @function FIXME
 */
static GtUword mrna_inspection(AgnASInspectVisitor *v, GtFeatureNode *gene,
                               GtArray *mrnas);

/**
 * @function FIXME
 */
static bool skipped_exon_event_equal(SkippedExonEvent *a, SkippedExonEvent *b);

/**
 * @function FIXME
 */
static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_as_inspect_stream_new(GtNodeStream *in, FILE *report)
{
  GtNodeVisitor *nv = agn_as_inspect_visitor_new(report);
  return gt_visitor_stream_new(in, nv);
}

GtNodeVisitor *agn_as_inspect_visitor_new(FILE *report)
{
  GtNodeVisitor *nv = gt_node_visitor_create(as_inspect_visitor_class());
  AgnASInspectVisitor *v = as_inspect_visitor_cast(nv);
  v->as_genes = gt_array_new( sizeof(GtFeatureNode *) );
  v->se_events = gt_array_new( sizeof(SkippedExonEvent) );
  v->report = report == NULL ? stdout : report;
  return nv;
}

static void as_inspect_free(GtNodeVisitor *nv)
{
  AgnASInspectVisitor *v = as_inspect_visitor_cast(nv);
  gt_array_delete(v->as_genes);
  gt_array_delete(v->se_events);
}

static const GtNodeVisitorClass *as_inspect_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnASInspectVisitor),
                                    as_inspect_free, NULL, visit_feature_node,
                                    NULL, NULL, NULL);
  }
  return nvc;
}

static GtUword check_skipped_exons(AgnASInspectVisitor *v, GtFeatureNode *gene,
                                   GtArray *mrnas)
{
  GtUword i,j,k,m,n, numevents = 0;
  for(i = 0; i < gt_array_size(mrnas); i++)
  {
    GtFeatureNode *mrna1 = *(GtFeatureNode **)gt_array_get(mrnas, i);
    for(j = 0; j < gt_array_size(mrnas); j++)
    {
      if(i == j)
        continue;

      GtFeatureNode *mrna2 = *(GtFeatureNode **)gt_array_get(mrnas, j);
      GtArray *m1_exons = gt_genome_node_get_user_data((GtGenomeNode *)mrna1,
                                                       "exons");
      GtIntervalTree *m2_exon_tree = gt_genome_node_get_user_data(
                                                           (GtGenomeNode*)mrna2,
                                                           "exon_tree");
      for(k = 1; k < gt_array_size(m1_exons); k++)
      {
        GtRange nullrange = {0,0};
        SkippedExonEvent se = {gene, nullrange, nullrange, nullrange};
        GtFeatureNode *leftexon  = *(GtFeatureNode **)gt_array_get(m1_exons, k-1);
        GtRange leftrange = gt_genome_node_get_range((GtGenomeNode *)leftexon);
        GtArray *overlap = gt_array_new( sizeof(GtFeatureNode *) );
        gt_interval_tree_find_all_overlapping(m2_exon_tree, leftrange.start,
                                              leftrange.end, overlap);
        for(m = 0; m < gt_array_size(overlap); m++)
        {
          GtFeatureNode *testexon = *(GtFeatureNode **)gt_array_get(overlap, m);
          GtRange testrange = gt_genome_node_get_range((GtGenomeNode *)testexon);
          if(gt_range_compare(&leftrange, &testrange) == 0)
          {
            se.left = leftrange;
            break;
          }
        }
        gt_array_delete(overlap);
        if(gt_range_compare(&se.left, &nullrange) == 0)
          continue;
  
        GtFeatureNode *rightexon = *(GtFeatureNode **)gt_array_get(m1_exons, k);
        GtRange rightrange = gt_genome_node_get_range((GtGenomeNode *)rightexon);
        overlap = gt_array_new( sizeof(GtFeatureNode *) );
        gt_interval_tree_find_all_overlapping(m2_exon_tree, rightrange.start,
                                              rightrange.end, overlap);
        for(m = 0; m < gt_array_size(overlap); m++)
        {
          GtFeatureNode *testexon = *(GtFeatureNode **)gt_array_get(overlap, m);
          GtRange testrange = gt_genome_node_get_range((GtGenomeNode *)testexon);
          if(gt_range_compare(&rightrange, &testrange) == 0)
          {
            se.right = rightrange;
            break;
          }
        }
        gt_array_delete(overlap);
        if(gt_range_compare(&se.right, &nullrange) == 0)
          continue;
  
        overlap = gt_array_new( sizeof(GtFeatureNode *) );
        gt_interval_tree_find_all_overlapping(m2_exon_tree, leftrange.end + 1,
                                              rightrange.start - 1, overlap);
        for(m = 0; m < gt_array_size(overlap); m++)
        {
          GtFeatureNode *skipped = *(GtFeatureNode **)gt_array_get(overlap, m);
          GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)skipped);
          se.skipped = gt_genome_node_get_range((GtGenomeNode *)skipped);
          bool newevent = true;
          for(n = 0; n < gt_array_size(v->se_events); n++)
          {
            SkippedExonEvent *testevent = gt_array_get(v->se_events, n);
            GtStr *testseqid =
                      gt_genome_node_get_seqid((GtGenomeNode *)testevent->gene);
            if(gt_str_cmp(seqid, testseqid) == 0 &&
               skipped_exon_event_equal(&se, testevent))
            {
              newevent = false;
              break;
            }
          }
          if(newevent)
          {
            GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)skipped);
            const char *geneid = gt_feature_node_get_attribute(gene, "ID");
            if(geneid == NULL) geneid = "";
            fprintf(v->report, "skipped exon: %s %s %lu-%lu %lu-%lu %lu-%lu\n",
                    gt_str_get(seqid), geneid, se.left.start, se.left.end,
                    se.skipped.start, se.skipped.end,
                    se.right.start, se.right.end);
            gt_array_add(v->se_events, se);
            numevents++;
          }
        }
        gt_array_delete(overlap);
      }
    }
  }
  return numevents;
}

static GtUword mrna_inspection(AgnASInspectVisitor *v, GtFeatureNode *gene,
                               GtArray *mrnas)
{
  GtUword i;
  agn_assert(v && mrnas && gt_array_size(mrnas) >= 2);
  for(i = 0; i < gt_array_size(mrnas); i++)
  {
    GtFeatureNode *mrna = *(GtFeatureNode **)gt_array_get(mrnas, i);
    GtArray *exons = agn_typecheck_select(mrna, agn_typecheck_exon);
    if(gt_array_size(exons) > 1)
      gt_array_sort(exons, (GtCompare)agn_genome_node_compare);
    GtIntervalTree *exon_tree = gt_interval_tree_new(NULL);
    GtUword j;
    for(j = 0; j < gt_array_size(exons); j++)
    {
      GtFeatureNode *exon = *(GtFeatureNode **)gt_array_get(exons, j);
      GtRange r = gt_genome_node_get_range((GtGenomeNode *)exon);
      GtIntervalTreeNode *itn = gt_interval_tree_node_new(exon, r.start, r.end);
      gt_interval_tree_insert(exon_tree, itn);
    }
    gt_genome_node_add_user_data((GtGenomeNode *)mrna, "exons", exons,
                                 (GtFree)gt_array_delete);
    gt_genome_node_add_user_data((GtGenomeNode *)mrna, "exon_tree", exon_tree,
                                 (GtFree)gt_interval_tree_delete);
  }

  GtUword skipped_exon_count = check_skipped_exons(v, gene, mrnas);
  return skipped_exon_count;
}

static bool skipped_exon_event_equal(SkippedExonEvent *a, SkippedExonEvent *b)
{
  if(gt_range_compare(&a->left, &b->left) != 0)
    return false;
  if(gt_range_compare(&a->skipped, &b->skipped) != 0)
    return false;
  if(gt_range_compare(&a->right, &b->right) != 0)
    return false;

  return true;
}

static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error)
{
  AgnASInspectVisitor *v = as_inspect_visitor_cast(nv);
  GtArray *genes = agn_typecheck_select(fn, agn_typecheck_gene);
  if(gt_array_size(genes) == 0)
  {
    gt_array_delete(genes);
    return 0;
  }

  GtUword i;
  for(i = 0; i < gt_array_size(genes); i++)
  {
    GtFeatureNode *gene = *(GtFeatureNode **)gt_array_get(genes, i);
    GtArray *mrnas = agn_typecheck_select(gene, agn_typecheck_mrna);
    if(gt_array_size(mrnas) <= 1)
    {
      gt_array_delete(mrnas);
      continue;
    }

    GtUword numevents = mrna_inspection(v, gene, mrnas);
    if(numevents > 0)
      gt_array_add(v->as_genes, gene);
    gt_array_delete(mrnas);
  }

  gt_array_delete(genes);
  return 0;
}
