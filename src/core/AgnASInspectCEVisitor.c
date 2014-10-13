/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "AgnASInspectCEVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define as_inspect_ce_visitor_cast(GV)\
        gt_node_visitor_cast(as_inspect_ce_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definitions
//----------------------------------------------------------------------------//

struct AgnASInspectCEVisitor
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
static void as_inspect_ce_free(GtNodeVisitor *nv);

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *as_inspect_ce_visitor_class();

/**
 * @function FIXME
 */
static GtUword check_skipped_exons(AgnASInspectCEVisitor *v,
                                   GtFeatureNode *gene, GtArray *mrnas);

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

GtNodeStream* agn_as_inspect_ce_stream_new(GtNodeStream *in, FILE *report)
{
  GtNodeVisitor *nv = agn_as_inspect_ce_visitor_new(report);
  return gt_visitor_stream_new(in, nv);
}

GtNodeVisitor *agn_as_inspect_ce_visitor_new(FILE *report)
{
  GtNodeVisitor *nv = gt_node_visitor_create(as_inspect_ce_visitor_class());
  AgnASInspectCEVisitor *v = as_inspect_ce_visitor_cast(nv);
  v->as_genes = gt_array_new( sizeof(GtFeatureNode *) );
  v->se_events = gt_array_new( sizeof(SkippedExonEvent) );
  v->report = report;
  return nv;
}

bool agn_as_inspect_ce_visitor_unit_test(AgnUnitTest *test)
{
  GtError *error = gt_error_new();
  const char *infile = "data/gff3/as-ce.gff3";
  GtNodeStream *annot = gt_gff3_in_stream_new_unsorted(1, &infile);
  GtNodeVisitor *nv = agn_as_inspect_ce_visitor_new(NULL);
  AgnASInspectCEVisitor *v = as_inspect_ce_visitor_cast(nv);
  GtNodeStream *as = gt_visitor_stream_new(annot, nv);
  int result = gt_node_stream_pull(as, error);
  if(result == -1)
  {
    fprintf(stderr, "[AgnASInspectCEVisitor::agn_as_inspect_ce_visitor_unit_"
            "test] error loading annotations into memory: %s\n",
            gt_error_get(error));
  }

  agn_assert(gt_array_size(v->as_genes) == 2);
  GtFeatureNode **gene1 = gt_array_get(v->as_genes, 0);
  const char *geneid = gt_feature_node_get_attribute(*gene1, "ID");
  bool test1 = strcmp(geneid, "gene2") == 0;
  GtFeatureNode **gene2 = gt_array_get(v->as_genes, 1);
  geneid = gt_feature_node_get_attribute(*gene2, "ID");
  test1 = test1 && strcmp(geneid, "gene3") == 0;
  agn_unit_test_result(test, "skipped exons: genes (GFF3)", test1);

  agn_assert(gt_array_size(v->se_events) == 3);
  SkippedExonEvent *event = gt_array_get(v->se_events, 0);
  SkippedExonEvent testevent1 = { *gene1, {615,838}, {952,1115}, {1198,1419} };
  bool test2 = skipped_exon_event_equal(event, &testevent1);
  event = gt_array_get(v->se_events, 1);
  SkippedExonEvent testevent2 = { *gene2, {615,838}, {952,1115}, {1198,1419} };
  test2 = test2 && skipped_exon_event_equal(event, &testevent2);
  event = gt_array_get(v->se_events, 2);
  SkippedExonEvent testevent3 = { *gene2, {952,1115}, {1198,1419}, {2067,2350}};
  test2 = test2 && skipped_exon_event_equal(event, &testevent3);
  agn_unit_test_result(test, "skipped exons: events (GFF3)", test2);

  gt_node_stream_delete(annot);
  gt_node_stream_delete(as);


  infile = "data/gff3/as-ce.gtf";
  annot = gt_gtf_in_stream_new(infile);
  GtNodeStream *sort = gt_sort_stream_new(annot);
  nv = agn_as_inspect_ce_visitor_new(NULL);
  v = as_inspect_ce_visitor_cast(nv);
  as = gt_visitor_stream_new(sort, nv);
  result = gt_node_stream_pull(as, error);
  if(result == -1)
  {
    fprintf(stderr, "[AgnASInspectCEVisitor::agn_as_inspect_ce_visitor_unit_"
            "test] error loading annotations into memory: %s\n",
            gt_error_get(error));
  }

  gt_assert(gt_array_size(v->as_genes) == 2);
  agn_assert(gt_array_size(v->se_events) == 3);
  gene1 = gt_array_get(v->as_genes, 0);
  gene2 = gt_array_get(v->as_genes, 1);

  event = gt_array_get(v->se_events, 0);
  SkippedExonEvent testevent4 = { *gene1, {488,528}, {952,1115}, {1198,1419} };
  bool test3 = skipped_exon_event_equal(event, &testevent4);
  event = gt_array_get(v->se_events, 1);
  SkippedExonEvent testevent5 = { *gene2, {488,528}, {952,1115}, {1198,1419} };
  test3 = test3 && skipped_exon_event_equal(event, &testevent5);
  event = gt_array_get(v->se_events, 2);
  SkippedExonEvent testevent6 = { *gene2, {952,1115}, {1198,1419}, {2067,2350}};
  test3 = test3 && skipped_exon_event_equal(event, &testevent6);
  agn_unit_test_result(test, "skipped exons: events (GTF)", test3);

  gt_error_delete(error);
  gt_node_stream_delete(annot);
  gt_node_stream_delete(sort);
  gt_node_stream_delete(as);

  return agn_unit_test_success(test);
}

static void as_inspect_ce_free(GtNodeVisitor *nv)
{
  AgnASInspectCEVisitor *v = as_inspect_ce_visitor_cast(nv);
  while(gt_array_size(v->as_genes) > 0)
  {
    GtGenomeNode **gene = gt_array_pop(v->as_genes);
    gt_genome_node_delete(*gene);
  }
  gt_array_delete(v->as_genes);
  gt_array_delete(v->se_events);
}

static const GtNodeVisitorClass *as_inspect_ce_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnASInspectCEVisitor),
                                    as_inspect_ce_free, NULL,
                                    visit_feature_node, NULL, NULL, NULL);
  }
  return nvc;
}

static GtUword check_exon_pair(AgnASInspectCEVisitor *v, GtFeatureNode *gene,
                               GtFeatureNode *leftexon,GtFeatureNode *rightexon,
                               GtFeatureNode *mrna2)
{
  agn_assert(v && gene && leftexon && rightexon && mrna2);
  GtIntervalTree *m2_exon_tree = agn_mrna_exon_tree(mrna2);
  GtRange nullrange = {0,0};
  SkippedExonEvent se = {gene, nullrange, nullrange, nullrange};
  GtUword i,j, numevents = 0;

  GtRange leftrange = gt_genome_node_get_range((GtGenomeNode *)leftexon);
  GtArray *overlap = gt_array_new( sizeof(GtGenomeNode *) );
  gt_interval_tree_find_all_overlapping(m2_exon_tree, leftrange.start,
                                        leftrange.end, overlap);
  for(i = 0; i < gt_array_size(overlap); i++)
  {
    GtGenomeNode *testexon = *(GtGenomeNode **)gt_array_get(overlap, i);
    GtRange testrange = gt_genome_node_get_range(testexon);
    if(gt_range_compare(&leftrange, &testrange) == 0)
    {
      se.left = leftrange;
      break;
    }
  }
  gt_array_delete(overlap);
  if(gt_range_compare(&se.left, &nullrange) == 0)
    return 0;

  GtRange rightrange = gt_genome_node_get_range((GtGenomeNode *)rightexon);
  overlap = gt_array_new( sizeof(GtGenomeNode *) );
  gt_interval_tree_find_all_overlapping(m2_exon_tree, rightrange.start,
                                        rightrange.end, overlap);
  for(i = 0; i < gt_array_size(overlap); i++)
  {
    GtGenomeNode *testexon = *(GtGenomeNode **)gt_array_get(overlap, i);
    GtRange testrange = gt_genome_node_get_range(testexon);
    if(gt_range_compare(&rightrange, &testrange) == 0)
    {
      se.right = rightrange;
      break;
    }
  }
  gt_array_delete(overlap);
  if(gt_range_compare(&se.right, &nullrange) == 0)
    return 0;

  overlap = gt_array_new( sizeof(GtFeatureNode *) );
  gt_interval_tree_find_all_overlapping(m2_exon_tree, leftrange.end + 1,
                                        rightrange.start - 1, overlap);
  for(i = 0; i < gt_array_size(overlap); i++)
  {
    GtGenomeNode *skipped = *(GtGenomeNode **)gt_array_get(overlap, i);
    se.skipped = gt_genome_node_get_range(skipped);
    bool newevent = true;
    for(j = 0; j < gt_array_size(v->se_events); j++)
    {
      SkippedExonEvent *testevent = gt_array_get(v->se_events, j);
      if(skipped_exon_event_equal(&se, testevent))
      {
        newevent = false;
        break;
      }
    }
    if(newevent)
    {
      GtStr *seqid = gt_genome_node_get_seqid(skipped);
      const char *geneid = gt_feature_node_get_attribute(gene, "ID");
      if(geneid == NULL)
      {
        geneid = gt_feature_node_get_attribute(gene, "gene_id");
        if(geneid == NULL)
          geneid = "";
      }
      if(v->report)
      {
        fprintf(v->report, "skipped exon: %s %s %lu-%lu %lu-%lu %lu-%lu\n",
                gt_str_get(seqid), geneid, se.left.start, se.left.end,
                se.skipped.start, se.skipped.end,
                se.right.start, se.right.end);
      }
      gt_array_add(v->se_events, se);
      numevents++;
    }
  }
  gt_array_delete(overlap);
  return numevents;
}

static GtUword check_skipped_exons(AgnASInspectCEVisitor *v,
                                   GtFeatureNode *gene, GtArray *mrnas)
{
  GtUword i,j,k, numevents = 0;
  agn_assert(v && gene && mrnas);
  agn_assert(gt_array_size(mrnas) > 1);
  for(i = 0; i < gt_array_size(mrnas); i++)
  {
    GtFeatureNode **mrna1 = gt_array_get(mrnas, i);
    GtArray *m1_exons = agn_mrna_exons(*mrna1);
    if(agn_typecheck_count(*mrna1, agn_typecheck_exon) < 2)
      continue;
    for(j = 0; j < gt_array_size(mrnas); j++)
    {
      GtFeatureNode **mrna2 = gt_array_get(mrnas, j);
      if(i == j || agn_typecheck_count(*mrna2, agn_typecheck_exon) < 3)
        continue;

      for(k = 1; k < gt_array_size(m1_exons); k++)
      {
        GtFeatureNode *lexon  = *(GtFeatureNode **)gt_array_get(m1_exons, k-1);
        GtFeatureNode *rexon  = *(GtFeatureNode **)gt_array_get(m1_exons, k);
        numevents += check_exon_pair(v, gene, lexon, rexon, *mrna2);
      }
    }
  }
  return numevents;
}

static bool skipped_exon_event_equal(SkippedExonEvent *a, SkippedExonEvent *b)
{
  agn_assert(a && b);
  GtStr *seqid_a = gt_genome_node_get_seqid((GtGenomeNode *)a->gene);
  GtStr *seqid_b = gt_genome_node_get_seqid((GtGenomeNode *)b->gene);
  if(gt_str_cmp(seqid_a, seqid_b) != 0)
    return false;
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
  AgnASInspectCEVisitor *v = as_inspect_ce_visitor_cast(nv);
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

    GtUword numevents = check_skipped_exons(v, gene, mrnas);
    if(numevents > 0)
    {
      gt_genome_node_ref((GtGenomeNode *)gene);
      gt_array_add(v->as_genes, gene);
    }
    gt_array_delete(mrnas);
  }

  gt_array_delete(genes);
  return 0;
}
