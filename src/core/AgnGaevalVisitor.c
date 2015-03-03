/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "core/array_api.h"
#include "core/queue_api.h"
#include "extended/feature_index_memory_api.h"
#include "AgnFilterStream.h"
#include "AgnGaevalVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define gaeval_visitor_cast(GV)\
        gt_node_visitor_cast(gaeval_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//

struct AgnGaevalVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureIndex *alignments;
  AgnGaevalParams params;
};


//----------------------------------------------------------------------------//
// Prototypes of private functions
//----------------------------------------------------------------------------//

/**
 * @function FIXME
 */
static double gaeval_visitor_calculate_coverage(AgnGaevalVisitor *v,
                                                GtFeatureNode *genemodel,
                                                GtError *error);

/**
 * @function Cast a node visitor object as a AgnGaevalVisitor.
 */
static const GtNodeVisitorClass* gaeval_visitor_class();

/**
 * @function FIXME
 */
static double gaeval_visitor_coverage_resolve(GtFeatureNode *genemodel,
                                              GtArray *exon_coverage);

/**
 * @function Destructor.
 */
static void gaeval_visitor_free(GtNodeVisitor *nv);

/**
 * @function FIXME
 */
static GtArray*
gaeval_visitor_intersect(GtGenomeNode *genemodel, GtGenomeNode *alignment);

/**
 * @function FIXME
 */
static GtRange gaeval_visitor_range_intersect(GtRange *r1, GtRange *r2);

/**
 * @function Generate data for unit testing.
 */
static void gaeval_visitor_test_data(GtQueue *queue);

/**
 * @function FIXME
 */
static void gaeval_visitor_union(GtArray *cov1, GtArray *cov2);

/**
 * @function Procedure for processing feature nodes (the only node of interest
 * for this node visitor).
 */
static int
gaeval_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                  GtError *error);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

GtNodeStream *agn_gaeval_stream_new(GtNodeStream *in, GtNodeStream *astream,
                                    AgnGaevalParams gparams)
{
  GtNodeVisitor *nv = agn_gaeval_visitor_new(astream, gparams);
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor*
agn_gaeval_visitor_new(GtNodeStream *astream, AgnGaevalParams gparams)
{
  agn_assert(astream);

  // Create the node visitor
  GtNodeVisitor *nv = gt_node_visitor_create(gaeval_visitor_class());
  AgnGaevalVisitor *v = gaeval_visitor_cast(nv);
  v->alignments = gt_feature_index_memory_new();
  v->params = gparams;

  // Set up node stream to load alignment features into memory
  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(typestokeep, "cDNA_match", "cDNA_match");
  gt_hashmap_add(typestokeep, "EST_match", "EST_match");
  gt_hashmap_add(typestokeep, "nucleotide_match", "nucleotide_match");
  GtNodeStream *fstream = agn_filter_stream_new(astream, typestokeep);
  GtNodeStream *align_stream = gt_feature_out_stream_new(fstream,v->alignments);

  // Process the node stream
  GtError *error = gt_error_new();
  int result = gt_node_stream_pull(align_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[AEGeAn::AgnGaevalStream] error parsing alignments: %s\n",
            gt_error_get(error));
    gt_node_visitor_delete(nv);
    return NULL;
  }
  gt_error_delete(error);
  gt_hashmap_delete(typestokeep);
  gt_node_stream_delete(fstream);
  gt_node_stream_delete(align_stream);

  return nv;
}

bool agn_gaeval_visitor_unit_test(AgnUnitTest *test)
{
  gaeval_visitor_test_data(NULL);
  return 0;
}

static const GtNodeVisitorClass* gaeval_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnGaevalVisitor),
                                    gaeval_visitor_free, NULL,
                                    gaeval_visitor_visit_feature_node,
                                    NULL, NULL, NULL);
  }
  return nvc;
}

static double gaeval_visitor_calculate_coverage(AgnGaevalVisitor *v,
                                                GtFeatureNode *genemodel,
                                                GtError *error)
{
  agn_assert(v && genemodel);

  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)genemodel);
  GtRange mrna_range = gt_genome_node_get_range((GtGenomeNode *)genemodel);
  GtArray *overlapping = gt_array_new( sizeof(GtFeatureNode *) );
  gt_feature_index_get_features_for_range(v->alignments, overlapping,
                                          gt_str_get(seqid), &mrna_range,
                                          error);

  GtArray *exon_coverage = gt_array_new( sizeof(GtRange) );
  GtUword i;
  for(i = 0; i < gt_array_size(overlapping); i++)
  {
    GtFeatureNode *alignment = *(GtFeatureNode **)gt_array_get(overlapping,i);
    GtArray *covered_parts = gaeval_visitor_intersect((GtGenomeNode*)genemodel,
                                                      (GtGenomeNode*)alignment);
    gaeval_visitor_union(exon_coverage, covered_parts);
    gt_array_delete(covered_parts);
  }
  double coverage = gaeval_visitor_coverage_resolve(genemodel, exon_coverage);
  gt_array_delete(exon_coverage);

  return coverage;
}

static double gaeval_visitor_coverage_resolve(GtFeatureNode *genemodel,
                                              GtArray *exon_coverage)
{
  agn_assert(genemodel && exon_coverage);
  agn_assert(gt_feature_node_has_type(genemodel, "mRNA"));

  GtUword cum_exon_length = 0;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(genemodel);
  GtFeatureNode *tempfeat;
  for(tempfeat  = gt_feature_node_iterator_next(iter);
      tempfeat != NULL;
      tempfeat  = gt_feature_node_iterator_next(iter))
  {
    if(gt_feature_node_has_type(tempfeat, "exon"))
      cum_exon_length += gt_genome_node_get_length((GtGenomeNode *)tempfeat);
  }
  gt_feature_node_iterator_delete(iter);

  GtUword i, covered = 0;
  for(i = 0; i < gt_array_size(exon_coverage); i++)
  {
    GtRange *range = gt_array_get(exon_coverage, i);
    covered += gt_range_length(range);
  }
  
  agn_assert(covered <= cum_exon_length);
  return (double)covered / (double)cum_exon_length;
}

static void gaeval_visitor_free(GtNodeVisitor *nv)
{
  AgnGaevalVisitor *v = gaeval_visitor_cast(nv);
  gt_feature_index_delete(v->alignments);
}

static GtArray*
gaeval_visitor_intersect(GtGenomeNode *genemodel, GtGenomeNode *alignment)
{
  agn_assert(genemodel && alignment);

  GtFeatureNode *genefn = gt_feature_node_cast(genemodel);
  GtFeatureNode *algnfn = gt_feature_node_cast(alignment);
  agn_assert(gt_feature_node_has_type(genefn, "mRNA"));

  GtArray *covered_parts = gt_array_new( sizeof(GtRange) );
  GtFeatureNodeIterator *gniter = gt_feature_node_iterator_new(genefn);
  GtFeatureNode *tempfeat;
  for(tempfeat  = gt_feature_node_iterator_next(gniter);
      tempfeat != NULL;
      tempfeat  = gt_feature_node_iterator_next(gniter))
  {
    if(gt_feature_node_has_type(tempfeat, "exon") == false)
      continue;

    GtRange featrange = gt_genome_node_get_range((GtGenomeNode *)tempfeat);
    GtFeatureNodeIterator *aniter = gt_feature_node_iterator_new(algnfn);
    GtFeatureNode *tempaln;
    GtRange nullrange = {0, 0};
    for(tempaln  = gt_feature_node_iterator_next(aniter);
        tempaln != NULL;
        tempaln  = gt_feature_node_iterator_next(aniter))
    {
      GtRange alnrange = gt_genome_node_get_range((GtGenomeNode *) tempaln);
      GtRange intr = gaeval_visitor_range_intersect(&featrange, &alnrange);
      if(gt_range_compare(&intr, &nullrange) != 0)
        gt_array_add(covered_parts, intr);
    }
    gt_feature_node_iterator_delete(aniter);
  }
  gt_feature_node_iterator_delete(gniter);

  return covered_parts;
}

static GtRange gaeval_visitor_range_intersect(GtRange *r1, GtRange *r2)
{
  agn_assert(r1 && r2);
  if(gt_range_overlap(r1, r2))
  {
    GtRange inter = *r1;
    if(r2->start < inter.start) inter.start = r2->start;
    if(r2->end   > inter.end)   inter.end   = r2->end;
    return inter;
  }
  
  GtRange nullrange = {0, 0};
  return nullrange;
}

static void gaeval_visitor_test_data(GtQueue *queue)
{
  //
}

static void gaeval_visitor_union(GtArray *cov1, GtArray *cov2)
{
  agn_assert(cov1 && cov2);

  GtUword i, j;
  for(i = 0; i < gt_array_size(cov2); i++)
  {
    GtRange *range2 = gt_array_get(cov2, i);
    bool overlaps = false;
    for(j = 0; j < gt_array_size(cov1); j++)
    {
      GtRange *range1 = gt_array_get(cov1, j);
      if(gt_range_overlap(range1, range2))
      {
        *range1 = gt_range_join(range1, range2);
        overlaps = true;
        break;
      }
    }
    if(overlaps == false)
      gt_array_add(cov1, *range2);
  }
}

static int
gaeval_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                  GtError *error)
{
  AgnGaevalVisitor *v = gaeval_visitor_cast(nv);
  gt_error_check(error);

  GtFeatureNodeIterator *feats = gt_feature_node_iterator_new(fn);
  GtFeatureNode *tempfeat;
  for(tempfeat  = gt_feature_node_iterator_next(feats);
      tempfeat != NULL;
      tempfeat  = gt_feature_node_iterator_next(feats))
  {
    if(agn_typecheck_mrna(tempfeat) == false)
      continue;

    double coverage = gaeval_visitor_calculate_coverage(v, tempfeat, error);
    char covstr[16];
    sprintf(covstr, "%.3lf", coverage);
    gt_feature_node_add_attribute(tempfeat, "gaeval_coverage", covstr);
  }
  gt_feature_node_iterator_delete(feats);

  return 0;
}
