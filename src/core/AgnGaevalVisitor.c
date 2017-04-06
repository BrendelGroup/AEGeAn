/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <math.h>
#include <string.h>
#include "core/array_api.h"
#include "core/queue_api.h"
#include "extended/feature_index_memory_api.h"
#include "AgnFilterStream.h"
#include "AgnGaevalVisitor.h"
#include "AgnInferCDSVisitor.h"
#include "AgnInferExonsVisitor.h"
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
  FILE *tsvout;
  AgnGaevalParams params;
};


//----------------------------------------------------------------------------//
// Prototypes of private functions
//----------------------------------------------------------------------------//

/**
 * @function Calculate coverage for the given gene model.
 */
static double gaeval_visitor_calculate_coverage(AgnGaevalVisitor *v,
                                                GtFeatureNode *genemodel,
                                                GtError *error);

/**
 * @function Calculate integrity for the given gene model.
 */
static double gaeval_visitor_calculate_integrity(AgnGaevalVisitor *v,
                                                 GtFeatureNode *genemodel,
                                                 double coverage,
                                                 GtError *error);

/**
 * @function Cast a node visitor object as a AgnGaevalVisitor.
 */
static const GtNodeVisitorClass* gaeval_visitor_class();

/**
 * @function Add up exon and match lengths to calculate coverage.
 */
static double gaeval_visitor_coverage_resolve(GtFeatureNode *genemodel,
                                              GtArray *exon_coverage);

/**
 * @function Destructor.
 */
static void gaeval_visitor_free(GtNodeVisitor *nv);

/**
 * @function Determine the ranges of overlap, if any, between the gene model and
 * the alignment. Returns NULL if there is no overlap.
 */
static GtArray*
gaeval_visitor_intersect(GtGenomeNode *genemodel, GtGenomeNode *alignment);

/**
 * @function Calculate the proportion of introns confirmed by gaps in
 * overlapping alignments.
 */
static double gaeval_visitor_introns_confirmed(GtArray *introns, GtArray *gaps);

/**
 * @function Determine the overlap, if any, between the two ranges. Returns the
 * null range {0,0} in case of no overlap.
 */
static GtRange gaeval_visitor_range_intersect(GtRange *r1, GtRange *r2);

/**
 * @function Typecheck select function for grabbing `match_gap` features.
 */
static bool gaeval_visitor_typecheck_gap(GtFeatureNode *fn);

/**
 * @function Used to bombine the coverage from individual alignments into a
 * single aggregate coverage.
 */
static GtArray *gaeval_visitor_union(GtArray *cov1, GtArray *cov2);

/**
 * @function Procedure for processing feature nodes (the only node of interest
 * for this node visitor).
 */
static int
gaeval_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                  GtError *error);

/**
 * @function Unit test for coverage calculations.
 */
static void gv_test_calc_coverage(AgnUnitTest *test);

/**
 * @function Unit test for integrity calculations.
 */
static void gv_test_calc_integrity(AgnUnitTest *test);

/**
 * @function Unit test for `gaeval_visitor_intersect` function.
 */
static void gv_test_intersect(AgnUnitTest *test);

/**
 * @function Unit test for `gaeval_visitor_introns_confirmed` function.
 */
static void gv_test_introns_confirmed(AgnUnitTest *test);

/**
 * @function Unit test for `gaeval_visitor_range_intersect` function.
 */
static void gv_test_range_intersect(AgnUnitTest *test);

/**
 * @function Unit test for `gaeval_visitor_union` function.
 */
static void gv_test_union(AgnUnitTest *test);


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
  v->tsvout = NULL;
  v->params = gparams;

  // Check that sum of weights is 1.0
  double weights_total = gparams.alpha + gparams.beta +
                         gparams.gamma + gparams.epsilon;
  if(fabs(weights_total - 1.0) > 0.0001)
  {
    fprintf(stderr, "[AgnGaevalVisitor::agn_gaeval_visitor_new] warning: "
            "sum of weights is not 1.0 %.3lf; integrity calculations will be "
            "incorrect\n", weights_total);
  }


  // Set up node stream to load alignment features into memory
  GtQueue *streams = gt_queue_new();
  GtNodeStream *stream, *last_stream;
  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(typestokeep, "cDNA_match", "cDNA_match");
  gt_hashmap_add(typestokeep, "EST_match", "EST_match");
  gt_hashmap_add(typestokeep, "nucleotide_match", "nucleotide_match");
  stream = agn_filter_stream_new(astream, typestokeep);
  gt_queue_add(streams, stream);
  last_stream = stream;

  stream = gt_feature_out_stream_new(last_stream, v->alignments);
  gt_queue_add(streams, stream);
  last_stream = stream;

  stream = gt_inter_feature_stream_new(last_stream, "cDNA_match", "match_gap");
  gt_queue_add(streams, stream);
  last_stream = stream;

  stream = gt_inter_feature_stream_new(last_stream, "EST_match", "match_gap");
  gt_queue_add(streams, stream);
  last_stream = stream;

  stream = gt_inter_feature_stream_new(last_stream, "nucleotide_match",
                                       "match_gap");
  gt_queue_add(streams, stream);
  last_stream = stream;

  // Process the node stream
  GtError *error = gt_error_new();
  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[AEGeAn::AgnGaevalStream] error parsing alignments: %s\n",
            gt_error_get(error));
    gt_node_visitor_delete(nv);
    return NULL;
  }
  gt_error_delete(error);
  gt_hashmap_delete(typestokeep);
  while(gt_queue_size(streams) > 0)
  {
    stream = gt_queue_get(streams);
    gt_node_stream_delete(stream);
  }
  gt_queue_delete(streams);

  return nv;
}

void agn_gaeval_visitor_tsv_out(AgnGaevalVisitor *v, GtStr *tsvfilename)
{
  v->tsvout = fopen(gt_str_get(tsvfilename), "w");
  if(v->tsvout == NULL)
  {
    fprintf(stderr, "error: unable to open output file %s",
            gt_str_get(tsvfilename));
    exit(1);
  }
  fprintf(v->tsvout, "ID\tLabel\tCoverage\tIntegrity\n");
}

bool agn_gaeval_visitor_unit_test(AgnUnitTest *test)
{
  gv_test_range_intersect(test);
  gv_test_union(test);
  gv_test_intersect(test);
  gv_test_calc_coverage(test);
  gv_test_introns_confirmed(test);
  gv_test_calc_integrity(test);
  return agn_unit_test_success(test);
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
  bool hasseqid;
  gt_feature_index_has_seqid(v->alignments, &hasseqid, gt_str_get(seqid),error);
  if(hasseqid)
  {
    gt_feature_index_get_features_for_range(v->alignments, overlapping,
                                            gt_str_get(seqid), &mrna_range,
                                            error);
  }

  GtArray *exon_coverage = gt_array_new( sizeof(GtRange) );
  GtUword i;
  for(i = 0; i < gt_array_size(overlapping); i++)
  {
    GtFeatureNode *alignment = *(GtFeatureNode **)gt_array_get(overlapping, i);
    GtArray *covered_parts = gaeval_visitor_intersect((GtGenomeNode*)genemodel,
                                                      (GtGenomeNode*)alignment);
    if(covered_parts != NULL)
    {
      GtArray *temp = gaeval_visitor_union(exon_coverage, covered_parts);
      gt_array_delete(covered_parts);
      gt_array_delete(exon_coverage);
      exon_coverage = temp;
    }
  }
  double coverage = gaeval_visitor_coverage_resolve(genemodel, exon_coverage);
  gt_array_delete(exon_coverage);
  gt_array_delete(overlapping);

  return coverage;
}

static double gaeval_visitor_calculate_integrity(AgnGaevalVisitor *v,
                                                 GtFeatureNode *genemodel,
                                                 double coverage,
                                                 GtError *error)
{
  agn_assert(v && genemodel);

  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)genemodel);
  GtRange mrna_range = gt_genome_node_get_range((GtGenomeNode *)genemodel);
  GtArray *overlapping = gt_array_new( sizeof(GtFeatureNode *) );
  bool hasseqid;
  gt_feature_index_has_seqid(v->alignments, &hasseqid, gt_str_get(seqid),error);
  if(hasseqid)
  {
    gt_feature_index_get_features_for_range(v->alignments, overlapping,
                                            gt_str_get(seqid), &mrna_range,
                                            error);
  }

  GtArray *gaps = gt_array_new( sizeof(GtFeatureNode *) );
  while(gt_array_size(overlapping) > 0)
  {
    GtFeatureNode *alignment = *(GtFeatureNode **)gt_array_pop(overlapping);
    GtArray *agaps = agn_typecheck_select(alignment,
                                          gaeval_visitor_typecheck_gap);
    gt_array_add_array(gaps, agaps);
    gt_array_delete(agaps);
  }
  gt_array_delete(overlapping);

  GtUword utr5p_len = agn_mrna_5putr_length(genemodel);
  double utr5p_score = 0.0;
  if(utr5p_len >= v->params.exp_5putr_len)
    utr5p_score = 1.0;
  else
    utr5p_score = (double)utr5p_len / (double)v->params.exp_5putr_len;

  GtUword utr3p_len = agn_mrna_3putr_length(genemodel);
  double utr3p_score = 0.0;
  if(utr3p_len >= v->params.exp_3putr_len)
    utr3p_score = 1.0;
  else
    utr3p_score = (double)utr3p_len / (double)v->params.exp_3putr_len;

  GtArray *introns = agn_typecheck_select(genemodel, agn_typecheck_intron);
  GtUword exoncount = agn_typecheck_count(genemodel, agn_typecheck_exon);
  agn_assert(gt_array_size(introns) == exoncount - 1);
  double structure_score = 0.0;
  if(gt_array_size(introns) == 0)
  {
    GtUword cdslen = agn_mrna_cds_length(genemodel);
    if(cdslen >= v->params.exp_cds_len)
      structure_score = 1.0;
    else
      structure_score = (double)cdslen / (double)v->params.exp_cds_len;
  }
  else
  {
    structure_score = gaeval_visitor_introns_confirmed(introns, gaps);
  }
  gt_array_delete(gaps);
  gt_array_delete(introns);

  double integrity = (v->params.alpha   * structure_score) +
                     (v->params.beta    * coverage)        +
                     (v->params.gamma   * utr5p_score)     +
                     (v->params.epsilon * utr3p_score);

  return integrity;
}

static double gaeval_visitor_coverage_resolve(GtFeatureNode *genemodel,
                                              GtArray *exon_coverage)
{
  agn_assert(genemodel && exon_coverage);
  agn_assert(gt_feature_node_has_type(genemodel, "mRNA"));

  GtUword cum_exon_length =
      agn_typecheck_feature_combined_length(genemodel, agn_typecheck_exon);

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
  GtStrand genestrand = gt_feature_node_get_strand(genefn);
  GtStrand algnstrand = gt_feature_node_get_strand(algnfn);
  if(genestrand != algnstrand)
    return NULL;

  GtArray *covered_parts = gt_array_new( sizeof(GtRange) );
  GtArray *exons = agn_typecheck_select(genefn, agn_typecheck_exon);
  GtWord i;
  for(i = 0; i < gt_array_size(exons); i++)
  {
    GtGenomeNode *exon = *(GtGenomeNode **)gt_array_get(exons, i);
    GtRange exonrange = gt_genome_node_get_range(exon);

    GtFeatureNodeIterator *aniter = gt_feature_node_iterator_new(algnfn);
    GtFeatureNode *tempaln;
    GtRange nullrange = {0, 0};
    for(tempaln  = gt_feature_node_iterator_next(aniter);
        tempaln != NULL;
        tempaln  = gt_feature_node_iterator_next(aniter))
    {
      if(gt_feature_node_has_type(tempaln, "match_gap"))
        continue;

      GtRange alnrange = gt_genome_node_get_range((GtGenomeNode *) tempaln);
      GtRange intr = gaeval_visitor_range_intersect(&exonrange, &alnrange);
      if(gt_range_compare(&intr, &nullrange) != 0)
        gt_array_add(covered_parts, intr);
    }
    gt_feature_node_iterator_delete(aniter);
  }
  gt_array_delete(exons);

  for(i = 0; i < gt_array_size(covered_parts); i++)
  {
    GtRange *r1 = gt_array_get(covered_parts, i);
    GtUword j;
    for(j = i+1; j < gt_array_size(covered_parts); j++)
    {
      GtRange *r2 = gt_array_get(covered_parts, j);
      agn_assert(gt_range_overlap(r1, r2) == false);
    }
  }

  return covered_parts;
}

static double gaeval_visitor_introns_confirmed(GtArray *introns, GtArray *gaps)
{
  agn_assert(introns && gaps);
  GtUword intron_count = gt_array_size(introns);
  GtUword gap_count = gt_array_size(gaps);
  agn_assert(intron_count > 0);

  if(gap_count == 0)
    return 0.0;

  GtUword i, j, num_confirmed = 0;
  for(i = 0; i < intron_count; i++)
  {
    GtGenomeNode *intron = *(GtGenomeNode **)gt_array_get(introns, i);
    GtRange intron_range = gt_genome_node_get_range(intron);
    for(j = 0; j < gap_count; j++)
    {
      GtGenomeNode *gap = *(GtGenomeNode **)gt_array_get(gaps, j);
      GtRange gap_range = gt_genome_node_get_range(gap);
      if(gt_range_compare(&intron_range, &gap_range) == 0)
      {
        num_confirmed++;
        break;
      }
    }
  }

  return (double)num_confirmed / (double)intron_count;
}

static GtRange gaeval_visitor_range_intersect(GtRange *r1, GtRange *r2)
{
  agn_assert(r1 && r2);
  if(gt_range_overlap(r1, r2))
  {
    GtRange inter = *r1;
    if(r2->start > inter.start) inter.start = r2->start;
    if(r2->end   < inter.end)   inter.end   = r2->end;
    return inter;
  }
  GtRange nullrange = {0, 0};
  return nullrange;
}

static bool gaeval_visitor_typecheck_gap(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "match_gap");
}

static GtArray *gaeval_visitor_union(GtArray *cov1, GtArray *cov2)
{
  agn_assert(cov1 && cov2);
  gt_array_add_array(cov1, cov2);
  if(gt_array_size(cov1) > 1)
    gt_array_sort(cov1, (GtCompare)gt_range_compare);

  GtArray *runion = gt_array_new(sizeof(GtRange));
  if(gt_array_size(cov1) == 0)
    return runion;
  GtRange *rng = gt_array_get(cov1, 0);
  gt_array_add(runion, *rng);
  GtRange *prev = gt_array_get(runion, 0);
  if(gt_array_size(cov1) == 1)
    return runion;

  GtUword i;
  for(i = 1; i < gt_array_size(cov1); i++)
  {
    rng = gt_array_get(cov1, i);
    if(gt_range_overlap(rng, prev))
      *prev = gt_range_join(rng, prev);
    else
    {
      gt_array_add(runion, *rng);
      prev = gt_array_get(runion, gt_array_size(runion) - 1);
    }
  }

  return runion;
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

    double integrity =
        gaeval_visitor_calculate_integrity(v, tempfeat, coverage, error);
    char intstr[16];
    sprintf(intstr, "%.3lf", integrity);
    gt_feature_node_add_attribute(tempfeat, "gaeval_integrity", intstr);

    if(v->tsvout)
    {
      const char *mrnaid = gt_feature_node_get_attribute(tempfeat, "ID");
      const char *mrnalabel = agn_feature_node_get_label(tempfeat);
      fprintf(v->tsvout, "%s\t%s\t%s\t%s\n", mrnaid, mrnalabel, covstr, intstr);
    }
  }
  gt_feature_node_iterator_delete(feats);

  return 0;
}

static void gv_test_calc_coverage(AgnUnitTest *test)
{
  const char *filename = "data/gff3/gaeval-stream-unit-test-1.gff3";
  GtNodeStream *align_in = gt_gff3_in_stream_new_unsorted(1, &filename);
  AgnGaevalParams params = { 0.6, 0.3, 0.05, 0.05, 400, 200, 100 };
  GtNodeVisitor *nv = agn_gaeval_visitor_new(align_in, params);
  AgnGaevalVisitor *gv = gaeval_visitor_cast(nv);
  gt_node_stream_delete(align_in);

  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &filename);
  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(typestokeep, "mRNA", "mRNA");
  GtNodeStream *filtstream = agn_filter_stream_new(gff3in, typestokeep);

  GtError *error = gt_error_new();
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *featstream = gt_array_out_stream_new(filtstream, feats, error);
  int result = gt_node_stream_pull(featstream, error);
  if(result == -1)
  {
    fprintf(stderr, "[AgnGaevalVisitor::gv_test_calc_coverage] error "
            "processing GFF3: %s\n", gt_error_get(error));
    return;
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(filtstream);
  gt_node_stream_delete(featstream);
  gt_hashmap_delete(typestokeep);

  agn_assert(gt_array_size(feats) == 3);
  GtFeatureNode *g1 = *(GtFeatureNode **)gt_array_get(feats, 0);
  GtFeatureNode *g2 = *(GtFeatureNode **)gt_array_get(feats, 1);
  GtFeatureNode *g3 = *(GtFeatureNode **)gt_array_get(feats, 2);

  double cov1 = gaeval_visitor_calculate_coverage(gv, g1, error);
  double cov2 = gaeval_visitor_calculate_coverage(gv, g2, error);
  double cov3 = gaeval_visitor_calculate_coverage(gv, g3, error);
  bool test1 = fabs(cov1 - 0.252) < 0.001 &&
               fabs(cov2 - 0.473) < 0.001 &&
               fabs(cov3 - 1.000) < 0.001;

  agn_unit_test_result(test, "calculate coverage", test1);

  gt_error_delete(error);
  gt_array_delete(feats);
  gt_genome_node_delete((GtGenomeNode *)g1);
  gt_genome_node_delete((GtGenomeNode *)g2);
  gt_genome_node_delete((GtGenomeNode *)g3);
  gt_node_visitor_delete(nv);
}

static void gv_test_calc_integrity(AgnUnitTest *test)
{
  const char *filename = "data/gff3/gaeval-stream-unit-test-2.gff3";
  GtNodeStream *align_in = gt_gff3_in_stream_new_unsorted(1, &filename);
  AgnGaevalParams params = { 0.6, 0.3, 0.05, 0.05, 400, 200, 100 };
  GtNodeVisitor *nv = agn_gaeval_visitor_new(align_in, params);
  AgnGaevalVisitor *gv = gaeval_visitor_cast(nv);
  gt_node_stream_delete(align_in);

  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &filename);
  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(typestokeep, "mRNA", "mRNA");
  GtNodeStream *filtstream = agn_filter_stream_new(gff3in, typestokeep);
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtNodeStream *ics = agn_infer_cds_stream_new(filtstream, NULL, logger);
  GtNodeStream *ies = agn_infer_exons_stream_new(ics, NULL, logger);

  GtError *error = gt_error_new();
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *featstream = gt_array_out_stream_new(ies, feats, error);
  int result = gt_node_stream_pull(featstream, error);
  if(result == -1)
  {
    fprintf(stderr, "[AgnGaevalVisitor::gv_test_calc_integrity] error "
            "processing GFF3: %s\n", gt_error_get(error));
    return;
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(filtstream);
  gt_node_stream_delete(featstream);
  gt_node_stream_delete(ics);
  gt_node_stream_delete(ies);
  gt_logger_delete(logger);
  gt_hashmap_delete(typestokeep);

  agn_assert(gt_array_size(feats) == 2);
  GtFeatureNode *g1 = *(GtFeatureNode **)gt_array_get(feats, 0);
  GtFeatureNode *g2 = *(GtFeatureNode **)gt_array_get(feats, 1);

  double cov1 = gaeval_visitor_calculate_coverage(gv,  g1, error);
  double cov2 = gaeval_visitor_calculate_coverage(gv,  g2, error);
  double int1 = gaeval_visitor_calculate_integrity(gv, g1, cov1, error);
  double int2 = gaeval_visitor_calculate_integrity(gv, g2, cov2, error);

  bool test1 = fabs(cov1 - 1.000) < 0.001 &&
               fabs(cov2 - 0.997) < 0.001 &&
               fabs(int1 - 0.850) < 0.001 &&
               fabs(int2 - 0.863) < 0.001;
  agn_unit_test_result(test, "calculate integrity", test1);

  gt_error_delete(error);
  gt_array_delete(feats);
  gt_genome_node_delete((GtGenomeNode *)g1);
  gt_genome_node_delete((GtGenomeNode *)g2);
  gt_node_visitor_delete(nv);
}

static void gv_test_intersect(AgnUnitTest *test)
{
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtError *error = gt_error_new();
  const char *filename = "data/gff3/gaeval-stream-unit-test-1.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &filename);
  GtNodeStream *fstream = gt_array_out_stream_new(gff3in, feats, error);
  int result = gt_node_stream_pull(fstream, error);
  if(result == -1)
  {
    fprintf(stderr, "[AgnGaevalVisitor::gv_test_intersect] error "
            "processing GFF3: %s\n", gt_error_get(error));
    return;
  }
  gt_error_delete(error);
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(fstream);

  agn_assert(gt_array_size(feats) == 9);
  GtGenomeNode *g1 = *(GtGenomeNode **)gt_array_get(feats, 1);
  GtGenomeNode *g2 = *(GtGenomeNode **)gt_array_get(feats, 3);
  GtGenomeNode *g3 = *(GtGenomeNode **)gt_array_get(feats, 7);

  GtGenomeNode *est1 = *(GtGenomeNode **)gt_array_get(feats, 0);
  GtGenomeNode *est2 = *(GtGenomeNode **)gt_array_get(feats, 2);
  GtGenomeNode *est3 = *(GtGenomeNode **)gt_array_get(feats, 4);
  GtGenomeNode *est4 = *(GtGenomeNode **)gt_array_get(feats, 5);
  GtGenomeNode *est5 = *(GtGenomeNode **)gt_array_get(feats, 6);
  GtGenomeNode *est6 = *(GtGenomeNode **)gt_array_get(feats, 8);

  GtArray *cov = gaeval_visitor_intersect(g1, est1);
  bool test1 = cov == NULL;
  cov = gaeval_visitor_intersect(g1, est2);
  test1 = gt_array_size(cov) == 1;
  if(test1)
  {
    GtRange *range01 = gt_array_pop(cov);
    GtRange testrange = { 400, 500 };
    test1 = gt_range_compare(range01, &testrange) == 0;
  }
  agn_unit_test_result(test, "intersect (1)", test1);
  gt_array_delete(cov);

  cov = gaeval_visitor_intersect(g2, est3);
  bool test2 = gt_array_size(cov) == 2;
  if(test2)
  {
    GtRange *range01 = gt_array_get(cov, 0);
    GtRange *range02 = gt_array_get(cov, 1);
    GtRange testrange1 = { 800, 900 };
    GtRange testrange2 = { 1050, 1075 };
    test2 = gt_range_compare(range01, &testrange1) == 0 &&
            gt_range_compare(range02, &testrange2) == 0;
  }
  agn_unit_test_result(test, "intersect (2)", test2);
  gt_array_delete(cov);

  cov = gaeval_visitor_intersect(g2, est4);
  bool test3 = gt_array_size(cov) == 2;
  if(test3)
  {
    GtRange *range01 = gt_array_get(cov, 0);
    GtRange *range02 = gt_array_get(cov, 1);
    GtRange testrange1 = { 1070, 1125 };
    GtRange testrange2 = { 1250, 1310 };
    test3 = gt_range_compare(range01, &testrange1) == 0 &&
            gt_range_compare(range02, &testrange2) == 0;
  }
  agn_unit_test_result(test, "intersect (3)", test3);
  gt_array_delete(cov);

  cov = gaeval_visitor_intersect(g3, est5);
  bool test4 = gt_array_size(cov) == 2;
  if(test4)
  {
    GtRange *range01 = gt_array_get(cov, 0);
    GtRange *range02 = gt_array_get(cov, 1);
    GtRange testrange1 = { 2000, 3000 };
    GtRange testrange2 = { 4000, 5000 };
    test4 = gt_range_compare(range01, &testrange1) == 0 &&
            gt_range_compare(range02, &testrange2) == 0;
  }
  agn_unit_test_result(test, "intersect (4)", test4);
  gt_array_delete(cov);

  cov = gaeval_visitor_intersect(g3, est6);
  bool test5 = gt_array_size(cov) == 2;
  if(test5)
  {
    GtRange *range01 = gt_array_get(cov, 0);
    GtRange *range02 = gt_array_get(cov, 1);
    GtRange testrange1 = { 2500, 3000 };
    GtRange testrange2 = { 4000, 5000 };
    test5 = gt_range_compare(range01, &testrange1) == 0 &&
            gt_range_compare(range02, &testrange2) == 0;
  }
  agn_unit_test_result(test, "intersect (5)", test5);
  gt_array_delete(cov);

  gt_array_delete(feats);
  gt_genome_node_delete(g1);
  gt_genome_node_delete(g2);
  gt_genome_node_delete(g3);
  gt_genome_node_delete(est1);
  gt_genome_node_delete(est2);
  gt_genome_node_delete(est3);
  gt_genome_node_delete(est4);
  gt_genome_node_delete(est5);
  gt_genome_node_delete(est6);
}

static void gv_test_introns_confirmed(AgnUnitTest *test)
{
  GtGenomeNode *intron, *gap;
  GtStr *seqid = gt_str_new_cstr("chr");
  GtArray *introns = gt_array_new( sizeof(GtGenomeNode *) );
  intron = gt_feature_node_new(seqid, "intron", 1000, 1170, GT_STRAND_REVERSE);
  gt_array_add(introns, intron);
  intron = gt_feature_node_new(seqid, "intron", 1225, 1305, GT_STRAND_REVERSE);
  gt_array_add(introns, intron);
  intron = gt_feature_node_new(seqid, "intron", 1950, 2110, GT_STRAND_REVERSE);
  gt_array_add(introns, intron);
  intron = gt_feature_node_new(seqid, "intron", 2545, 2655, GT_STRAND_REVERSE);
  gt_array_add(introns, intron);
  intron = gt_feature_node_new(seqid, "intron", 2800, 2950, GT_STRAND_REVERSE);
  gt_array_add(introns, intron);

  GtArray *gaps = gt_array_new( sizeof(GtGenomeNode *) );

  double intcon = gaeval_visitor_introns_confirmed(introns, gaps);
  bool test1 = fabs(intcon - 0.0) < 0.0001;
  agn_unit_test_result(test, "introns confirmed (no gaps)", test1);

  gap = gt_feature_node_new(seqid, "match_gap", 1000, 1170, GT_STRAND_REVERSE);
  gt_array_add(gaps, gap);
  gap = gt_feature_node_new(seqid, "match_gap", 1225, 1302, GT_STRAND_REVERSE);
  gt_array_add(gaps, gap);
  gap = gt_feature_node_new(seqid, "match_gap", 1950, 2110, GT_STRAND_REVERSE);
  gt_array_add(gaps, gap);
  gap = gt_feature_node_new(seqid, "match_gap", 2575, 2655, GT_STRAND_REVERSE);
  gt_array_add(gaps, gap);
  gap = gt_feature_node_new(seqid, "match_gap", 2800, 2950, GT_STRAND_REVERSE);
  gt_array_add(gaps, gap);

  intcon = gaeval_visitor_introns_confirmed(introns, gaps);
  bool test2 = fabs(intcon - 0.6) < 0.0001;
  agn_unit_test_result(test, "introns confirmed (gaps)", test2);

  while(gt_array_size(introns) > 0)
  {
    intron = *(GtGenomeNode **)gt_array_pop(introns);
    gt_genome_node_delete(intron);
  }
  gt_array_delete(introns);

  while(gt_array_size(gaps) > 0)
  {
    gap = *(GtGenomeNode **)gt_array_pop(gaps);
    gt_genome_node_delete(gap);
  }
  gt_array_delete(gaps);
  gt_str_delete(seqid);
}

static void gv_test_range_intersect(AgnUnitTest *test)
{
  GtRange nullrange = {0,0};

  GtRange r1 = { 500, 750 };
  GtRange r2 = { 645, 900 };
  GtRange inter1 = { 645, 750 };
  GtRange testrange1 = gaeval_visitor_range_intersect(&r1, &r2);
  agn_unit_test_result(test, "range intersect (1)",
                       gt_range_compare(&inter1, &testrange1) == 0);

  GtRange r3 = { 500, 750 };
  GtRange r4 = { 750, 900 };
  GtRange inter2 = { 750, 750 };
  GtRange testrange2 = gaeval_visitor_range_intersect(&r3, &r4);
  agn_unit_test_result(test, "range intersect (2)",
                       gt_range_compare(&inter2, &testrange2) == 0);

  GtRange r5 = { 500, 750 };
  GtRange r6 = { 751, 900 };
  GtRange testrange3 = gaeval_visitor_range_intersect(&r5, &r6);
  agn_unit_test_result(test, "range intersect (3)",
                       gt_range_compare(&nullrange, &testrange3) == 0);
}

static void gv_test_union(AgnUnitTest *test)
{
  GtArray *r1 = gt_array_new( sizeof(GtRange) );
  GtArray *r2 = gt_array_new( sizeof(GtRange) );
  GtRange rng01 = {1050,  9005};
  GtRange rng02 = {11525, 14070};
  gt_array_add(r2, rng01);
  gt_array_add(r2, rng02);
  GtArray *ru = gaeval_visitor_union(r1, r2);
  bool test1 = gt_array_size(ru) == 2;
  if(test1)
  {
    GtRange *temp1 = gt_array_get(ru, 0);
    GtRange *temp2 = gt_array_get(ru, 1);
    test1 = gt_range_compare(temp1, &rng01) == 0 &&
            gt_range_compare(temp2, &rng02) == 0;
  }
  agn_unit_test_result(test, "union (1)", test1);
  gt_array_delete(r1);
  gt_array_delete(r2);
  gt_array_delete(ru);

  r1 = gt_array_new( sizeof(GtRange) );
  r2 = gt_array_new( sizeof(GtRange) );
  GtRange rng03 = { 300, 500 };
  GtRange rng04 = { 700, 800 };
  GtRange rng05 = { 200, 400 };
  GtRange rng06 = { 700, 900 };
  gt_array_add(r1, rng03);
  gt_array_add(r1, rng04);
  gt_array_add(r2, rng05);
  gt_array_add(r2, rng06);
  ru = gaeval_visitor_union(r1, r2);
  bool test2 = gt_array_size(ru) == 2;
  if(test2)
  {
    GtRange *temp1 = gt_array_get(ru, 0);
    GtRange *temp2 = gt_array_get(ru, 1);
    GtRange testr1 = { 200, 500 };
    GtRange testr2 = { 700, 900 };
    test2 = gt_range_compare(temp1, &testr1) == 0 &&
            gt_range_compare(temp2, &testr2) == 0;
  }
  agn_unit_test_result(test, "union (2)", test2);
  gt_array_delete(r1);
  gt_array_delete(r2);
  gt_array_delete(ru);

  r1 = gt_array_new( sizeof(GtRange) );
  r2 = gt_array_new( sizeof(GtRange) );
  GtRange rng07 = { 300, 500 };
  GtRange rng08 = { 700, 800 };
  GtRange rng09 = { 200, 400 };
  GtRange rng10 = { 700, 900 };
  GtRange rng11 = { 100, 150 };
  gt_array_add(r1, rng07);
  gt_array_add(r1, rng08);
  gt_array_add(r2, rng09);
  gt_array_add(r2, rng10);
  gt_array_add(r2, rng11);
  ru = gaeval_visitor_union(r1, r2);
  bool test3 = gt_array_size(ru) == 3;

  if(test3)
  {
    GtRange *temp1 = gt_array_get(ru, 0);
    GtRange *temp2 = gt_array_get(ru, 1);
    GtRange *temp3 = gt_array_get(ru, 2);
    GtRange testr1 = { 100, 150 };
    GtRange testr2 = { 200, 500 };
    GtRange testr3 = { 700, 900 };
    test3 = gt_range_compare(temp1, &testr1) == 0 &&
            gt_range_compare(temp2, &testr2) == 0 &&
            gt_range_compare(temp3, &testr3) == 0;
  }
  agn_unit_test_result(test, "union (3)", test3);
  gt_array_delete(r1);
  gt_array_delete(r2);
  gt_array_delete(ru);
}
