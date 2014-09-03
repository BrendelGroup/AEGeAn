/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include "core/array_api.h"
#include "AgnFilterStream.h"
#include "AgnInferCDSVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define infer_cds_visitor_cast(GV)\
        gt_node_visitor_cast(infer_cds_visitor_class(), GV)

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

struct AgnInferCDSVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureNode *mrna;
  GtArray *cds;
  GtArray *utrs;
  GtArray *exons;
  GtArray *starts;
  GtArray *stops;
  GtUword cdscounter;
  GtLogger *logger;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function If the mRNA's CDS is discontinuous, ensure each CDS feature is
 * labeled as a multifeature.
 */
static void infer_cds_visitor_check_cds_multi(AgnInferCDSVisitor *v);

/**
 * @function Check and correct phase attributes for each CDS feature.
 */
static void infer_cds_visitor_check_cds_phase(AgnInferCDSVisitor *v);

/**
 * @function If start codon is provided explicitly, ensure it agrees with CDS,
 * whether the CDS is provided explicitly or implicitly inferred. If start codon
 * is not provided explicitly, infer it from CDS if possible.
 */
static void infer_cds_visitor_check_start(AgnInferCDSVisitor *v);

/**
 * @function If stop codon is provided explicitly, ensure it agrees with CDS,
 * whether the CDS is provided explicitly or implicitly inferred. If stop codon
 * is not provided explicitly, infer it from CDS if possible.
 */
static void infer_cds_visitor_check_stop(AgnInferCDSVisitor *v);

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *infer_cds_visitor_class();

/**
 * @function Infer CDS for any mRNAs that have none specified but have exons and
 * start/stop codons explicitly specified.
 */
static void infer_cds_visitor_infer_cds(AgnInferCDSVisitor *v);

/**
 * @function Given an exon and the start/stop codons associated with its
 * corresponding mRNA, determine which parts of the exon (if any) correspond to
 * coding sequence. If the exon contains coding sequence, the range of that
 * coding sequence will be stored in ``cds_range`` and the function will return
 * true. Otherwise, the function will return false. If the mRNA is on the
 * forward strand, ``left_codon_range`` should contain the coordinates for the
 * start codon and ``right_codon_range`` should contain coordinates for the stop
 * codon. If the mRNA is on the reverse strand, these should be swapped.
 */
static bool infer_cds_visitor_infer_range(GtRange *exon_range,
                                          GtRange *leftcodon_range,
                                          GtRange *rightcodon_range,
                                          GtRange *cds_range);

/**
 * @function Infer UTRs for any mRNAs that have none specified but do have exons
 * and start/stop codons and/or CDS explicitly specified.
 */
static void infer_cds_visitor_infer_utrs(AgnInferCDSVisitor *v);

/**
 * @function Ensure that UTR types are correctly encoded.
 */
static void infer_cds_visitor_set_utrs(AgnInferCDSVisitor *v);

/**
 * @function Generate data for unit testing.
 */
static void infer_cds_visitor_test_data(GtQueue *queue);

/**
 * @function Identify any mRNA subfeatures associated with this top-level
 * feature and apply the CDS inference procedure.
 */
static int infer_cds_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                GtFeatureNode *fn,
                                                GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_infer_cds_stream_new(GtNodeStream *in, GtLogger *logger)
{
  GtNodeVisitor *nv = agn_infer_cds_visitor_new(logger);
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_infer_cds_visitor_new(GtLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(infer_cds_visitor_class());
  AgnInferCDSVisitor *v = infer_cds_visitor_cast(nv);
  v->logger = logger;
  v->cdscounter = 0;
  return nv;
}

bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  infer_cds_visitor_test_data(queue);
  agn_assert(gt_queue_size(queue) == 4);

  GtFeatureNode *fn = gt_queue_get(queue);
  GtArray *cds = agn_typecheck_select(fn, agn_typecheck_cds);
  bool grape1 = (gt_array_size(cds) == 4);
  if(grape1)
  {
    GtGenomeNode *cds2 = *(GtGenomeNode **)gt_array_get(cds, 1);
    GtRange range = gt_genome_node_get_range(cds2);
    grape1 = (range.start == 349 && range.end == 522);
  }
  agn_unit_test_result(test, "grape test sans UTRs", grape1);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(cds);

  fn = gt_queue_get(queue);
  cds = agn_typecheck_select(fn, agn_typecheck_cds);
  bool grape2 = (gt_array_size(cds) == 1);
  if(grape2)
  {
    GtGenomeNode *cds1 = *(GtGenomeNode **)gt_array_get(cds, 0);
    GtRange range = gt_genome_node_get_range(cds1);
    GtStrand strand = gt_feature_node_get_strand((GtFeatureNode *)cds1);
    grape2 = (range.start == 10747 && range.end == 11577 &&
              strand == GT_STRAND_REVERSE);
  }
  agn_unit_test_result(test, "grape test with UTRs, strand check", grape2);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(cds);

  fn = gt_queue_get(queue);
  cds = agn_typecheck_select(fn, agn_typecheck_cds);
  bool grape3 = (gt_array_size(cds) == 2);
  if(grape3)
  {
    GtGenomeNode *cds2 = *(GtGenomeNode **)gt_array_get(cds, 1);
    GtRange range = gt_genome_node_get_range(cds2);
    grape3 = (range.start == 22651 && range.end == 23022);
  }
  agn_unit_test_result(test, "grape test 3", grape3);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(cds);

  fn = gt_queue_get(queue);
  cds = agn_typecheck_select(fn, agn_typecheck_cds);
  bool grape4 = (gt_array_size(cds) == 12);
  if(grape4)
  {
    GtGenomeNode *cds7 = *(GtGenomeNode **)gt_array_get(cds, 6);
    GtRange range = gt_genome_node_get_range(cds7);
    grape4 = (range.start == 27956 && range.end == 27996);
  }
  agn_unit_test_result(test, "grape test 4", grape4);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(cds);

  while(gt_queue_size(queue) > 0)
  {
    GtGenomeNode *cds_n = gt_queue_get(queue);
    gt_genome_node_delete(cds_n);
  }
  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static void infer_cds_visitor_check_cds_multi(AgnInferCDSVisitor *v)
{
  if(gt_array_size(v->cds) <= 1)
  {
    return;
  }

  GtFeatureNode **firstsegment = gt_array_get(v->cds, 0);
  const char *id = gt_feature_node_get_attribute(*firstsegment, "ID");
  if(id == NULL)
  {
    char newid[64];
    sprintf(newid, "CDS%lu", v->cdscounter++);
    gt_feature_node_add_attribute(*firstsegment, "ID", newid);
  }
  gt_feature_node_make_multi_representative(*firstsegment);
  GtUword i;
  for(i = 0; i < gt_array_size(v->cds); i++)
  {
    GtFeatureNode **segment = gt_array_get(v->cds, i);
    if(!gt_feature_node_is_multi(*segment))
    {
      gt_feature_node_set_multi_representative(*segment, *firstsegment);
    }
  }
}

static void infer_cds_visitor_check_cds_phase(AgnInferCDSVisitor *v)
{
  unsigned long num_cds_feats = gt_array_size(v->cds);
  if(num_cds_feats == 0)
    return;

  GtFeatureNode *cdsf1 = *(GtFeatureNode **)gt_array_get(v->cds, 0);
  gt_feature_node_set_phase(cdsf1, GT_PHASE_ZERO);
  if(num_cds_feats == 1)
    return;

  GtStrand strand = gt_feature_node_get_strand(cdsf1);
  unsigned long cds_length = gt_genome_node_get_length((GtGenomeNode *)cdsf1);
  unsigned long i;
  if(strand == GT_STRAND_REVERSE)
  {
    for(i = num_cds_feats - 1; i > 0; i--)
    {
      GtFeatureNode *cds = *(GtFeatureNode **)gt_array_get(v->cds, i);
      int phasenum = cds_length % 3;
      GtPhase phase = GT_PHASE_ZERO;
      if(phasenum == 1)
        phase = GT_PHASE_TWO;
      else if(phasenum == 2)
        phase = GT_PHASE_ONE;
      gt_feature_node_set_phase(cds, phase);
      cds_length += gt_genome_node_get_length((GtGenomeNode *)cds);
    }
  }
  else
  {
    for(i = 1; i < num_cds_feats; i++)
    {
      GtFeatureNode *cds = *(GtFeatureNode **)gt_array_get(v->cds, i);
      int phasenum = cds_length % 3;
      GtPhase phase = GT_PHASE_ZERO;
      if(phasenum == 1)
        phase = GT_PHASE_TWO;
      else if(phasenum == 2)
        phase = GT_PHASE_ONE;
      gt_feature_node_set_phase(cds, phase);
      cds_length += gt_genome_node_get_length((GtGenomeNode *)cds);
    }
  }
}

static void infer_cds_visitor_check_start(AgnInferCDSVisitor *v)
{
  if(gt_array_size(v->cds) == 0)
    return;

  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
  GtStrand strand = gt_feature_node_get_strand(v->mrna);

  GtRange startrange;
  GtGenomeNode **fiveprimesegment = gt_array_get(v->cds, 0);
  startrange = gt_genome_node_get_range(*fiveprimesegment);
  startrange.end = startrange.start + 2;
  if(strand == GT_STRAND_REVERSE)
  {
    GtUword fiveprimeindex = gt_array_size(v->cds) - 1;
    fiveprimesegment = gt_array_get(v->cds, fiveprimeindex);
    startrange = gt_genome_node_get_range(*fiveprimesegment);
    startrange.start = startrange.end - 2;
  }

  if(gt_array_size(v->starts) > 1)
  {
    gt_logger_log(v->logger, "mRNA '%s' (line %u) has %lu start codons", mrnaid,
                  ln, gt_array_size(v->starts));
  }
  else if(gt_array_size(v->starts) == 1)
  {
    GtGenomeNode **codon = gt_array_get(v->starts, 0);
    GtRange testrange = gt_genome_node_get_range(*codon);
    if(gt_range_compare(&startrange, &testrange) != 0)
    {
      gt_logger_log(v->logger, "start codon inferred from CDS [%lu, %lu] does "
                    "not match explicitly provided start codon [%lu, %lu] for "
                    "mRNA '%s'", startrange.start, startrange.end,
                    testrange.start, testrange.end, mrnaid);
    }
  }
  else // agn_assert(gt_array_size(v->starts) == 0)
  {
    GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)v->mrna);
    GtGenomeNode *codonfeature = gt_feature_node_new(seqid, "start_codon",
                                                     startrange.start,
                                                     startrange.end,
                                                     strand);
    GtFeatureNode *cf = (GtFeatureNode *)codonfeature;
    gt_feature_node_add_child(v->mrna, cf);
    gt_array_add(v->starts, cf);
  }
}

static void infer_cds_visitor_check_stop(AgnInferCDSVisitor *v)
{
  if(gt_array_size(v->cds) == 0)
    return;

  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
  GtStrand strand = gt_feature_node_get_strand(v->mrna);

  GtRange stoprange;
  GtUword threeprimeindex = gt_array_size(v->cds) - 1;
  GtGenomeNode **threeprimesegment = gt_array_get(v->cds, threeprimeindex);
  stoprange = gt_genome_node_get_range(*threeprimesegment);
  stoprange.start = stoprange.end - 2;
  if(strand == GT_STRAND_REVERSE)
  {
    threeprimesegment = gt_array_get(v->cds, 0);
    stoprange = gt_genome_node_get_range(*threeprimesegment);
    stoprange.end = stoprange.start + 2;
  }

  if(gt_array_size(v->stops) > 1)
  {
    gt_logger_log(v->logger, "mRNA '%s' (line %u) has %lu stop codons", mrnaid,
                  ln, gt_array_size(v->starts));
  }
  else if(gt_array_size(v->stops) == 1)
  {
    GtGenomeNode **codon = gt_array_get(v->stops, 0);
    GtRange testrange = gt_genome_node_get_range(*codon);
    if(gt_range_compare(&stoprange, &testrange) != 0)
    {
      gt_logger_log(v->logger, "stop codon inferred from CDS [%lu, %lu] does "
                    "not match explicitly provided stop codon [%lu, %lu] for "
                    "mRNA '%s'", stoprange.start, stoprange.end,
                    testrange.start, testrange.end, mrnaid);
    }
  }
  else // agn_assert(gt_array_size(v->stops) == 0)
  {
    GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)v->mrna);
    GtGenomeNode *codonfeature = gt_feature_node_new(seqid, "stop_codon",
                                                     stoprange.start,
                                                     stoprange.end,
                                                     strand);
    GtFeatureNode *cf = (GtFeatureNode *)codonfeature;
    gt_feature_node_add_child(v->mrna, cf);
    gt_array_add(v->stops, cf);
  }
}

static const GtNodeVisitorClass *infer_cds_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInferCDSVisitor), NULL, NULL,
                                    infer_cds_visitor_visit_feature_node, NULL,
                                    NULL, NULL);
  }
  return nvc;
}

static void infer_cds_visitor_infer_cds(AgnInferCDSVisitor *v)
{
  GtFeatureNode **start_codon, **stop_codon;

  bool exonsexplicit    = gt_array_size(v->exons) > 0;
  bool startcodon_check = gt_array_size(v->starts) == 1 &&
                          (start_codon = gt_array_get(v->starts, 0)) != NULL;
  bool stopcodon_check  = gt_array_size(v->stops)  == 1 &&
                          (stop_codon  = gt_array_get(v->stops,  0)) != NULL;

  if(gt_array_size(v->cds) > 0)
  {
    return;
  }
  else if(!exonsexplicit || !startcodon_check || !stopcodon_check)
  {
    return;
  }

  GtRange left_codon_range, right_codon_range;
  left_codon_range  = gt_genome_node_get_range(*(GtGenomeNode **)start_codon);
  right_codon_range = gt_genome_node_get_range(*(GtGenomeNode **)stop_codon);
  if(gt_feature_node_get_strand(v->mrna) == GT_STRAND_REVERSE)
  {
    left_codon_range  = gt_genome_node_get_range(*(GtGenomeNode **)stop_codon);
    right_codon_range = gt_genome_node_get_range(*(GtGenomeNode **)start_codon);
  }
  GtUword i;
  for(i = 0; i < gt_array_size(v->exons); i++)
  {
    GtFeatureNode *exon = *(GtFeatureNode **)gt_array_get(v->exons, i);
    GtGenomeNode *exon_gn = (GtGenomeNode *)exon;
    GtRange exon_range = gt_genome_node_get_range(exon_gn);
    GtStrand exon_strand = gt_feature_node_get_strand(exon);

    GtRange cdsrange;
    bool exon_includes_cds = infer_cds_visitor_infer_range(&exon_range,
                                                           &left_codon_range,
                                                           &right_codon_range,
                                                           &cdsrange);
    if(exon_includes_cds)
    {
      GtGenomeNode *cdsfeat;
      cdsfeat = gt_feature_node_new(gt_genome_node_get_seqid(exon_gn), "CDS",
                                    cdsrange.start, cdsrange.end, exon_strand);
      gt_feature_node_add_child(v->mrna, (GtFeatureNode *)cdsfeat);
      gt_array_add(v->cds, cdsfeat);
    }
  }
}

static bool infer_cds_visitor_infer_range(GtRange *exon_range,
                                          GtRange *leftcodon_range,
                                          GtRange *rightcodon_range,
                                          GtRange *cds_range)
{
  cds_range->start = 0;
  cds_range->end   = 0;

  // UTR
  if(exon_range->end < leftcodon_range->start ||
     exon_range->start > rightcodon_range->end)
    return false;

  bool overlap_left  = gt_range_overlap(exon_range, leftcodon_range);
  bool overlap_right = gt_range_overlap(exon_range, rightcodon_range);
  if(overlap_left && overlap_right)
  {
    cds_range->start = leftcodon_range->start;
    cds_range->end   = rightcodon_range->end;
  }
  else if(overlap_left)
  {
    cds_range->start = leftcodon_range->start;
    cds_range->end   = exon_range->end;
  }
  else if(overlap_right)
  {
    cds_range->start = exon_range->start;
    cds_range->end   = rightcodon_range->end;
  }
  else
  {
    cds_range->start = exon_range->start;
    cds_range->end   = exon_range->end;
  }

  return true;
}

static void infer_cds_visitor_infer_utrs(AgnInferCDSVisitor *v)
{
  GtFeatureNode *start_codon, *stop_codon;

  bool exonsexplicit    = gt_array_size(v->exons) > 0;
  bool cdsexplicit      = gt_array_size(v->cds) > 0;
  bool startcodon_check = gt_array_size(v->starts) == 1 &&
                          (start_codon = gt_array_get(v->starts, 0)) != NULL;
  bool stopcodon_check  = gt_array_size(v->stops)  == 1 &&
                          (stop_codon  = gt_array_get(v->stops,  0)) != NULL;
  bool caninferutrs     = exonsexplicit && startcodon_check && stopcodon_check;

  if(gt_array_size(v->utrs) > 0)
  {
    return;
  }
  else if(!cdsexplicit && !caninferutrs)
  {
    return;
  }

  GtGenomeNode **leftcodon = gt_array_get(v->starts, 0);
  GtGenomeNode **rightcodon = gt_array_get(v->stops, 0);
  GtStrand strand = gt_feature_node_get_strand(v->mrna);
  const char *lefttype  = "five_prime_UTR";
  const char *righttype = "three_prime_UTR";
  if(strand == GT_STRAND_REVERSE)
  {
    lefttype   = "three_prime_UTR";
    righttype  = "five_prime_UTR";
    void *temp = leftcodon;
    leftcodon  = rightcodon;
    rightcodon = temp;
  }
  GtRange leftrange  = gt_genome_node_get_range(*leftcodon);
  GtRange rightrange = gt_genome_node_get_range(*rightcodon);

  GtUword i;
  for(i = 0; i < gt_array_size(v->exons); i++)
  {
    GtGenomeNode **exon = gt_array_get(v->exons, i);
    GtRange exonrange = gt_genome_node_get_range(*exon);
    if(exonrange.start < leftrange.start)
    {
      GtRange utrrange;
      if(gt_range_overlap(&exonrange, &leftrange))
      {
        utrrange.start = exonrange.start;
        utrrange.end   = leftrange.start - 1;
      }
      else
      {
        utrrange = exonrange;
      }
      GtGenomeNode *utr = gt_feature_node_new(gt_genome_node_get_seqid(*exon),
                                              lefttype, utrrange.start,
                                              utrrange.end, strand);
      gt_feature_node_add_child(v->mrna, (GtFeatureNode *)utr);
      gt_array_add(v->utrs, utr);
    }

    if(exonrange.end > rightrange.end)
    {
      GtRange utrrange;
      if(gt_range_overlap(&exonrange, &rightrange))
      {
        utrrange.start = rightrange.end + 1;
        utrrange.end   = exonrange.end;
      }
      else
      {
        utrrange = exonrange;
      }
      GtGenomeNode *utr = gt_feature_node_new(gt_genome_node_get_seqid(*exon),
                                              righttype, utrrange.start,
                                              utrrange.end, strand);
      gt_feature_node_add_child(v->mrna, (GtFeatureNode *)utr);
      gt_array_add(v->utrs, utr);
    }
  }
}

static void infer_cds_visitor_set_utrs(AgnInferCDSVisitor *v)
{
  GtGenomeNode **start;
  GtUword i, cds_start;

  if(!v->starts || gt_array_size(v->starts) != 1)
    return;
  start = gt_array_get(v->starts, 0);
  cds_start = gt_genome_node_get_start(*start);

  for(i = 0; i < gt_array_size(v->utrs); i++)
  {
    GtFeatureNode *utr = *(GtFeatureNode **)gt_array_get(v->utrs, i);
    GtStrand strand = gt_feature_node_get_strand(utr);
    GtUword utr_start = gt_genome_node_get_start((GtGenomeNode *)utr);

    if(!gt_feature_node_has_type(utr, "five_prime_UTR") &&
       !gt_feature_node_has_type(utr, "three_prime_UTR"))
    {
      if(strand == GT_STRAND_FORWARD)
      {
        if(utr_start < cds_start)
          gt_feature_node_set_type(utr, "five_prime_UTR");
        else
          gt_feature_node_set_type(utr, "three_prime_UTR");
      }
      else
      {
        if(utr_start < cds_start)
          gt_feature_node_set_type(utr, "three_prime_UTR");
        else
          gt_feature_node_set_type(utr, "five_prime_UTR");
      }
    }
  }
}

static void infer_cds_visitor_test_data(GtQueue *queue)
{
  GtError *error = gt_error_new();
  const char *file = "data/gff3/grape-codons.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtNodeStream *icv_stream = agn_infer_cds_stream_new(gff3in, logger);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(icv_stream, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnInferCDSVisitor::infer_cds_visitor_test_data] error "
            "processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(icv_stream);
  gt_node_stream_delete(arraystream);
  gt_logger_delete(logger);
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

static int infer_cds_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                GtFeatureNode *fn,
                                                GtError *error)
{
  AgnInferCDSVisitor *v = infer_cds_visitor_cast(nv);
  gt_error_check(error);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_typecheck_mrna(current))
      continue;

    v->cds    = agn_typecheck_select(current, agn_typecheck_cds);
    v->utrs   = agn_typecheck_select(current, agn_typecheck_utr);
    v->exons  = agn_typecheck_select(current, agn_typecheck_exon);
    v->starts = agn_typecheck_select(current, agn_typecheck_start_codon);
    v->stops  = agn_typecheck_select(current, agn_typecheck_stop_codon);
    v->mrna   = current;

    infer_cds_visitor_infer_cds(v);
    infer_cds_visitor_check_start(v);
    infer_cds_visitor_check_stop(v);
    infer_cds_visitor_infer_utrs(v);
    infer_cds_visitor_check_cds_multi(v);
    infer_cds_visitor_check_cds_phase(v);
    infer_cds_visitor_set_utrs(v);

    v->mrna = NULL;
    gt_array_delete(v->cds);
    gt_array_delete(v->utrs);
    gt_array_delete(v->exons);
    gt_array_delete(v->starts);
    gt_array_delete(v->stops);
  }
  gt_feature_node_iterator_delete(iter);

  return 0;
}
