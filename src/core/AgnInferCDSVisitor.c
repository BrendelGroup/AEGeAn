#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "AgnInferCDSVisitor.h"
#include "AgnGtExtensions.h"
#include "AgnTestData.h"
#include "AgnUtils.h"

#define agn_infer_cds_visitor_cast(GV)\
        gt_node_visitor_cast(agn_infer_cds_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
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
  AgnLogger *logger;
};


//----------------------------------------------------------------------------//
// Prototypes of private functions
//----------------------------------------------------------------------------//

/**
 * FIXME check phase!
 * FIXME check UTR type!
 */

/**
 * Cast a node visitor object as a AgnInferCDSVisitor
 *
 * @returns    a node visitor object cast as a AgnInferCDSVisitor
 */
static const GtNodeVisitorClass* agn_infer_cds_visitor_class();

/**
 * Identify any mRNA subfeatures associated with this top-level feature and
 * apply the CDS inference procedure.
 *
 * @param[in]  nv       a node visitor
 * @param[in]  fn       node representing a top-level GFF3 feature entry
 * @param[out] error    stores error messages
 * @returns             0 for success
 */
static int visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                              GtError *error);

/**
 * If the mRNA's CDS is discontinuous, ensure each CDS feature is labeled as a
 * multifeature.
 *
 * @param[in]  nv   the node visitor
 */
static void visit_mrna_check_cds_multi(AgnInferCDSVisitor *v);

/**
 * If start codon is provided explicitly, ensure it agrees with CDS, whether the
 * CDS is provided explicitly or implicitly inferred. If start codon is not
 * provided explicitly, infer it from CDS if possible.
 *
 * @param[in]  nv   the node visitor
 */
static void visit_mrna_check_start(AgnInferCDSVisitor *v);

/**
 * If stop codon is provided explicitly, ensure it agrees with CDS, whether the
 * CDS is provided explicitly or implicitly inferred. If stop codon is not
 * provided explicitly, infer it from CDS if possible.
 *
 * @param[in]  nv   the node visitor
 */
static void visit_mrna_check_stop(AgnInferCDSVisitor *v);

/**
 * Infer CDS for any mRNAs that have none specified but have exons and
 * start/stop codons explicitly specified.
 *
 * @param[in]  nv   the node visitor
 */
static void visit_mrna_infer_cds(AgnInferCDSVisitor *v);

/**
 * Infer UTRs for any mRNAs that have none specified but do have exons and
 * start/stop codons and/or CDS explicitly specified.
 *
 * @param[in]  nv   the node visitor
 */
static void visit_mrna_infer_utrs(AgnInferCDSVisitor *v);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

static const GtNodeVisitorClass* agn_infer_cds_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInferCDSVisitor), NULL, NULL,
                                    visit_feature_node, NULL, NULL, NULL);
  }
  return nvc;
}

GtNodeVisitor* agn_infer_cds_visitor_new(AgnLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(agn_infer_cds_visitor_class());
  AgnInferCDSVisitor *v = agn_infer_cds_visitor_cast(nv);
  v->logger = logger;
  v->cdscounter = 0;
  return nv;
}

bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test)
{
  GtArray *genes = agn_test_data_genes_codons();
  GtError *error = gt_error_new();
  AgnLogger *logger = agn_logger_new();
  GtNodeStream *genestream = gt_array_in_stream_new(genes, NULL, error);
  GtNodeVisitor *icv = agn_infer_cds_visitor_new(logger);
  GtNodeStream *icvstream = gt_visitor_stream_new(genestream, icv);
  int result;
  GtGenomeNode *gn;
  GtFeatureNode *fn;

  result = gt_node_stream_next(icvstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  GtArray *cds = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_cds_feature);
  GtArray *utrs = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_utr_feature);
  bool cds1correct;
  if(gt_array_size(cds) != 3)
  {
    cds1correct = false;
  }
  else
  {
    GtGenomeNode **cds1 = gt_array_get(cds, 0);
    GtGenomeNode **cds2 = gt_array_get(cds, 1);
    GtGenomeNode **cds3 = gt_array_get(cds, 2);
    GtRange range1 = gt_genome_node_get_range(*cds1);
    GtRange range2 = gt_genome_node_get_range(*cds2);
    GtRange range3 = gt_genome_node_get_range(*cds3);
    cds1correct = (range1.start == 22167 && range1.end == 22382 &&
                   range2.start == 22497 && range2.end == 22550 &&
                   range3.start == 22651 && range3.end == 23022);
  }
  agn_unit_test_result(test, "CDS correct for mRNA 1", cds1correct);

  bool utrs1correct;
  if(gt_array_size(utrs) != 2)
  {
    utrs1correct = false;
  }
  else
  {
    GtGenomeNode **utr1 = gt_array_get(utrs, 0);
    GtGenomeNode **utr2 = gt_array_get(utrs, 1);
    GtFeatureNode *fn1 = *(GtFeatureNode **)utr1;
    GtFeatureNode *fn2 = *(GtFeatureNode **)utr2;
    GtRange range1 = gt_genome_node_get_range(*utr1);
    GtRange range2 = gt_genome_node_get_range(*utr2);
    utrs1correct = (range1.start == 22057 && range1.end == 22166 &&
                    gt_feature_node_has_type(fn1, "five_prime_UTR") &&
                    range2.start == 23023 && range2.end == 23119 &&
                    gt_feature_node_has_type(fn2, "three_prime_UTR"));
  }
  agn_unit_test_result(test, "UTRs correct for mRNA 1", utrs1correct);
  gt_array_delete(cds);
  gt_array_delete(utrs);
  gt_genome_node_delete(gn);

  result = gt_node_stream_next(icvstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  cds = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_cds_feature);
  utrs = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_utr_feature);
  bool cds2correct;
  if(gt_array_size(cds) != 3)
  {
    cds2correct = false;
  }
  else
  {
    GtGenomeNode **cds1 = gt_array_get(cds, 0);
    GtGenomeNode **cds2 = gt_array_get(cds, 1);
    GtGenomeNode **cds3 = gt_array_get(cds, 2);
    GtRange range1 = gt_genome_node_get_range(*cds1);
    GtRange range2 = gt_genome_node_get_range(*cds2);
    GtRange range3 = gt_genome_node_get_range(*cds3);
    cds2correct = (range1.start == 48411 && range1.end == 48537 &&
                   range2.start == 48637 && range2.end == 48766 &&
                   range3.start == 48870 && range3.end == 48984);
  }
  agn_unit_test_result(test, "CDS correct for mRNA 2", cds2correct);

  bool utrs2correct;
  if(gt_array_size(utrs) != 1)
  {
    utrs2correct = false;
  }
  else
  {
    GtGenomeNode **utr = gt_array_get(utrs, 0);
    GtFeatureNode *fn = *(GtFeatureNode **)utr;
    GtRange range = gt_genome_node_get_range(*utr);
    utrs2correct = (range.start == 48012 && range.end == 48410 &&
                    gt_feature_node_has_type(fn, "three_prime_UTR"));
  }
  agn_unit_test_result(test, "UTRs correct for mRNA 2", utrs2correct);
  gt_array_delete(cds);
  gt_array_delete(utrs);
  gt_genome_node_delete(gn);

  result = gt_node_stream_next(icvstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  cds = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_cds_feature);
  utrs = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_utr_feature);
  bool cds3correct;
  if(gt_array_size(cds) != 6)
  {
    cds3correct = false;
  }
  else
  {
    GtGenomeNode **cds1 = gt_array_get(cds, 0);
    GtGenomeNode **cds2 = gt_array_get(cds, 1);
    GtGenomeNode **cds3 = gt_array_get(cds, 2);
    GtGenomeNode **cds4 = gt_array_get(cds, 3);
    GtGenomeNode **cds5 = gt_array_get(cds, 4);
    GtGenomeNode **cds6 = gt_array_get(cds, 5);
    GtRange range1 = gt_genome_node_get_range(*cds1);
    GtRange range2 = gt_genome_node_get_range(*cds2);
    GtRange range3 = gt_genome_node_get_range(*cds3);
    GtRange range4 = gt_genome_node_get_range(*cds4);
    GtRange range5 = gt_genome_node_get_range(*cds5);
    GtRange range6 = gt_genome_node_get_range(*cds6);
    cds3correct = (range1.start == 88892 && range1.end == 89029 &&
                   range2.start == 89265 && range2.end == 89549 &&
                   range3.start == 90074 && range3.end == 90413 &&
                   range4.start == 90728 && range4.end == 90833 &&
                   range5.start == 91150 && range5.end == 91362 &&
                   range6.start == 91810 && range6.end == 91963);
  }
  agn_unit_test_result(test, "CDS correct for mRNA 3", cds3correct);
  
  bool utrs3correct;
  if(gt_array_size(utrs) != 2)
  {
    utrs3correct = false;
  }
  else
  {
    GtGenomeNode **utr1 = gt_array_get(utrs, 0);
    GtGenomeNode **utr2 = gt_array_get(utrs, 1);
    GtFeatureNode *fn1 = *(GtFeatureNode **)utr1;
    GtFeatureNode *fn2 = *(GtFeatureNode **)utr2;
    GtRange range1 = gt_genome_node_get_range(*utr1);
    GtRange range2 = gt_genome_node_get_range(*utr2);
    utrs3correct = (range1.start == 88551 && range1.end == 88891 &&
                    gt_feature_node_has_type(fn1, "three_prime_UTR") &&
                    range2.start == 91964 && range2.end == 92176 &&
                    gt_feature_node_has_type(fn2, "five_prime_UTR"));
  }
  agn_unit_test_result(test, "UTRs correct for mRNA 3", utrs3correct);
  gt_array_delete(cds);
  gt_array_delete(utrs);
  gt_genome_node_delete(gn);

  agn_logger_delete(logger);
  gt_node_stream_delete(icvstream);
  gt_array_delete(genes);
  gt_node_stream_delete(genestream);
  gt_error_delete(error);
  return cds1correct && utrs1correct && cds2correct && utrs2correct &&
         cds3correct && utrs3correct;
}

static int visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                              GtError *error)
{
  AgnInferCDSVisitor *v = agn_infer_cds_visitor_cast(nv);
  gt_error_check(error);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_gt_feature_node_is_mrna_feature(current))
      continue;

    v->cds    = agn_gt_feature_node_children_of_type(current,
                                            agn_gt_feature_node_is_cds_feature);
    v->utrs   = agn_gt_feature_node_children_of_type(current,
                                            agn_gt_feature_node_is_utr_feature);
    v->exons  = agn_gt_feature_node_children_of_type(current,
                                           agn_gt_feature_node_is_exon_feature);
    v->starts = agn_gt_feature_node_children_of_type(current,
                                    agn_gt_feature_node_is_start_codon_feature);
    v->stops  = agn_gt_feature_node_children_of_type(current,
                                     agn_gt_feature_node_is_stop_codon_feature);
    v->mrna   = current;

    visit_mrna_infer_cds(v);
    visit_mrna_check_start(v);
    visit_mrna_check_stop(v);
    visit_mrna_infer_utrs(v);
    visit_mrna_check_cds_multi(v);

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

static void visit_mrna_check_cds_multi(AgnInferCDSVisitor *v)
{
  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
  if(gt_array_size(v->cds) <= 1)
  {
    if(gt_array_size(v->cds) == 0)
    {
      agn_logger_log_error(v->logger, "error inferring CDS from codons and "
                           "exons for mRNA '%s' (line %u)", mrnaid, ln);
    }
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

static void visit_mrna_check_start(AgnInferCDSVisitor *v)
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
    agn_logger_log_error(v->logger, "mRNA '%s' (line %u) has %lu start codons",
                         mrnaid, ln, gt_array_size(v->starts));
  }
  else if(gt_array_size(v->starts) == 1)
  {
    GtGenomeNode **codon = gt_array_get(v->starts, 0);
    GtRange testrange = gt_genome_node_get_range(*codon);
    if(gt_range_compare(&startrange, &testrange) != 0)
    {
      agn_logger_log_error(v->logger, "start codon inferred from CDS [%lu, "
                           "%lu] does not match explicitly provided start "
                           "codon [%lu, %lu] for mRNA '%s'", startrange.start,
                           startrange.end, testrange.start, testrange.end,
                           mrnaid);
    }
  }
  else // gt_assert(gt_array_size(v->starts) == 0)
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

static void visit_mrna_check_stop(AgnInferCDSVisitor *v)
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
    agn_logger_log_error(v->logger, "mRNA '%s' (line %u) has %lu stop codons",
                         mrnaid, ln, gt_array_size(v->starts));
  }
  else if(gt_array_size(v->stops) == 1)
  {
    GtGenomeNode **codon = gt_array_get(v->stops, 0);
    GtRange testrange = gt_genome_node_get_range(*codon);
    if(gt_range_compare(&stoprange, &testrange) != 0)
    {
      agn_logger_log_error(v->logger, "stop codon inferred from CDS [%lu, "
                           "%lu] does not match explicitly provided stop "
                           "codon [%lu, %lu] for mRNA '%s'", stoprange.start,
                           stoprange.end, testrange.start, testrange.end,
                           mrnaid);
    }
  }
  else // gt_assert(gt_array_size(v->stops) == 0)
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

static void visit_mrna_infer_cds(AgnInferCDSVisitor *v)
{
  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
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
    agn_logger_log_error(v->logger, "cannot infer missing CDS for mRNA '%s'"
                         "(line %u) without exons and start/stop codons",
                         mrnaid, ln);
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
    bool exon_includes_cds = agn_infer_cds_range_from_exon_and_codons(
                                 &exon_range, &left_codon_range,
                                 &right_codon_range, &cdsrange);
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

static void visit_mrna_infer_utrs(AgnInferCDSVisitor *v)
{
  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
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
    agn_logger_log_error(v->logger, "cannot infer missing UTRs for mRNA '%s'"
                         "(line %u) without exons and start/stop codons or CDS",
                         mrnaid, ln);
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
