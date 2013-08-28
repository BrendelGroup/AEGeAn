#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "AgnInferCDSVisitor.h"
#include "AgnGtExtensions.h"
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
const GtNodeVisitorClass* agn_infer_cds_visitor_class();

/**
 * Destructor for the AgnInferCDSVisitor class
 *
 * @param[in] nv    the node visitor object
 */
static void infer_cds_visitor_free(GtNodeVisitor *nv);

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
void visit_mrna_check_cds_multi(AgnInferCDSVisitor *v);

/**
 * If start codon is provided explicitly, ensure it agrees with CDS, whether the
 * CDS is provided explicitly or implicitly inferred. If start codon is not
 * provided explicitly, infer it from CDS if possible.
 *
 * @param[in]  nv   the node visitor
 */
void visit_mrna_check_start(AgnInferCDSVisitor *v);

/**
 * If stop codon is provided explicitly, ensure it agrees with CDS, whether the
 * CDS is provided explicitly or implicitly inferred. If stop codon is not
 * provided explicitly, infer it from CDS if possible.
 *
 * @param[in]  nv   the node visitor
 */
void visit_mrna_check_stop(AgnInferCDSVisitor *v);

/**
 * Infer CDS for any mRNAs that have none specified but have exons and
 * start/stop codons explicitly specified.
 *
 * @param[in]  nv   the node visitor
 */
void visit_mrna_infer_cds(AgnInferCDSVisitor *v);

/**
 * Infer UTRs for any mRNAs that have none specified but do have exons and
 * start/stop codons and/or CDS explicitly specified.
 *
 * @param[in]  nv   the node visitor
 */
void visit_mrna_infer_utrs(AgnInferCDSVisitor *v);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

const GtNodeVisitorClass* agn_infer_cds_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInferCDSVisitor),
                                    infer_cds_visitor_free,
                                    NULL,
                                    visit_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

static void infer_cds_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED AgnInferCDSVisitor *v = agn_infer_cds_visitor_cast(nv);
}

GtNodeVisitor* agn_infer_cds_visitor_new(AgnLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(agn_infer_cds_visitor_class());
  AgnInferCDSVisitor *v = agn_infer_cds_visitor_cast(nv);
  v->logger = logger;
  return nv;
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

void visit_mrna_check_cds_multi(AgnInferCDSVisitor *v)
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
  gt_feature_node_make_multi_representative(*firstsegment);
  unsigned long i;
  for(i = 0; i < gt_array_size(v->cds); i++)
  {
    GtFeatureNode **segment = gt_array_get(v->cds, i);
    if(!gt_feature_node_is_multi(*segment))
    {
      gt_feature_node_set_multi_representative(*segment, *firstsegment);
    }
  }
}

void visit_mrna_check_start(AgnInferCDSVisitor *v)
{
  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
  GtStrand strand = gt_feature_node_get_strand(v->mrna);

  GtRange startrange;
  GtGenomeNode **fiveprimesegment = gt_array_get(v->cds, 0);
  startrange = gt_genome_node_get_range(*fiveprimesegment);
  startrange.end = startrange.start + 2;
  if(strand == GT_STRAND_REVERSE)
  {
    unsigned long fiveprimeindex = gt_array_size(v->cds) - 1;
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
    gt_feature_node_add_attribute(cf, "Parent", mrnaid);
    gt_array_add(v->starts, cf);
  }
}

void visit_mrna_check_stop(AgnInferCDSVisitor *v)
{
  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
  GtStrand strand = gt_feature_node_get_strand(v->mrna);

  GtRange stoprange;
  unsigned long threeprimeindex = gt_array_size(v->cds) - 1;
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
    gt_feature_node_add_attribute(cf, "Parent", mrnaid);
    gt_array_add(v->stops, cf);
  }
}

void visit_mrna_infer_cds(AgnInferCDSVisitor *v)
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
  unsigned long i;
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
      GtGenomeNode *cdsfeat = gt_feature_node_new
      (
        gt_genome_node_get_seqid(exon_gn), "CDS", cdsrange.start,
        cdsrange.end, exon_strand
      );
      gt_feature_node_add_child(v->mrna, (GtFeatureNode *)cdsfeat);
      gt_feature_node_add_attribute((GtFeatureNode *)cdsfeat, "Parent", mrnaid);
      gt_array_add(v->cds, cdsfeat);
    }
  }
}

void visit_mrna_infer_utrs(AgnInferCDSVisitor *v)
{
  const char *mrnaid = gt_feature_node_get_attribute(v->mrna, "ID");
  unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)v->mrna);
  GtFeatureNode *start_codon, *stop_codon;

  bool exonsexplicit    = gt_array_size(v->exons) > 0;
  bool startcodon_check = gt_array_size(v->starts) == 1 &&
                          (start_codon = gt_array_get(v->starts, 0)) != NULL;
  bool stopcodon_check  = gt_array_size(v->stops)  == 1 &&
                          (stop_codon  = gt_array_get(v->stops,  0)) != NULL;

  if(gt_array_size(v->utrs) > 0)
  {
    return;
  }
  else if(!exonsexplicit || !startcodon_check || !stopcodon_check)
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

  unsigned long i;
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
      gt_feature_node_add_attribute((GtFeatureNode *)utr, "Parent", mrnaid);
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
      gt_feature_node_add_attribute((GtFeatureNode *)utr, "Parent", mrnaid);
      gt_array_add(v->utrs, utr);
    }
  }
}
