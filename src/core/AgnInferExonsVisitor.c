#include <string.h>
#include "core/array_api.h"
#include "core/queue_api.h"
#include "AgnInferExonsVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define infer_exons_visitor_cast(GV)\
        gt_node_visitor_cast(infer_exons_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//

struct AgnInferExonsVisitor
{
  const GtNodeVisitor parent_instance;
  GtFeatureNode *gene;
  GtIntervalTree *exonsbyrange;
  GtIntervalTree *intronsbyrange;
  GtArray *exons;
  GtArray *introns;
  GtLogger *logger;
};


//----------------------------------------------------------------------------//
// Prototypes of private functions
//----------------------------------------------------------------------------//

/**
 * @function Cast a node visitor object as a AgnInferExonsVisitor.
 */
static const GtNodeVisitorClass* infer_exons_visitor_class();

/**
 * @function Generate data for unit testing.
 */
static void infer_exons_visitor_test_data(GtQueue *queue);

/**
 * @function Procedure for processing feature nodes (the only node of interest
 * for this node visitor).
 */
static int
infer_exons_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                       GtError *error);

/**
 * @function If an exon with the same coordinates already exists and belongs to
 * another mRNA, associate it with this mRNA as well instead of creating a
 * duplicate feature.
 */
static bool
infer_exons_visitor_visit_gene_collapse_feature(AgnInferExonsVisitor *v,
                                                GtFeatureNode *mrna,
                                                GtRange *range,
                                                GtIntervalTree *featsbyrange);

/**
 * @function Infer exons from CDS and UTR segments if possible.
 */
static void
infer_exons_visitor_visit_gene_infer_exons(AgnInferExonsVisitor *v);

/**
 * @function Infer introns from exons if necessary.
 */
static void
infer_exons_visitor_visit_gene_infer_introns(AgnInferExonsVisitor *v);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

GtNodeStream* agn_infer_exons_stream_new(GtNodeStream *in, GtLogger *logger)
{
  GtNodeVisitor *nv = agn_infer_exons_visitor_new(logger);
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor* agn_infer_exons_visitor_new(GtLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(infer_exons_visitor_class());
  AgnInferExonsVisitor *v = infer_exons_visitor_cast(nv);
  v->logger = logger;
  return nv;
}

bool agn_infer_exons_visitor_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  infer_exons_visitor_test_data(queue);
  gt_assert(gt_queue_size(queue) == 4);

  GtFeatureNode *fn = gt_queue_get(queue);
  GtArray *exons = agn_typecheck_select(fn, agn_typecheck_exon);
  bool grape1 = (gt_array_size(exons) == 4);
  if(grape1)
  {
    GtGenomeNode *exon2 = *(GtGenomeNode **)gt_array_get(exons, 1);
    GtRange range = gt_genome_node_get_range(exon2);
    grape1 = (range.start == 349 && range.end == 522);
  }
  agn_unit_test_result(test, "grape test sans UTRs", grape1);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(exons);

  fn = gt_queue_get(queue);
  exons = agn_typecheck_select(fn, agn_typecheck_exon);
  bool grape2 = (gt_array_size(exons) == 1);
  if(grape2)
  {
    GtGenomeNode *exon1 = *(GtGenomeNode **)gt_array_get(exons, 0);
    GtRange range = gt_genome_node_get_range(exon1);
    GtStrand strand = gt_feature_node_get_strand((GtFeatureNode *)exon1);
    grape2 = (range.start == 10538 && range.end == 11678 &&
              strand == GT_STRAND_REVERSE);
  }
  agn_unit_test_result(test, "grape test with UTRs, strand check", grape2);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(exons);

  fn = gt_queue_get(queue);
  exons = agn_typecheck_select(fn, agn_typecheck_exon);
  bool grape3 = (gt_array_size(exons) == 2);
  if(grape3)
  {
    GtGenomeNode *exon2 = *(GtGenomeNode **)gt_array_get(exons, 1);
    GtRange range = gt_genome_node_get_range(exon2);
    grape3 = (range.start == 22651 && range.end == 23448);
  }
  agn_unit_test_result(test, "grape test 3", grape3);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(exons);

  fn = gt_queue_get(queue);
  exons = agn_typecheck_select(fn, agn_typecheck_exon);
  bool grape4 = (gt_array_size(exons) == 12);
  if(grape4)
  {
    GtGenomeNode *cds7 = *(GtGenomeNode **)gt_array_get(exons, 6);
    GtRange range = gt_genome_node_get_range(cds7);
    grape4 = (range.start == 27956 && range.end == 27996);
  }
  agn_unit_test_result(test, "grape test 4", grape4);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(exons);

  while(gt_queue_size(queue) > 0)
  {
    GtGenomeNode *cds_n = gt_queue_get(queue);
    gt_genome_node_delete(cds_n);
  }
  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static const GtNodeVisitorClass* infer_exons_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInferExonsVisitor), NULL, NULL,
                                    infer_exons_visitor_visit_feature_node,
                                    NULL, NULL, NULL);
  }
  return nvc;
}


static void infer_exons_visitor_test_data(GtQueue *queue)
{
  GtError *error = gt_error_new();
  const char *file = "data/gff3/grape-utrs.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtNodeStream *iev_stream = agn_infer_exons_stream_new(gff3in, logger);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(iev_stream, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnInferExonsVisitor::infer_exons_visitor_test_data] "
            "error processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(iev_stream);
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

static int
infer_exons_visitor_visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                       GtError *error)
{
  AgnInferExonsVisitor *v = infer_exons_visitor_cast(nv);
  gt_error_check(error);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    if(!agn_typecheck_gene(current) && !agn_typecheck_transcript(current))
      continue;
    
    GtUword i;
    v->gene = current;

    v->exonsbyrange = gt_interval_tree_new(NULL);
    v->exons = agn_typecheck_select(current, agn_typecheck_exon);
    for(i = 0; i < gt_array_size(v->exons); i++)
    {
      GtGenomeNode **exon = gt_array_get(v->exons, i);
      GtRange range = gt_genome_node_get_range(*exon);
      GtIntervalTreeNode *itn = gt_interval_tree_node_new(*exon, range.start,
                                                          range.end);
      gt_interval_tree_insert(v->exonsbyrange, itn);
    }
    if(gt_array_size(v->exons) == 0)
      infer_exons_visitor_visit_gene_infer_exons(v);

    v->intronsbyrange = gt_interval_tree_new(NULL);
    v->introns = agn_typecheck_select(current, agn_typecheck_intron);
    if(gt_array_size(v->introns) == 0 && gt_array_size(v->exons) > 1)
      infer_exons_visitor_visit_gene_infer_introns(v);

    gt_array_delete(v->exons);
    gt_array_delete(v->introns);
    gt_interval_tree_delete(v->exonsbyrange);
    gt_interval_tree_delete(v->intronsbyrange);
  }
  gt_feature_node_iterator_delete(iter);

  return 0;
}

static bool
infer_exons_visitor_visit_gene_collapse_feature(AgnInferExonsVisitor *v,
                                                GtFeatureNode *mrna,
                                                GtRange *range,
                                                GtIntervalTree *featsbyrange)
{
  GtArray *overlapping = gt_array_new( sizeof(GtFeatureNode *) );
  gt_interval_tree_find_all_overlapping(featsbyrange, range->start,
                                        range->end, overlapping);
  bool collapsed = false;
  while(gt_array_size(overlapping) > 0)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(overlapping);
    GtRange fnrange = gt_genome_node_get_range((GtGenomeNode *)fn);
    if(gt_range_compare(range, &fnrange) == 0)
    {
      gt_feature_node_add_child(mrna, fn);
      gt_genome_node_ref((GtGenomeNode *)fn);
      const char *parentattr = gt_feature_node_get_attribute(fn, "Parent");
      const char *tid = gt_feature_node_get_attribute(mrna, "ID");
      char parentstr[1024];
      strcpy(parentstr, parentattr);
      sprintf(parentstr + strlen(parentstr), ",%s", tid);
      gt_feature_node_set_attribute(fn, "Parent", parentstr);

      collapsed = true;
      break;
    }
  }
  gt_array_delete(overlapping);
  return collapsed;
}

static void
infer_exons_visitor_visit_gene_infer_exons(AgnInferExonsVisitor *v)
{
  GtFeatureNode *fn;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(v->gene);
  for(fn  = gt_feature_node_iterator_next(iter);
      fn != NULL;
      fn  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_typecheck_mrna(fn))
     continue;

    const char *mrnaid = gt_feature_node_get_attribute(fn, "ID");
    unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)fn);
    GtArray *cds  = agn_typecheck_select(fn, agn_typecheck_cds);
    GtArray *utrs = agn_typecheck_select(fn, agn_typecheck_utr);

    bool cds_explicit = gt_array_size(cds) > 0;
    if(!cds_explicit)
    {
      gt_logger_log(v->logger, "cannot infer missing exons for mRNA '%s' "
                    "(line %u) without CDS feature(s)", mrnaid, ln);
      continue;
    }

    GtUword i,j;
    GtHashmap *adjacent_utrs = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    GtArray *exons_to_add = gt_array_new( sizeof(GtRange) );
    for(i = 0; i < gt_array_size(cds); i++)
    {
      GtGenomeNode **cdssegment = gt_array_get(cds, i);
      GtRange crange = gt_genome_node_get_range(*cdssegment);
      GtRange erange = crange;
      for(j = 0; j < gt_array_size(utrs); j++)
      {
        GtGenomeNode **utrsegment = gt_array_get(utrs, j);
        GtRange urange = gt_genome_node_get_range(*utrsegment);

        // If the UTR segment is adjacent to the CDS, merge the ranges
        if(urange.end+1 == crange.start || crange.end+1 == urange.start)
        {
          erange = gt_range_join(&erange, &urange);
          gt_hashmap_add(adjacent_utrs, utrsegment, utrsegment);
        }
      }
      gt_array_add(exons_to_add, erange);
    }

    // Now create UTR-only exons
    for(i = 0; i < gt_array_size(utrs); i++)
    {
      GtGenomeNode **utrsegment = gt_array_get(utrs, i);
      GtRange urange = gt_genome_node_get_range(*utrsegment);
      if(gt_hashmap_get(adjacent_utrs, utrsegment) == NULL)
      {
        gt_array_add(exons_to_add, urange);
      }
    }

    for(i = 0; i < gt_array_size(exons_to_add); i++)
    {
      GtRange *erange = gt_array_get(exons_to_add, i);
      if(infer_exons_visitor_visit_gene_collapse_feature(v, fn, erange,
                                                         v->exonsbyrange))
      {
        continue;
      }

      GtGenomeNode **firstcds = gt_array_get(cds, 0);
      GtGenomeNode *exon = gt_feature_node_new
      (
        gt_genome_node_get_seqid(*firstcds), "exon", erange->start, erange->end,
        gt_feature_node_get_strand(*(GtFeatureNode **)firstcds)
      );
      GtFeatureNode *fn_exon = (GtFeatureNode *)exon;
      gt_feature_node_add_child(fn, fn_exon);
      if(mrnaid)
        gt_feature_node_add_attribute(fn_exon, "Parent", mrnaid);
      gt_array_add(v->exons, exon);
      GtIntervalTreeNode *node = gt_interval_tree_node_new(exon, erange->start,
                                                           erange->end);
      gt_interval_tree_insert(v->exonsbyrange, node);
    }
    gt_array_delete(exons_to_add);

    if(gt_array_size(v->exons) == 0)
    {
      gt_logger_log(v->logger, "unable to infer exons for mRNA '%s' (line %u)",
                    mrnaid, ln);
    }
    gt_array_delete(cds);
    gt_array_delete(utrs);
    gt_hashmap_delete(adjacent_utrs);
  }
  gt_feature_node_iterator_delete(iter);
}

static void
infer_exons_visitor_visit_gene_infer_introns(AgnInferExonsVisitor *v)
{
  GtFeatureNode *fn;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(v->gene);
  for(fn  = gt_feature_node_iterator_next(iter);
      fn != NULL;
      fn  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_typecheck_mrna(fn))
     continue;

    const char *mrnaid = gt_feature_node_get_attribute(fn, "ID");
    unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)fn);
    GtArray *exons = agn_typecheck_select(fn, agn_typecheck_exon);
    if(gt_array_size(exons) < 2)
    {
      gt_array_delete(exons);
      continue;
    }

    GtUword i;
    GtArray *introns_to_add = gt_array_new( sizeof(GtRange) );
    for(i = 1; i < gt_array_size(exons); i++)
    {
      GtGenomeNode **exon1 = gt_array_get(exons, i-1);
      GtGenomeNode **exon2 = gt_array_get(exons, i);
      GtRange first_range  = gt_genome_node_get_range(*exon1);
      GtRange second_range = gt_genome_node_get_range(*exon2);
      
      if(first_range.end == second_range.start - 1)
      {
        gt_logger_log(v->logger, "mRNA '%s' (line %u) has directly adjacent "
                      "exons", mrnaid, ln);
        return;
      }
      else
      {
        GtRange irange = { first_range.end + 1, second_range.start - 1 };
        gt_array_add(introns_to_add, irange);
      }
    }
    
    for(i = 0; i < gt_array_size(introns_to_add); i++)
    {
      GtRange *irange = gt_array_get(introns_to_add, i);
      if(infer_exons_visitor_visit_gene_collapse_feature(v, fn, irange,
                                                         v->intronsbyrange))
      {
        continue;
      }

      GtGenomeNode **firstexon = gt_array_get(exons, 0);
      GtGenomeNode *intron = gt_feature_node_new
      (
        gt_genome_node_get_seqid(*firstexon), "intron", irange->start,
        irange->end, gt_feature_node_get_strand(*(GtFeatureNode **)firstexon)
      );
      GtFeatureNode *fn_intron = (GtFeatureNode *)intron;
      gt_feature_node_add_child(fn, fn_intron);
      if(mrnaid)
        gt_feature_node_add_attribute(fn_intron, "Parent", mrnaid);
      gt_array_add(v->introns, fn_intron);
      GtIntervalTreeNode *node = gt_interval_tree_node_new(intron,irange->start,
                                                           irange->end);
      gt_interval_tree_insert(v->intronsbyrange, node);
    }
    gt_array_delete(introns_to_add);
    gt_array_delete(exons);
  }
  gt_feature_node_iterator_delete(iter);
}
