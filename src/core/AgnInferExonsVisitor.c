#include <string.h>
#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "AgnInferExonsVisitor.h"
#include "AgnGtExtensions.h"
#include "AgnTestData.h"
#include "AgnUtils.h"

#define agn_infer_exons_visitor_cast(GV)\
        gt_node_visitor_cast(agn_infer_exons_visitor_class(), GV)

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
  AgnLogger *logger;
};


//----------------------------------------------------------------------------//
// Prototypes of private functions
//----------------------------------------------------------------------------//

/**
 * @function Cast a node visitor object as a AgnInferExonsVisitor.
 */
static const GtNodeVisitorClass* agn_infer_exons_visitor_class();

/**
 * @function Run unit tests using the basic grape example data.
 */
static bool unit_test_grape(AgnUnitTest *test);

/**
 * @function Run unit tests using the grape example data with CDSs and UTRs (no
 * exons explicitly defined).
 */
static bool unit_test_grape_sansexons(AgnUnitTest *test);

/**
 * @function Procedure for processing feature nodes (the only node of interest
 * for this node visitor).
 */
static int visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                              GtError *error);

/**
 * @function If an exon with the same coordinates already exists and belongs to
 * another mRNA, associate it with this mRNA as well instead of creating a
 * duplicate feature.
 */
static bool visit_gene_collapse_feature(AgnInferExonsVisitor *v,
                                        GtFeatureNode *mrna, GtRange *range,
                                        GtIntervalTree *featsbyrange);

/**
 * @function Infer exons from CDS and UTR segments if possible.
 */
static void visit_gene_infer_exons(AgnInferExonsVisitor *v);

/**
 * @function Infer introns from exons if necessary.
 */
static void visit_gene_infer_introns(AgnInferExonsVisitor *v);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

static const GtNodeVisitorClass* agn_infer_exons_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInferExonsVisitor),
                                    NULL,
                                    NULL,
                                    visit_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* agn_infer_exons_visitor_new(AgnLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(agn_infer_exons_visitor_class());
  AgnInferExonsVisitor *v = agn_infer_exons_visitor_cast(nv);
  v->logger = logger;
  return nv;
}

bool agn_infer_exons_visitor_unit_test(AgnUnitTest *test)
{
  bool grape = unit_test_grape(test);
  bool grape_sansexons = unit_test_grape_sansexons(test);
  return grape && grape_sansexons;
}

static bool unit_test_grape(AgnUnitTest *test)
{
  GtArray *genes = agn_test_data_grape();
  GtError *error = gt_error_new();
  AgnLogger *logger = agn_logger_new();
  GtNodeStream *genestream = gt_array_in_stream_new(genes, NULL, error);
  GtNodeVisitor *iev = agn_infer_exons_visitor_new(logger);
  GtNodeStream *ievstream = gt_visitor_stream_new(genestream, iev);
  int result;
  GtGenomeNode *gn;
  GtFeatureNode *fn;

  result = gt_node_stream_next(ievstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  GtArray *exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
  bool exons1correct = gt_array_size(exons) == 3;
  agn_unit_test_result(test, "grape: exon check 1", exons1correct);
  gt_array_delete(exons);
  gt_genome_node_delete(gn);

  result = gt_node_stream_next(ievstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
  bool exons2correct = gt_array_size(exons) == 3;
  agn_unit_test_result(test, "grape: exon check 2", exons2correct);
  gt_array_delete(exons);
  gt_genome_node_delete(gn);

  result = gt_node_stream_next(ievstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
  bool exons3correct = gt_array_size(exons) == 6;
  agn_unit_test_result(test, "grape: exon check 3", exons3correct);
  gt_array_delete(exons);
  gt_genome_node_delete(gn);
  
  return exons1correct && exons2correct && exons3correct;
}

static bool unit_test_grape_sansexons(AgnUnitTest *test)
{
  GtArray *genes = agn_test_data_grape_sansexons();
  GtError *error = gt_error_new();
  AgnLogger *logger = agn_logger_new();
  GtNodeStream *genestream = gt_array_in_stream_new(genes, NULL, error);
  GtNodeVisitor *iev = agn_infer_exons_visitor_new(logger);
  GtNodeStream *ievstream = gt_visitor_stream_new(genestream, iev);
  int result;
  GtGenomeNode *gn;
  GtFeatureNode *fn;

  result = gt_node_stream_next(ievstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  GtArray *exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
  bool exons1correct = gt_array_size(exons) == 3;
  if(exons1correct)
  {
    GtGenomeNode **exon1 = gt_array_get(exons, 0);
    GtGenomeNode **exon2 = gt_array_get(exons, 1);
    GtGenomeNode **exon3 = gt_array_get(exons, 2);
    GtRange range1 = gt_genome_node_get_range(*exon1);
    GtRange range2 = gt_genome_node_get_range(*exon2);
    GtRange range3 = gt_genome_node_get_range(*exon3);
    exons1correct = (range1.start == 22057 && range1.end == 22382 &&
                     range2.start == 22497 && range2.end == 22550 &&
                     range3.start == 22651 && range3.end == 23119);
  }
  agn_unit_test_result(test, "grape::sansexons: exon check 1", exons1correct);
  gt_array_delete(exons);
  gt_genome_node_delete(gn);

  result = gt_node_stream_next(ievstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
  bool exons2correct = gt_array_size(exons) == 3;
  if(exons2correct)
  {
    GtGenomeNode **exon1 = gt_array_get(exons, 0);
    GtGenomeNode **exon2 = gt_array_get(exons, 1);
    GtGenomeNode **exon3 = gt_array_get(exons, 2);
    GtRange range1 = gt_genome_node_get_range(*exon1);
    GtRange range2 = gt_genome_node_get_range(*exon2);
    GtRange range3 = gt_genome_node_get_range(*exon3);
    exons2correct = (range1.start == 48012 && range1.end == 48537 &&
                     range2.start == 48637 && range2.end == 48766 &&
                     range3.start == 48870 && range3.end == 48984);
  }
  agn_unit_test_result(test, "grape::sansexons: exon check 2", exons2correct);
  gt_array_delete(exons);
  gt_genome_node_delete(gn);

  result = gt_node_stream_next(ievstream, &gn, error);
  if(result == -1)
  {
    fprintf(stderr, "node stream error: %s\n", gt_error_get(error));
    return false;
  }
  fn = (GtFeatureNode *)gn;
  exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
  bool exons3correct = gt_array_size(exons) == 6;
  if(exons3correct)
  {
    GtGenomeNode **exon1 = gt_array_get(exons, 0);
    GtGenomeNode **exon2 = gt_array_get(exons, 1);
    GtGenomeNode **exon3 = gt_array_get(exons, 2);
    GtGenomeNode **exon4 = gt_array_get(exons, 3);
    GtGenomeNode **exon5 = gt_array_get(exons, 4);
    GtGenomeNode **exon6 = gt_array_get(exons, 5);
    GtRange range1 = gt_genome_node_get_range(*exon1);
    GtRange range2 = gt_genome_node_get_range(*exon2);
    GtRange range3 = gt_genome_node_get_range(*exon3);
    GtRange range4 = gt_genome_node_get_range(*exon4);
    GtRange range5 = gt_genome_node_get_range(*exon5);
    GtRange range6 = gt_genome_node_get_range(*exon6);
    exons3correct = (range1.start == 88551 && range1.end == 89029 &&
                     range2.start == 89265 && range2.end == 89549 &&
                     range3.start == 90074 && range3.end == 90413 &&
                     range4.start == 90728 && range4.end == 90833 &&
                     range5.start == 91150 && range5.end == 91362 &&
                     range6.start == 91810 && range6.end == 92176);
  }
  agn_unit_test_result(test, "grape::sansexons: exon check 3", exons3correct);
  gt_array_delete(exons);
  gt_genome_node_delete(gn);
  
  return exons1correct && exons2correct && exons3correct;
}

static int visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                              GtError *error)
{
  AgnInferExonsVisitor *v = agn_infer_exons_visitor_cast(nv);
  gt_error_check(error);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    if(!agn_gt_feature_node_is_gene_feature(current))
      continue;
    
    GtUword i;
    v->gene = current;

    v->exonsbyrange = gt_interval_tree_new(NULL);
    v->exons = agn_gt_feature_node_children_of_type(current,
                                           agn_gt_feature_node_is_exon_feature);
    for(i = 0; i < gt_array_size(v->exons); i++)
    {
      GtGenomeNode **exon = gt_array_get(v->exons, i);
      GtRange range = gt_genome_node_get_range(*exon);
      GtIntervalTreeNode *itn = gt_interval_tree_node_new(*exon, range.start,
                                                          range.end);
      gt_interval_tree_insert(v->exonsbyrange, itn);
    }
    if(gt_array_size(v->exons) == 0)
      visit_gene_infer_exons(v);

    v->intronsbyrange = gt_interval_tree_new(NULL);
    v->introns = agn_gt_feature_node_children_of_type(current,
                                         agn_gt_feature_node_is_intron_feature);
    if(gt_array_size(v->introns) == 0 && gt_array_size(v->exons) > 1)
      visit_gene_infer_introns(v);

    gt_array_delete(v->exons);
    gt_array_delete(v->introns);
    gt_interval_tree_delete(v->exonsbyrange);
    gt_interval_tree_delete(v->intronsbyrange);
  }
  gt_feature_node_iterator_delete(iter);

  return 0;
}

static bool visit_gene_collapse_feature(AgnInferExonsVisitor *v,
                                        GtFeatureNode *mrna, GtRange *range,
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

static void visit_gene_infer_exons(AgnInferExonsVisitor *v)
{
  GtFeatureNode *fn;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(v->gene);
  for(fn  = gt_feature_node_iterator_next(iter);
      fn != NULL;
      fn  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_gt_feature_node_is_mrna_feature(fn))
     continue;

    const char *mrnaid = gt_feature_node_get_attribute(fn, "ID");
    unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)fn);
    GtArray *cds  = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_cds_feature);
    GtArray *utrs = agn_gt_feature_node_children_of_type(fn,
                                            agn_gt_feature_node_is_utr_feature);

    bool cds_explicit = gt_array_size(cds) > 0;
    if(!cds_explicit)
    {
      agn_logger_log_error(v->logger, "cannot infer missing exons for mRNA "
                           "'%s' (line %u) without CDS feature(s)", mrnaid, ln);
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
      if(visit_gene_collapse_feature(v, fn, erange, v->exonsbyrange))
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
      agn_logger_log_error(v->logger, "unable to infer exons for mRNA '%s'"
                           "(line %u)", mrnaid, ln);
    }
    gt_array_delete(cds);
    gt_array_delete(utrs);
    gt_hashmap_delete(adjacent_utrs);
  }
  gt_feature_node_iterator_delete(iter);
}

void visit_gene_infer_introns(AgnInferExonsVisitor *v)
{
  GtFeatureNode *fn;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(v->gene);
  for(fn  = gt_feature_node_iterator_next(iter);
      fn != NULL;
      fn  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_gt_feature_node_is_mrna_feature(fn))
     continue;

    const char *mrnaid = gt_feature_node_get_attribute(fn, "ID");
    unsigned int ln = gt_genome_node_get_line_number((GtGenomeNode *)fn);
    GtArray *exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
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
        agn_logger_log_error(v->logger, "mRNA '%s' has directly adjacent exons",
                             mrnaid);
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
      if(visit_gene_collapse_feature(v, fn, irange, v->intronsbyrange))
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

    if(gt_array_size(v->introns) == 0)
    {
      agn_logger_log_error(v->logger, "unable to infer introns for mRNA '%s'"
                           "(line %u)", mrnaid, ln);
    }
    gt_array_delete(exons);
  }
  gt_feature_node_iterator_delete(iter);
}
