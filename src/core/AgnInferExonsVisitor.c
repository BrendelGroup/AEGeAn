#include <string.h>
#include "extended/feature_index_memory_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "AgnInferExonsVisitor.h"
#include "AgnGtExtensions.h"
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
  GtArray *exons;
  AgnLogger *logger;
};


//----------------------------------------------------------------------------//
// Prototypes of private functions
//----------------------------------------------------------------------------//

/**
 * Cast a node visitor object as a AgnInferExonsVisitor
 *
 * @returns    a node visitor object cast as a AgnInferExonsVisitor
 */
static const GtNodeVisitorClass* agn_infer_exons_visitor_class();

/**
 * Procedure for processing feature nodes (the only node of interest for this
 * node visitor).
 *
 * @param[in]  nv       a node visitor
 * @param[in]  fn       node representing a top-level GFF3 feature entry
 * @param[out] error    stores error messages
 * @returns             0 for success
 */
static int visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                              GtError *error);

/**
 * If an exon with the same coordinates already exists and belongs to another
 * mRNA, associate it with this mRNA as well instead of creating a duplicate
 * feature.
 *
 * @param[in] v        visitor object
 * @param[in] mrna     mRNA for which an exon is to be added
 * @param[in] range    range of the exon to (potentially) be created
 * @returns            true if the exon already exists and has been associated
 *                     with the mRNA, false if the exon feature needs to be
 *                     created
 */
static bool visit_gene_collapse_feature(AgnInferExonsVisitor *v,
                                        GtFeatureNode *mrna, GtRange *range);

/**
 * Infer exons from CDS and UTR segments if possible.
 *
 * @param[in] v    visitor object
 */
static void visit_gene_infer_exons(AgnInferExonsVisitor *v);


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

static int visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                              GtError *error)
{
  AgnInferExonsVisitor *v = agn_infer_exons_visitor_cast(nv);
  gt_error_check(error);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(!agn_gt_feature_node_is_gene_feature(current))
      continue;
    
    GtUword i;
    v->exonsbyrange = gt_interval_tree_new(NULL);
    v->exons   = agn_gt_feature_node_children_of_type(current,
                                         agn_gt_feature_node_is_exon_feature);
    v->gene = current;
    for(i = 0; i < gt_array_size(v->exons); i++)
    {
      GtGenomeNode **exon = gt_array_get(v->exons, i);
      GtRange range = gt_genome_node_get_range(*exon);
      GtIntervalTreeNode *itn = gt_interval_tree_node_new(*exon, range.start,
                                                          range.end);
      gt_interval_tree_insert(v->exonsbyrange, itn);
    }

    visit_gene_infer_exons(v);

    gt_array_delete(v->exons);
    gt_interval_tree_delete(v->exonsbyrange);
    v->gene = NULL;
  }
  gt_feature_node_iterator_delete(iter);

  return 0;
}

static bool visit_gene_collapse_feature(AgnInferExonsVisitor *v,
                                        GtFeatureNode *mrna, GtRange *range)
{
  GtArray *overlapping = gt_array_new( sizeof(GtFeatureNode *) );
  gt_interval_tree_find_all_overlapping(v->exonsbyrange, range->start,
                                        range->end, overlapping);
  bool collapsed = false;
  while(gt_array_size(overlapping) > 0)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(overlapping);
    GtRange fnrange = gt_genome_node_get_range((GtGenomeNode *)fn);
    if(gt_range_compare(range, &fnrange) == 0)
    {
      gt_feature_node_add_child(mrna, fn);
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
    GtArray *exons = agn_gt_feature_node_children_of_type(fn,
                                           agn_gt_feature_node_is_exon_feature);
    bool exons_explicit = gt_array_size(exons) > 0;
    gt_array_delete(exons);
    bool cds_explicit = gt_array_size(cds) > 0;
    if(exons_explicit)
      continue;

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
      if(visit_gene_collapse_feature(v, fn, erange))
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
  // FIXME: poor reference handling somewhere; if I delete this iterator, I get
  //        a segfault later on; if I don't delete it, I leak memory
  // gt_feature_node_iterator_delete(iter);
}
