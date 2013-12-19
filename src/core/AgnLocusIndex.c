#include "AgnLocusIndex.h"
#include "AgnGeneLocus.h"
#include "AgnTestData.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnLocusIndex
{
  GtStrArray *seqids;
  GtHashmap *seqranges;
  GtHashmap *locus_trees;
  GtFree locusfreefunc;
};


//------------------------------------------------------------------------------
// Prototypes for private methods
//------------------------------------------------------------------------------

/**
 * Function used to traverse the given interval tree and add all its loci to the
 * given array.
 *
 * @param[in]  itn    node in the interval tree
 * @param[out] lp     array to which loci are added
 * @returns           0 for success, 1 otherwise
 */
static int agn_locus_index_it_traverse(GtIntervalTreeNode *itn, void *lp);

/**
 * Given two sets of annotations for the same sequence (a reference set and a
 * prediction set), this function associates each gene annotation with the
 * appropriate locus.
 *
 * @param[in] idx        the locus index
 * @param[in] seqid      the sequence for which loci are requested
 * @param[in] refr       reference annotations for the sequence
 * @param[in] pred       prediction annotations for the sequence
 * @param[in] filters    filtering criteria
 * @param[in] logger     object to which warning/error messages are written
 * @returns              an interval tree containing the sequence's gene loci
 */
static GtIntervalTree *agn_locus_index_parse_pairwise(AgnLocusIndex *idx,
                                                     const char *seqid,
                                                     GtFeatureIndex *refr,
                                                     GtFeatureIndex *pred,
                                                     AgnCompareFilters *filters,
                                                     AgnLogger *logger);

/**
 * Given a locus and a collection of gene features, search for genes that
 * overlap with the locus and assign them to the locus.
 *
 * @param[in]  idx              locus index object
 * @param[in]  features         collection of gene features
 * @param[in]  visited_genes    list of genes that have already been assigned to
 *                              a locus and which should be ignored
 * @param[out] locus            the locus of interest, to which overlapping
 *                              genes will be assigned
 * @param[in]  add_func         the function that should be used to add genes to
 *                              the locus (distinguishes between reference and
 *                              prediction)
 * @param[in]  logger           object to which warning/error messages are
 *                              written
 * @returns                     the number of genes assigned to the locus
 */
static int agn_locus_index_pairwise_test_overlap(AgnLocusIndex *idx,
                                                 GtFeatureIndex *features,
                                                 GtHashmap *visited_genes,
                                                 AgnGeneLocus *locus,
                                                 AgnComparisonSource source,
                                                 AgnLogger *logger);

/**
 * Given a locus and a collection of genomic features, search for genes that
 * overlap with the locus and assign them to the locus.
 *
 * @param[in]  idx              locus index object
 * @param[in]  features         collection of genomic features
 * @param[in]  visited_genes    list of genes that have already been assigned to
 *                              a locus and which should be ignored
 * @param[out] locus            the locus of interest, to which overlapping
 *                              genes will be assigned
 * @param[in]  logger           object to which warning/error messages are
 *                              written
 * @returns                     the number of genes assigned to the locus
 */
static int agn_locus_index_test_overlap(AgnLocusIndex *idx,
                                        GtFeatureIndex *features,
                                        GtHashmap *visited_genes,
                                        AgnGeneLocus *locus,
                                        AgnLogger *logger);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_locus_index_comparative_analysis(AgnLocusIndex *idx, const char *seqid,
                                          AgnLocusIndexVisitFunc preanalyfunc,
                                          AgnLocusIndexVisitFunc postanalyfunc,
                                          void *analyfuncdata,
                                          AgnLogger *logger)
{
  GtArray *seqloci = agn_locus_index_get(idx, seqid);
  GtUword nloci = gt_array_size(seqloci);

  int i;
  for(i = 0; i < nloci; i++)
  {
    GtTimer *timer = gt_timer_new();
    gt_timer_start(timer);
    AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(seqloci, i);
    preanalyfunc(locus, analyfuncdata);
    agn_gene_locus_comparative_analysis(locus);
    postanalyfunc(locus, analyfuncdata);
    gt_timer_stop(timer);
    //agn_logger_log_status(logger, "proc=%d; locus=%s[%lu, %lu]; length=%lu; "
    //    "trans=%lu,%lu; pairs=%lu,%lu", rank, agn_gene_locus_get_seqid(locus),
    //    agn_gene_locus_get_start(locus), agn_gene_locus_get_end(locus),
    //    agn_gene_locus_get_length(locus),
    //    agn_gene_locus_num_refr_transcripts(locus),
    //    agn_gene_locus_num_pred_transcripts(locus), npairs,
    //    gt_array_size(reportedpairs));
    gt_timer_delete(timer);
  }
}

void agn_locus_index_delete(AgnLocusIndex *idx)
{
  gt_hashmap_delete(idx->locus_trees);
  gt_str_array_delete(idx->seqids);
  gt_hashmap_delete(idx->seqranges);
  gt_free(idx);
  idx = NULL;
}

void agn_locus_index_find(AgnLocusIndex *idx, const char *seqid, GtRange *range,
                          GtArray *loci)
{
  gt_assert(loci != NULL);
  GtIntervalTree *it = gt_hashmap_get(idx->locus_trees, seqid);
  if(it == NULL)
    return;

  gt_interval_tree_find_all_overlapping(it, range->start, range->end, loci);
}

GtArray *agn_locus_index_get(AgnLocusIndex *idx, const char *seqid)
{
  GtIntervalTree *it = gt_hashmap_get(idx->locus_trees, seqid);
  if(it == NULL)
    return NULL;

  GtArray *loci = gt_array_new( sizeof(AgnGeneLocus *) );
  gt_interval_tree_traverse(it, agn_locus_index_it_traverse, loci);
  return loci;
}

GtArray *agn_locus_index_interval_loci(AgnLocusIndex *idx, const char *seqid,
                                       GtUword delta, bool skipterminal)
{
  GtArray *loci = agn_locus_index_get(idx, seqid);
  if(loci == NULL)
    return NULL;
  gt_array_sort(loci, (GtCompare)agn_gene_locus_array_compare);

  GtUword nloci = gt_array_size(loci);
  GtArray *iloci = gt_array_new( sizeof(AgnGeneLocus *) );
  GtRange *seqrange = gt_hashmap_get(idx->seqranges, seqid);

  GtUword i;
  GtArray *clonedloci = gt_array_new( sizeof(AgnGeneLocus *) );
  for(i = 0; i < nloci; i++)
  {
    AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(loci, i);
    AgnGeneLocus *clone = agn_gene_locus_clone(locus);
    gt_array_add(clonedloci, clone);
  }
  gt_array_delete(loci);
  loci = clonedloci;
  
  // Handle trivial case
  if(nloci == 0)
  {
    AgnGeneLocus *ilocus = agn_gene_locus_new(seqid);
    agn_gene_locus_set_range(ilocus, seqrange->start, seqrange->end);
    gt_array_add(iloci, ilocus);
    gt_array_delete(loci);
    return iloci;
  }

  // Handle initial ilocus
  AgnGeneLocus *l1 = *(AgnGeneLocus **)gt_array_get(loci, 0);
  GtRange r1 = agn_gene_locus_range(l1);
  if(r1.start >= seqrange->start + (2*delta))
  {
    if(nloci == 1)
      agn_gene_locus_set_range(l1, r1.start - delta, r1.end);

    if(!skipterminal)
    {
      AgnGeneLocus *ilocus = agn_gene_locus_new(seqid);
      agn_gene_locus_set_range(ilocus, seqrange->start, r1.start - delta - 1);
      gt_array_add(iloci, ilocus);
    }
  }
  else
  {
    agn_gene_locus_set_range(l1, seqrange->start, r1.end);
  }

  // Handle internal iloci
  for(i = 0; i < nloci - 1; i++)
  {
    AgnGeneLocus *llocus = *(AgnGeneLocus **)gt_array_get(loci, i);
    AgnGeneLocus *rlocus = *(AgnGeneLocus **)gt_array_get(loci, i+1);
    GtRange lrange = agn_gene_locus_range(llocus);
    GtRange rrange = agn_gene_locus_range(rlocus);

    gt_array_add(iloci, llocus);
    if(lrange.end + delta >= rrange.start)
    {
      agn_gene_locus_set_range(llocus, lrange.start, rrange.start - 1);
      agn_gene_locus_set_range(rlocus, lrange.end + 1, rrange.end);
    }
    else if(lrange.end + (2*delta) >= rrange.start)
    {
      agn_gene_locus_set_range(llocus, lrange.start, lrange.end + delta);
      agn_gene_locus_set_range(rlocus, rrange.start - delta, rrange.end);
    }
    else if(lrange.end + (3*delta) >= rrange.start)
    {
      GtUword midpoint = (lrange.end + rrange.start) / 2;
      agn_gene_locus_set_range(llocus, lrange.start, midpoint);
      agn_gene_locus_set_range(rlocus, midpoint + 1, rrange.end);
    }
    else
    {
      agn_gene_locus_set_range(llocus, lrange.start, lrange.end + delta);
      agn_gene_locus_set_range(rlocus, rrange.start - delta, rrange.end);
      AgnGeneLocus *ilocus = agn_gene_locus_new(seqid);
      agn_gene_locus_set_range(ilocus, lrange.end + delta + 1,
                               rrange.start - delta - 1);
      gt_array_add(iloci, ilocus);
    }

    if(i == nloci - 2)
      gt_array_add(iloci, rlocus);
  }

  // Handle terminal ilocus
  l1 = *(AgnGeneLocus **)gt_array_get(loci, nloci - 1);
  r1 = agn_gene_locus_range(l1);
  if(r1.end <= seqrange->end - (2*delta))
  {
    agn_gene_locus_set_range(l1, r1.start, r1.end + delta);

    if(!skipterminal)
    {
      AgnGeneLocus *ilocus = agn_gene_locus_new(seqid);
      agn_gene_locus_set_range(ilocus, r1.end + delta + 1, seqrange->end);
      gt_array_add(iloci, ilocus);
    }
  }
  else
  {
    agn_gene_locus_set_range(l1, r1.start, seqrange->end);
  }
  if(nloci == 1)
    gt_array_add(iloci, l1);

  gt_array_delete(loci);
  gt_array_sort(iloci, (GtCompare)agn_gene_locus_array_compare);
  return iloci;
}

static int agn_locus_index_it_traverse(GtIntervalTreeNode *itn, void *lp)
{
  GtArray *loci = (GtArray *)lp;
  AgnGeneLocus *locus = gt_interval_tree_node_get_data(itn);
  gt_array_add(loci, locus);
  return 0;
}

AgnLocusIndex *agn_locus_index_new(bool freeondelete)
{
  AgnLocusIndex *idx = gt_malloc( sizeof(AgnLocusIndex) );
  idx->seqids = gt_str_array_new();
  idx->seqranges = gt_hashmap_new(GT_HASH_STRING, NULL, (GtFree)gt_free_func);
  idx->locus_trees = gt_hashmap_new(GT_HASH_STRING, NULL,
                                    (GtFree)gt_interval_tree_delete);
  idx->locusfreefunc = NULL;
  if(freeondelete)
    idx->locusfreefunc = (GtFree)agn_gene_locus_delete;

  return idx;
}

GtUword agn_locus_index_parse_disk(AgnLocusIndex * idx, int numfiles,
                                   const char **filenames, AgnLogger *logger)
{
  gt_assert(idx != NULL);
  GtUword nloci;
  GtFeatureIndex *features = agn_import_simple(numfiles, filenames, "gene",
                                               logger);
  if(agn_logger_has_error(logger))
  {
    gt_feature_index_delete(features);
    return 0;
  }

  nloci = agn_locus_index_parse_memory(idx, features, logger);
  gt_feature_index_delete(features);
  return nloci;
}

static int agn_locus_index_pairwise_test_overlap(AgnLocusIndex *idx,
                                                 GtFeatureIndex *features,
                                                 GtHashmap *visited_genes,
                                                 AgnGeneLocus *locus,
                                                 AgnComparisonSource source,
                                                 AgnLogger *logger)
{
  GtError *error = gt_error_new();
  bool has_seqid;
  gt_feature_index_has_seqid(features, &has_seqid,
                             agn_gene_locus_get_seqid(locus),
                             error);
  if(!has_seqid)
  {
    gt_error_delete(error);
    return 0;
  }

  GtRange locusrange;
  GtArray *genes_to_add = gt_array_new( sizeof(GtFeatureNode *) );
  GtUword new_gene_count = 0;

  locusrange.start = agn_gene_locus_get_start(locus);
  locusrange.end = agn_gene_locus_get_end(locus);
  gt_feature_index_get_features_for_range(features, genes_to_add,
          agn_gene_locus_get_seqid(locus), &locusrange, error);
  if(gt_error_is_set(error))
  {
    agn_logger_log_error(logger, "error fetching features for range %s[%lu, "
                         "%lu]: %s",agn_gene_locus_get_seqid(locus),
                         locusrange.start, locusrange.end, gt_error_get(error));
    gt_error_unset(error);
  }

  while(gt_array_size(genes_to_add) > 0)
  {
    GtFeatureNode *gene_to_add = *(GtFeatureNode **)gt_array_pop(genes_to_add);
    if(gt_hashmap_get(visited_genes, gene_to_add) == NULL)
    {
      gt_hashmap_add(visited_genes, gene_to_add, gene_to_add);
      agn_gene_locus_add(locus, gene_to_add, source);
      new_gene_count++;
    }
  }

  gt_array_delete(genes_to_add);
  gt_error_delete(error);

  return new_gene_count;
}

GtIntervalTree *agn_locus_index_parse(AgnLocusIndex *idx, const char *seqid,
                                      GtFeatureIndex *features,
                                      AgnLogger *logger)
{
  GtUword i;
  GtError *error = gt_error_new();
  GtHashmap *visited_genes = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  GtIntervalTree *loci = gt_interval_tree_new(idx->locusfreefunc);

  GtArray *seqfeatures = gt_feature_index_get_features_for_seqid(features,
                                 seqid, error);
  if(gt_error_is_set(error))
  {
    agn_logger_log_error(logger, "error fetching features for sequence '%s': "
                         "%s", seqid, gt_error_get(error));
    gt_error_delete(error);
    gt_array_delete(seqfeatures);
    gt_hashmap_delete(visited_genes);
    gt_interval_tree_delete(loci);
    return NULL;
  }
  for(i = 0; i < gt_array_size(seqfeatures); i++)
  {
    GtFeatureNode **feature = gt_array_get(seqfeatures, i);
    if(gt_hashmap_get(visited_genes, *feature) != NULL)
      continue;

    gt_hashmap_add(visited_genes, *feature, *feature);
    AgnGeneLocus *locus = agn_gene_locus_new(seqid);
    agn_gene_locus_add_gene(locus, *feature);

    int new_gene_count = 0;
    do
    {
      int temp_new_gene_count = agn_locus_index_test_overlap(idx,
                                    features, visited_genes, locus,
                                    logger);
      new_gene_count = temp_new_gene_count;
    } while(new_gene_count > 0);

    GtIntervalTreeNode *itn = gt_interval_tree_node_new(locus,
                                            agn_gene_locus_get_start(locus),
                                            agn_gene_locus_get_end(locus));
    gt_interval_tree_insert(loci, itn);
  }

  gt_error_delete(error);
  gt_hashmap_delete(visited_genes);
  gt_array_delete(seqfeatures);
  return loci;
}

static GtIntervalTree *agn_locus_index_parse_pairwise(AgnLocusIndex *idx,
                                                     const char *seqid,
                                                     GtFeatureIndex *refr,
                                                     GtFeatureIndex *pred,
                                                     AgnCompareFilters *filters,
                                                     AgnLogger *logger)
{
  GtUword i;
  GtError *error = gt_error_new();
  GtHashmap *visited_genes = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  GtIntervalTree *loci = gt_interval_tree_new(idx->locusfreefunc);

  // Seed new loci with reference genes
  GtArray *refr_list = gt_feature_index_get_features_for_seqid(refr, seqid,
                                                               error);
  if(gt_error_is_set(error))
  {
    agn_logger_log_error(logger, "error fetching reference features for "
                         "sequence '%s': %s", seqid, gt_error_get(error));
    gt_error_delete(error);
    gt_array_delete(refr_list);
    gt_hashmap_delete(visited_genes);
    gt_interval_tree_delete(loci);
    return NULL;
  }
  for(i = 0; i < gt_array_size(refr_list); i++)
  {
    GtFeatureNode *refr_gene = *(GtFeatureNode**)gt_array_get(refr_list, i);
    if(gt_hashmap_get(visited_genes, refr_gene) != NULL)
      continue;

    gt_hashmap_add(visited_genes, refr_gene, refr_gene);
    AgnGeneLocus *locus = agn_gene_locus_new(seqid);
    agn_gene_locus_add_refr_gene(locus, refr_gene);

    int new_gene_count = 0;
    do
    {
      int new_refr_gene_count = agn_locus_index_pairwise_test_overlap(idx, refr,
                                    visited_genes, locus, REFERENCESOURCE,
                                    logger);
      int new_pred_gene_count = agn_locus_index_pairwise_test_overlap(idx, pred,
                                    visited_genes, locus, PREDICTIONSOURCE,
                                    logger);
      if(agn_logger_has_error(logger))
      {
        gt_error_delete(error);
        gt_array_delete(refr_list);
        gt_hashmap_delete(visited_genes);
        gt_interval_tree_delete(loci);
        return NULL;
      }
      new_gene_count = new_refr_gene_count + new_pred_gene_count;
    } while(new_gene_count > 0);

    bool failfilter = agn_gene_locus_filter(locus, filters);
    if(!failfilter)
    {
      GtIntervalTreeNode *itn = gt_interval_tree_node_new(locus,
                                    agn_gene_locus_get_start(locus),
                                    agn_gene_locus_get_end(locus));
      gt_interval_tree_insert(loci, itn);
    }
    else
    {
      agn_logger_log_status(logger, "locus %s[%lu, %lu] did not pass filtering "
                            "criteria; moving on", seqid,
                            agn_gene_locus_get_start(locus),
                            agn_gene_locus_get_end(locus));
      agn_gene_locus_delete(locus);
    }
  }
  gt_array_delete(refr_list);

  // All reference genes, and some prediction genes, have been assigned to loci.
  // Now seed new loci with prediction genes, ignoring reference genes.
  GtArray *pred_list = gt_feature_index_get_features_for_seqid(pred, seqid,
                                                               error);
  if(gt_error_is_set(error))
  {
    agn_logger_log_error(logger, "error fetching prediction features for "
                         "sequence '%s': %s", seqid, gt_error_get(error));
    gt_error_delete(error);
    gt_array_delete(refr_list);
    gt_hashmap_delete(visited_genes);
    gt_interval_tree_delete(loci);
    return NULL;
  }
  for(i = 0; i < gt_array_size(pred_list); i++)
  {
    GtFeatureNode *pred_gene = *(GtFeatureNode**)gt_array_get(pred_list,i);
    if(gt_hashmap_get(visited_genes, pred_gene) != NULL)
      continue;

    gt_hashmap_add(visited_genes, pred_gene, pred_gene);
    AgnGeneLocus *locus = agn_gene_locus_new(seqid);
    agn_gene_locus_add_pred_gene(locus, pred_gene);

    int new_gene_count = 0;
    do
    {
      int new_pred_gene_count = agn_locus_index_pairwise_test_overlap(idx, pred,
                                    visited_genes, locus, PREDICTIONSOURCE,
                                    logger);
      if(agn_logger_has_error(logger))
      {
        gt_error_delete(error);
        gt_array_delete(refr_list);
        gt_hashmap_delete(visited_genes);
        gt_interval_tree_delete(loci);
        return NULL;
      }
      new_gene_count = new_pred_gene_count;
    } while(new_gene_count > 0);

    bool failfilter = agn_gene_locus_filter(locus, filters);
    if(!failfilter)
    {
      GtIntervalTreeNode *itn = gt_interval_tree_node_new(locus,
                                    agn_gene_locus_get_start(locus),
                                    agn_gene_locus_get_end(locus));
      gt_interval_tree_insert(loci, itn);
    }
    else
    {
      agn_logger_log_status(logger, "locus %s[%lu, %lu] did not pass filtering "
                            "criteria; moving on", seqid,
                            agn_gene_locus_get_start(locus),
                            agn_gene_locus_get_end(locus));
      agn_gene_locus_delete(locus);
    }
  }

  gt_error_delete(error);
  gt_hashmap_delete(visited_genes);
  gt_array_delete(pred_list);

  return loci;
}

GtUword agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx,
                                              GtFeatureIndex *refrfeats,
                                              GtFeatureIndex *predfeats,
                                              AgnCompareFilters *filters,
                                              AgnLogger *logger)
{
  int i;
  gt_assert(idx != NULL);
  gt_assert(refrfeats != NULL && predfeats != NULL);
  GtStrArray *seqids = agn_seq_union(refrfeats, predfeats, logger);
  if(agn_logger_has_error(logger))
  {
    gt_str_array_delete(seqids);
    return 0;
  }
  gt_str_array_delete(idx->seqids);
  idx->seqids = seqids;

  GtUword totalloci = 0;

  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtIntervalTree *loci = agn_locus_index_parse_pairwise(idx, seqid,
                                                          refrfeats,
                                                          predfeats, filters,
                                                          logger);
    GtRange *seqrange = gt_malloc( sizeof(GtRange) );
    GtRange refrrange, predrange;
    GtError *error = gt_error_new();
    bool refrhasseq, predhasseq;
    gt_feature_index_has_seqid(refrfeats, &refrhasseq, seqid, error);
    gt_feature_index_has_seqid(predfeats, &predhasseq, seqid, error);
    if(refrhasseq && predhasseq)
    {
      gt_feature_index_get_orig_range_for_seqid(refrfeats, &refrrange, seqid,
                                                error);
    }
    if(predhasseq)
    {
      gt_feature_index_get_orig_range_for_seqid(predfeats, &predrange, seqid,
                                                error);
    }
    GtRange trange;
    if(refrhasseq && predhasseq) gt_range_join(&refrrange, &predrange);
    else if(refrhasseq) trange = refrrange;
    else if(predhasseq) trange = predrange;
    seqrange->start = trange.start;
    seqrange->end   = trange.end;

    totalloci += gt_interval_tree_size(loci);
    gt_hashmap_add(idx->locus_trees, (char *)seqid, loci);
    gt_hashmap_add(idx->seqranges, (char *)seqid, seqrange);
    agn_logger_log_status(logger, "computed loci for sequence '%s'", seqid);

    gt_error_delete(error);
  }

  return totalloci;
}

GtUword agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx,
                                            const char *refrfile,
                                            const char *predfile,
                                            AgnCompareFilters *filters,
                                            AgnLogger *logger)
{
  gt_assert(idx != NULL);
  GtUword nloci;
  GtFeatureIndex *refrfeats = agn_import_canonical(1, &refrfile, logger);
  GtFeatureIndex *predfeats = agn_import_canonical(1, &predfile, logger);
  if(agn_logger_has_error(logger))
  {
    gt_feature_index_delete(refrfeats);
    gt_feature_index_delete(predfeats);
    return 0;
  }
  nloci = agn_locus_index_parse_pairwise_memory(idx, refrfeats, predfeats,
                                                filters, logger);
  gt_feature_index_delete(refrfeats);
  gt_feature_index_delete(predfeats);
  return nloci;
}

GtUword agn_locus_index_parse_memory(AgnLocusIndex *idx,
                                     GtFeatureIndex *features,
                                     AgnLogger *logger)
{
  int i;
  gt_assert(idx != NULL && features != NULL);
  GtError *error = gt_error_new();
  GtStrArray *seqids = gt_feature_index_get_seqids(features, error);
  if(gt_error_is_set(error))
  {
    agn_logger_log_error(logger, "error fetching seqids: %s",
                         gt_error_get(error));
    gt_error_delete(error);
    gt_str_array_delete(seqids);
    return 0;
  }
  gt_str_array_delete(idx->seqids);
  idx->seqids = seqids;

  GtUword totalloci = 0;

  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtIntervalTree *loci = agn_locus_index_parse(idx,seqid,features,logger);

    GtRange *seqrange = gt_malloc( sizeof(GtRange) );
    GtRange trange;
    GtError *err = gt_error_new();
    gt_feature_index_get_orig_range_for_seqid(features, &trange, seqid, err);
    seqrange->start = trange.start;
    seqrange->end   = trange.end;

    totalloci += gt_interval_tree_size(loci);
    gt_hashmap_add(idx->locus_trees, (char *)seqid, loci);
    gt_hashmap_add(idx->seqranges, (char *)seqid, seqrange);
    agn_logger_log_status(logger, "computed loci for sequence '%s'", seqid);

    gt_error_delete(err);
  }
  gt_error_delete(error);

  return totalloci;
}

int agn_locus_index_test_overlap(AgnLocusIndex *idx, GtFeatureIndex *features,
                                 GtHashmap *visited_genes, AgnGeneLocus *locus,
                                 AgnLogger *logger)
{
  GtError *error = gt_error_new();
  const char *seqid = agn_gene_locus_get_seqid(locus);
  GtRange locusrange = {agn_gene_locus_get_start(locus),
                        agn_gene_locus_get_end(locus)};
  bool has_seqid;
  gt_feature_index_has_seqid(features, &has_seqid, seqid, error);
  if(!has_seqid)
  {
    gt_error_delete(error);
    return 0;
  }

  int new_gene_count = 0;
  GtArray *overlapping_features = gt_array_new( sizeof(GtFeatureNode *) );

  gt_feature_index_get_features_for_range(features, overlapping_features, seqid,
                                          &locusrange, error);
  if(gt_error_is_set(error))
  {
    agn_logger_log_error(logger, "error fetching features for range "
                         "%s[%lu, %lu]: %s", seqid, locusrange.start,
                         locusrange.end, gt_error_get(error));
    gt_error_delete(error);
    gt_array_delete(overlapping_features);
    return 0;
  }
  gt_error_delete(error);

  while(gt_array_size(overlapping_features) > 0)
  {
    GtFeatureNode **fn = gt_array_pop(overlapping_features);
    if(gt_hashmap_get(visited_genes, *fn) == NULL)
    {
      gt_hashmap_add(visited_genes, *fn, *fn);
      agn_gene_locus_add_gene(locus, *fn);
      new_gene_count++;
    }
  }
  gt_array_delete(overlapping_features);

  return new_gene_count;
}

GtStrArray *agn_locus_index_seqids(AgnLocusIndex *idx)
{
  return idx->seqids;
}

bool agn_locus_index_unit_test(AgnUnitTest *test)
{
  AgnLocusIndex *index = agn_locus_index_new(true);
  AgnLogger *logger = agn_logger_new();
  GtFeatureIndex *features = agn_test_data_ilocus_data();
  agn_locus_index_parse_memory(index, features, logger);

  GtArray *seqids = gt_array_new( sizeof(GtStr *) );
  GtUword i;
  for(i = 1; i <= 25; i++)
  {
    char buffer[16];
    sprintf(buffer, "seq%02lu", i);
    GtStr *seqid = gt_str_new_cstr(buffer);
    gt_array_add(seqids, seqid);
  }

  GtArray *iloci;
  AgnGeneLocus *ilocus;
  GtStr *seqid;
  GtUword delta = 200;

  bool test1pass = true;
  seqid = *(GtStr **)gt_array_get(seqids, 0);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange = agn_gene_locus_range(ilocus);
    GtRange testrange = {1, 900};
    if(gt_range_compare(&locusrange, &testrange) != 0)
      test1pass = false;
  }
  else
    test1pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 1);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 2)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange = agn_gene_locus_range(ilocus);
    GtRange testrange = {1, 200};
    if(gt_range_compare(&locusrange, &testrange) != 0)
      test1pass = false;
    else
    {
      ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
      locusrange = agn_gene_locus_range(ilocus);
      GtRange newtestrange = {201, 900};
      if(gt_range_compare(&locusrange, &newtestrange) != 0)
        test1pass = false;
    }
  }
  else
    test1pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 2);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 2)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange = agn_gene_locus_range(ilocus);
    GtRange testrange = {1, 201};
    if(gt_range_compare(&locusrange, &testrange) != 0)
      test1pass = false;
    else
    {
      ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
      locusrange = agn_gene_locus_range(ilocus);
      GtRange newtestrange = {202, 900};
      if(gt_range_compare(&locusrange, &newtestrange) != 0)
        test1pass = false;
    }
  }
  else
    test1pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 3);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange = agn_gene_locus_range(ilocus);
    GtRange testrange = {1, 900};
    if(gt_range_compare(&locusrange, &testrange) != 0)
      test1pass = false;
  }
  else
    test1pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 4);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange = agn_gene_locus_range(ilocus);
    GtRange testrange = {1, 900};
    if(gt_range_compare(&locusrange, &testrange) != 0)
      test1pass = false;
  }
  else
    test1pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 5);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange = agn_gene_locus_range(ilocus);
    GtRange testrange = {1, 900};
    if(gt_range_compare(&locusrange, &testrange) != 0)
      test1pass = false;
  }
  else
    test1pass = false;
  gt_array_delete(iloci);

  agn_unit_test_result(test, "iLocus parsing: initial iLoci", test1pass);


  bool test2pass = true;
  seqid = *(GtStr **)gt_array_get(seqids, 6);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 4)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {801, 1001};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1002, 1600};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 3);
    GtRange locusrange4 = agn_gene_locus_range(ilocus);
    GtRange testrange4 = {1601, 2000};
    if(gt_range_compare(&locusrange4, &testrange4) != 0)
      test2pass = false;
  }
  else
    test2pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 7);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 4)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {801, 1000};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1001, 1600};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 3);
    GtRange locusrange4 = agn_gene_locus_range(ilocus);
    GtRange testrange4 = {1601, 2000};
    if(gt_range_compare(&locusrange4, &testrange4) != 0)
      test2pass = false;
  }
  else
    test2pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 8);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 900};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {901, 1600};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test2pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1601, 2000};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test2pass = false;
  }
  else
    test2pass = false;
  gt_array_delete(iloci);

  agn_unit_test_result(test, "iLocus parsing: internal 3delta", test2pass);


  bool test3pass = true;
  seqid = *(GtStr **)gt_array_get(seqids, 9);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 801};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test3pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {802, 1500};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test3pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1501, 2000};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test3pass = false;
  }
  else
    test3pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 10);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test3pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {801, 1500};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test3pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1501, 2000};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test3pass = false;
  }
  else
    test3pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 11);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test3pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {800, 1500};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test3pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1501, 2000};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test3pass = false;
  }
  else
    test3pass = false;
  gt_array_delete(iloci);

  agn_unit_test_result(test, "iLocus parsing: internal 2delta", test3pass);


  bool test4pass = true;
  seqid = *(GtStr **)gt_array_get(seqids, 12);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {603, 1300};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1301, 1500};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test4pass = false;
  }
  else
    test4pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 13);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {602, 1300};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1301, 1500};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test4pass = false;
  }
  else
    test4pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 14);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {601, 1300};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1301, 1500};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test4pass = false;
  }
  else
    test4pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 15);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 799};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {601, 1300};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test4pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1301, 1500};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test4pass = false;
  }
  else
    test4pass = false;
  gt_array_delete(iloci);

  agn_unit_test_result(test, "iLocus parsing: internal delta", test4pass);


  bool test5pass = true;
  seqid = *(GtStr **)gt_array_get(seqids, 16);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 605};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test5pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {601, 1200};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test5pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1201, 1500};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test5pass = false;
  }
  else
    test5pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 17);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 601};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test5pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {601, 1200};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test5pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1201, 1500};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test5pass = false;
  }
  else
    test5pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 18);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 3)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 600};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test5pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {601, 1200};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test5pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 2);
    GtRange locusrange3 = agn_gene_locus_range(ilocus);
    GtRange testrange3 = {1201, 1500};
    if(gt_range_compare(&locusrange3, &testrange3) != 0)
      test5pass = false;
  }
  else
    test5pass = false;
  gt_array_delete(iloci);

  agn_unit_test_result(test, "iLocus parsing: internal adjacent", test5pass);


  bool test6pass = true;
  seqid = *(GtStr **)gt_array_get(seqids, 19);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 2)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test6pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {801, 1001};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test6pass = false;
  }
  else
    test6pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 20);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 2)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test6pass = false;
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 1);
    GtRange locusrange2 = agn_gene_locus_range(ilocus);
    GtRange testrange2 = {801, 1000};
    if(gt_range_compare(&locusrange2, &testrange2) != 0)
      test6pass = false;
  }
  else
    test6pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 21);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 999};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test6pass = false;
  }
  else
    test6pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 22);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 801};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test6pass = false;
  }
  else
    test6pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 23);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 800};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test6pass = false;
  }
  else
    test6pass = false;
  gt_array_delete(iloci);

  seqid = *(GtStr **)gt_array_get(seqids, 24);
  iloci = agn_locus_index_interval_loci(index, gt_str_get(seqid), delta, false);
  if(gt_array_size(iloci) == 1)
  {
    ilocus = *(AgnGeneLocus **)gt_array_get(iloci, 0);
    GtRange locusrange1 = agn_gene_locus_range(ilocus);
    GtRange testrange1 = {1, 799};
    if(gt_range_compare(&locusrange1, &testrange1) != 0)
      test6pass = false;
  }
  else
    test6pass = false;
  gt_array_delete(iloci);

  agn_unit_test_result(test, "iLocus parsing: terminal iLoci", test6pass);


  while(gt_array_size(seqids) > 0)
  {
    GtStr **seqid = gt_array_pop(seqids);
    gt_str_delete(*seqid);
  }
  gt_array_delete(seqids);
  agn_locus_index_delete(index);

  return test1pass && test2pass && test3pass && test4pass && test5pass &&
         test6pass;
}
