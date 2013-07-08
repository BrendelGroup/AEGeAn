#include <omp.h>
#include "AgnLocusIndex.h"
#include "AgnGeneLocus.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnLocusIndex
{
  GtStrArray *seqids;
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
int agn_locus_index_it_traverse(GtIntervalTreeNode *itn, void *lp);

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
GtIntervalTree *agn_locus_index_parse_pairwise(AgnLocusIndex *idx,
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
int agn_locus_index_pairwise_test_overlap(AgnLocusIndex *idx,
                  GtFeatureIndex *features, GtHashmap *visited_genes,
                  AgnGeneLocus *locus, AgnComparisonSource source,
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
int agn_locus_index_test_overlap(AgnLocusIndex *idx, GtFeatureIndex *features,
                                 GtHashmap *visited_genes, AgnGeneLocus *locus,
                                 AgnLogger *logger);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_locus_index_delete(AgnLocusIndex *idx)
{
  gt_hashmap_delete(idx->locus_trees);
  gt_str_array_delete(idx->seqids);
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
                                       unsigned long delta)
{
  GtArray *loci = agn_locus_index_get(idx, seqid);
  if(loci == NULL)
    return NULL;
  gt_array_sort(loci, (GtCompare)agn_gene_locus_array_compare);

  unsigned long nloci = gt_array_size(loci);
  gt_assert(nloci > 0);
  GtArray *iloci = gt_array_new( sizeof(AgnGeneLocus *) );
  
  // Handle intial iloci (unless there is only 1 gene locus)
  AgnGeneLocus *l1 = *(AgnGeneLocus **)gt_array_get(loci, 0);
  GtRange r1 = agn_gene_locus_range(l1);
  if(r1.start > 1 && r1.start > delta + 1)
  {
    AgnGeneLocus *ilocus = agn_gene_locus_new(seqid);
    agn_gene_locus_set_range(ilocus, 1, r1.start - delta - 1);
    gt_array_add(iloci, ilocus);
  }
  if(nloci > 1)
  {
    AgnGeneLocus *l2 = *(AgnGeneLocus **)gt_array_get(loci, 1);
    GtRange r2 = agn_gene_locus_range(l2);

    unsigned long newstart, newend;
    newstart = (r1.start > 1 && r1.start > delta + 1) ? r1.start - delta : 1;
    newend   = (r1.end + delta < r2.start) ? r1.end + delta : r2.start - 1;
    l1 = agn_gene_locus_new(seqid);
    agn_gene_locus_set_range(l1, newstart, newend);
    gt_array_add(iloci, l1);

    newstart = newend + 1;
    newend   = (r2.start - delta >= newend) ? r2.start - delta :
                                             r2.start - 1;
    AgnGeneLocus *ilocus = agn_gene_locus_new(seqid);
    agn_gene_locus_set_range(ilocus, newstart, newend);
    gt_array_add(iloci, ilocus);
  }

  // Handle internal (not initial or terminal) iloci
  unsigned long i;
  for(i = 1; i < nloci - 1; i++)
  {
    AgnGeneLocus *llocus = *(AgnGeneLocus **)gt_array_get(loci, i-1);
    AgnGeneLocus *locus  = *(AgnGeneLocus **)gt_array_get(loci, i);
    AgnGeneLocus *rlocus = *(AgnGeneLocus **)gt_array_get(loci, i+1);
    GtRange lrange = agn_gene_locus_range(llocus);
    GtRange range = agn_gene_locus_range(locus);
    GtRange rrange = agn_gene_locus_range(rlocus);

    unsigned long newstart, newend;
    newstart = (range.start - delta > lrange.end) ? range.start - delta :
                                                    lrange.end + 1;
    newend   = (range.end + delta < rrange.start) ? range.end + delta :
                                                    rrange.start - 1;
    AgnGeneLocus *ilocus = agn_gene_locus_new(seqid);
    agn_gene_locus_set_range(ilocus, newstart, newend);
    gt_array_add(iloci, ilocus);
    
    if(newend < rrange.start - 1)
    {
      newstart = newend + 1;
      newend   = (rrange.start - delta > newend) ? rrange.start - delta :
                                                   rrange.start - 1;
      ilocus = agn_gene_locus_new(seqid);
      agn_gene_locus_set_range(ilocus, newstart, newend);
      gt_array_add(iloci, ilocus);
    }
  }

  // FIXME Handle terminal ilocus
  l1 = *(AgnGeneLocus **)gt_array_get(loci, nloci-1);
  r1 = agn_gene_locus_range(l1);

  gt_array_delete(loci);
  return iloci;
}

int agn_locus_index_it_traverse(GtIntervalTreeNode *itn, void *lp)
{
  GtArray *loci = (GtArray *)lp;
  // Stored as a void* since we don't know whether these are AgnLocus objects or
  // AgnGeneLocus objects.
  void *locus = gt_interval_tree_node_get_data(itn);
  gt_array_add(loci, locus);
  return 0;
}

AgnLocusIndex *agn_locus_index_new(bool freeondelete)
{
  AgnLocusIndex *idx = gt_malloc( sizeof(AgnLocusIndex) );
  idx->seqids = gt_str_array_new();
  idx->locus_trees = gt_hashmap_new(GT_HASH_STRING, NULL,
                                    (GtFree)gt_interval_tree_delete);
  idx->locusfreefunc = NULL;
  if(freeondelete)
    idx->locusfreefunc = (GtFree)agn_gene_locus_delete;

  return idx;
}

unsigned long agn_locus_index_parse_disk(AgnLocusIndex * idx, int numfiles,
                  const char **filenames, int numprocs, AgnLogger *logger)
{
  gt_assert(idx != NULL);
  unsigned long nloci;
  // Do I want to use the import canonical function or the add_gff3 function?
  GtFeatureIndex *features = agn_import_canonical(numfiles, filenames, logger);
  if(agn_logger_has_error(logger))
  {
    gt_feature_index_delete(features);
    return 0;
  }

  nloci = agn_locus_index_parse_memory(idx, features, numprocs, logger);
  gt_feature_index_delete(features);
  return nloci;
}

int agn_locus_index_pairwise_test_overlap(AgnLocusIndex *idx,
                  GtFeatureIndex *features, GtHashmap *visited_genes,
                  AgnGeneLocus *locus, AgnComparisonSource source,
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
  unsigned long new_gene_count = 0;

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
  unsigned long i;
  GtError *error = gt_error_new();
  GtHashmap *visited_genes = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  //GtIntervalTree *loci = gt_interval_tree_new((GtFree)agn_locus_delete);
  GtIntervalTree *loci = gt_interval_tree_new(NULL);

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
    GtFeatureNode *feature = *(GtFeatureNode **)gt_array_get(seqfeatures, i);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(feature);
    GtFeatureNode *fn;
    for(fn = gt_feature_node_iterator_next(iter);
        fn != NULL;
        fn = gt_feature_node_iterator_next(iter))
    {
      if(!gt_feature_node_has_type(fn, "gene") ||
         gt_hashmap_get(visited_genes, fn) != NULL)
      {
        continue;
      }

      gt_hashmap_add(visited_genes, fn, fn);
      AgnGeneLocus *locus = agn_gene_locus_new(seqid);
      agn_gene_locus_add_gene(locus, fn);

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
    gt_feature_node_iterator_delete(iter);
  }

  gt_error_delete(error);
  gt_hashmap_delete(visited_genes);
  gt_array_delete(seqfeatures);
  return loci;
}

GtIntervalTree *agn_locus_index_parse_pairwise(AgnLocusIndex *idx,
                                               const char *seqid,
                                               GtFeatureIndex *refr,
                                               GtFeatureIndex *pred,
                                               AgnCompareFilters *filters,
                                               AgnLogger *logger)
{
  unsigned long i;
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

unsigned long agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx,
                  GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats,
                  int numprocs, AgnCompareFilters *filters, AgnLogger *logger)
{
  int i, rank;
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

  int orig_numprocs = omp_get_num_threads();
  omp_set_num_threads(numprocs);
  unsigned long totalloci = 0;
  #pragma omp parallel private(i, rank)
  {
    rank = omp_get_thread_num();

    #pragma omp for schedule(static)
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      GtIntervalTree *loci = agn_locus_index_parse_pairwise(idx, seqid,
                                                            refrfeats,
                                                            predfeats, filters,
                                                            logger);
      #pragma omp critical
      {
        totalloci += gt_interval_tree_size(loci);
        gt_hashmap_add(idx->locus_trees, (char *)seqid, loci);
        agn_logger_log_status(logger, "loci for sequence '%s' identified by "
                              "processor %d", seqid, rank);
      }
    }
  } // End parallelize
  omp_set_num_threads(orig_numprocs);

  return totalloci;
}

unsigned long agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx,
                  const char *refrfile, const char *predfile, int numprocs,
                  AgnCompareFilters *filters, AgnLogger *logger)
{
  gt_assert(idx != NULL);
  unsigned long nloci;
  GtFeatureIndex *refrfeats = agn_import_canonical(1, &refrfile, logger);
  GtFeatureIndex *predfeats = agn_import_canonical(1, &predfile, logger);
  if(agn_logger_has_error(logger))
  {
    gt_feature_index_delete(refrfeats);
    gt_feature_index_delete(predfeats);
    return 0;
  }
  nloci = agn_locus_index_parse_pairwise_memory(idx, refrfeats, predfeats,
                                                numprocs, filters, logger);
  gt_feature_index_delete(refrfeats);
  gt_feature_index_delete(predfeats);
  return nloci;
}

unsigned long agn_locus_index_parse_memory(AgnLocusIndex * idx,
                  GtFeatureIndex *features, int numprocs, AgnLogger *logger)
{
  int i, rank;
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

  int orig_numprocs = omp_get_num_threads();
  omp_set_num_threads(numprocs);
  unsigned long totalloci = 0;
  #pragma omp parallel private(i, rank)
  {
    rank = omp_get_thread_num();

    #pragma omp for schedule(static)
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      GtIntervalTree *loci = agn_locus_index_parse(idx,seqid,features,logger);
      #pragma omp critical
      {
        totalloci += gt_interval_tree_size(loci);
        gt_hashmap_add(idx->locus_trees, (char *)seqid, loci);
        agn_logger_log_status(logger, "loci for sequence '%s' identified by "
                              "processor %d", seqid, rank);
      }
    } // End parallelize
  }
  omp_set_num_threads(orig_numprocs);

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
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(overlapping_features);
    if(gt_feature_node_has_type(fn, "gene") &&
       gt_hashmap_get(visited_genes, fn) == NULL)
    {
      gt_hashmap_add(visited_genes, fn, fn);
      agn_gene_locus_add_gene(locus, fn);
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
