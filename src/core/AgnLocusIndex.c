#include <omp.h>
#include "AgnLocus.h"
#include "AgnLocusIndex.h"
#include "AgnUtils.h"

struct AgnLocusIndex
{
  GtStrArray *seqids;
  GtHashmap *locus_trees;
};

void agn_locus_index_delete(AgnLocusIndex *idx)
{
  gt_hashmap_delete(idx->locus_trees);
  gt_str_array_delete(idx->seqids);
  gt_free(idx);
  idx = NULL;
}

AgnLocusIndex *agn_locus_index_new()
{
  AgnLocusIndex *idx = gt_malloc( sizeof(AgnLocusIndex *) );
  idx->seqids = gt_str_array_new();
  //idx->locus_trees = gt_hashmap_new(GT_HASH_STRING, NULL,
  //                                  (GtFree)gt_interval_tree_delete);
  idx->locus_trees = gt_hashmap_new(GT_HASH_STRING, NULL,
                                    (GtFree)gt_array_delete);
  return idx;
}

unsigned long agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx,
                  GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats,
                  int numprocs, AgnLogger *logger)
{
  int i, rank;
  gt_assert(idx != NULL);
  gt_assert(refrfeats != NULL && predfeats != NULL);
  GtStrArray *seqids = agn_seq_intersection(refrfeats, predfeats, logger);
  if(agn_logger_has_error(logger))
  {
    gt_str_array_delete(seqids);
    return 0;
  }
  gt_str_array_delete(idx->seqids);
  idx->seqids = seqids;
  
  unsigned long totalloci = 0;
  #pragma omp parallel private(i, rank)
  {
    rank = omp_get_thread_num();

    #pragma omp for schedule(static)
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      GtArray *loci = agn_parse_loci(seqid, refrfeats, predfeats);
      #pragma omp critical
      {
        totalloci += gt_array_size(loci);
        gt_hashmap_add(idx->locus_trees, (char *)seqid, loci);
        agn_logger_log_status(logger, "loci for sequence '%s' identified by "
                              "processor %d", seqid, rank);
      }
    }
  } // End parallelize
  
  return totalloci;
}

unsigned long agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx,
                  const char *refrfile, const char *predfile, int numprocs,
                  AgnLogger *logger)
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
                                                numprocs, logger);
  gt_feature_index_delete(refrfeats);
  gt_feature_index_delete(predfeats);
  return nloci;
}

unsigned long agn_locus_index_parse_memory(AgnLocusIndex * idx,
                  GtFeatureIndex *features, int numprocs, AgnLogger *logger)
{
  int i, j, rank;
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
  
  unsigned long totalloci = 0;
  #pragma omp parallel private(i, rank)
  {
    rank = omp_get_thread_num();
    
    #pragma omp for schedule(static)
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      GtArray *seqfeatures = gt_feature_index_get_features_for_seqid(features,
                                 seqid, error);
      if(gt_error_is_set(error))
      {
        agn_logger_log_error(logger, "error fetching features for seqid '%s': "
                             "%s", seqid, gt_error_get(error));
        gt_error_unset(error);
        //break; // Break not supported in OpenMP for loops
      }
      
      GtArray *loci = gt_array_new( sizeof(AgnLocus *) );
      GtHashmap *genes_visited = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
      for(j = 0; j < gt_array_size(seqfeatures); j++)
      {
        GtFeatureNode *feature = *(GtFeatureNode **)gt_array_get(seqfeatures,j);
        GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(feature);
        GtFeatureNode *fn;
        for(fn = gt_feature_node_iterator_next(iter);
            fn != NULL;
            fn = gt_feature_node_iterator_next(iter))
        {
          if(!gt_feature_node_has_type(fn, "gene") ||
             gt_hashmap_get(genes_visited, fn) != NULL)
          {
            continue;
          }
          
          AgnLocus *locus = agn_locus_new(seqid);
          agn_locus_add(locus, fn);
          
          int new_gene_count = 0;
          do
          {
            int temp_new_gene_count = 0;
            GtArray *overlapping_feats = gt_array_new(sizeof(GtFeatureNode *));
            gt_feature_index_get_features_for_range(features,overlapping_feats,
                                                    seqid,&locus->range,error);
            if(gt_error_is_set(error))
            {
              agn_logger_log_error(logger, "error fetching features for range "
                                   "%s[%lu, %lu]: %s", seqid,
                                   locus->range.start, locus->range.end,
                                   gt_error_get(error));
              gt_error_unset(error);
              //break; // Break not supported in OpenMP for loops; plus, it's
              // the outer loop that I want to break out of anyway
            }
            
            while(gt_array_size(overlapping_feats) > 0)
            {
              GtFeatureNode *fn2add = *(GtFeatureNode **)
                                          gt_array_pop(overlapping_feats);
              if(gt_feature_node_has_type(fn2add, "gene") &&
                 gt_hashmap_get(genes_visited, fn2add) == NULL)
              {
                gt_hashmap_add(genes_visited, fn2add, fn2add);
                agn_locus_add(locus, fn2add);
                temp_new_gene_count++;
              }
            }
            gt_array_delete(overlapping_feats);
            new_gene_count = temp_new_gene_count;
          } while(new_gene_count > 0);
          
          gt_array_add(loci, locus);
        }
      }
      gt_hashmap_delete(genes_visited);
      
      #pragma omp critical
      {
        totalloci += gt_array_size(loci);
        gt_hashmap_add(idx->locus_trees, (char *)seqid, loci);
        agn_logger_log_status(logger, "loci for sequence '%s' identified by "
                              "processor %d", seqid, rank);
      }
    } // End parallelize
  }
  
  gt_error_delete(error);
  return totalloci;
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
