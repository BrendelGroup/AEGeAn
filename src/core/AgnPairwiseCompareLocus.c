#include <math.h>
#include <string.h>
#include "extended/feature_node.h"
#include "AgnPairwiseCompareLocus.h"
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnPairwiseCompareLocus
{
  char *seqid;
  GtDlist *refr_genes;
  GtDlist *pred_genes;
  GtRange range;
  GtArray *refr_cliques;
  GtArray *pred_cliques;
  GtArray *clique_pairs;
  bool     clique_pairs_formed;
  GtArray *reported_pairs;
  GtArray *unique_refr_cliques;
  GtArray *unique_pred_cliques;
  AgnComparisonStats stats;
  AgnComparisonCounts counts;
  AgnCompareClassAggregateDescription results;
  double refr_splice_complexity;
  double pred_splice_complexity;
};


//----------------------------------------------------------------------------//
// Prototypes for private method(s)
//----------------------------------------------------------------------------//

#ifndef WITHOUT_CAIRO
/**
 * Track selector function for PNG graphics
 *
 * @param[in]  block    the object to be printed in the graphic
 * @param[out] track    string to which the appropriate track name will be
 *                      written
 * @param[in]  data     any auxilliary data needed
 */
void agn_pairwise_compare_locus_png_track_selector(GtBlock *block, GtStr *track, void *data);
#endif

/**
 * Update this locus' start and end coordinates based on the gene being merged.
 *
 * @param[out] locus           the locus
 * @param[in]  gene_feature    gene that was just merged with the locus
 */
void agn_pairwise_compare_locus_update_range( AgnPairwiseCompareLocus *locus,
                                 GtFeatureNode *gene_feature );


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
void agn_pairwise_compare_locus_add_pred_gene( AgnPairwiseCompareLocus *locus,
                                  GtFeatureNode *gene_feature )
{
  if(locus->pred_genes == NULL)
  {
    locus->pred_genes = gt_dlist_new( (GtCompare)gt_genome_node_cmp );
  }

  gt_dlist_add(locus->pred_genes, gene_feature);
  agn_pairwise_compare_locus_update_range(locus, gene_feature);
}

void agn_pairwise_compare_locus_add_refr_gene( AgnPairwiseCompareLocus *locus,
                                  GtFeatureNode *gene_feature )
{
  if(locus->refr_genes == NULL)
  {
    locus->refr_genes = gt_dlist_new( (GtCompare)gt_genome_node_cmp );
  }

  gt_dlist_add(locus->refr_genes, gene_feature);
  agn_pairwise_compare_locus_update_range(locus, gene_feature);
}

void agn_pairwise_compare_locus_aggregate_results( AgnPairwiseCompareLocus *locus,
                                      AgnSummaryData *data )
{
  unsigned long i;
  GtArray *reported_pairs = agn_pairwise_compare_locus_find_best_pairs(locus);
  unsigned long pairs_to_report = gt_array_size(reported_pairs);
  gt_assert(pairs_to_report > 0 && reported_pairs != NULL);

  for(i = 0; i < pairs_to_report; i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(reported_pairs, i);
    gt_assert(agn_clique_pair_needs_comparison(pair));

    data->counts.num_comparisons++;

    // Classify this comparison: perfect match, CDS match, etc.
    unsigned int compareclass = agn_clique_pair_classify(pair);
    switch(compareclass)
    {
      case PE_CLIQUE_PAIR_PERFECT_MATCH:
        data->counts.num_perfect++;
        agn_clique_pair_record_characteristics
        (
          pair,
          &data->results.perfect_matches
        );
        break;

      case PE_CLIQUE_PAIR_MISLABELED:
        data->counts.num_mislabeled++;
        agn_clique_pair_record_characteristics
        (
          pair,
          &data->results.perfect_mislabeled
        );
        break;

      case PE_CLIQUE_PAIR_CDS_MATCH:
        data->counts.num_cds_match++;
        agn_clique_pair_record_characteristics
        (
          pair,
          &data->results.cds_matches
        );
        break;

      case PE_CLIQUE_PAIR_EXON_MATCH:
        data->counts.num_exon_match++;
        agn_clique_pair_record_characteristics
        (
          pair,
          &data->results.exon_matches
        );
        break;

      case PE_CLIQUE_PAIR_UTR_MATCH:
        data->counts.num_utr_match++;
        agn_clique_pair_record_characteristics
        (
          pair,
          &data->results.utr_matches
        );
        break;

      case PE_CLIQUE_PAIR_NON_MATCH:
        data->counts.non_match++;
        agn_clique_pair_record_characteristics
        (
          pair,
          &data->results.non_matches
        );
        break;

      default:
        fprintf(stderr, "error: unknown classification %d\n", compareclass);
        exit(1);
        break;
    }

    // Record structure-level counts
    AgnComparisonStats *pairstats = agn_clique_pair_get_stats(pair);
    data->stats.cds_struc_stats.correct  += pairstats->cds_struc_stats.correct;
    data->stats.cds_struc_stats.missing  += pairstats->cds_struc_stats.missing;
    data->stats.cds_struc_stats.wrong    += pairstats->cds_struc_stats.wrong;
    data->stats.exon_struc_stats.correct += pairstats->exon_struc_stats.correct;
    data->stats.exon_struc_stats.missing += pairstats->exon_struc_stats.missing;
    data->stats.exon_struc_stats.wrong   += pairstats->exon_struc_stats.wrong;
    data->stats.utr_struc_stats.correct  += pairstats->utr_struc_stats.correct;
    data->stats.utr_struc_stats.missing  += pairstats->utr_struc_stats.missing;
    data->stats.utr_struc_stats.wrong    += pairstats->utr_struc_stats.wrong;

    // Record nucleotide-level counts
    data->stats.cds_nuc_stats.tp += pairstats->cds_nuc_stats.tp;
    data->stats.cds_nuc_stats.fn += pairstats->cds_nuc_stats.fn;
    data->stats.cds_nuc_stats.fp += pairstats->cds_nuc_stats.fp;
    data->stats.cds_nuc_stats.tn += pairstats->cds_nuc_stats.tn;
    data->stats.utr_nuc_stats.tp += pairstats->utr_nuc_stats.tp;
    data->stats.utr_nuc_stats.fn += pairstats->utr_nuc_stats.fn;
    data->stats.utr_nuc_stats.fp += pairstats->utr_nuc_stats.fp;
    data->stats.utr_nuc_stats.tn += pairstats->utr_nuc_stats.tn;
    data->stats.overall_matches  += pairstats->overall_matches;
    data->stats.overall_length   += agn_pairwise_compare_locus_get_length(locus);
  }
}

void agn_pairwise_compare_locus_calc_splice_complexity_pred(AgnPairwiseCompareLocus *locus)
{
  GtArray *pred_trans = agn_pairwise_compare_locus_get_pred_transcripts(locus);
// fprintf(stderr, "DELETEME pred sc\n");
  locus->pred_splice_complexity = agn_calc_splice_complexity(pred_trans);
  gt_array_delete(pred_trans);
}

void agn_pairwise_compare_locus_calc_splice_complexity_refr(AgnPairwiseCompareLocus *locus)
{
  GtArray *refr_trans = agn_pairwise_compare_locus_get_refr_transcripts(locus);
// fprintf(stderr, "DELETEME refr sc\n");
  locus->refr_splice_complexity = agn_calc_splice_complexity(refr_trans);
  gt_array_delete(refr_trans);
}

int agn_pairwise_compare_locus_array_compare(const void *p1, const void *p2)
{
  AgnPairwiseCompareLocus *l1 = *(AgnPairwiseCompareLocus **)p1;
  AgnPairwiseCompareLocus *l2 = *(AgnPairwiseCompareLocus **)p2;

  if(l1->range.start == l2->range.start && l1->range.end == l2->range.end)
    return 0;
  else if
  (
    (l1->range.start < l2->range.start) ||
    (l1->range.start == l2->range.start && l1->range.end < l2->range.end)
  )
  {
    return -1;
  }
  
  return 1;
}

void agn_pairwise_compare_locus_delete(AgnPairwiseCompareLocus *locus)
{
  if(locus->refr_genes != NULL)
  {
    GtDlistelem *current;
    for( current = gt_dlist_first(locus->refr_genes);
         current != NULL;
         current = gt_dlistelem_next(current) )
    {
      GtGenomeNode *gn = gt_dlistelem_get_data(current);
      gt_genome_node_delete(gn);
    }
    gt_dlist_delete(locus->refr_genes);
  }
  if(locus->pred_genes != NULL)
  {
    GtDlistelem *current;
    for( current = gt_dlist_first(locus->pred_genes);
         current != NULL;
         current = gt_dlistelem_next(current) )
    {
      GtGenomeNode *gn = gt_dlistelem_get_data(current);
      gt_genome_node_delete(gn);
    }
    gt_dlist_delete(locus->pred_genes);
  }
  if(locus->refr_cliques != NULL)
  {
    while(gt_array_size(locus->refr_cliques) > 0)
    {
      AgnTranscriptClique *clique = *(AgnTranscriptClique **)
                                              gt_array_pop(locus->refr_cliques);
      agn_transcript_clique_delete(clique);
    }
    gt_array_delete(locus->refr_cliques);
  }
  if(locus->pred_cliques != NULL)
  {
    while(gt_array_size(locus->pred_cliques) > 0)
    {
      AgnTranscriptClique *clique = *(AgnTranscriptClique **)
                                              gt_array_pop(locus->pred_cliques);
      agn_transcript_clique_delete(clique);
    }
    gt_array_delete(locus->pred_cliques);
  }
  if(locus->clique_pairs != NULL)
  {
    while(gt_array_size(locus->clique_pairs) > 0)
    {
      AgnCliquePair *pair = *(AgnCliquePair **)gt_array_pop(locus->clique_pairs);
      agn_clique_pair_delete(pair);
    }
    gt_array_delete(locus->clique_pairs);
  }

  if(locus->reported_pairs != NULL)
    gt_array_delete(locus->reported_pairs);
  if(locus->unique_refr_cliques != NULL)
    gt_array_delete(locus->unique_refr_cliques);
  if(locus->unique_pred_cliques != NULL)
    gt_array_delete(locus->unique_pred_cliques);

  gt_free(locus->seqid);
  gt_free(locus);
  locus = NULL;
}

bool agn_pairwise_compare_locus_filter( AgnPairwiseCompareLocus *locus,
                                        AgnCompareFilters *filters )
{
  if(filters == NULL)
    return false;

  // Ignore all filter configuration settings set to 0

  // Filter by locus length
  unsigned long length = agn_pairwise_compare_locus_get_length(locus);
  if(filters->LocusLengthUpperLimit > 0)
  {
    if(length > filters->LocusLengthUpperLimit)
      return true;
  }
  if(filters->LocusLengthLowerLimit > 0)
  {
    if(length < filters->LocusLengthLowerLimit)
      return true;
  }

  // Filter by reference gene models
  unsigned long num_refr_genes = agn_pairwise_compare_locus_num_refr_genes(locus);
  if(filters->MinReferenceGeneModels > 0)
  {
    if(num_refr_genes < filters->MinReferenceGeneModels)
      return true;
  }
  if(filters->MaxReferenceGeneModels > 0)
  {
    if(num_refr_genes > filters->MaxReferenceGeneModels)
      return true;
  }

  // Filter by prediction gene models
  unsigned long num_pred_genes = agn_pairwise_compare_locus_num_pred_genes(locus);
  if(filters->MinPredictionGeneModels > 0)
  {
    if(num_pred_genes < filters->MinPredictionGeneModels)
      return true;
  }
  if(filters->MaxPredictionGeneModels > 0)
  {
    if(num_pred_genes > filters->MaxPredictionGeneModels)
      return true;
  }

  // Filter by reference transcript models
  unsigned long num_refr_transcripts = agn_pairwise_compare_locus_num_refr_transcripts(locus);
  if(filters->MinReferenceTranscriptModels > 0)
  {
    if(num_refr_transcripts < filters->MinReferenceTranscriptModels)
      return true;
  }
  if(filters->MaxReferenceTranscriptModels > 0)
  {
    if(num_refr_transcripts > filters->MaxReferenceTranscriptModels)
      return true;
  }

  // Filter by prediction transcript models
  unsigned long num_pred_transcripts = agn_pairwise_compare_locus_num_pred_transcripts(locus);
  if(filters->MinPredictionTranscriptModels > 0)
  {
    if(num_pred_transcripts < filters->MinPredictionTranscriptModels)
      return true;
  }
  if(filters->MaxPredictionTranscriptModels > 0)
  {
    if(num_pred_transcripts > filters->MaxPredictionTranscriptModels)
      return true;
  }

  // Filter by number of transcripts per gene model
  GtArray *refr_genes = agn_pairwise_compare_locus_get_refr_genes(locus);
  if(filters->MinTranscriptsPerReferenceGeneModel > 0)
  {
    unsigned long i;
    bool success = false;
    for(i = 0; i < gt_array_size(refr_genes); i++)
    {
      GtFeatureNode *gene = gt_array_get(refr_genes, i);
      if( agn_gt_feature_node_num_transcripts(gene) >=
          filters->MinTranscriptsPerReferenceGeneModel )
      {
        success = true;
        break;
      }
    }
    if(!success)
    {
      gt_array_delete(refr_genes);
      return true;
    }
  }
  if(filters->MaxTranscriptsPerReferenceGeneModel > 0)
  {
    unsigned long i;
    bool success = false;
    for(i = 0; i < gt_array_size(refr_genes); i++)
    {
      GtFeatureNode *gene = gt_array_get(refr_genes, i);
      if( agn_gt_feature_node_num_transcripts(gene) <=
          filters->MaxTranscriptsPerReferenceGeneModel )
      {
        success = true;
        break;
      }
    }
    if(!success)
    {
      gt_array_delete(refr_genes);
      return true;
    }
  }
  gt_array_delete(refr_genes);
  GtArray *pred_genes = agn_pairwise_compare_locus_get_pred_genes(locus);
  if(filters->MinTranscriptsPerPredictionGeneModel > 0)
  {
    unsigned long i;
    bool success = false;
    for(i = 0; i < gt_array_size(pred_genes); i++)
    {
      GtFeatureNode *gene = gt_array_get(pred_genes, i);
      if( agn_gt_feature_node_num_transcripts(gene) >=
          filters->MinTranscriptsPerPredictionGeneModel )
      {
        success = true;
        break;
      }
    }
    if(!success)
    {
      gt_array_delete(pred_genes);
      return true;
    }
  }
  if(filters->MaxTranscriptsPerPredictionGeneModel > 0)
  {
    unsigned long i;
    bool success = false;
    for(i = 0; i < gt_array_size(pred_genes); i++)
    {
      GtFeatureNode *gene = gt_array_get(pred_genes, i);
      if( agn_gt_feature_node_num_transcripts(gene) <=
          filters->MaxTranscriptsPerPredictionGeneModel )
      {
        success = true;
        break;
      }
    }
    if(!success)
    {
      gt_array_delete(refr_genes);
      return true;
    }
  }
  gt_array_delete(pred_genes);

  // Filter by reference exons
  unsigned long num_refr_exons = agn_pairwise_compare_locus_num_refr_exons(locus);
  if(filters->MinReferenceExons > 0)
  {
    if(num_refr_exons < filters->MinReferenceExons)
      return true;
  }
  if(filters->MaxReferenceExons > 0)
  {
    if(num_refr_exons > filters->MaxReferenceExons)
      return true;
  }

  // Filter by prediction exons
  unsigned long num_pred_exons = agn_pairwise_compare_locus_num_pred_exons(locus);
  if(filters->MinPredictionExons > 0)
  {
    if(num_pred_exons < filters->MinPredictionExons)
      return true;
  }
  if(filters->MaxPredictionExons > 0)
  {
    if(num_pred_exons > filters->MaxPredictionExons)
      return true;
  }

  // Filter by reference CDS length
  unsigned long refr_cds_length = agn_pairwise_compare_locus_refr_cds_length(locus);
  if(filters->MinReferenceCDSLength > 0)
  {
    if(refr_cds_length < filters->MinReferenceCDSLength)
      return true;
  }
  if(filters->MaxReferenceCDSLength > 0)
  {
    if(refr_cds_length > filters->MaxReferenceCDSLength)
      return true;
  }

  // Filter by prediction CDS length
  unsigned long pred_cds_length = agn_pairwise_compare_locus_pred_cds_length(locus);
  if(filters->MinPredictionCDSLength > 0)
  {
    if(pred_cds_length < filters->MinPredictionCDSLength)
      return true;
  }
  if(filters->MaxPredictionCDSLength > 0)
  {
    if(pred_cds_length > filters->MaxPredictionCDSLength)
      return true;
  }

  return false;
}

GtArray *agn_pairwise_compare_locus_find_best_pairs(AgnPairwiseCompareLocus *locus)
{
  if( locus->refr_cliques == NULL || gt_array_size(locus->refr_cliques) == 0 ||
      locus->pred_cliques == NULL || gt_array_size(locus->pred_cliques) == 0 )
    return NULL;

  if(locus->reported_pairs != NULL)
    return locus->reported_pairs;

  gt_array_sort(locus->clique_pairs, (GtCompare)agn_clique_pair_compare_reverse);
  GtHashmap *refr_cliques_acctd = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  GtHashmap *pred_cliques_acctd = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  locus->reported_pairs = gt_array_new( sizeof(AgnCliquePair *) );

  unsigned long i;
  for(i = 0; i < gt_array_size(locus->clique_pairs); i++)
  {
    AgnCliquePair *pair = *(AgnCliquePair **)
                                           gt_array_get(locus->clique_pairs, i);
    AgnTranscriptClique *rclique = agn_clique_pair_get_refr_clique(pair);
    AgnTranscriptClique *pclique = agn_clique_pair_get_pred_clique(pair);
    if
    (
      agn_transcript_clique_has_id_in_hash(rclique, refr_cliques_acctd) == false &&
      agn_transcript_clique_has_id_in_hash(pclique, pred_cliques_acctd) == false
    )
    {
      gt_array_add(locus->reported_pairs, pair);
      agn_transcript_clique_put_ids_in_hash(rclique, refr_cliques_acctd);
      agn_transcript_clique_put_ids_in_hash(pclique, pred_cliques_acctd);
    }
  }

  locus->unique_refr_cliques = gt_array_new( sizeof(AgnTranscriptClique *) );
  for(i = 0; i < gt_array_size(locus->refr_cliques); i++)
  {
    AgnTranscriptClique *refr_clique = *(AgnTranscriptClique **)
                                           gt_array_get(locus->refr_cliques, i);
    if(agn_transcript_clique_has_id_in_hash(refr_clique, refr_cliques_acctd) == false)
    {
      gt_array_add(locus->unique_refr_cliques, refr_clique);
      agn_transcript_clique_put_ids_in_hash(refr_clique, refr_cliques_acctd);
    }
  }

  locus->unique_pred_cliques = gt_array_new( sizeof(AgnTranscriptClique *) );
  for(i = 0; i < gt_array_size(locus->pred_cliques); i++)
  {
    AgnTranscriptClique *pred_clique = *(AgnTranscriptClique **)
                                           gt_array_get(locus->pred_cliques, i);
    if(agn_transcript_clique_has_id_in_hash(pred_clique, pred_cliques_acctd) == false)
    {
      gt_array_add(locus->unique_pred_cliques, pred_clique);
      agn_transcript_clique_put_ids_in_hash(pred_clique, pred_cliques_acctd);
    }
  }

  gt_hashmap_delete(refr_cliques_acctd);
  gt_hashmap_delete(pred_cliques_acctd);
  return locus->reported_pairs;
}

GtArray* agn_pairwise_compare_locus_get_clique_pairs( AgnPairwiseCompareLocus *locus,
                                          unsigned int trans_per_locus )
{
  GtArray *refr_trans, *pred_trans;
  GtDlistelem *elem;

  // No need to do this multiple times
  if(locus->clique_pairs_formed)
    return locus->clique_pairs;

  // Grab reference transcripts
  refr_trans = gt_array_new( sizeof(GtFeatureNode *) );
  if(locus->refr_genes)
  {
    for( elem = gt_dlist_first(locus->refr_genes);
         elem != NULL;
         elem = gt_dlistelem_next(elem) )
    {
      GtFeatureNode *gene = (GtFeatureNode *)gt_dlistelem_get_data(elem);
      GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
      GtFeatureNode *feature;

      for( feature = gt_feature_node_iterator_next(iter);
           feature != NULL;
           feature = gt_feature_node_iterator_next(iter) )
      {
        if(agn_gt_feature_node_is_mrna_feature(feature))
          gt_array_add(refr_trans, feature);
      }
      gt_feature_node_iterator_delete(iter);
    }
  }

  // Grab prediction transcripts
  pred_trans = gt_array_new( sizeof(GtFeatureNode *) );
  if(locus->pred_genes)
  {
    for( elem = gt_dlist_first(locus->pred_genes);
         elem != NULL;
         elem = gt_dlistelem_next(elem) )
    {
      GtFeatureNode *gene = (GtFeatureNode *)gt_dlistelem_get_data(elem);
      GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
      GtFeatureNode *feature;

      for( feature = gt_feature_node_iterator_next(iter);
           feature != NULL;
           feature = gt_feature_node_iterator_next(iter) )
      {
        if(agn_gt_feature_node_is_mrna_feature(feature))
          gt_array_add(pred_trans, feature);
      }
      gt_feature_node_iterator_delete(iter);
    }
  }

  // Compute maximal transcript cliques
  bool have_refr_trans = gt_array_size(refr_trans) > 0;
  bool refr_trans_reasonable =
  (
    gt_array_size(refr_trans) <= trans_per_locus ||
    trans_per_locus == 0
  );
  if(have_refr_trans && refr_trans_reasonable)
  {
    locus->refr_cliques = agn_enumerate_feature_cliques(refr_trans);
  }
  bool have_pred_trans = gt_array_size(pred_trans) > 0;
  bool pred_trans_reasonable =
  (
    gt_array_size(pred_trans) <= trans_per_locus ||
    trans_per_locus == 0
  );
  if(have_pred_trans && pred_trans_reasonable)
  {
    locus->pred_cliques = agn_enumerate_feature_cliques(pred_trans);
  }

  // Pairs require both reference and prediction cliques
  if(!have_refr_trans || !have_pred_trans)
  {
    locus->clique_pairs_formed = true;
    //if(options->debug)
    if(false)
    {
      fprintf
      (
        stderr,
        "debug: skipping locus %s[%lu, %lu] with %lu reference transcripts and"
        " %lu prediction transcripts (must have at least 1 transcript for"
        " both)\n",
        locus->seqid,
        locus->range.start,
        locus->range.end,
        gt_array_size(refr_trans),
        gt_array_size(pred_trans)
      );
    }
    gt_array_delete(refr_trans);
    gt_array_delete(pred_trans);
    return NULL;
  }

  if(!refr_trans_reasonable || !pred_trans_reasonable)
  {
    locus->clique_pairs_formed = true;
    //if(options->debug)
    if(false)
    {
      fprintf
      (
        stderr,
        "debug: skipping locus %s[%lu, %lu] with %lu reference transcripts and"
        " %lu prediction transcripts (exceeds reasonable limit of %u)\n",
        locus->seqid,
        locus->range.start,
        locus->range.end,
        gt_array_size(refr_trans),
        gt_array_size(pred_trans),
        trans_per_locus
      );
    }
    gt_array_delete(refr_trans);
    gt_array_delete(pred_trans);
    return NULL;
  }

  gt_array_delete(refr_trans);
  gt_array_delete(pred_trans);

  // Form all possible pairs of reference and prediction cliques
  GtArray *clique_pairs = gt_array_new( sizeof(AgnCliquePair *) );
  unsigned long i,j;
  for(i = 0; i < gt_array_size(locus->refr_cliques); i++)
  {
    AgnTranscriptClique *refr_clique = *(AgnTranscriptClique **)
                                           gt_array_get(locus->refr_cliques, i);
    for(j = 0; j < gt_array_size(locus->pred_cliques); j++)
    {
      AgnTranscriptClique *pred_clique = *(AgnTranscriptClique **)
                                            gt_array_get(locus->pred_cliques,j);
      AgnCliquePair *pair = agn_clique_pair_new( locus->seqid, refr_clique,
                                                 pred_clique, &locus->range );
      gt_array_add(clique_pairs, pair);
    }
  }

  locus->clique_pairs = clique_pairs;
  locus->clique_pairs_formed = true;
  return locus->clique_pairs;
}

unsigned long agn_pairwise_compare_locus_get_end(AgnPairwiseCompareLocus *locus)
{
  return locus->range.end;
}

unsigned long agn_pairwise_compare_locus_get_length(AgnPairwiseCompareLocus *locus)
{
  return gt_range_length(&locus->range);
}

AgnCliquePair* agn_pairwise_compare_locus_get_optimal_clique_pair
(
  AgnPairwiseCompareLocus *locus,
  AgnTranscriptClique *refr_clique
)
{
  AgnCliquePair *bestpair = NULL;
  unsigned long i;

  for(i = 0; i < gt_array_size(locus->clique_pairs); i++)
  {
    AgnCliquePair *currentpair = *(AgnCliquePair **)
                                    gt_array_get(locus->clique_pairs, i);
    if(refr_clique == agn_clique_pair_get_refr_clique(currentpair))
    {
      if( bestpair == NULL ||
          agn_clique_pair_compare_direct(bestpair, currentpair) == -1 )
      {
        bestpair = currentpair;
      }
    }
  }

  return bestpair;
}

GtArray *agn_pairwise_compare_locus_get_pred_genes(AgnPairwiseCompareLocus *locus)
{
  if(locus->pred_genes == NULL)
    return NULL;

  GtArray *genes = gt_array_new( sizeof(GtFeatureNode *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->pred_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    gt_array_add(genes, gene);
  }
  gt_array_sort(genes, (GtCompare)agn_gt_genome_node_compare);
  return genes;
}

GtArray *agn_pairwise_compare_locus_get_pred_gene_ids(AgnPairwiseCompareLocus *locus)
{
  if(locus->pred_genes == NULL)
    return NULL;

  GtArray *ids = gt_array_new( sizeof(char *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->pred_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    const char *id = gt_feature_node_get_attribute(gene, "ID");
    gt_array_add(ids, id);
  }
  gt_array_sort(ids, (GtCompare)agn_string_compare);
  return ids;
}

double agn_pairwise_compare_locus_get_pred_splice_complexity(AgnPairwiseCompareLocus *locus)
{
  return locus->pred_splice_complexity;
}

GtArray *agn_pairwise_compare_locus_get_pred_transcripts(AgnPairwiseCompareLocus *locus)
{
  if(locus->pred_genes == NULL)
    return NULL;

  GtArray *transcripts = gt_array_new( sizeof(GtFeatureNode *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->pred_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
        gt_array_add(transcripts, feature);
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_sort(transcripts, (GtCompare)agn_gt_genome_node_compare);
  return transcripts;
}

GtArray *agn_pairwise_compare_locus_get_pred_transcript_ids(AgnPairwiseCompareLocus *locus)
{
  if(locus->pred_genes == NULL)
    return NULL;

  GtArray *ids = gt_array_new( sizeof(char *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->pred_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
      {
        const char *id = gt_feature_node_get_attribute(feature, "ID");
        gt_array_add(ids, id);
      }
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_sort(ids, (GtCompare)agn_string_compare);
  return ids;
}

GtArray *agn_pairwise_compare_locus_get_refr_genes(AgnPairwiseCompareLocus *locus)
{
  if(locus->refr_genes == NULL)
    return NULL;

  GtArray *genes = gt_array_new( sizeof(GtFeatureNode *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->refr_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    gt_array_add(genes, gene);
  }
  gt_array_sort(genes, (GtCompare)gt_genome_node_compare);
  return genes;
}

GtArray *agn_pairwise_compare_locus_get_refr_gene_ids(AgnPairwiseCompareLocus *locus)
{
  if(locus->refr_genes == NULL)
    return NULL;

  GtArray *ids = gt_array_new( sizeof(char *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->refr_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    const char *id = gt_feature_node_get_attribute(gene, "ID");
    gt_array_add(ids, id);
  }
  gt_array_sort(ids, (GtCompare)agn_string_compare);
  return ids;
}

double agn_pairwise_compare_locus_get_refr_splice_complexity(AgnPairwiseCompareLocus *locus)
{
  return locus->refr_splice_complexity;
}

GtArray *agn_pairwise_compare_locus_get_refr_transcripts(AgnPairwiseCompareLocus *locus)
{
  if(locus->refr_genes == NULL)
    return NULL;

  GtArray *transcripts = gt_array_new( sizeof(GtFeatureNode *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->refr_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
        gt_array_add(transcripts, feature);
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_sort(transcripts, (GtCompare)gt_genome_node_compare);
  return transcripts;
}

GtArray *agn_pairwise_compare_locus_get_refr_transcript_ids(AgnPairwiseCompareLocus *locus)
{
  if(locus->refr_genes == NULL)
    return NULL;

  GtArray *ids = gt_array_new( sizeof(char *) );
  GtDlistelem *elem;
  for( elem = gt_dlist_first(locus->refr_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
      {
        const char *id = gt_feature_node_get_attribute(feature, "ID");
        gt_array_add(ids, id);
      }
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_sort(ids, (GtCompare)agn_string_compare);
  return ids;
}

const char* agn_pairwise_compare_locus_get_seqid(AgnPairwiseCompareLocus *locus)
{
  return locus->seqid;
}

unsigned long agn_pairwise_compare_locus_get_start(AgnPairwiseCompareLocus *locus)
{
  return locus->range.start;
}

void agn_pairwise_compare_locus_get_summary_data( AgnPairwiseCompareLocus *locus,
                                     AgnSummaryData *data )
{
  data->counts  = locus->counts;
  data->stats   = locus->stats;
  data->results = locus->results;
}

GtArray *agn_pairwise_compare_locus_get_unique_pred_cliques(AgnPairwiseCompareLocus *locus)
{
  return locus->unique_pred_cliques;
}

GtArray *agn_pairwise_compare_locus_get_unique_refr_cliques(AgnPairwiseCompareLocus *locus)
{
  return locus->unique_refr_cliques;
}

bool agn_pairwise_compare_locus_is_complex(AgnPairwiseCompareLocus *locus)
{
  return gt_array_size(locus->refr_cliques) > 1 ||
         gt_array_size(locus->pred_cliques) > 1;
}

AgnPairwiseCompareLocus* agn_pairwise_compare_locus_new(const char *seqid)
{
  AgnPairwiseCompareLocus *locus = (AgnPairwiseCompareLocus *)
                                     gt_malloc(sizeof(AgnPairwiseCompareLocus));
  locus->refr_genes = NULL;
  locus->pred_genes = NULL;
  locus->range.start = 0;
  locus->range.end = 0;
  locus->clique_pairs = NULL;
  locus->clique_pairs_formed = false;
  locus->refr_cliques = NULL;
  locus->pred_cliques = NULL;
  locus->reported_pairs = NULL;
  locus->unique_refr_cliques = NULL;
  locus->unique_pred_cliques = NULL;

  locus->seqid = (char *)gt_malloc(sizeof(char)*(strlen(seqid) + 1));
  strcpy(locus->seqid, seqid);

  agn_comparison_counts_init(&locus->counts);
  agn_comparison_stats_init(&locus->stats);
  agn_compare_class_agg_desc_init(&locus->results);
  locus->refr_splice_complexity = 0.0;
  locus->pred_splice_complexity = 0.0;

  return locus;
}

unsigned long agn_pairwise_compare_locus_num_pred_exons(AgnPairwiseCompareLocus *locus)
{
  if(locus->pred_genes == NULL)
    return 0;

  GtArray *transcripts = agn_pairwise_compare_locus_get_pred_transcripts(locus);
  unsigned long exon_count = 0, i;
  for(i = 0; i < gt_array_size(transcripts); i++)
  {
    GtFeatureNode *transcript = *(GtFeatureNode **)gt_array_get(transcripts, i);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_exon_feature(feature))
        exon_count++;
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_delete(transcripts);
  return exon_count;
}

unsigned long agn_pairwise_compare_locus_num_pred_genes(AgnPairwiseCompareLocus *locus)
{
  if(locus->pred_genes == NULL)
    return 0;

  return gt_dlist_size(locus->pred_genes);
}

unsigned long agn_pairwise_compare_locus_num_pred_transcripts(AgnPairwiseCompareLocus *locus)
{
  GtDlistelem *elem;

  if(locus->pred_genes == NULL)
    return 0;

  unsigned long transcript_count = 0;
  for( elem = gt_dlist_first(locus->pred_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = (GtFeatureNode *)gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
        transcript_count++;
    }
    gt_feature_node_iterator_delete(iter);
  }
  return transcript_count;
}

unsigned long agn_pairwise_compare_locus_num_refr_exons(AgnPairwiseCompareLocus *locus)
{
  if(locus->refr_genes == NULL)
    return 0;

  GtArray *transcripts = agn_pairwise_compare_locus_get_refr_transcripts(locus);
  unsigned long exon_count = 0, i;
  for(i = 0; i < gt_array_size(transcripts); i++)
  {
    GtFeatureNode *transcript = *(GtFeatureNode **)gt_array_get(transcripts, i);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_exon_feature(feature))
        exon_count++;
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_delete(transcripts);
  return exon_count;
}

unsigned long agn_pairwise_compare_locus_num_refr_genes(AgnPairwiseCompareLocus *locus)
{
  if(locus->refr_genes == NULL)
    return 0;

  return gt_dlist_size(locus->refr_genes);
}

unsigned long agn_pairwise_compare_locus_num_refr_transcripts(AgnPairwiseCompareLocus *locus)
{
  GtDlistelem *elem;

  if(locus->refr_genes == NULL)
    return 0;

  unsigned long transcript_count = 0;
  for( elem = gt_dlist_first(locus->refr_genes);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *gene = (GtFeatureNode *)gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for( feature = gt_feature_node_iterator_next(iter);
         feature != NULL;
         feature = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
        transcript_count++;
    }
    gt_feature_node_iterator_delete(iter);
  }
  return transcript_count;
}

#ifndef WITHOUT_CAIRO
void agn_pairwise_compare_locus_png_track_selector(GtBlock *block, GtStr *track, void *data)
{
  GtFeatureNode *fn = gt_block_get_top_level_feature(block);
  GtGenomeNode *gn = (GtGenomeNode *)fn;
  const char *filename = gt_genome_node_get_filename(gn);
  AgnPairwiseCompareLocusPngMetadata *metadata =
                                 (AgnPairwiseCompareLocusPngMetadata *)data;
  char trackname[512];

  if(strcmp(filename, metadata->refrfile) == 0)
  {
    if(strcmp(metadata->refrlabel, "") == 0)
    {
      sprintf(trackname, "Reference annotations (%s)", metadata->refrfile);
      gt_str_set(track, trackname);
    }
    else
    {
      sprintf(trackname, "%s (Reference)", metadata->refrlabel);
      gt_str_set(track, trackname);
    }
    return;
  }
  else if(strcmp(filename, metadata->predfile) == 0)
  {
    if(strcmp(metadata->predlabel, "") == 0)
    {
      sprintf(trackname, "Prediction annotations (%s)", metadata->predfile);
      gt_str_set(track, trackname);
    }
    else
    {
      sprintf(trackname, "%s (Prediction)", metadata->predlabel);
      gt_str_set(track, trackname);
    }
    return;
  }
  else
  {
    fprintf(stderr, "Error: unknown filename '%s'\n", filename);
    exit(1);
  }
}

void agn_pairwise_compare_locus_print_png
(
  AgnPairwiseCompareLocus *locus,
  AgnPairwiseCompareLocusPngMetadata *metadata
)
{
  // Build a gene index
  GtError *error = gt_error_new();
  GtFeatureIndex *index = gt_feature_index_memory_new();
  GtDlistelem *elem;
  if(locus->refr_genes != NULL)
  {
    for( elem = gt_dlist_first(locus->refr_genes);
         elem != NULL;
         elem = gt_dlistelem_next(elem) )
    {
      GtFeatureNode *gene = gt_dlistelem_get_data(elem);
      gt_feature_index_add_feature_node(index, gene, error);
      if(gt_error_is_set(error))
      {
        fprintf(stderr, "error: %s\n", gt_error_get(error));
        exit(1);
      }
    }
  }
  if(locus->pred_genes != NULL)
  {
    for( elem = gt_dlist_first(locus->pred_genes);
         elem != NULL;
         elem = gt_dlistelem_next(elem) )
    {
      GtFeatureNode *gene = gt_dlistelem_get_data(elem);
      gt_feature_index_add_feature_node(index, gene, error);
      if(gt_error_is_set(error))
      {
        fprintf(stderr, "error: %s\n", gt_error_get(error));
        exit(1);
      }
    }
  }

  // Generate the graphic...this is going to get a bit hairy
  GtStyle *style;
  if(!(style = gt_style_new(error)))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  if(gt_style_load_file(style, metadata->stylefile, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  GtDiagram *diagram = gt_diagram_new( index, locus->seqid, &locus->range,
                                       style, error );
  gt_diagram_set_track_selector_func
  (
    diagram,
    (GtTrackSelectorFunc)agn_pairwise_compare_locus_png_track_selector,
    metadata
  );
  GtLayout *layout = gt_layout_new(diagram, metadata->graphic_width, style, error);
  if(!layout)
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  gt_layout_set_track_ordering_func( layout,
                                     (GtTrackOrderingFunc)metadata->track_order_func,
                                     NULL );
  unsigned long image_height;
  if(gt_layout_get_height(layout, &image_height, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  GtCanvas *canvas = gt_canvas_cairo_file_new( style, GT_GRAPHICS_PNG,
                                               metadata->graphic_width, image_height, NULL,
                                               error );
  if(!canvas)
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  if(gt_layout_sketch(layout, canvas, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  if(gt_canvas_cairo_file_to_file((GtCanvasCairoFile*) canvas, metadata->filename, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }

  gt_feature_index_delete(index);
  gt_canvas_delete(canvas);
  gt_layout_delete(layout);
  gt_diagram_delete(diagram);
  gt_style_delete(style);
  gt_error_delete(error);
}
#endif

unsigned long agn_pairwise_compare_locus_pred_cds_length(AgnPairwiseCompareLocus *locus)
{
  unsigned long length = 0, i;
  GtArray *transcripts = agn_pairwise_compare_locus_get_pred_transcripts(locus);
  for(i = 0; i < gt_array_size(transcripts); i++)
  {
    GtFeatureNode *transcript = *(GtFeatureNode **)gt_array_get(transcripts, i);
    length = agn_gt_feature_node_cds_length(transcript);
  }
  gt_array_delete(transcripts);
  return length;
}

unsigned long agn_pairwise_compare_locus_refr_cds_length(AgnPairwiseCompareLocus *locus)
{
  unsigned long length = 0, i;
  GtArray *transcripts = agn_pairwise_compare_locus_get_refr_transcripts(locus);
  for(i = 0; i < gt_array_size(transcripts); i++)
  {
    GtFeatureNode *transcript = *(GtFeatureNode **)gt_array_get(transcripts, i);
    length = agn_gt_feature_node_cds_length(transcript);
  }
  gt_array_delete(transcripts);
  return length;
}

void agn_pairwise_compare_locus_to_gff3(AgnPairwiseCompareLocus *locus, FILE *outstream)
{
  float score;
  if(locus->counts.non_match > 0)
    score = 0.0;
  else if(locus->counts.num_utr_match > 0)
    score = 0.25;
  else if(locus->counts.num_exon_match > 0)
    score = 0.5;
  else if(locus->counts.num_cds_match > 0)
    score = 0.75;
  else if(locus->counts.num_perfect + locus->counts.num_mislabeled > 0)
    score = 1.0;
  else
    score = 0.0;

  fprintf
  (
    outstream,
    "%s\t%s\t%s\t%lu\t%lu\t%.2f\t%c\t.\tperfect_matches=%d;cds_matches=%d;"
    "exon_matches=%d;utr_matches=%d;non_matches=%d\n",
    agn_pairwise_compare_locus_get_seqid(locus),
    "AEGeAn",
    "locus",
    agn_pairwise_compare_locus_get_start(locus),
    agn_pairwise_compare_locus_get_end(locus),
    score,
    '.',
    locus->counts.num_perfect + locus->counts.num_mislabeled,
    locus->counts.num_cds_match,
    locus->counts.num_exon_match,
    locus->counts.num_utr_match,
    locus->counts.non_match
  );
}

void agn_pairwise_compare_locus_update_range( AgnPairwiseCompareLocus *locus,
                                 GtFeatureNode *gene_feature )
{
  GtGenomeNode *gn = (GtGenomeNode *)gene_feature;
  GtRange gene_range = gt_genome_node_get_range(gn);
  if(locus->range.start == 0 && locus->range.end == 0)
    locus->range = gene_range;
  else
    locus->range = gt_range_join(&locus->range, &gene_range);
}
