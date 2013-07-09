#include <math.h>
#include <string.h>
#include "extended/feature_node.h"
#include "AgnGeneLocus.h"
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnGeneLocus
{
  AgnLoc locus;
  GtDlist *genes;
  GtHashmap *refr_genes;
  GtHashmap *pred_genes;
  GtArray *refr_cliques;
  GtArray *pred_cliques;
  GtArray *clique_pairs;
  GtArray *reported_pairs;
  GtArray *unique_refr_cliques;
  GtArray *unique_pred_cliques;
};


//----------------------------------------------------------------------------//
// Prototypes for private method(s)
//----------------------------------------------------------------------------//

/**
 * Update this locus' start and end coordinates based on the gene being merged.
 *
 * @param[out] locus    the locus
 * @param[in]  gene     gene that was just merged with the locus
 */
void agn_gene_locus_update_range(AgnGeneLocus *locus, GtFeatureNode *gene);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
void agn_gene_locus_add(AgnGeneLocus *locus, GtFeatureNode *gene,
                             AgnComparisonSource source)
{
  gt_dlist_add(locus->genes, gene);
  agn_gene_locus_update_range(locus, gene);
  if(source == REFERENCESOURCE)
    gt_hashmap_add(locus->refr_genes, gene, gene);
  else if(source == PREDICTIONSOURCE)
    gt_hashmap_add(locus->pred_genes, gene, gene);
}

AgnGeneLocus *agn_gene_locus_clone(AgnGeneLocus *locus)
{
  AgnGeneLocus *newlocus = (AgnGeneLocus *)gt_malloc(sizeof(AgnGeneLocus));
  newlocus->locus = locus->locus;
  newlocus->genes = locus->genes;
  newlocus->refr_genes = locus->refr_genes;
  newlocus->pred_genes = locus->pred_genes;
  newlocus->refr_cliques = locus->refr_cliques;
  newlocus->pred_cliques = locus->pred_cliques;
  newlocus->clique_pairs = locus->clique_pairs;
  newlocus->reported_pairs = locus->reported_pairs;
  newlocus->unique_refr_cliques = locus->unique_refr_cliques;
  newlocus->unique_pred_cliques = locus->unique_pred_cliques;
  return newlocus;
}

int agn_gene_locus_array_compare(const void *p1, const void *p2)
{
  AgnGeneLocus *l1 = *(AgnGeneLocus **)p1;
  AgnGeneLocus *l2 = *(AgnGeneLocus **)p2;
  GtRange l1r = l1->locus.range;
  GtRange l2r = l2->locus.range;
  return gt_range_compare(&l1r, &l2r);
}

unsigned long agn_gene_locus_cds_length(AgnGeneLocus *locus,
                                        AgnComparisonSource src)
{
  unsigned long length = 0, i;
  GtArray *transcripts = agn_gene_locus_transcripts(locus, src);
  for(i = 0; i < gt_array_size(transcripts); i++)
  {
    GtFeatureNode *transcript = *(GtFeatureNode **)gt_array_get(transcripts, i);
    length = agn_gt_feature_node_cds_length(transcript);
  }
  gt_array_delete(transcripts);
  return length;
}

GtArray *agn_gene_locus_comparative_analysis(AgnGeneLocus *locus)
{
  if(locus->reported_pairs != NULL)
    return locus->reported_pairs;

  unsigned long i;
  for(i = 0; i < gt_array_size(locus->clique_pairs); i++)
  {
    AgnCliquePair *p = *(AgnCliquePair **)gt_array_get(locus->clique_pairs, i);
    agn_clique_pair_build_model_vectors(p);
    agn_clique_pair_comparative_analysis(p);
  }
  gt_array_sort(locus->clique_pairs,(GtCompare)agn_clique_pair_compare_reverse);

  GtHashmap *refr_cliques_acctd = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  GtHashmap *pred_cliques_acctd = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  locus->reported_pairs = gt_array_new( sizeof(AgnCliquePair *) );
  for(i = 0; i < gt_array_size(locus->clique_pairs); i++)
  {
    AgnCliquePair *pair;
    pair = *(AgnCliquePair **)gt_array_get(locus->clique_pairs, i);
    AgnTranscriptClique *rclique = agn_clique_pair_get_refr_clique(pair);
    AgnTranscriptClique *pclique = agn_clique_pair_get_pred_clique(pair);
    if(!agn_transcript_clique_has_id_in_hash(rclique, refr_cliques_acctd) &&
       !agn_transcript_clique_has_id_in_hash(pclique, pred_cliques_acctd))
    {
      gt_array_add(locus->reported_pairs, pair);
      agn_transcript_clique_put_ids_in_hash(rclique, refr_cliques_acctd);
      agn_transcript_clique_put_ids_in_hash(pclique, pred_cliques_acctd);
    }
  }

  locus->unique_refr_cliques = gt_array_new( sizeof(AgnTranscriptClique *) );
  for(i = 0; i < gt_array_size(locus->refr_cliques); i++)
  {
    AgnTranscriptClique *refr_clique;
    refr_clique = *(AgnTranscriptClique **)gt_array_get(locus->refr_cliques, i);
    if(!agn_transcript_clique_has_id_in_hash(refr_clique, refr_cliques_acctd))
    {
      gt_array_add(locus->unique_refr_cliques, refr_clique);
      agn_transcript_clique_put_ids_in_hash(refr_clique, refr_cliques_acctd);
    }
  }

  locus->unique_pred_cliques = gt_array_new( sizeof(AgnTranscriptClique *) );
  for(i = 0; i < gt_array_size(locus->pred_cliques); i++)
  {
    AgnTranscriptClique *pred_clique;
    pred_clique= *(AgnTranscriptClique **)gt_array_get(locus->pred_cliques, i);
    if(!agn_transcript_clique_has_id_in_hash(pred_clique, pred_cliques_acctd))
    {
      gt_array_add(locus->unique_pred_cliques, pred_clique);
      agn_transcript_clique_put_ids_in_hash(pred_clique, pred_cliques_acctd);
    }
  }

  gt_hashmap_delete(refr_cliques_acctd);
  gt_hashmap_delete(pred_cliques_acctd);
  return locus->reported_pairs;
}

void agn_gene_locus_delete(AgnGeneLocus *locus)
{
  GtDlistelem *current;
  for(current = gt_dlist_first(locus->genes);
      current != NULL;
      current = gt_dlistelem_next(current))
  {
    GtGenomeNode *gn = gt_dlistelem_get_data(current);
    gt_genome_node_delete(gn);
  }
  gt_dlist_delete(locus->genes);
  gt_hashmap_delete(locus->refr_genes);
  gt_hashmap_delete(locus->pred_genes);

  AgnTranscriptClique *clique;
  if(locus->refr_cliques != NULL)
  {
    while(gt_array_size(locus->refr_cliques) > 0)
    {
      clique = *(AgnTranscriptClique **)gt_array_pop(locus->refr_cliques);
      agn_transcript_clique_delete(clique);
    }
    gt_array_delete(locus->refr_cliques);
  }
  if(locus->pred_cliques != NULL)
  {
    while(gt_array_size(locus->pred_cliques) > 0)
    {
      clique = *(AgnTranscriptClique **)gt_array_pop(locus->pred_cliques);
      agn_transcript_clique_delete(clique);
    }
    gt_array_delete(locus->pred_cliques);
  }

  if(locus->clique_pairs != NULL)
  {
    while(gt_array_size(locus->clique_pairs) > 0)
    {
      AgnCliquePair *pair = *(AgnCliquePair**)gt_array_pop(locus->clique_pairs);
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

  gt_free(locus->locus.seqid);
  gt_free(locus);
  locus = NULL;
}

unsigned long agn_gene_locus_enumerate_clique_pairs(AgnGeneLocus *locus)
{
  if(locus->clique_pairs != NULL)
    return gt_array_size(locus->clique_pairs);

  if(agn_gene_locus_num_refr_transcripts(locus) > 0)
  {
    GtArray *refr_trans = agn_gene_locus_refr_transcripts(locus);
    locus->refr_cliques = agn_enumerate_feature_cliques(refr_trans);
    gt_array_delete(refr_trans);
  }
  if(agn_gene_locus_num_pred_transcripts(locus) > 0)
  {
    GtArray *pred_trans = agn_gene_locus_pred_transcripts(locus);
    locus->pred_cliques = agn_enumerate_feature_cliques(pred_trans);
    gt_array_delete(pred_trans);
  }

  locus->clique_pairs = gt_array_new( sizeof(AgnCliquePair *) );

  unsigned long i,j;
  for(i = 0; i < gt_array_size(locus->refr_cliques); i++)
  {
    AgnTranscriptClique *refr_clique, *pred_clique;
    refr_clique = *(AgnTranscriptClique **)gt_array_get(locus->refr_cliques, i);
    for(j = 0; j < gt_array_size(locus->pred_cliques); j++)
    {
      pred_clique = *(AgnTranscriptClique**)gt_array_get(locus->pred_cliques,j);
      AgnCliquePair *pair = agn_clique_pair_new(locus->locus.seqid,
                                                refr_clique, pred_clique,
                                                &locus->locus.range);
      gt_array_add(locus->clique_pairs, pair);
    }
  }

  return gt_array_size(locus->clique_pairs);
}

unsigned long agn_gene_locus_exon_num(AgnGeneLocus *locus,
                                      AgnComparisonSource src)
{
  GtDlistelem *elem;
  unsigned long exon_count = 0;
  for(elem = gt_dlist_first(locus->genes);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *gene = (GtFeatureNode *)gt_dlistelem_get_data(elem);
    bool isrefr = gt_hashmap_get(locus->refr_genes, gene) != NULL;
    bool ispred = gt_hashmap_get(locus->pred_genes, gene) != NULL;

    if(src == DEFAULTSOURCE ||
       (src == REFERENCESOURCE  && isrefr) ||
       (src == PREDICTIONSOURCE && ispred))
    {
      GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(gene);
      GtFeatureNode *feature;
      for(feature = gt_feature_node_iterator_next(iter);
          feature != NULL;
          feature = gt_feature_node_iterator_next(iter))
      {
        if(agn_gt_feature_node_is_exon_feature(feature))
          exon_count++;
      }
      gt_feature_node_iterator_delete(iter);
    }
  }
  return exon_count;
}

bool agn_gene_locus_filter(AgnGeneLocus *locus, AgnCompareFilters *filters)
{
  if(filters == NULL)
    return false;

  // Ignore all filter configuration settings set to 0

  // Filter by locus length
  unsigned long length = agn_gene_locus_get_length(locus);
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
  unsigned long num_refr_genes = agn_gene_locus_num_refr_genes(locus);
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
  unsigned long num_pred_genes = agn_gene_locus_num_pred_genes(locus);
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
  unsigned long nrefr_transcripts = agn_gene_locus_num_refr_transcripts(locus);
  if(filters->MinReferenceTranscriptModels > 0)
  {
    if(nrefr_transcripts < filters->MinReferenceTranscriptModels)
      return true;
  }
  if(filters->MaxReferenceTranscriptModels > 0)
  {
    if(nrefr_transcripts > filters->MaxReferenceTranscriptModels)
      return true;
  }

  // Filter by prediction transcript models
  unsigned long npred_transcripts = agn_gene_locus_num_pred_transcripts(locus);
  if(filters->MinPredictionTranscriptModels > 0)
  {
    if(npred_transcripts < filters->MinPredictionTranscriptModels)
      return true;
  }
  if(filters->MaxPredictionTranscriptModels > 0)
  {
    if(npred_transcripts > filters->MaxPredictionTranscriptModels)
      return true;
  }

  // Filter by number of transcripts per gene model
  GtArray *refr_genes = agn_gene_locus_refr_genes(locus);
  if(filters->MinTranscriptsPerReferenceGeneModel > 0)
  {
    unsigned long i;
    bool success = false;
    for(i = 0; i < gt_array_size(refr_genes); i++)
    {
      GtFeatureNode *gene = gt_array_get(refr_genes, i);
      if(agn_gt_feature_node_num_transcripts(gene) >=
         filters->MinTranscriptsPerReferenceGeneModel)
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
      if(agn_gt_feature_node_num_transcripts(gene) <=
         filters->MaxTranscriptsPerReferenceGeneModel)
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
  GtArray *pred_genes = agn_gene_locus_pred_genes(locus);
  if(filters->MinTranscriptsPerPredictionGeneModel > 0)
  {
    unsigned long i;
    bool success = false;
    for(i = 0; i < gt_array_size(pred_genes); i++)
    {
      GtFeatureNode *gene = gt_array_get(pred_genes, i);
      if(agn_gt_feature_node_num_transcripts(gene) >=
         filters->MinTranscriptsPerPredictionGeneModel)
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
      if(agn_gt_feature_node_num_transcripts(gene) <=
         filters->MaxTranscriptsPerPredictionGeneModel)
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
  unsigned long num_refr_exons = agn_gene_locus_num_refr_exons(locus);
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
  unsigned long num_pred_exons = agn_gene_locus_num_pred_exons(locus);
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
  unsigned long refr_cds_length = agn_gene_locus_refr_cds_length(locus);
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
  unsigned long pred_cds_length = agn_gene_locus_pred_cds_length(locus);
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

GtArray *agn_gene_locus_genes(AgnGeneLocus *locus, AgnComparisonSource src)
{
  GtArray *genes = gt_array_new( sizeof(GtFeatureNode *) );
  GtDlistelem *elem;
  for(elem = gt_dlist_first(locus->genes);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    bool isrefr = gt_hashmap_get(locus->refr_genes, gene) != NULL;
    bool ispred = gt_hashmap_get(locus->pred_genes, gene) != NULL;

    if(src == DEFAULTSOURCE)
      gt_array_add(genes, gene);
    else if(src == REFERENCESOURCE && isrefr)
      gt_array_add(genes, gene);
    else if(src == PREDICTIONSOURCE && ispred)
      gt_array_add(genes, gene);
  }
  gt_array_sort(genes, (GtCompare)agn_gt_genome_node_compare);
  return genes;
}

GtArray *agn_gene_locus_gene_ids(AgnGeneLocus *locus, AgnComparisonSource src)
{
  GtArray *ids = gt_array_new( sizeof(char *) );
  GtDlistelem *elem;
  for(elem = gt_dlist_first(locus->genes);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    bool isrefr = gt_hashmap_get(locus->refr_genes, gene) != NULL;
    bool ispred = gt_hashmap_get(locus->pred_genes, gene) != NULL;
    const char *id = gt_feature_node_get_attribute(gene, "ID");

    if(src == DEFAULTSOURCE)
      gt_array_add(ids, id);
    else if(src == REFERENCESOURCE && isrefr)
      gt_array_add(ids, id);
    else if(src == PREDICTIONSOURCE && ispred)
      gt_array_add(ids, id);
  }
  gt_array_sort(ids, (GtCompare)agn_string_compare);
  return ids;
}

unsigned long agn_gene_locus_gene_num(AgnGeneLocus *locus,
                                      AgnComparisonSource src)
{
  GtDlistelem *elem;
  unsigned long count = 0;
  for(elem = gt_dlist_first(locus->genes);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *gene = (GtFeatureNode *)gt_dlistelem_get_data(elem);
    bool isrefr = gt_hashmap_get(locus->refr_genes, gene) != NULL;
    bool ispred = gt_hashmap_get(locus->pred_genes, gene) != NULL;

    if(src == DEFAULTSOURCE)
      count++;
    else if(src == REFERENCESOURCE && isrefr)
      count++;
    else if(src == PREDICTIONSOURCE && ispred)
      count++;
  }
  return count;
}

unsigned long agn_gene_locus_get_end(AgnGeneLocus *locus)
{
  return locus->locus.range.end;
}

unsigned long agn_gene_locus_get_length(AgnGeneLocus *locus)
{
  return gt_range_length(&locus->locus.range);
}

AgnCliquePair* agn_gene_locus_get_optimal_clique_pair(AgnGeneLocus *locus,
                                              AgnTranscriptClique *refr_clique)
{
  AgnCliquePair *bestpair = NULL;
  unsigned long i;

  for(i = 0; i < gt_array_size(locus->clique_pairs); i++)
  {
    AgnCliquePair *currentpair;
    currentpair = *(AgnCliquePair **)gt_array_get(locus->clique_pairs, i);
    if(refr_clique == agn_clique_pair_get_refr_clique(currentpair))
    {
      if(bestpair == NULL ||
         agn_clique_pair_compare_direct(bestpair, currentpair) == -1)
      {
        bestpair = currentpair;
      }
    }
  }

  return bestpair;
}

const char* agn_gene_locus_get_seqid(AgnGeneLocus *locus)
{
  return locus->locus.seqid;
}

unsigned long agn_gene_locus_get_start(AgnGeneLocus *locus)
{
  return locus->locus.range.start;
}

GtArray *agn_gene_locus_get_unique_pred_cliques(AgnGeneLocus *locus)
{
  return locus->unique_pred_cliques;
}

GtArray *agn_gene_locus_get_unique_refr_cliques(AgnGeneLocus *locus)
{
  return locus->unique_refr_cliques;
}

bool agn_gene_locus_is_complex(AgnGeneLocus *locus)
{
  return gt_array_size(locus->refr_cliques) > 1 ||
         gt_array_size(locus->pred_cliques) > 1;
}

AgnGeneLocus* agn_gene_locus_new(const char *seqid)
{
  AgnGeneLocus *locus = (AgnGeneLocus *)gt_malloc(sizeof(AgnGeneLocus));

  locus->locus.seqid = gt_cstr_dup(seqid);
  locus->locus.range.start = 0;
  locus->locus.range.end = 0;
  locus->genes = gt_dlist_new( (GtCompare)gt_genome_node_cmp );
  locus->refr_genes = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  locus->pred_genes = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  locus->clique_pairs = NULL;
  locus->refr_cliques = NULL;
  locus->pred_cliques = NULL;
  locus->reported_pairs = NULL;
  locus->unique_refr_cliques = NULL;
  locus->unique_pred_cliques = NULL;

  return locus;
}

GtRange agn_gene_locus_range(AgnGeneLocus *locus)
{
  return locus->locus.range;
}

void agn_gene_locus_set_range(AgnGeneLocus *locus, unsigned long start,
                              unsigned long end)
{
  locus->locus.range.start = start;
  locus->locus.range.end   = end;
}

double agn_gene_locus_splice_complexity(AgnGeneLocus *locus,
                                        AgnComparisonSource src)
{
  GtArray *trans = agn_gene_locus_transcripts(locus, src);
  double sc = agn_calc_splice_complexity(trans);
  gt_array_delete(trans);
  return sc;
}

void agn_gene_locus_to_gff3(AgnGeneLocus *locus, FILE *outstream,
                            const char *source)
{
  const char *src = "AEGeAn";
  if(source != NULL)
    src = source;
  fprintf(outstream,
          "%s\t%s\tlocus\t%lu\t%lu\t.\t.\t.\tnum_genes=%lu\n",
          locus->locus.seqid, src, locus->locus.range.start,
          locus->locus.range.end, gt_dlist_size(locus->genes));
}

GtArray *agn_gene_locus_transcripts(AgnGeneLocus *locus,
                                    AgnComparisonSource src)
{
  GtArray *transcripts = gt_array_new( sizeof(GtFeatureNode *) );
  GtDlistelem *elem;
  for(elem = gt_dlist_first(locus->genes);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for(feature = gt_feature_node_iterator_next(iter);
        feature != NULL;
        feature = gt_feature_node_iterator_next(iter))
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
      {
        bool isrefr = gt_hashmap_get(locus->refr_genes, gene) != NULL;
        bool ispred = gt_hashmap_get(locus->pred_genes, gene) != NULL;

        if(src == DEFAULTSOURCE)
          gt_array_add(transcripts, feature);
        else if(src == REFERENCESOURCE && isrefr)
          gt_array_add(transcripts, feature);
        else if(src == PREDICTIONSOURCE && ispred)
          gt_array_add(transcripts, feature);
      }
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_sort(transcripts, (GtCompare)agn_gt_genome_node_compare);
  return transcripts;
}

GtArray *agn_gene_locus_transcript_ids(AgnGeneLocus *locus,
                                    AgnComparisonSource src)
{
  GtArray *ids = gt_array_new( sizeof(GtFeatureNode *) );
  GtDlistelem *elem;
  for(elem = gt_dlist_first(locus->genes);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *gene = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
    GtFeatureNode *feature;

    for(feature = gt_feature_node_iterator_next(iter);
        feature != NULL;
        feature = gt_feature_node_iterator_next(iter))
    {
      if(agn_gt_feature_node_is_mrna_feature(feature))
      {
        bool isrefr = gt_hashmap_get(locus->refr_genes, gene) != NULL;
        bool ispred = gt_hashmap_get(locus->pred_genes, gene) != NULL;
        const char *id = gt_feature_node_get_attribute(feature, "ID");

        if(src == DEFAULTSOURCE)
          gt_array_add(ids, id);
        else if(src == REFERENCESOURCE && isrefr)
          gt_array_add(ids, id);
        else if(src == PREDICTIONSOURCE && ispred)
          gt_array_add(ids, id);
      }
    }
    gt_feature_node_iterator_delete(iter);
  }
  gt_array_sort(ids, (GtCompare)agn_string_compare);
  return ids;
}

unsigned long agn_gene_locus_transcript_num(AgnGeneLocus *locus,
                                            AgnComparisonSource src)
{
  GtDlistelem *elem;
  unsigned long transcript_count = 0;
  for(elem = gt_dlist_first(locus->genes);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *gene = (GtFeatureNode *)gt_dlistelem_get_data(elem);
    bool isrefr = gt_hashmap_get(locus->refr_genes, gene) != NULL;
    bool ispred = gt_hashmap_get(locus->pred_genes, gene) != NULL;

    if(src == DEFAULTSOURCE ||
       (src == REFERENCESOURCE  && isrefr) ||
       (src == PREDICTIONSOURCE && ispred))
    {
      GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(gene);
      GtFeatureNode *feature;
      for(feature = gt_feature_node_iterator_next(iter);
          feature != NULL;
          feature = gt_feature_node_iterator_next(iter))
      {
        if(agn_gt_feature_node_is_mrna_feature(feature))
          transcript_count++;
      }
      gt_feature_node_iterator_delete(iter);
    }
  }
  return transcript_count;
}

void agn_gene_locus_update_range(AgnGeneLocus *locus, GtFeatureNode *gene)
{
  GtRange gene_range = gt_genome_node_get_range((GtGenomeNode *)gene);
  if(locus->locus.range.start == 0 && locus->locus.range.end == 0)
    locus->locus.range = gene_range;
  else
    locus->locus.range = gt_range_join(&locus->locus.range, &gene_range);
}
