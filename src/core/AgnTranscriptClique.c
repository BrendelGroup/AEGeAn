#include "AgnGtExtensions.h"
#include "AgnTranscriptClique.h"
#include "AgnUtils.h"

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnTranscriptClique
{
  GtDlist *transcripts;
  GtDlistelem *current;
};


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
void agn_transcript_clique_add( AgnTranscriptClique *clique,
                                GtFeatureNode *transcript )
{
  if(clique->transcripts == NULL)
    clique->transcripts = gt_dlist_new( (GtCompare)gt_genome_node_cmp );

  gt_dlist_add(clique->transcripts, transcript);
  agn_transcript_clique_reset(clique);
}

unsigned long agn_transcript_clique_cds_length(AgnTranscriptClique *clique)
{
  unsigned long length = 0;
  GtDlistelem *elem;
  for( elem = gt_dlist_first(clique->transcripts);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    length += agn_gt_feature_node_cds_length(transcript);
  }
  return length;
}

AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique)
{
  AgnTranscriptClique *new = agn_transcript_clique_new();
  GtFeatureNode *transcript;
  while((transcript = agn_transcript_clique_next(clique)))
    agn_transcript_clique_add(new, transcript);

  return new;
}

void agn_transcript_clique_delete(AgnTranscriptClique *clique)
{
  if(clique->transcripts != NULL)
    gt_dlist_delete(clique->transcripts);
  gt_free(clique);
  clique = NULL;
}

bool agn_transcript_clique_has_id_in_hash( AgnTranscriptClique *clique,
                                           GtHashmap *map )
{
  GtDlistelem *elem;
  for( elem = gt_dlist_first(clique->transcripts);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    const char *tid = gt_feature_node_get_attribute(transcript, "ID");
    if(gt_hashmap_get(map, tid) != NULL)
      return true;
  }
  
  return false;
}

AgnTranscriptClique* agn_transcript_clique_new()
{
  AgnTranscriptClique *clique = (AgnTranscriptClique *)gt_malloc
  (
    sizeof(AgnTranscriptClique)
  );
  clique->transcripts = NULL;
  agn_transcript_clique_reset(clique);
  return clique;
}

GtFeatureNode* agn_transcript_clique_next(AgnTranscriptClique *clique)
{
  if(clique->current == NULL)
  {
    agn_transcript_clique_reset(clique);
    return NULL;
  }

  GtFeatureNode *transcript = gt_dlistelem_get_data(clique->current);
  clique->current = gt_dlistelem_next(clique->current);

  return transcript;
}

unsigned long agn_transcript_clique_num_exons(AgnTranscriptClique *clique)
{
  unsigned long exon_count = 0;
  GtDlistelem *elem;
  for( elem = gt_dlist_first(clique->transcripts);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
    GtFeatureNode *current;
    for( current = gt_feature_node_iterator_next(iter);
         current != NULL;
         current = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_exon_feature(current))
        exon_count++;
    }
    gt_feature_node_iterator_delete(iter);
  }
  return exon_count;
}

unsigned long agn_transcript_clique_num_utrs(AgnTranscriptClique *clique)
{
  unsigned long utr_count = 0;
  GtDlistelem *elem;
  for( elem = gt_dlist_first(clique->transcripts);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
    GtFeatureNode *current;
    for( current = gt_feature_node_iterator_next(iter);
         current != NULL;
         current = gt_feature_node_iterator_next(iter) )
    {
      if(agn_gt_feature_node_is_utr_feature(current))
        utr_count++;
    }
    gt_feature_node_iterator_delete(iter);
  }
  return utr_count;
}

void agn_transcript_clique_print_ids( AgnTranscriptClique *clique,
                                      FILE *outstream )
{
  if(clique->transcripts == NULL || gt_dlist_size(clique->transcripts) == 0)
  {
    fprintf(outstream, "[]");
  }
  else if(gt_dlist_size(clique->transcripts) == 1)
  {
    GtDlistelem *elem = gt_dlist_first(clique->transcripts);
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    fprintf(outstream, "%s", gt_feature_node_get_attribute(transcript, "ID"));
  }
  else
  {
    fprintf(outstream, "[");
    GtDlistelem *elem;
    for( elem = gt_dlist_first(clique->transcripts);
         elem != NULL;
         elem = gt_dlistelem_next(elem) )
    {
      if(elem != gt_dlist_first(clique->transcripts))
        fprintf(outstream, ", ");

      GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
      fprintf(outstream, "%s", gt_feature_node_get_attribute(transcript, "ID"));
    }
    fprintf(outstream, "]");
  }
}

void agn_transcript_clique_put_ids_in_hash( AgnTranscriptClique *clique,
                                            GtHashmap *map )
{
  GtDlistelem *elem;
  for( elem = gt_dlist_first(clique->transcripts);
       elem != NULL;
       elem = gt_dlistelem_next(elem) )
  {
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    const char *tid = gt_feature_node_get_attribute(transcript, "ID");
    gt_assert(gt_hashmap_get(map, tid) == NULL);
    gt_hashmap_add(map, (char *)tid, (char *)tid);
  }
}

void agn_transcript_clique_reset(AgnTranscriptClique *clique)
{
  if(clique->transcripts == NULL)
    clique->current = NULL;
  else
    clique->current = gt_dlist_first(clique->transcripts);
}

unsigned long agn_transcript_clique_size(AgnTranscriptClique *clique)
{
  if(clique->transcripts == NULL)
    return 0;

  return gt_dlist_size(clique->transcripts);
}

GtArray* agn_transcript_clique_to_array(AgnTranscriptClique *clique)
{
  GtArray *new = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNode *transcript;
  while((transcript = agn_transcript_clique_next(clique)))
    gt_array_add(new, transcript);

  return new;
}
