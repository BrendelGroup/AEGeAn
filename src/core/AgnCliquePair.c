#include <math.h>
#include "AgnCliquePair.h"
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

//----------------------------------------------------------------------------//
// Data structure definitions
//----------------------------------------------------------------------------//
struct AgnCliquePair
{
  AgnSequenceRegion region;
  AgnTranscriptClique *refr_clique;
  AgnTranscriptClique *pred_clique;
  char *refr_vector;
  char *pred_vector;
  AgnComparison stats;
};

typedef struct
{
  char *modelvector;
  GtRange *locusrange;
} ModelVectorData;


//----------------------------------------------------------------------------//
// Prototypes for private method(s)
//----------------------------------------------------------------------------//

/**
 * Add a transcript and all its features to a model vector.
 *
 * @param[in]  transcript     a transcript in the clique
 * @param[out] modelvector    ModelVectorData for the vector being built
 */
void clique_pair_add_transcript_to_vector(GtFeatureNode *transcript,
                                          void *data);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
void clique_pair_add_transcript_to_vector(GtFeatureNode *transcript,
                                          void *data)
{
  ModelVectorData *mvd = data;
  GtFeatureNode *fn;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
  for(fn = gt_feature_node_iterator_next(iter);
      fn != NULL;
      fn = gt_feature_node_iterator_next(iter))
  {
    char c;
    if(agn_gt_feature_node_is_cds_feature(fn))
      c = 'C';
    else if(agn_gt_feature_node_is_utr_feature(fn))
    {
      gt_assert(gt_feature_node_has_type(fn, "five_prime_UTR") ||
                gt_feature_node_has_type(fn, "three_prime_UTR"));
      if(gt_feature_node_has_type(fn, "five_prime_UTR"))
        c = 'F';
      else
        c = 'T';
    }
    else if(agn_gt_feature_node_is_intron_feature(fn))
      c = 'I';
    else
      c = 'G';

    if(c != 'G')
    {
      unsigned long fn_start = gt_genome_node_get_start((GtGenomeNode *)fn);
      unsigned long fn_end = gt_genome_node_get_end((GtGenomeNode *)fn);
      unsigned long i;

      for(i = fn_start - mvd->locusrange->start;
          i < fn_end - mvd->locusrange->start + 1;
          i++)
      {
        mvd->modelvector[i] = c;
      }
    }
  }
  gt_feature_node_iterator_delete(iter);
}

void agn_clique_pair_build_model_vectors(AgnCliquePair *pair)
{
  int vector_length = agn_clique_pair_length(pair) + 1;
  pair->refr_vector = (char *)gt_malloc(sizeof(char) * (vector_length));
  pair->pred_vector = (char *)gt_malloc(sizeof(char) * (vector_length));

  int i;
  pair->refr_vector[vector_length - 1] = '\0';
  pair->pred_vector[vector_length - 1] = '\0';
  for(i = 0; i < vector_length - 1; i++)
  {
    pair->refr_vector[i] = 'G';
    pair->pred_vector[i] = 'G';
  }

  ModelVectorData data = { pair->refr_vector, &pair->region.range };
  agn_transcript_clique_traverse(pair->refr_clique,
      (AgnCliqueVisitFunc)clique_pair_add_transcript_to_vector, &data);
  data.modelvector = pair->pred_vector;
  agn_transcript_clique_traverse(pair->pred_clique,
      (AgnCliqueVisitFunc)clique_pair_add_transcript_to_vector, &data);
}

void agn_clique_pair_comparative_analysis(AgnCliquePair *pair)
{
  unsigned long i, j;
  unsigned long locus_length = agn_clique_pair_length(pair);

  // Counts for exon structure
  int num_refr_exons = 0;
  int num_pred_exons = 0;
  unsigned long refr_exon_starts[MAX_EXONS];
  unsigned long pred_exon_starts[MAX_EXONS];
  unsigned long refr_exon_ends[MAX_EXONS];
  unsigned long pred_exon_ends[MAX_EXONS];

  // Counts for CDS structure
  int num_refr_cdss = 0;
  int num_pred_cdss = 0;
  unsigned long refr_cds_starts[MAX_EXONS];
  unsigned long pred_cds_starts[MAX_EXONS];
  unsigned long refr_cds_ends[MAX_EXONS];
  unsigned long pred_cds_ends[MAX_EXONS];

  // Counts for UTR structure
  int num_refr_utrs = 0;
  int num_pred_utrs = 0;
  unsigned long refr_utr_starts[MAX_UTRS];
  unsigned long pred_utr_starts[MAX_UTRS];
  unsigned long refr_utr_ends[MAX_UTRS];
  unsigned long pred_utr_ends[MAX_UTRS];

  // Collect counts
  for(i = 0; i < locus_length; i++)
  {
    // Coding nucleotide counts
    if(pair->refr_vector[i] == 'C' && pair->pred_vector[i] == 'C')
      pair->stats.cds_nuc_stats.tp++;
    else if(pair->refr_vector[i] == 'C' && pair->pred_vector[i] != 'C')
      pair->stats.cds_nuc_stats.fn++;
    else if(pair->refr_vector[i] != 'C' && pair->pred_vector[i] == 'C')
      pair->stats.cds_nuc_stats.fp++;
    else if(pair->refr_vector[i] != 'C' && pair->pred_vector[i] != 'C')
      pair->stats.cds_nuc_stats.tn++;

    // UTR nucleotide counts
    bool refr_utr = (pair->refr_vector[i] == 'F' || pair->refr_vector[i] == 'T');
    bool pred_utr = (pair->pred_vector[i] == 'F' || pair->pred_vector[i] == 'T');
    if(refr_utr && pred_utr)        pair->stats.utr_nuc_stats.tp++;
    else if(refr_utr && !pred_utr)  pair->stats.utr_nuc_stats.fn++;
    else if(!refr_utr && pred_utr)  pair->stats.utr_nuc_stats.fp++;
    else if(!refr_utr && !pred_utr) pair->stats.utr_nuc_stats.tn++;

    // Overall matches
    if(pair->refr_vector[i] == pair->pred_vector[i])
      pair->stats.overall_matches++;

    // CDS structure counts
    if(pair->refr_vector[i] == 'C')
    {
      if(i == 0 || pair->refr_vector[i-1] != 'C')
      {
        refr_cds_starts[num_refr_cdss] = i;
      }

      if(i == locus_length - 1 || pair->refr_vector[i+1] != 'C')
      {
        refr_cds_ends[num_refr_cdss] = i;
        num_refr_cdss++;
      }
    }
    if(pair->pred_vector[i] == 'C')
    {
      if(i == 0 || pair->pred_vector[i-1] != 'C')
      {
        pred_cds_starts[num_pred_cdss] = i;
      }

      if(i == locus_length - 1 || pair->pred_vector[i+1] != 'C')
      {
        pred_cds_ends[num_pred_cdss] = i;
        num_pred_cdss++;
      }
    }

    // Exon structure counts
    if
    (
      pair->refr_vector[i] == 'C' ||
      pair->refr_vector[i] == 'F' ||
      pair->refr_vector[i] == 'T'
    )
    {
      if( i == 0 || ( pair->refr_vector[i-1] != 'C' &&
                      pair->refr_vector[i-1] != 'F' &&
                      pair->refr_vector[i-1] != 'T'    ))
      {
        refr_exon_starts[num_refr_exons] = i;
      }
      else if( i == locus_length - 1 || ( pair->refr_vector[i+1] != 'C' &&
                                          pair->refr_vector[i+1] != 'F' &&
                                          pair->refr_vector[i+1] != 'T'    ))
      {
        refr_exon_ends[num_refr_exons] = i;
        num_refr_exons++;
      }
    }
    if
    (
      pair->pred_vector[i] == 'C' ||
      pair->pred_vector[i] == 'F' ||
      pair->pred_vector[i] == 'T'
    )
    {
      if( i == 0 || ( pair->pred_vector[i-1] != 'C' &&
                      pair->pred_vector[i-1] != 'F' &&
                      pair->pred_vector[i-1] != 'T'    ))
      {
        pred_exon_starts[num_pred_exons] = i;
      }

      if( i == locus_length - 1 || ( pair->pred_vector[i+1] != 'C' &&
                                     pair->pred_vector[i+1] != 'F' &&
                                     pair->pred_vector[i+1] != 'T'    ))
      {
        pred_exon_ends[num_pred_exons] = i;
        num_pred_exons++;
      }
    }

    // UTR structure counts
    if(pair->refr_vector[i] == 'F' || pair->refr_vector[i] == 'T')
    {
      if( i == 0 || (pair->refr_vector[i-1] != 'F' && pair->refr_vector[i-1] != 'T' ) )
      {
        refr_utr_starts[num_refr_utrs] = i;
      }

      if(i == locus_length - 1 || ( pair->refr_vector[i+1] != 'F' &&
                                    pair->refr_vector[i+1] != 'T' ))
      {
        refr_utr_ends[num_refr_utrs] = i;
        num_refr_utrs++;
      }
    }
    if(pair->pred_vector[i] == 'F' || pair->pred_vector[i] == 'T')
    {
      if( i == 0 || ( pair->pred_vector[i-1] != 'F' &&
                      pair->pred_vector[i-1] != 'T' ))
      {
        pred_utr_starts[num_pred_utrs] = i;
      }

      if(i == locus_length - 1 || ( pair->pred_vector[i+1] != 'F' &&
                                    pair->pred_vector[i+1] != 'T' ))
      {
        pred_utr_ends[num_pred_utrs] = i;
        num_pred_utrs++;
      }
    }
  }

  // Calculate nucleotide-level statistics from counts
  agn_comp_stats_scaled_resolve(&pair->stats.cds_nuc_stats);
  agn_comp_stats_scaled_resolve(&pair->stats.utr_nuc_stats);
  pair->stats.overall_identity = pair->stats.overall_matches / (double)locus_length;

  // Calculate statistics for CDS structure from counts
  for(i = 0; i < num_refr_cdss; i++)
  {
    bool found_match = 0;
    for(j = 0; j < num_pred_cdss; j++)
    {
      if(refr_cds_starts[i] == pred_cds_starts[j] && refr_cds_ends[i] == pred_cds_ends[j])
      {
        pair->stats.cds_struc_stats.correct++;
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      pair->stats.cds_struc_stats.missing++;
  }
  for(i = 0; i < num_pred_cdss; i++)
  {
    bool found_match = 0;
    for(j = 0; j < num_refr_cdss; j++)
    {
      if(refr_cds_starts[j] == pred_cds_starts[i] && refr_cds_ends[j] == pred_cds_ends[i])
      {
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      pair->stats.cds_struc_stats.wrong++;
  }
  agn_comp_stats_binary_resolve(&pair->stats.cds_struc_stats);

  // Calculate statistics for exon structure from counts
  for(i = 0; i < num_refr_exons; i++)
  {
    bool found_match = 0;
    for(j = 0; j < num_pred_exons; j++)
    {
      if(refr_exon_starts[i] == pred_exon_starts[j] && refr_exon_ends[i] == pred_exon_ends[j])
      {
        pair->stats.exon_struc_stats.correct++;
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      pair->stats.exon_struc_stats.missing++;
  }
  for(i = 0; i < num_pred_exons; i++)
  {
    bool found_match = 0;
    for(j = 0; j < num_refr_exons; j++)
    {
      if(refr_exon_starts[j] == pred_exon_starts[i] && refr_exon_ends[j] == pred_exon_ends[i])
      {
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      pair->stats.exon_struc_stats.wrong++;
  }
  agn_comp_stats_binary_resolve(&pair->stats.exon_struc_stats);

  // Calculate statistics for UTR structure from counts
  for(i = 0; i < num_refr_utrs; i++)
  {
    bool found_match = 0;
    for(j = 0; j < num_pred_utrs; j++)
    {
      if(refr_utr_starts[i] == pred_utr_starts[j] && refr_utr_ends[i] == pred_utr_ends[j])
      {
        pair->stats.utr_struc_stats.correct++;
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      pair->stats.utr_struc_stats.missing++;
  }
  for(i = 0; i < num_pred_utrs; i++)
  {
    bool found_match = 0;
    for(j = 0; j < num_refr_utrs; j++)
    {
      if(refr_utr_starts[j] == pred_utr_starts[i] && refr_utr_ends[j] == pred_utr_ends[i])
      {
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      pair->stats.utr_struc_stats.wrong++;
  }
  agn_comp_stats_binary_resolve(&pair->stats.utr_struc_stats);
}

AgnCliquePairClassification agn_clique_pair_classify(AgnCliquePair *pair)
{
  if(pair->stats.overall_identity == 1.0 ||
     fabs(pair->stats.overall_identity - 1.0) < pair->stats.tolerance)
  {
    return AGN_CLIQUE_PAIR_PERFECT_MATCH;
  }
  else if(pair->stats.cds_struc_stats.missing == 0 &&
          pair->stats.cds_struc_stats.wrong   == 0)
  {
    if(pair->stats.exon_struc_stats.missing == 0 &&
       pair->stats.exon_struc_stats.wrong   == 0)
    {
      return AGN_CLIQUE_PAIR_MISLABELED;
    }
    else
    {
      return AGN_CLIQUE_PAIR_CDS_MATCH;
    }
  }
  else
  {
    if(pair->stats.exon_struc_stats.missing == 0 &&
       pair->stats.exon_struc_stats.wrong   == 0)
    {
      return AGN_CLIQUE_PAIR_EXON_MATCH;
    }
    else if(pair->stats.utr_struc_stats.missing == 0 &&
            pair->stats.utr_struc_stats.wrong   == 0 &&
            agn_clique_pair_has_utrs(pair))
    {
      return AGN_CLIQUE_PAIR_UTR_MATCH;
    }
  }

  return AGN_CLIQUE_PAIR_NON_MATCH;
}

int agn_clique_pair_compare(void *p1, void *p2)
{
  AgnCliquePair *pair1 = *(AgnCliquePair **)p1;
  AgnCliquePair *pair2 = *(AgnCliquePair **)p2;
  return agn_clique_pair_compare_direct(pair1, pair2);
}

int agn_clique_pair_compare_direct(AgnCliquePair *p1, AgnCliquePair *p2)
{
  // First, check if either pair is a perfect match
  bool p1_perfect = p1->stats.overall_identity == 1.0 ||
                    fabs(p1->stats.overall_identity - 1.0) < p1->stats.tolerance;
  bool p2_perfect = p2->stats.overall_identity == 1.0 ||
                    fabs(p2->stats.overall_identity - 1.0) < p2->stats.tolerance;
  if(p1_perfect && p2_perfect)
    return 0;
  else if(p1_perfect)
    return 1;
  else if(p2_perfect)
    return -1;

  // Now check whether there is a CDS structure match
  bool p1_cds = p1->stats.cds_struc_stats.missing == 0 &&
                p1->stats.cds_struc_stats.wrong   == 0;
  bool p2_cds = p2->stats.cds_struc_stats.missing == 0 &&
                p2->stats.cds_struc_stats.wrong   == 0;
  if(p1_cds && !p2_cds)
    return 1;
  else if(!p1_cds && p2_cds)
    return -1;

  // Now check whether there is an exon structure match
  bool p1_exon = p1->stats.exon_struc_stats.missing == 0 &&
                 p1->stats.exon_struc_stats.wrong   == 0;
  bool p2_exon = p2->stats.exon_struc_stats.missing == 0 &&
                 p2->stats.exon_struc_stats.wrong   == 0;
  if(p1_exon && !p2_exon)
    return 1;
  else if(!p1_exon && p2_exon)
    return -1;

  // At this point, the coding nucleotide correlation coefficient is probably
  // the best measure of similarity
  if(fabs(p1->stats.cds_nuc_stats.cc - p2->stats.cds_nuc_stats.cc) < p1->stats.tolerance)
  {
    if(fabs(p1->stats.utr_nuc_stats.cc - p2->stats.utr_nuc_stats.cc) < p1->stats.tolerance)
    {
      if(p1->stats.overall_identity > p2->stats.overall_identity)
        return 1;
      else if(p1->stats.overall_identity < p2->stats.overall_identity)
        return -1;
      else
        return 0;
    }
    else if(p1->stats.utr_nuc_stats.cc > p2->stats.utr_nuc_stats.cc)
      return 1;
    else
      return -1;
  }
  else if(p1->stats.cds_nuc_stats.cc > p2->stats.cds_nuc_stats.cc)
    return 1;
  else
    return -1;
}

int agn_clique_pair_compare_reverse(void *p1, void *p2)
{
  int result = agn_clique_pair_compare(p1, p2);
  if(result == 1)
    return -1;
  else if(result == -1)
    return 1;

  return 0;
}

void agn_clique_pair_delete(AgnCliquePair *pair)
{
  if(pair->refr_vector != NULL)
    gt_free(pair->refr_vector);
  if(pair->pred_vector != NULL)
    gt_free(pair->pred_vector);
  gt_free(pair);
  pair = NULL;
}

double agn_clique_pair_get_edit_distance(AgnCliquePair *pair)
{
  return pair->stats.exon_struc_stats.ed;
}

AgnTranscriptClique *agn_clique_pair_get_pred_clique(AgnCliquePair *pair)
{
  return pair->pred_clique;
}

const char *agn_clique_pair_get_pred_vector(AgnCliquePair *pair)
{
  return pair->pred_vector;
}

AgnTranscriptClique *agn_clique_pair_get_refr_clique(AgnCliquePair *pair)
{
  return pair->refr_clique;
}

const char *agn_clique_pair_get_refr_vector(AgnCliquePair *pair)
{
  return pair->refr_vector;
}

AgnComparison *agn_clique_pair_get_stats(AgnCliquePair *pair)
{
  return &pair->stats;
}

bool agn_clique_pair_has_utrs(AgnCliquePair *pair)
{
  unsigned long refr_utrs = agn_transcript_clique_num_utrs(pair->refr_clique);
  unsigned long pred_utrs = agn_transcript_clique_num_utrs(pair->pred_clique);
  return(refr_utrs + pred_utrs > 0);
}

bool agn_clique_pair_is_simple(AgnCliquePair *pair)
{
  return( agn_transcript_clique_size(pair->refr_clique) == 1 &&
          agn_transcript_clique_size(pair->pred_clique) == 1 );
}

bool agn_clique_pair_needs_comparison(AgnCliquePair *pair)
{
  return ( agn_transcript_clique_size(pair->refr_clique) > 0 &&
           agn_transcript_clique_size(pair->pred_clique) > 0 );
}

unsigned long agn_clique_pair_length(AgnCliquePair *pair)
{
  return gt_range_length(&pair->region.range);
}

AgnCliquePair* agn_clique_pair_new(const char *seqid,
                                   AgnTranscriptClique *refr_clique,
                                   AgnTranscriptClique *pred_clique,
                                   GtRange *locus_range)
{
  gt_assert(refr_clique != NULL && pred_clique != NULL);

  AgnCliquePair *pair = (AgnCliquePair *)gt_malloc(sizeof(AgnCliquePair));
  pair->region.seqid = (char *)seqid;
  pair->region.range = *locus_range;
  pair->refr_clique = refr_clique;
  pair->pred_clique = pred_clique;

  agn_comparison_init(&pair->stats);
  double perc = 1.0 / (double)gt_range_length(locus_range);
  pair->stats.tolerance = 1.0;
  while(pair->stats.tolerance > perc)
    pair->stats.tolerance /= 10;

  pair->refr_vector = NULL;
  pair->pred_vector = NULL;

  return pair;
}

void agn_clique_pair_record_characteristics(AgnCliquePair *pair,
                                            AgnCompResultDesc *desc)
{
  AgnTranscriptClique *refr = agn_clique_pair_get_refr_clique(pair);
  AgnTranscriptClique *pred = agn_clique_pair_get_pred_clique(pair);

  desc->transcript_count += 1;
  desc->total_length += agn_clique_pair_length(pair);
  desc->refr_cds_length += agn_transcript_clique_cds_length(refr);
  desc->pred_cds_length += agn_transcript_clique_cds_length(pred);
  desc->refr_exon_count += agn_transcript_clique_num_exons(refr);
  desc->pred_exon_count += agn_transcript_clique_num_exons(pred);
}
