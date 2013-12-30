#include <math.h>
#include "core/queue_api.h"
#include "AgnCliquePair.h"

#define char_is_exonic(C) (C == 'F' || C == 'T' || C == 'C')
#define char_is_utric(C)  (C == 'F' || C == 'T')
#define clique_pair_has_utrs(CP) \
        (agn_transcript_clique_num_utrs(CP->refr_clique) + \
         agn_transcript_clique_num_utrs(CP->pred_clique) > 0)

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

struct AgnCliquePair
{
  AgnTranscriptClique *refr_clique;
  AgnTranscriptClique *pred_clique;
  AgnComparison stats;
  double tolerance;
};

typedef struct
{
  GtArray *refrstarts;
  GtArray *refrends;
  GtArray *predstarts;
  GtArray *predends;
  AgnCompStatsBinary *stats;
} StructuralData;


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Given a set of of start and end coordinates for reference and
 * prediction structures (exons, CDS segments, or UTR segments), determine the
 * number of congruent and incongruent structures.
 */
static void clique_pair_calc_struct_stats(StructuralData *dat);

/**
 * @function FIXME
 */
static void clique_pair_comparative_analysis(AgnCliquePair *pair);

/**
 * @function Initialize the data structure used to store start and end
 * coordinates for reference and prediction structures (exons, CDS segments, or
 * UTR segments) and associated statistics.
 */
static void clique_pair_init_struct_dat(StructuralData *dat,
                                        AgnCompStatsBinary *stats);

/**
 * @function Free the memory previously occupied by the data structure.
 */
static void clique_pair_term_struct_dat(StructuralData *dat);

/**
 * @function Generate data for unit testing.
 */
static void clique_pair_test_data(GtQueue *queue);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

AgnCompClassification agn_clique_pair_classify(AgnCliquePair *pair)
{
  double identity = (double)pair->stats.overall_matches /
                    (double)pair->stats.overall_length;
  if(identity == 1.0 || fabs(identity - 1.0) < pair->tolerance)
  {
    return AGN_COMP_CLASS_PERFECT_MATCH;
  }
  else if(pair->stats.cds_struc_stats.missing == 0 &&
          pair->stats.cds_struc_stats.wrong   == 0)
  {
    if(pair->stats.exon_struc_stats.missing == 0 &&
       pair->stats.exon_struc_stats.wrong   == 0)
    {
      return AGN_COMP_CLASS_MISLABELED;
    }
    else
    {
      return AGN_COMP_CLASS_CDS_MATCH;
    }
  }
  else
  {
    if(pair->stats.exon_struc_stats.missing == 0 &&
       pair->stats.exon_struc_stats.wrong   == 0)
    {
      return AGN_COMP_CLASS_EXON_MATCH;
    }
    else if(pair->stats.utr_struc_stats.missing == 0 &&
            pair->stats.utr_struc_stats.wrong   == 0 &&
            clique_pair_has_utrs(pair))
    {
      return AGN_COMP_CLASS_UTR_MATCH;
    }
  }

  return AGN_COMP_CLASS_NON_MATCH;
}

void agn_clique_pair_comparison_aggregate(AgnCliquePair *pair,
                                          AgnComparison *comp)
{
  agn_comparison_aggregate(comp, &pair->stats);
}

int agn_clique_pair_compare(void *p1, void *p2)
{
  AgnCliquePair *pair1 = *(AgnCliquePair **)p1;
  AgnCliquePair *pair2 = *(AgnCliquePair **)p2;
  return agn_clique_pair_compare_direct(pair1, pair2);
}

int agn_clique_pair_compare_direct(AgnCliquePair *p1, AgnCliquePair *p2)
{
  double p1identity = (double)p1->stats.overall_matches /
                      (double)p1->stats.overall_length;
  double p2identity = (double)p1->stats.overall_matches /
                      (double)p1->stats.overall_length;

  // First, check if either pair is a perfect match
  bool p1_perfect = p1identity == 1.0 || fabs(p1identity - 1.0) < p1->tolerance;
  bool p2_perfect = p2identity == 1.0 || fabs(p2identity - 1.0) < p2->tolerance;
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
  double cdiff = fabs(p1->stats.cds_nuc_stats.cc - p2->stats.cds_nuc_stats.cc);
  double udiff = fabs(p1->stats.utr_nuc_stats.cc - p2->stats.utr_nuc_stats.cc);
  if(cdiff < p1->tolerance)
  {
    if(udiff < p1->tolerance)
    {
      if(p1identity > p2identity)
        return 1;
      else if(p1identity < p2identity)
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
  gt_free(pair);
}

AgnCliquePair* agn_clique_pair_new(AgnTranscriptClique *refr,
                                   AgnTranscriptClique *pred)
{
  GtStr *seqidrefr = gt_genome_node_get_seqid(refr);
  GtStr *seqidpred = gt_genome_node_get_seqid(pred);
  gt_assert(gt_genome_node_get_start(refr) == gt_genome_node_get_start(pred) &&
            gt_genome_node_get_end(refr) == gt_genome_node_get_end(pred) &&
            gt_str_cmp(seqidrefr, seqidpred) == 0);

  AgnCliquePair *pair = (AgnCliquePair *)gt_malloc( sizeof(AgnCliquePair) );
  pair->refr_clique = refr;
  pair->pred_clique = pred;

  agn_comparison_init(&pair->stats);
  double perc = 1.0 / (double)gt_genome_node_get_length(refr);
  pair->tolerance = 1.0;
  while(pair->tolerance > perc)
    pair->tolerance /= 10;

  clique_pair_comparative_analysis(pair);
  return pair;
}

bool agn_clique_pair_unit_test(AgnUnitTest *test)
{
  clique_pair_test_data(NULL);
  return false;
}

static void clique_pair_calc_struct_stats(StructuralData *dat)
{
  GtUword num_refr = gt_array_size(dat->refrstarts);
  GtUword num_pred = gt_array_size(dat->predstarts);
  gt_assert(num_refr == gt_array_size(dat->refrends));
  gt_assert(num_pred == gt_array_size(dat->predends));
  GtUword i, j;
  for(i = 0; i < num_refr; i++)
  {
    bool found_match = 0;
    GtUword *refrstart = gt_array_get(dat->refrstarts, i);
    GtUword *refrend   = gt_array_get(dat->refrends,   i);
    for(j = 0; j < num_pred; j++)
    {
      GtUword *predstart = gt_array_get(dat->predstarts, j);
      GtUword *predend   = gt_array_get(dat->predends,   j);
      if(*refrstart == *predstart && *refrend == *predend)
      {
        dat->stats->correct++;
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      dat->stats->missing++;
  }
  for(i = 0; i < num_pred; i++)
  {
    bool found_match = 0;
    GtUword *predstart = gt_array_get(dat->predstarts, i);
    GtUword *predend   = gt_array_get(dat->predends,   i);
    for(j = 0; j < num_refr; j++)
    {
      GtUword *refrstart = gt_array_get(dat->refrstarts, j);
      GtUword *refrend   = gt_array_get(dat->refrends,   j);
      if(*refrstart == *predstart && *refrend == *predend)
      {
        found_match = 1;
        break;
      }
    }
    if(!found_match)
      dat->stats->wrong++;
  }
  agn_comp_stats_binary_resolve(dat->stats);
  clique_pair_term_struct_dat(dat);
}

static void clique_pair_comparative_analysis(AgnCliquePair *pair)
{
  GtUword locus_length = gt_genome_node_get_length(pair->refr_clique);
  char *refr_vector = gt_genome_node_get_user_data(pair->refr_clique,
                                                   "modelvector");
  char *pred_vector = gt_genome_node_get_user_data(pair->pred_clique,
                                                   "modelvector");

  StructuralData cdsstruct;
  clique_pair_init_struct_dat(&cdsstruct, &pair->stats.cds_struc_stats);
  StructuralData exonstruct;
  clique_pair_init_struct_dat(&exonstruct, &pair->stats.exon_struc_stats);
  StructuralData utrstruct;
  clique_pair_init_struct_dat(&utrstruct, &pair->stats.utr_struc_stats);

  // Collect counts
  GtUword i;
  for(i = 0; i < locus_length; i++)
  {
    // Coding nucleotide counts
    if(refr_vector[i] == 'C' && pred_vector[i] == 'C')
      pair->stats.cds_nuc_stats.tp++;
    else if(refr_vector[i] == 'C' && pred_vector[i] != 'C')
      pair->stats.cds_nuc_stats.fn++;
    else if(refr_vector[i] != 'C' && pred_vector[i] == 'C')
      pair->stats.cds_nuc_stats.fp++;
    else if(refr_vector[i] != 'C' && pred_vector[i] != 'C')
      pair->stats.cds_nuc_stats.tn++;

    // UTR nucleotide counts
    bool refr_utr = char_is_utric(refr_vector[i]);
    bool pred_utr = char_is_utric(pred_vector[i]);
    if(refr_utr && pred_utr)        pair->stats.utr_nuc_stats.tp++;
    else if(refr_utr && !pred_utr)  pair->stats.utr_nuc_stats.fn++;
    else if(!refr_utr && pred_utr)  pair->stats.utr_nuc_stats.fp++;
    else if(!refr_utr && !pred_utr) pair->stats.utr_nuc_stats.tn++;

    // Overall matches
    if(refr_vector[i] == pred_vector[i])
      pair->stats.overall_matches++;

    // CDS structure counts
    if(refr_vector[i] == 'C')
    {
      if(i == 0 || refr_vector[i-1] != 'C')
        gt_array_add(cdsstruct.refrstarts, i);

      if(i == locus_length - 1 || refr_vector[i+1] != 'C')
        gt_array_add(cdsstruct.refrends, i);
    }
    if(pred_vector[i] == 'C')
    {
      if(i == 0 || pred_vector[i-1] != 'C')
        gt_array_add(cdsstruct.predstarts, i);

      if(i == locus_length - 1 || pred_vector[i+1] != 'C')
        gt_array_add(cdsstruct.predends, i);
    }

    // Exon structure counts
    if(char_is_exonic(refr_vector[i]))
    {
      if(i == 0 || !char_is_exonic(refr_vector[i-1]))
        gt_array_add(exonstruct.refrstarts, i);

      if(i == locus_length - 1 || !char_is_exonic(refr_vector[i+1]))
        gt_array_add(exonstruct.refrends, i);
    }
    if(char_is_exonic(pred_vector[i]))
    {
      if(i == 0 || !char_is_exonic(pred_vector[i-1]))
        gt_array_add(exonstruct.predstarts, i);

      if(i == locus_length - 1 || !char_is_exonic(pred_vector[i+1]))
        gt_array_add(exonstruct.predends, i);
    }

    // UTR structure counts
    if(char_is_utric(refr_vector[i]))
    {
      if(i == 0 || !char_is_utric(refr_vector[i-1]))
        gt_array_add(utrstruct.refrstarts, i);

      if(i == locus_length - 1 || !char_is_utric(refr_vector[i+1]))
        gt_array_add(utrstruct.refrends, i);
    }
    if(char_is_utric(pred_vector[i]))
    {
      if(i == 0 || !char_is_utric(pred_vector[i-1]))
        gt_array_add(utrstruct.predstarts, i);

      if(i == locus_length - 1 || !char_is_utric(pred_vector[i+1]))
        gt_array_add(utrstruct.predends, i);
    }
  }

  // Calculate nucleotide-level statistics from counts
  agn_comp_stats_scaled_resolve(&pair->stats.cds_nuc_stats);
  agn_comp_stats_scaled_resolve(&pair->stats.utr_nuc_stats);

  // Calculate statistics for structure from counts
  clique_pair_calc_struct_stats(&cdsstruct);
  clique_pair_calc_struct_stats(&exonstruct);
  clique_pair_calc_struct_stats(&utrstruct);
}

static void clique_pair_init_struct_dat(StructuralData *dat,
                                        AgnCompStatsBinary *stats)
{
  dat->refrstarts = gt_array_new( sizeof(GtUword) );
  dat->refrends   = gt_array_new( sizeof(GtUword) );
  dat->predstarts = gt_array_new( sizeof(GtUword) );
  dat->predends   = gt_array_new( sizeof(GtUword) );
  dat->stats      = stats;
}

static void clique_pair_term_struct_dat(StructuralData *dat)
{
  gt_array_delete(dat->refrstarts);
  gt_array_delete(dat->refrends);
  gt_array_delete(dat->predstarts);
  gt_array_delete(dat->predends);
}

static void clique_pair_test_data(GtQueue *queue)
{
  
}