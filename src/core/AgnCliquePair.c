#include <math.h>
#include <string.h>
#include "core/queue_api.h"
#include "AgnCliquePair.h"
#include "AgnUtils.h"

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
  GtQueue *pairs = gt_queue_new();
  clique_pair_test_data(pairs);
  gt_assert(gt_queue_size(pairs) == 2);

  AgnCliquePair *pair = gt_queue_get(pairs);
  AgnCompClassification result = agn_clique_pair_classify(pair);
  bool simplecheck = (result == AGN_COMP_CLASS_PERFECT_MATCH);
  agn_unit_test_result(test, "perfect match vs. self", simplecheck);
  agn_transcript_clique_delete(pair->refr_clique);
  agn_transcript_clique_delete(pair->pred_clique);
  agn_clique_pair_delete(pair);

  pair = gt_queue_get(pairs);
  result = agn_clique_pair_classify(pair);
  bool cdscheck = result == (AGN_COMP_CLASS_CDS_MATCH);
  agn_unit_test_result(test, "CDS match", cdscheck);
  agn_transcript_clique_delete(pair->refr_clique);
  agn_transcript_clique_delete(pair->pred_clique);
  agn_clique_pair_delete(pair);

  gt_queue_delete(pairs);
  return simplecheck && cdscheck;
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
  gt_assert(
      strlen(refr_vector) == gt_genome_node_get_length(pair->refr_clique) &&
      strlen(pred_vector) == gt_genome_node_get_length(pair->refr_clique)
  );
  pair->stats.overall_length = locus_length;

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
  gt_assert(queue != NULL);
  GtArray *refrfeats, *predfeats;

  GtError *error = gt_error_new();
  const char *refrfile = "data/gff3/grape-refr-mrnas.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &refrfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  refrfeats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(gff3in, refrfeats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnCliquePair::clique_pair_test_data] error processing "
            "reference features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(arraystream);
  gt_array_sort(refrfeats, (GtCompare)agn_genome_node_compare);

  const char *predfile = "data/gff3/grape-pred-mrnas.gff3";
  gff3in = gt_gff3_in_stream_new_unsorted(1, &predfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  predfeats = gt_array_new( sizeof(GtFeatureNode *) );
  arraystream = gt_array_out_stream_new(gff3in, predfeats, error);
  pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnCliquePair::clique_pair_test_data] error processing "
            "prediction features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(arraystream);
  gt_array_sort(predfeats, (GtCompare)agn_genome_node_compare);

  gt_assert(gt_array_size(refrfeats) == 12 && gt_array_size(predfeats) == 13);

  GtRange range = { 26493, 29591 };
  GtStr *seqid = gt_str_new_cstr("chr8");
  AgnSequenceRegion region = { gt_str_get(seqid), range };
  GtFeatureNode *pred = *(GtFeatureNode **)gt_array_get(predfeats, 3);
  AgnTranscriptClique *refrclique = agn_transcript_clique_new(&region);
  agn_transcript_clique_add(refrclique, pred);
  AgnTranscriptClique *predclique = agn_transcript_clique_new(&region);
  agn_transcript_clique_add(predclique, pred);
  AgnCliquePair *pair = agn_clique_pair_new(refrclique, predclique);
  gt_queue_add(queue, pair);

  region.range.end = 29602;
  GtFeatureNode *refr = *(GtFeatureNode **)gt_array_get(refrfeats, 2);
  refrclique = agn_transcript_clique_new(&region);
  agn_transcript_clique_add(refrclique, refr);
  predclique = agn_transcript_clique_new(&region);
  agn_transcript_clique_add(predclique, pred);
  pair = agn_clique_pair_new(refrclique, predclique);
  gt_queue_add(queue, pair);

  gt_str_delete(seqid);
  while(gt_array_size(refrfeats) > 0)
  {
    GtGenomeNode **gn = gt_array_pop(refrfeats);
    gt_genome_node_delete(*gn);
  }
  gt_array_delete(refrfeats);
  while(gt_array_size(predfeats) > 0)
  {
    GtGenomeNode **gn = gt_array_pop(predfeats);
    gt_genome_node_delete(*gn);
  }
  gt_array_delete(predfeats);
  gt_error_delete(error);
}