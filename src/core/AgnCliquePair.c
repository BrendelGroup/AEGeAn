#include <string.h>
#include <math.h>
#include "AgnCanonGeneStream.h"
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

typedef struct
{
  GtArray *refrstarts;
  GtArray *refrends;
  GtArray *predstarts;
  GtArray *predends;
  AgnCompStatsBinary *stats;
} StructuralData;


//----------------------------------------------------------------------------//
// Prototypes for private method(s)
//----------------------------------------------------------------------------//

/**
 * Add a transcript and all its features to a model vector.
 *
 * @param[in]  transcript     a transcript in the clique
 * @param[out] modelvector    ModelVectorData for the vector being built
 */
static void clique_pair_add_transcript_to_vector(GtFeatureNode *transcript,
                                                 void *data);

/**
 * Given a set of of start and end coordinates for reference and prediction
 * structures (exons, CDS segments, or UTR segments), determine the number of
 * congruent and incongruent structures.
 *
 * @param[in] dat    the data
 */
static void clique_pair_calc_struct_stats(StructuralData *dat);

/**
 * Initialize the data structure used to store start and end coordinates for
 * reference and prediction structures (exons, CDS segments, or UTR segments)
 * and associated statistics.
 *
 * @param[out] dat      the data structure
 * @param[in]  stats    stats to be stored in the structure
 */
static void clique_pair_init_struct_dat(StructuralData *dat,
                                        AgnCompStatsBinary *stats);

/**
 * Free the memory previously occupied by the data structure.
 *
 * @param[out] dat      the data structure
 */
static void clique_pair_term_struct_dat(StructuralData *dat);

// Macros for convenience when checking how a given nucleotide is annotated
#define char_is_exonic(C) (C == 'F' || C == 'T' || C == 'C')
#define char_is_utric(C)  (C == 'F' || C == 'T')


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

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
  GtUword locus_length = agn_clique_pair_length(pair);

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
    if(pair->refr_vector[i] == 'C' && pair->pred_vector[i] == 'C')
      pair->stats.cds_nuc_stats.tp++;
    else if(pair->refr_vector[i] == 'C' && pair->pred_vector[i] != 'C')
      pair->stats.cds_nuc_stats.fn++;
    else if(pair->refr_vector[i] != 'C' && pair->pred_vector[i] == 'C')
      pair->stats.cds_nuc_stats.fp++;
    else if(pair->refr_vector[i] != 'C' && pair->pred_vector[i] != 'C')
      pair->stats.cds_nuc_stats.tn++;

    // UTR nucleotide counts
    bool refr_utr = char_is_utric(pair->refr_vector[i]);
    bool pred_utr = char_is_utric(pair->pred_vector[i]);
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
        gt_array_add(cdsstruct.refrstarts, i);

      if(i == locus_length - 1 || pair->refr_vector[i+1] != 'C')
        gt_array_add(cdsstruct.refrends, i);
    }
    if(pair->pred_vector[i] == 'C')
    {
      if(i == 0 || pair->pred_vector[i-1] != 'C')
        gt_array_add(cdsstruct.predstarts, i);

      if(i == locus_length - 1 || pair->pred_vector[i+1] != 'C')
        gt_array_add(cdsstruct.predends, i);
    }

    // Exon structure counts
    if(char_is_exonic(pair->refr_vector[i]))
    {
      if(i == 0 || !char_is_exonic(pair->refr_vector[i-1]))
        gt_array_add(exonstruct.refrstarts, i);

      if(i == locus_length - 1 || !char_is_exonic(pair->refr_vector[i+1]))
        gt_array_add(exonstruct.refrends, i);
    }
    if(char_is_exonic(pair->pred_vector[i]))
    {
      if(i == 0 || !char_is_exonic(pair->pred_vector[i-1]))
        gt_array_add(exonstruct.predstarts, i);

      if(i == locus_length - 1 || !char_is_exonic(pair->pred_vector[i+1]))
        gt_array_add(exonstruct.predends, i);
    }

    // UTR structure counts
    if(char_is_utric(pair->refr_vector[i]))
    {
      if(i == 0 || !char_is_utric(pair->refr_vector[i-1]))
        gt_array_add(utrstruct.refrstarts, i);

      if(i == locus_length - 1 || !char_is_utric(pair->refr_vector[i+1]))
        gt_array_add(utrstruct.refrends, i);
    }
    if(char_is_utric(pair->pred_vector[i]))
    {
      if(i == 0 || !char_is_utric(pair->pred_vector[i-1]))
        gt_array_add(utrstruct.predstarts, i);

      if(i == locus_length - 1 || !char_is_utric(pair->pred_vector[i+1]))
        gt_array_add(utrstruct.predends, i);
    }
  }

  // Calculate nucleotide-level statistics from counts
  agn_comp_stats_scaled_resolve(&pair->stats.cds_nuc_stats);
  agn_comp_stats_scaled_resolve(&pair->stats.utr_nuc_stats);
  pair->stats.overall_identity = pair->stats.overall_matches /
                                 (double)locus_length;

  // Calculate statistics for structure from counts
  clique_pair_calc_struct_stats(&cdsstruct);
  clique_pair_calc_struct_stats(&exonstruct);
  clique_pair_calc_struct_stats(&utrstruct);
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
  GtUword refr_utrs = agn_transcript_clique_num_utrs(pair->refr_clique);
  GtUword pred_utrs = agn_transcript_clique_num_utrs(pair->pred_clique);
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

GtUword agn_clique_pair_length(AgnCliquePair *pair)
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

bool agn_clique_pair_unit_test(AgnUnitTest *test)
{
  GtError *error = gt_error_new();
  AgnLogger *logger = agn_logger_new();
  int pullresult;

  const char *refrfile = "data/gff3/grape-refr.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &refrfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtNodeStream *cgstream = agn_canon_gene_stream_new(gff3in, logger);
  GtArray *refrfeats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(cgstream,refrfeats,error);
  pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    agn_logger_log_error(logger, "error processing refr node stream: %s",
                         gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(cgstream);
  gt_node_stream_delete(arraystream);
  gt_array_sort(refrfeats, (GtCompare)agn_gt_genome_node_compare);

  const char *predfile = "data/gff3/grape-pred.gff3";
  gff3in = gt_gff3_in_stream_new_unsorted(1, &predfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  cgstream = agn_canon_gene_stream_new(gff3in, logger);
  GtArray *predfeats = gt_array_new( sizeof(GtFeatureNode *) );
  arraystream = gt_array_out_stream_new(cgstream, predfeats, error);
  pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    agn_logger_log_error(logger, "error processing refr node stream: %s",
                         gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(cgstream);
  gt_node_stream_delete(arraystream);
  gt_array_sort(predfeats, (GtCompare)agn_gt_genome_node_compare);

  bool genenumpass = (gt_array_size(refrfeats) == 12 &&
                      gt_array_size(predfeats) == 13);
  agn_unit_test_result(test, "parse grape example", genenumpass);

  GtFeatureNode **r2 = gt_array_get(refrfeats, 0);
  GtFeatureNode **p1 = gt_array_get(predfeats, 0);
  GtFeatureNode **p2 = gt_array_get(predfeats, 1);
  AgnTranscriptClique *tcr1 = agn_transcript_clique_new();
  AgnTranscriptClique *tcr2 = agn_transcript_clique_new();
  AgnTranscriptClique *tcp1 = agn_transcript_clique_new();
  AgnTranscriptClique *tcp2 = agn_transcript_clique_new();
  agn_transcript_clique_add(tcr2, *r2);
  agn_transcript_clique_add(tcp1, *p1);
  agn_transcript_clique_add(tcp2, *p2);
  GtRange lr1 = {72, 5081};
  GtRange lr2 = {10503, 11678};
  AgnCliquePair *pair1 = agn_clique_pair_new("chr8", tcr1, tcp1, &lr1);
  AgnCliquePair *pair2 = agn_clique_pair_new("chr8", tcr2, tcp2, &lr2);

  bool issimplepass = (!agn_clique_pair_is_simple(pair1) &&
                       agn_clique_pair_is_simple(pair2));
  agn_unit_test_result(test, "simple clique pair", issimplepass);
  bool needscomppass = (!agn_clique_pair_needs_comparison(pair1) &&
                        agn_clique_pair_needs_comparison(pair2));
  agn_unit_test_result(test, "comparison test", needscomppass);
  bool hasutrspass = (!agn_clique_pair_has_utrs(pair1) &&
                      agn_clique_pair_has_utrs(pair2));
  agn_unit_test_result(test, "UTR test", hasutrspass);

  agn_clique_pair_build_model_vectors(pair1);
  agn_clique_pair_build_model_vectors(pair2);
  const char *v1 = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTT"
"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
"FFFFFFFFFFFFFFFFFFFFFFFFF";
  const char *v2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
"IIIIIIIIIIIIIIIIIIIIIIIIIIICCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCIIIIIIIIIIIIIIIIIIIIIII"
"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIICCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
"CCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
"GGGGGGGGGGGGGGGGGGGGGGGGG";
  const char *mv1test = agn_clique_pair_get_pred_vector(pair2);
  const char *mv2test = agn_clique_pair_get_refr_vector(pair2);
  bool modelvectorpass = strcmp(v1, mv1test) == 0 && strcmp(v2, mv2test) == 0;
  agn_unit_test_result(test, "model vectors", modelvectorpass);

  agn_clique_pair_delete(pair1);
  agn_clique_pair_delete(pair2);
  agn_transcript_clique_delete(tcr1);
  agn_transcript_clique_delete(tcr2);
  agn_transcript_clique_delete(tcp1);
  agn_transcript_clique_delete(tcp2);

  GtFeatureNode **r3 = gt_array_get(refrfeats, 2);
  GtFeatureNode **p3 = gt_array_get(predfeats, 3);
  AgnTranscriptClique *tcr3 = agn_transcript_clique_new();
  AgnTranscriptClique *tcp3 = agn_transcript_clique_new();
  agn_transcript_clique_add(tcr3, *r3);
  agn_transcript_clique_add(tcp3, *p3);
  GtRange lr3 = {26493, 29602};
  AgnCliquePair *pair3 = agn_clique_pair_new("chr8", tcr3, tcp3, &lr3);
  agn_clique_pair_build_model_vectors(pair3);
  agn_clique_pair_comparative_analysis(pair3);
  bool companalypass = agn_clique_pair_classify(pair3) ==
                       AGN_CLIQUE_PAIR_CDS_MATCH;
  agn_unit_test_result(test, "analysis and classification", companalypass);
  agn_clique_pair_delete(pair3);
  agn_transcript_clique_delete(tcr3);
  agn_transcript_clique_delete(tcp3);

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
  agn_logger_delete(logger);
  return true;
}

static void clique_pair_add_transcript_to_vector(GtFeatureNode *transcript,
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
      GtUword fn_start = gt_genome_node_get_start((GtGenomeNode *)fn);
      GtUword fn_end = gt_genome_node_get_end((GtGenomeNode *)fn);
      GtUword i;

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
