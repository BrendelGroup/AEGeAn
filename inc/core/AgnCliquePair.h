#ifndef AEGEAN_CLIQUE_PAIR
#define AEGEAN_CLIQUE_PAIR

#include "genometools.h"
#include "AgnComparEval.h"
#include "AgnTranscriptClique.h"

/**
 * The purpose of the AgnCliquePair class is to associate all of the data needed
 * to compare transcripts from two alternative sources of annotation for the
 * same sequence.
 */
typedef struct AgnCliquePair AgnCliquePair;

/**
 * This enumerated type refers to all the possible outcomes when transcript
 * cliques from two different sources are compared.
 */
typedef enum
{
  AGN_CLIQUE_PAIR_UNCLASSIFIED,
  AGN_CLIQUE_PAIR_PERFECT_MATCH,
  AGN_CLIQUE_PAIR_MISLABELED,
  AGN_CLIQUE_PAIR_CDS_MATCH,
  AGN_CLIQUE_PAIR_EXON_MATCH,
  AGN_CLIQUE_PAIR_UTR_MATCH,
  AGN_CLIQUE_PAIR_NON_MATCH
} AgnCliquePairClassification;

/**
 * Build a pair of model vectors to represent this pair of maximal transcript
 * cliques.
 *
 * @param[in]  pair     pointer to a pair object
 */
void agn_clique_pair_build_model_vectors(AgnCliquePair *pair);

/**
 * Based on the already-computed comparison statistics, classify this clique
 * pair as a perfect match, a CDS match, etc.
 *
 * @param[in] pair    a clique pair
 * @returns           a classification code (see codes above)
 */
AgnCliquePairClassification agn_clique_pair_classify(AgnCliquePair *pair);

/**
 * Compare the annotations for this pair, reference vs prediction.
 *
 * @param[in] pair    the clique pair to compare
 */
void agn_clique_pair_comparative_analysis(AgnCliquePair *pair);

/**
 * Determine which pair has higher comparison scores.
 *
 * @param[in] p1    a clique pair
 * @param[in] p2    another clique pair
 * @returns         1 if the first pair has better scores, -1 if the second pair
 *                  has better scores, 0 if they are equal
 */
int agn_clique_pair_compare(void *p1, void *p2);

/**
 * Same as agn_clique_pair_compare, but without pointer dereferencing.
 *
 * @param[in] p1    a clique pair
 * @param[in] p2    another clique pair
 * @returns         1 if the first pair has better scores, -1 if the second pair
 *                  has better scores, 0 if they are equal
 */
int agn_clique_pair_compare_direct(AgnCliquePair *p1, AgnCliquePair *p2);

/**
 * Same as agn_clique_pair_compare, but in reverse.
 *
 * @param[in] p1    a clique pair
 * @param[in] p2    another clique pair
 * @returns         1 if the first pair has better scores, -1 if the second pair
 *                  has better scores, 0 if they are equal
 */
int agn_clique_pair_compare_reverse(void *p1, void *p2);

/**
 * Free the memory previously occupied by this pair object
 *
 * @param[out] pair    object to be deleted
 */
void agn_clique_pair_delete(AgnCliquePair *pair);

/**
 * Return the calculated annotation edit distance between this pair of
 * transcript annotation cliques.
 *
 * @param[in] pair    a clique pair
 * @returns           the associated annotation edit distance
 */
double agn_clique_pair_get_edit_distance(AgnCliquePair *pair);

/**
 * Return a pointer to this pair's prediction transcript clique.
 *
 * @param[in] pair    a clique pair
 * @returns           a pointer to the pred clique
 */
AgnTranscriptClique *agn_clique_pair_get_pred_clique(AgnCliquePair *pair);

/**
 * Get the model vector associated with this pair's prediction transcript clique.
 *
 * @param[in] pair    a clique pair
 * @returns           the prediction model vector
 */
const char *agn_clique_pair_get_pred_vector(AgnCliquePair *pair);

/**
 * Return a pointer to this pair's reference transcript clique.
 *
 * @param[in] pair    a clique pair
 * @returns           a pointer to the refr clique
 */
AgnTranscriptClique *agn_clique_pair_get_refr_clique(AgnCliquePair *pair);

/**
 * Get the model vector associated with this pair's reference transcript clique.
 *
 * @param[in] pair    a clique pair
 * @returns           the reference model vector
 */
const char *agn_clique_pair_get_refr_vector(AgnCliquePair *pair);

/**
 * Return a pointer to this clique pair's comparison statistics.
 *
 * @param[in] pair    a clique pair
 * @returns           a pointer to the pair's stats
 */
AgnComparison *agn_clique_pair_get_stats(AgnCliquePair *pair);

/**
 * Determine whether there are UTRs in this clique pair
 *
 * @param[in] pair    the clique pair in question
 * @returns           boolean indicating whether the clique pair contains any
 *                    UTRs
 */
bool agn_clique_pair_has_utrs(AgnCliquePair *pair);

/**
 * Does this clique pair contain a single reference transcript and a single
 * prediction transcript?
 *
 * @param[in] pair    the clique pair
 * @returns           boolean indicating whether the clique pair is simple
 */
bool agn_clique_pair_is_simple(AgnCliquePair *pair);

/**
 * Get the length of the locus to which this clique pair belongs.
 *
 * @param[in] pair    the clique pair
 * @returns           the length (in bp) of the corresponding locus
 */
unsigned long agn_clique_pair_length(AgnCliquePair *pair);

/**
 * Determine whether this clique pair needs comparison (i.e., whether there are
 * both reference and prediction transcripts).
 *
 * @param[in] pair    the clique pair
 * @returns           1 indicating comparison needed, 0 indicating no comparison
 *                    needed
 */
bool agn_clique_pair_needs_comparison(AgnCliquePair *pair);

/**
 * Allocate some memory for a pair object.
 *
 * @param[in] refr_clique    reference transcript(s)
 * @param[in] pred_clique    prediction transcript(s)
 * @param[in] locus_range    the coordinates of the locus to which this pair
 *                           belongs
 * @returns                  pointer to a new pair object
 */
AgnCliquePair* agn_clique_pair_new(const char *seqid,
                                   AgnTranscriptClique *refr_clique,
                                   AgnTranscriptClique *pred_clique,
                                   GtRange *locus_range);

/**
 * Add information about a clique pair to a set of aggregate characteristics.
 *
 * @param[in]  pair               the clique pair
 * @param[out] characteristics    a set of aggregate characteristics
 */
void agn_clique_pair_record_characteristics(AgnCliquePair *pair,
                                            AgnCompResultDesc *desc);

/**
 * Run unit tests for this class.
 *
 * @param[out] test    object for storing test results
 * @returns            true for success, false for failure
 */
bool agn_clique_pair_unit_test(AgnUnitTest *test);

#endif
