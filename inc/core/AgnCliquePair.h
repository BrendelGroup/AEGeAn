#ifndef AEGEAN_CLIQUE_PAIR
#define AEGEAN_CLIQUE_PAIR

#include "AgnComparison.h"
#include "AgnTranscriptClique.h"

/**
 * @class AgnCliquePair
 *
 * The AgnCliquePair class facilitates comparison of two alternative sources of
 * annotation for the same sequence.
 */
typedef struct AgnCliquePair AgnCliquePair;

/**
 * @function Based on the already-computed comparison statistics, classify this
 * clique pair as a perfect match, a CDS match, etc. See
 * :c:type:`AgnCompClassification`.
 */
AgnCompClassification agn_clique_pair_classify(AgnCliquePair *pair);

/**
 * @function Add this clique pair's internal comparison stats to a larger set of
 * aggregate stats.
 */
void agn_clique_pair_comparison_aggregate(AgnCliquePair *pair,
                                          AgnComparison *comp);

/**
 * @function Same as c:func:`agn_clique_pair_compare_direct`, but with pointer
 * dereferencing.
 */
int agn_clique_pair_compare(void *p1, void *p2);

/**
 * @function Determine which pair has higher comparison scores. Returns 1 if the
 * first pair has better scores, -1 if the second pair has better scores, 0 if
 * they are equal.
 */
int agn_clique_pair_compare_direct(AgnCliquePair *p1, AgnCliquePair *p2);

/**
 * @function Negation of c:func:`agn_clique_pair_compare`.
 */
int agn_clique_pair_compare_reverse(void *p1, void *p2);

/**
 * @function Class destructor.
 */
void agn_clique_pair_delete(AgnCliquePair *pair);

/**
 * @function Return a pointer to the prediction annotation from this pair.
 */
AgnTranscriptClique *agn_clique_pair_get_pred_clique(AgnCliquePair *pair);

/**
 * @function Return a pointer to the reference annotation from this pair.
 */
AgnTranscriptClique *agn_clique_pair_get_refr_clique(AgnCliquePair *pair);

/**
 * @function Class constructor.
 */
AgnCliquePair* agn_clique_pair_new(AgnTranscriptClique *refr,
                                   AgnTranscriptClique *pred);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_clique_pair_unit_test(AgnUnitTest *test);

#endif
