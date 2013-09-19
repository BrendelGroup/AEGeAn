#ifndef AEGEAN_CLIQUE_PAIR
#define AEGEAN_CLIQUE_PAIR

#include "genometools.h"
#include "AgnComparEval.h"
#include "AgnTranscriptClique.h"

/**
 * @class AgnCliquePair
 *
 * The purpose of the AgnCliquePair class is to associate all of the data needed
 * to compare transcripts from two alternative sources of annotation for the
 * same sequence.
 */
typedef struct AgnCliquePair AgnCliquePair;

/**
 * @type This enumerated type refers to all the possible outcomes when
 * transcript cliques from two different sources are compared:
 * ``AGN_CLIQUE_PAIR_UNCLASSIFIED``, ``AGN_CLIQUE_PAIR_PERFECT_MATCH``,
 * ``AGN_CLIQUE_PAIR_MISLABELED``,   ``AGN_CLIQUE_PAIR_CDS_MATCH``,
 * ``AGN_CLIQUE_PAIR_EXON_MATCH``,   ``AGN_CLIQUE_PAIR_UTR_MATCH``, and
 * ``AGN_CLIQUE_PAIR_NON_MATCH``.
 */
enum AgnCliquePairClassification
{
  AGN_CLIQUE_PAIR_UNCLASSIFIED,
  AGN_CLIQUE_PAIR_PERFECT_MATCH,
  AGN_CLIQUE_PAIR_MISLABELED,
  AGN_CLIQUE_PAIR_CDS_MATCH,
  AGN_CLIQUE_PAIR_EXON_MATCH,
  AGN_CLIQUE_PAIR_UTR_MATCH,
  AGN_CLIQUE_PAIR_NON_MATCH
};
typedef enum AgnCliquePairClassification AgnCliquePairClassification;

/**
 * @function Build a pair of model vectors to represent this pair of maximal
 * transcripts or transcript cliques.
 */
void agn_clique_pair_build_model_vectors(AgnCliquePair *pair);

/**
 * @function Based on the already-computed comparison statistics, classify this
 * clique pair as a perfect match, a CDS match, etc. See
 * :c:type:`AgnCliquePairClassification`.
 */
AgnCliquePairClassification agn_clique_pair_classify(AgnCliquePair *pair);

/**
 * @function Compare the annotations for this pair, reference vs prediction.
 */
void agn_clique_pair_comparative_analysis(AgnCliquePair *pair);

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
 * @function Return the calculated annotation edit distance between this pair of
 * transcripts or transcript cliques.
 */
double agn_clique_pair_get_edit_distance(AgnCliquePair *pair);

/**
 * @function Return a pointer to this pair's prediction transcript or transcript
 * clique.
 */
AgnTranscriptClique *agn_clique_pair_get_pred_clique(AgnCliquePair *pair);

/**
 * @function Get the model vector associated with this pair's prediction
 * transcript clique.
 */
const char *agn_clique_pair_get_pred_vector(AgnCliquePair *pair);

/**
 * @function Return a pointer to this pair's reference transcript or transcript
 * clique.
 */
AgnTranscriptClique *agn_clique_pair_get_refr_clique(AgnCliquePair *pair);

/**
 * @function Get the model vector associated with this pair's reference
 * transcript clique.
 */
const char *agn_clique_pair_get_refr_vector(AgnCliquePair *pair);

/**
 * @function Return a pointer to this clique pair's comparison statistics.
 */
AgnComparison *agn_clique_pair_get_stats(AgnCliquePair *pair);

/**
 * @function Determine whether there are UTRs in this clique pair
 */
bool agn_clique_pair_has_utrs(AgnCliquePair *pair);

/**
 * @function Does this clique pair contain a single reference transcript and a
 * single prediction transcript?
 */
bool agn_clique_pair_is_simple(AgnCliquePair *pair);

/**
 * @function Get the length of the locus to which this clique pair belongs.
 */
GtUword agn_clique_pair_length(AgnCliquePair *pair);

/**
 * @function Determine whether this clique pair needs comparison (i.e., whether
 * there are both reference and prediction transcripts).
 */
bool agn_clique_pair_needs_comparison(AgnCliquePair *pair);

/**
 * @function Class constructor.
 */
AgnCliquePair* agn_clique_pair_new(const char *seqid,
                                   AgnTranscriptClique *refr_clique,
                                   AgnTranscriptClique *pred_clique,
                                   GtRange *locus_range);

/**
 * @function Add information about a clique pair to a set of aggregate
 * characteristics.
 */
void agn_clique_pair_record_characteristics(AgnCliquePair *pair,
                                            AgnCompResultDesc *desc);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_clique_pair_unit_test(AgnUnitTest *test);

#endif
