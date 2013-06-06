#ifndef AEGEAN_PAIRWISE_COMPARE_LOCUS
#define AEGEAN_PAIRWISE_COMPARE_LOCUS

#include "genometools.h"
#include "AgnCliquePair.h"

/**
 * The purpose of the AgnPairwiseCompareLocus class is to store all of the data
 * associated with a distinct locus to facilitate the comparison of two sets of
 * gene structure annotations for that locus.
 */
typedef struct AgnPairwiseCompareLocus AgnPairwiseCompareLocus;

/**
 * Associate the given gene annotation with this gene locus.
 *
 * @param[out] locus           locus to which the gene annotation will be added
 * @param[in]  gene_feature    annotation to associate with this locus
 */
void agn_pairwise_compare_locus_add_pred_gene( AgnPairwiseCompareLocus *locus,
                                               GtFeatureNode *gene_feature );

/**
 * Associate the given gene annotation with this gene locus.
 *
 * @param[out] locus           locus to which the gene annotation will be added
 * @param[in]  gene_feature    annotation to associate with this locus
 */
void agn_pairwise_compare_locus_add_refr_gene( AgnPairwiseCompareLocus *locus,
                                               GtFeatureNode *gene_feature );

/**
 * Calculate and store the splice complexity of this locus' reference and
 * prediction annotations, respectively.
 *
 * @param[in] locus    this locus
 */
void agn_pairwise_compare_locus_calc_splice_complexity_pred
(
  AgnPairwiseCompareLocus *locus
);
void agn_pairwise_compare_locus_calc_splice_complexity_refr
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Array comparison function for locus objects
 *
 * @param[in] p1    pointer to one locus object
 * @param[in] p2    pointer to another locus object
 * @returns         -1, 0, or 1 depending on their relative position
 */
int agn_pairwise_compare_locus_array_compare(const void *p1, const void *p2);

/**
 * Free the memory previously occupied by this locus object
 *
 * @param[out] locus    object to be deleted
 */
void agn_pairwise_compare_locus_delete(AgnPairwiseCompareLocus *locus);

/**
 * Determine whether a locus should be filtered (excluded) from reports based on
 * user-provided filter settings.
 *
 * @param[in] locus    the locus in question
 * @returns            true if the locus should be excluded, false otherwise
 */
bool agn_pairwise_compare_locus_filter( AgnPairwiseCompareLocus *locus,
                                        AgnCompareFilters *filters );

/**
 * For gene loci with multiple transcripts or transcript cliques, we do not want
 * to report every pairwise comparison of every reference clique with every
 * prediction clique. Instead, we report each reference clique along with its
 * best matching prediction clique. If there are any prediction cliques that are
 * not included in these matches, they are reported separately.
 *
 * @param[in]  locus    the locus
 * @returns             the pairs to be reported
 */
GtArray *agn_pairwise_compare_locus_find_best_pairs
(
  AgnPairwiseCompareLocus *locus
);

/**
 * We use the Bron-Kerbosch algorithm to separate reference transcripts and
 * prediction transcripts into a set of maximal cliques, or subsets of
 * transcripts that do not overlap within the group. We then compare each
 * reference clique with each prediction clique.
 *
 * @param[in] locus              the locus
 * @param[in] trans_per_locus    the maximum number of reference or prediction
 *                               mRNAs allowed for a given locus; since the
 *                               maximal clique enumeration algorithm used is of
 *                               exponential complexity, there needs to be a
 *                               reasonable limit to the number of mRNAs
 *                               considered when identifying transcript cliques
 * @returns                      an array of locus objects representing every
 *                               pair of a reference transcript with a
 *                               prediction transcript
 */
GtArray* agn_pairwise_compare_locus_get_clique_pairs
(
  AgnPairwiseCompareLocus *locus,
  unsigned int trans_per_locus
);

/**
 * Get this locus' end coordinate.
 *
 * @param[in] locus    the locus
 * @returns            the end coordinate of the locus
 */
unsigned long agn_pairwise_compare_locus_get_end(AgnPairwiseCompareLocus *locus);

/**
 * Get this locus' length.
 *
 * @param[in] locus    the locus
 * @returns            the length of the locus
 */
unsigned long agn_pairwise_compare_locus_get_length
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Given a reference transcript, find the clique pair with the best scores of
 * which the transcript is a part.
 *
 * @param[in] locus          the locus
 * @param[in] refr_clique    a reference transcript clique belonging to the
 * @returns                  locus the pair containing the transcript that has
 *                           the highest comparison scores
 */
AgnCliquePair* agn_pairwise_compare_locus_get_optimal_clique_pair
(
  AgnPairwiseCompareLocus *locus,
  AgnTranscriptClique *refr_clique
);

/**
 * Get the prediction genes for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the prediction genes
 */
GtArray *agn_pairwise_compare_locus_get_pred_genes
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the prediction gene IDs for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the prediction gene IDs
 */
GtArray *agn_pairwise_compare_locus_get_pred_gene_ids
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the splice complexity of prediction annotations associated with this
 * locus.
 *
 * @param[in] locus    the locus
 * @returns            the prediction splice complexity
 */
double agn_pairwise_compare_locus_get_pred_splice_complexity
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the prediction transcripts for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the prediction transcripts
 */
GtArray *agn_pairwise_compare_locus_get_pred_transcripts
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the prediction transcript IDs for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the prediction transcript IDs
 */
GtArray *agn_pairwise_compare_locus_get_pred_transcript_ids
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the reference genes for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the reference genes
 */
GtArray *agn_pairwise_compare_locus_get_refr_genes
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the reference gene IDs for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the reference gene IDs
 */
GtArray *agn_pairwise_compare_locus_get_refr_gene_ids
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the splice complexity of reference annotations associated with this
 * locus.
 *
 * @param[in] locus    the locus
 * @returns            the reference splice complexity
 */
double agn_pairwise_compare_locus_get_refr_splice_complexity
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the reference transcripts for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the reference transcripts
 */
GtArray *agn_pairwise_compare_locus_get_refr_transcripts
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the reference transcript IDs for this locus.
 *
 * @param[in] locus    the locus
 * @returns            the reference transcript IDs
 */
GtArray *agn_pairwise_compare_locus_get_refr_transcript_ids
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get this locus' sequence ID.
 *
 * @param[in] locus    the locus
 * @returns            the sequence ID of the locus
 */
const char* agn_pairwise_compare_locus_get_seqid(AgnPairwiseCompareLocus *locus);

/**
 * Get this locus' start coordinate.
 *
 * @param[in] locus    the locus
 * @returns            the start coordinate of the locus
 */
unsigned long agn_pairwise_compare_locus_get_start
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get a list of all the prediction transcript cliques that have no
 * corresponding reference transcript clique.
 *
 * @param[in] locus    the locus
 * @returns            a list of all unique prediction cliques
 */
GtArray *agn_pairwise_compare_locus_get_unique_pred_cliques
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get a list of all the reference transcript cliques that have no
 * corresponding prediction transcript clique.
 *
 * @param[in] locus    the locus
 * @returns            a list of all unique reference cliques
 */
GtArray *agn_pairwise_compare_locus_get_unique_refr_cliques
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Loci with more than one reference transcript clique and/or more than one
 * prediction transcript clique are complex and require special handling.
 *
 * @param[in] locus    a gene locus
 * @returns            boolean indicating whether the locus is complex
 */
bool agn_pairwise_compare_locus_is_complex(AgnPairwiseCompareLocus *locus);

/**
 * Allocate some memory for a locus object.
 *
 * @param[in] seqid       ID of the locus' sequence
 * @returns               pointer to a new locus object
 */
AgnPairwiseCompareLocus* agn_pairwise_compare_locus_new(const char *seqid);

/**
 * Get the total number of prediction exons for the locus.
 *
 * @param[in]    the locus
 * @returns      the number of prediction exons
 */
unsigned long agn_pairwise_compare_locus_num_pred_exons
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the number of prediction gene annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of prediction gene annotations for that locus
 */
unsigned long agn_pairwise_compare_locus_num_pred_genes
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the number of prediction mRNA annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of prediction mRNA annotations for that locus
 */
unsigned long agn_pairwise_compare_locus_num_pred_transcripts
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the total number of reference exons for the locus.
 *
 * @param[in]    the locus
 * @returns      the number of reference exons
 */
unsigned long agn_pairwise_compare_locus_num_refr_exons
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the number of reference gene annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of reference gene annotations for that locus
 */
unsigned long agn_pairwise_compare_locus_num_refr_genes
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Get the number of reference mRNA annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of reference mRNA annotations for that locus
 */
unsigned long agn_pairwise_compare_locus_num_refr_transcripts
(
  AgnPairwiseCompareLocus *locus
);

/**
 * The combined length of all prediction coding sequences associated with this
 * locus.
 *
 * @param[in] locus    the locus
 * @returns            the combined CDS length
 */
unsigned long agn_pairwise_compare_locus_pred_cds_length
(
  AgnPairwiseCompareLocus *locus
);

/**
 * Return the coordinates of this locus.
 *
 * @param[in] locus    the locus
 * @returns            its coordinates
 */
GtRange agn_pairwise_compare_locus_range(AgnPairwiseCompareLocus *locus);

/**
 * The combined length of all reference coding sequences associated with this
 * locus.
 *
 * @param[in] locus    the locus
 * @returns            the combined CDS length
 */
unsigned long agn_pairwise_compare_locus_refr_cds_length
(
  AgnPairwiseCompareLocus *locus
);

#endif
