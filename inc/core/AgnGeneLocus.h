#ifndef AEGEAN_GENE_LOCUS
#define AEGEAN_GENE_LOCUS

#include "genometools.h"
#include "AgnCliquePair.h"

/**
 * The purpose of the AgnGeneLocus class is to store all of the data
 * associated with a distinct locus to facilitate the comparison of two sets of
 * gene structure annotations for that locus.
 */
typedef struct AgnGeneLocus AgnGeneLocus;
enum AgnComparisonSource { REFERENCESOURCE, PREDICTIONSOURCE, DEFAULTSOURCE };
typedef enum AgnComparisonSource AgnComparisonSource;

/**
 * Associate the given gene annotation with this gene locus. Rather than calling
 * this function directly, users are recommended to use one of the following
 * macros: agn_gene_locus_add_pred_gene(locus, gene)' and
 * `agn_gene_locus_add_refr_gene(locus, gene)', to be used when keeping track of
 * an annotation's source is important (i.e. for pairwise comparison); and
 * `agn_gene_locus_add_gene(locus, gene)' otherwise.
 *
 * @param[out] locus     locus to which the gene annotation will be added
 * @param[in]  gene      annotation to associate with this locus
 * @param[in]  source    indication as to whether this gene is a reference gene
 *                       (REFERENCESOURCE) or prediction gene (PREDICTIONSOURCE)
 *                       if this locus is intended for pairwise comparison; use
 *                       DEFAULTSOURCE otherwise
 */
void agn_gene_locus_add(AgnGeneLocus *locus, GtFeatureNode *gene,
                        AgnComparisonSource source);
#define agn_gene_locus_add_pred_gene(LC, GN)\
        agn_gene_locus_add(LC, GN, PREDICTIONSOURCE)
#define agn_gene_locus_add_refr_gene(LC, GN)\
        agn_gene_locus_add(LC, GN, REFERENCESOURCE)
#define agn_gene_locus_add_gene(LC, GN)\
        agn_gene_locus_add(LC, GN, DEFAULTSOURCE)

/**
 * Analog of strcmp for comparing AgnGeneLocus objects, used for sorting GtArray
 * objects containing AgnGeneLocus objects.
 *
 * @param[in] p1    pointer to a pointer to one locus object
 * @param[in] p2    pointer to a pointer to another locus object
 * @returns         -1, 0, or 1 depending on their relative position
 */
int agn_gene_locus_array_compare(const void *p1, const void *p2);

/**
 * Free the memory previously occupied by this locus object.
 *
 * @param[out] locus    locus to be deleted
 */
void agn_gene_locus_delete(AgnGeneLocus *locus);

/**
 * Given a set of filtering criteria, determine whether a locus meets those
 * criteria.
 *
 * @param[in] locus    the locus in question
 * @returns            true if the locus should be filtered (if it does not meet
 *                     the criteria), false otherwise
 */
bool agn_gene_locus_filter(AgnGeneLocus *locus, AgnCompareFilters *filters);

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
GtArray *agn_gene_locus_find_best_pairs(AgnGeneLocus *locus);

/**
 * Get the genes associated with this locus. Rather than calling this function
 * directly, users are encouraged to use one of the following macros:
 * `agn_gene_locus_pred_genes(locus)' to retrieve prediction genes,
 * `agn_gene_locus_refr_genes(locus)' to retrieve reference genes, or
 * `agn_gene_locus_get_genes(locus)' if the source of annotation is undesignated
 * or irrelevant.
 *
 * @param[in] locus    a gene locus
 * @param[in] src      REFERENCESOURCE will return only reference genes,
 *                     PREDICTIONSOURCE will return only prediction genes,
 *                     DEFAULTSOURCE will return all genes
 * @returns            array containing the genes
 */
GtArray *agn_gene_locus_genes(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_pred_genes(LC)\
        agn_gene_locus_genes(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_genes(LC)\
        agn_gene_locus_genes(LC, REFERENCESOURCE)
#define agn_gene_locus_get_genes(LC)\
        agn_gene_locus_genes(LC, DEFAULTSOURCE)

/**
 * Get IDs of the genes associated with this locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * `agn_gene_locus_pred_gene_ids(locus)' to retrieve prediction genes IDs,
 * `agn_gene_locus_refr_gene_ids(locus)' to retrieve reference genes IDs, or
 * `agn_gene_locus_get_gene_ids(locus)' if the source of annotation is
 * undesignated or irrelevant.
 *
 * @param[in] locus    a gene locus
 * @param[in] src      REFERENCESOURCE will return only reference gene IDs,
 *                     PREDICTIONSOURCE will return only prediction gene IDs,
 *                     DEFAULTSOURCE will return all gene IDs
 * @returns            array containing the genes
 */
GtArray *agn_gene_locus_gene_ids(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_pred_gene_ids(LC)\
        agn_gene_locus_genes(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_gene_ids(LC)\
        agn_gene_locus_genes(LC, REFERENCESOURCE)
#define agn_gene_locus_get_gene_ids(LC)\
        agn_gene_locus_genes(LC, DEFAULTSOURCE)

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
GtArray* agn_gene_locus_get_clique_pairs(AgnGeneLocus *locus,
                                         unsigned int trans_per_locus);

/**
 * Get this locus' end coordinate.
 *
 * @param[in] locus    the locus
 * @returns            the end coordinate of the locus
 */
unsigned long agn_gene_locus_get_end(AgnGeneLocus *locus);

/**
 * Get this locus' length.
 *
 * @param[in] locus    the locus
 * @returns            the length of the locus
 */
unsigned long agn_gene_locus_get_length(AgnGeneLocus *locus);

/**
 * Given a reference transcript, find the clique pair with the best scores of
 * which the transcript is a part.
 *
 * @param[in] locus          the locus
 * @param[in] refr_clique    a reference transcript clique belonging to the
 * @returns                  locus the pair containing the transcript that has
 *                           the highest comparison scores
 */
AgnCliquePair* agn_gene_locus_get_optimal_clique_pair(AgnGeneLocus *locus,
                                              AgnTranscriptClique *refr_clique);

/**
 * Get this locus' sequence ID.
 *
 * @param[in] locus    the locus
 * @returns            the sequence ID of the locus
 */
const char* agn_gene_locus_get_seqid(AgnGeneLocus *locus);

/**
 * Get this locus' start coordinate.
 *
 * @param[in] locus    the locus
 * @returns            the start coordinate of the locus
 */
unsigned long agn_gene_locus_get_start(AgnGeneLocus *locus);

/**
 * Get a list of all the prediction transcript cliques that have no
 * corresponding reference transcript clique.
 *
 * @param[in] locus    the locus
 * @returns            a list of all unique prediction cliques
 */
GtArray *agn_gene_locus_get_unique_pred_cliques(AgnGeneLocus *locus);

/**
 * Get a list of all the reference transcript cliques that have no
 * corresponding prediction transcript clique.
 *
 * @param[in] locus    the locus
 * @returns            a list of all unique reference cliques
 */
GtArray *agn_gene_locus_get_unique_refr_cliques(AgnGeneLocus *locus);

/**
 * Loci with more than one reference transcript clique and/or more than one
 * prediction transcript clique are complex and require special handling.
 *
 * @param[in] locus    a gene locus
 * @returns            boolean indicating whether the locus is complex
 */
bool agn_gene_locus_is_complex(AgnGeneLocus *locus);

/**
 * Allocate some memory for a locus object.
 *
 * @param[in] seqid       ID of the locus' sequence
 * @returns               pointer to a new locus object
 */
AgnGeneLocus* agn_gene_locus_new(const char *seqid);

/**
 * Get the total number of prediction exons for the locus.
 *
 * @param[in]    the locus
 * @returns      the number of prediction exons
 */
unsigned long agn_gene_locus_num_pred_exons(AgnGeneLocus *locus);

/**
 * Get the number of prediction gene annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of prediction gene annotations for that locus
 */
unsigned long agn_gene_locus_num_pred_genes(AgnGeneLocus *locus);

/**
 * Get the number of prediction mRNA annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of prediction mRNA annotations for that locus
 */
unsigned long agn_gene_locus_num_pred_transcripts(AgnGeneLocus *locus);

/**
 * Get the total number of reference exons for the locus.
 *
 * @param[in]    the locus
 * @returns      the number of reference exons
 */
unsigned long agn_gene_locus_num_refr_exons(AgnGeneLocus *locus);

/**
 * Get the number of reference gene annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of reference gene annotations for that locus
 */
unsigned long agn_gene_locus_num_refr_genes(AgnGeneLocus *locus);

/**
 * Get the number of reference mRNA annotations associated with this locus.
 *
 * @param[in] locus    the locus
 * @returns            the number of reference mRNA annotations for that locus
 */
unsigned long agn_gene_locus_num_refr_transcripts(AgnGeneLocus *locus);

/**
 * The combined length of all prediction coding sequences associated with this
 * locus.
 *
 * @param[in] locus    the locus
 * @returns            the combined CDS length
 */
unsigned long agn_gene_locus_pred_cds_length(AgnGeneLocus *locus);

/**
 * Return the coordinates of this locus.
 *
 * @param[in] locus    the locus
 * @returns            its coordinates
 */
GtRange agn_gene_locus_range(AgnGeneLocus *locus);

/**
 * The combined length of all reference coding sequences associated with this
 * locus.
 *
 * @param[in] locus    the locus
 * @returns            the combined CDS length
 */
unsigned long agn_gene_locus_refr_cds_length(AgnGeneLocus *locus);

/**
 * Calculate the splice complexity of this gene locus. Rather than calling this
 * method directly, users are recommended to use one of the following macros:
 * `agn_gene_locus_prep_splice_complexity(locus)' to calculate the splice
 * complexity of just the prediction transcripts,
 * `agn_gene_locus_refr_splice_complexity(locus)' to calculate the splice
 * complexity of just the reference transcripts, and
 * `agn_gene_locus_calc_splice_complexity(locus)' to calculate the splice
 * complexity taking into account all transcripts.
 *
 * @param[in] locus    the locus
 * @param[in] src      indication as to whether to calculate splice complexity
 *                     for just reference transcripts (REFERENCESOURCE),
 *                     just prediction transcripts (PREDICTIONSOURCE),
 *                     or all transcripts (DEFAULTSOURCE)
 * @returns            the splice complexity
 */
double agn_gene_locus_splice_complexity(AgnGeneLocus *locus,
                                        AgnComparisonSource src);
#define agn_gene_locus_pred_splice_complexity(LC)\
        agn_gene_locus_splice_complexity(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_splice_complexity(LC)\
        agn_gene_locus_splice_complexity(LC, REFERENCESOURCE)
#define agn_gene_locus_calc_splice_complexity(LC)\
        agn_gene_locus_splice_complexity(LC, DEFAULTSOURCE)


/**
 * Print the locus in GFF3 format
 *
 * @param[in] locus        the locus
 * @param[in] outstream    the file to which the locus will be printed
 * @param[in] source       source string to use as second column of GFF3 output;
 *                         if NULL is provided, "AEGeAn" will be used
 */
void agn_gene_locus_to_gff3(AgnGeneLocus *locus, FILE *outstream,
                            const char *source);

/**
 * Get the transcripts associated with this locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * `agn_gene_locus_pred_transcripts(locus)' to retrieve prediction transcripts,
 * `agn_gene_locus_refr_transcripts(locus)' to retrieve reference transcripts,
 * or `agn_gene_locus_get_genes(locus)' if the source of annotation is
 * undesignated or irrelevant.
 *
 * @param[in] locus    the locus
 * @param[in] src      REFERENCESOURCE will return only reference transcripts,
 *                     PREDICTIONSOURCE will return only prediction transcripts,
 *                     DEFAULTSOURCE will return all transcripts
 * @returns            the prediction transcripts
 */
GtArray *agn_gene_locus_transcripts(AgnGeneLocus *locus,
                                    AgnComparisonSource src);
#define agn_gene_locus_pred_transcripts(LC)\
        agn_gene_locus_transcripts(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_transcripts(LC)\
        agn_gene_locus_transcripts(LC, REFERENCESOURCE)
#define agn_gene_locus_get_transcripts(LC)\
        agn_gene_locus_transcripts(LC, DEFAULTSOURCE)

/**
 * Get the transcript IDs associated with this locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * `agn_gene_locus_pred_transcripts(locus)' to retrieve prediction IDs,
 * `agn_gene_locus_refr_transcripts(locus)' to retrieve reference IDs,
 * or `agn_gene_locus_get_genes(locus)' if the source of annotation is
 * undesignated or irrelevant.
 *
 * @param[in] locus    the locus
 * @param[in] src      REFERENCESOURCE will return only reference IDs,
 *                     PREDICTIONSOURCE will return only prediction IDs,
 *                     DEFAULTSOURCE will return all transcript IDs
 * @returns            the prediction transcripts
 */
GtArray *agn_gene_locus_transcript_ids(AgnGeneLocus *locus,
                                       AgnComparisonSource src);
#define agn_gene_locus_pred_transcript_ids(LC)\
        agn_gene_locus_transcript_ids(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_transcript_ids(LC)\
        agn_gene_locus_transcript_ids(LC, REFERENCESOURCE)
#define agn_gene_locus_get_transcript_ids(LC)\
        agn_gene_locus_transcript_ids(LC, DEFAULTSOURCE)

#endif
