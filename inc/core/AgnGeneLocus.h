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
 * Do a shallow copy of this data structure.
 *
 * @param[in] locus    a locus object
 * @returns            a clone of the locus object, all of whose internal data
 *                     point to the same objects
 */
AgnGeneLocus *agn_gene_locus_clone(AgnGeneLocus *locus);

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
 * The combined length of all coding sequences associated with this locus.
 * Rather than calling this function directly, users are encouraged to use one
 * of the following macros: `agn_gene_locus_refr_cds_length(locus)' for the
 * combined length of all reference CDSs,
 * `agn_gene_locus_pred_cds_length(locus)' for the combined length of all
 * prediction CDSs, and `agn_gene_locus_get_cds_length(locus)' for the combined
 * length of all CDSs.
 *
 * @param[in] locus    the locus
 * @param[in] src      REFERENCESOURCE will consider only reference CDSs,
 *                     PREDICTIONSOURCE will consider only prediction CDSs,
 *                     DEFAULTSOURCE will consider all CDSs
 * @returns            the combined CDS length
 */
unsigned long agn_gene_locus_cds_length(AgnGeneLocus *locus,
                                        AgnComparisonSource src);
#define agn_gene_locus_pred_cds_length(LC)\
        agn_gene_locus_cds_length(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_cds_length(LC)\
        agn_gene_locus_cds_length(LC, REFERENCESOURCE)
#define agn_gene_locus_get_cds_length(LC)\
        agn_gene_locus_cds_length(LC, DEFAULTSOURCE)

/**
 * Compare every reference transcript clique with every prediction transcript
 * clique. For gene loci with multiple transcript cliques, each comparison is
 * not necessarily reported. Instead, we report the set of clique pairs that
 * provides the optimal pairing of reference and prediction transcripts. If
 * there are more reference transcript cliques than prediction cliques (or vice
 * versa), these unmatched cliques are reported separately.
 *
 * @param[in]  locus    the locus
 * @returns             the pairs to be reported; the first time this method is
 *                      called, the comparative analysis is performed and the
 *                      clique pairs are stored and returned; subsequence method
 *                      calls simply return the previously stored results;
 *                      the macro `agn_gene_locus_pairs_to_report(locus)' is
 *                      provided for this use case
 */
GtArray *agn_gene_locus_comparative_analysis(AgnGeneLocus *locus);
#define agn_gene_locus_pairs_to_report(LOC)\
        agn_gene_locus_comparative_analysis(LOC)

/**
 * Free the memory previously occupied by this locus object.
 *
 * @param[out] locus    locus to be deleted
 */
void agn_gene_locus_delete(AgnGeneLocus *locus);

/**
 * We use the Bron-Kerbosch algorithm to enumerate maximal cliques of non-
 * overlapping transcripts. This is done separately for the reference
 * transcripts and the prediction transcripts. Every possible pairing of
 * reference cliques and prediction cliques is enumerated, to enable subsequent
 * pairwise comparison.
 *
 * @param[in] locus    the locus
 * @returns            the number of clique pairs formed; the first
 *                     time this method is called, the clique pairs are
 *                     enumerated and stored and their number is returned;
 *                     subsequent method calls simply return the number of
 *                     previously enumerated cliques; the macro
 *                     `agn_gene_locus_num_clique_pairs(locus)' has been
 *                     provided for this use case
 */
unsigned long agn_gene_locus_enumerate_clique_pairs(AgnGeneLocus *locus);
#define agn_gene_locus_num_clique_pairs(LOC)\
        agn_gene_locus_enumerate_clique_pairs(LOC)

/**
 * Get the number of exons for the locus. Rather than calling this function
 * directly, users are encouraged to use one of the following macros:
 * `agn_gene_locus_num_pred_exons(locus)' for the number of prediction exons,
 * `agn_gene_locus_num_refr_exons(locus)' for the number of reference exons, or
 * `agn_gene_locus_num_exons(locus)' if the source of annotation is undesignated
 * or irrelevant.
 *
 * @param[in] locus   the locus
 * @param[in] src     REFERENCESOURCE will consider only reference exons,
 *                    PREDICTIONSOURCE will consider only prediction exons,
 *                    DEFAULTSOURCE will consider all exons
 * @returns           the number of exons associated with the locus
 */
unsigned long agn_gene_locus_exon_num(AgnGeneLocus *locus,
                                      AgnComparisonSource src);
#define agn_gene_locus_num_pred_exons(LC)\
        agn_gene_locus_exon_num(LC, PREDICTIONSOURCE)
#define agn_gene_locus_num_refr_exons(LC)\
        agn_gene_locus_exon_num(LC, REFERENCESOURCE)
#define agn_gene_locus_num_exons(LC)\
        agn_gene_locus_exon_num(LC, DEFAULTSOURCE)

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
 * @returns            array containing the gene IDs
 */
GtArray *agn_gene_locus_gene_ids(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_pred_gene_ids(LC)\
        agn_gene_locus_gene_ids(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_gene_ids(LC)\
        agn_gene_locus_gene_ids(LC, REFERENCESOURCE)
#define agn_gene_locus_get_gene_ids(LC)\
        agn_gene_locus_gene_ids(LC, DEFAULTSOURCE)

/**
 * Get the number of genes for the locus. Rather than calling this function
 * directly, users are encouraged to use one of the following macros:
 * `agn_gene_locus_num_pred_genes(locus)' for the number of prediction genes,
 * `agn_gene_locus_num_refr_genes(locus)' for the number of reference genes, or
 * `agn_gene_locus_num_genes(locus)' if the source of annotation is undesignated
 * or irrelevant.
 *
 * @param[in] locus   the locus
 * @param[in] src     REFERENCESOURCE will consider only reference genes,
 *                    PREDICTIONSOURCE will consider only prediction genes,
 *                    DEFAULTSOURCE will consider all genes
 * @returns           the number of genes associated with the locus
 */
unsigned long agn_gene_locus_gene_num(AgnGeneLocus *locus,
                                      AgnComparisonSource src);
#define agn_gene_locus_num_pred_genes(LC)\
        agn_gene_locus_gene_num(LC, PREDICTIONSOURCE)
#define agn_gene_locus_num_refr_genes(LC)\
        agn_gene_locus_gene_num(LC, REFERENCESOURCE)
#define agn_gene_locus_num_genes(LC)\
        agn_gene_locus_gene_num(LC, DEFAULTSOURCE)

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
 * Return the coordinates of this locus.
 *
 * @param[in] locus    the locus
 * @returns            its coordinates
 */
GtRange agn_gene_locus_range(AgnGeneLocus *locus);

/**
 * Set the range of this locus, no questions asked.
 *
 * @param[out] locus    locus object
 * @param[in]  start    new start coordinate
 * @param[in]  end      new end coordinate
 */
void agn_gene_locus_set_range(AgnGeneLocus *locus, unsigned long start,
                              unsigned long end);

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

/**
 * Get the number of transcripts for the locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * `agn_transcript_locus_num_pred_transcripts(locus)' for the number of
 * prediction transcripts, `agn_transcript_locus_num_refr_transcripts(locus)'
 * for the number of reference transcripts, or
 * `agn_transcript_locus_num_transcripts(locus)' if the source of annotation is
 * undesignated or irrelevant.
 *
 * @param[in] locus   the locus
 * @param[in] src     REFERENCESOURCE will consider only reference transcripts,
 *                    PREDICTIONSOURCE will consider only prediction
 *                    transcripts, DEFAULTSOURCE will consider all transcripts
 * @returns           the number of transcripts associated with the locus
 */
unsigned long agn_gene_locus_transcript_num(AgnGeneLocus *locus,
                                                  AgnComparisonSource src);
#define agn_gene_locus_num_pred_transcripts(LC)\
        agn_gene_locus_transcript_num(LC, PREDICTIONSOURCE)
#define agn_gene_locus_num_refr_transcripts(LC)\
        agn_gene_locus_transcript_num(LC, REFERENCESOURCE)
#define agn_gene_locus_num_transcripts(LC)\
        agn_gene_locus_transcript_num(LC, DEFAULTSOURCE)

#endif
