#ifndef AEGEAN_TRANSCRIPT_CLIQUE
#define AEGEAN_TRANSCRIPT_CLIQUE

#include "genometools.h"
#include "AgnUnitTest.h"

/**
 * The purpose of the AgnTranscriptClique class is to store data pertaining to
 * an individual maximal transcript clique. This clique may only contain a
 * single transcript, or it may contain many. The only stipulation is that the
 * transcripts do not overlap.
 */
typedef struct AgnTranscriptClique AgnTranscriptClique;

/**
 * Prototype functions much match to be applied to each transcript in the given
 * clique.
 */
typedef void (*AgnCliqueVisitFunc)(GtFeatureNode*, void*);

/**
 * Add a transcript to this clique.
 *
 * @param[out] clique        the clique
 * @param[in]  transcript    the transcript feature
 */
void agn_transcript_clique_add(AgnTranscriptClique *clique,
                               GtFeatureNode *transcript);

/**
 * Get the CDS length for this transcript clique.
 *
 * @param[in] clique    the clique
 * @returns             the length of the CDS in codons/amino acids
 */
unsigned long agn_transcript_clique_cds_length(AgnTranscriptClique *clique);

/**
 * Make a copy of this transcript clique.
 *
 * @param[in] clique    the clique
 * @returns             a copy of the clique
 */
AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique);

/**
 * Free the memory previously held by this clique.
 *
 * @param[in] clique    the clique
 */
void agn_transcript_clique_delete(AgnTranscriptClique *clique);

/**
 * Determine whether any of the transcript IDs associated with this clique are
 * contained in the given hash map.
 *
 * @param[in] clique    the clique
 * @param[in] map       a hashmap containing a set of IDs
 * @returns             true if the hashmap contains any of the clique's IDs,
 *                      false otherwise
 */
bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique,
                                          GtHashmap *map);

/**
 * FIXME
 */
const char *agn_transcript_clique_id(AgnTranscriptClique *clique);

/**
 * Allocate some memory for a new transcript clique.
 *
 * @returns    a pointer to the newly allocated memory
 */
AgnTranscriptClique* agn_transcript_clique_new();

/**
 * Get the number of exons in this clique.
 *
 * @param[in] clique    the clique
 * @returns             the number of exons in the clique
 */
unsigned long agn_transcript_clique_num_exons(AgnTranscriptClique *clique);

/**
 * Get the number of UTRs in this clique.
 *
 * @param[in] clique    the clique
 * @returns             the number of UTRs in the clique
 */
unsigned long agn_transcript_clique_num_utrs(AgnTranscriptClique *clique);

/**
 * Print this clique's IDs to the given output stream.
 *
 * @param[in]  clique       the clique
 * @param[out] outstream    the output stream to which the IDs will be printed
 */
void agn_transcript_clique_print_ids(AgnTranscriptClique *clique,
                                     FILE *outstream);

/**
 * Add all of the IDs associated with this clique to the given hash map.
 *
 * @param[in]  clique    the clique
 * @param[out] map       the hash map to which the IDs will be written
 */
void agn_transcript_clique_put_ids_in_hash(AgnTranscriptClique *clique,
                                           GtHashmap *map);

/**
 * Get the number of transcripts in this clique.
 *
 * @param[in] clique    the clique
 * @returns             the number of transcripts in the clique
 */
unsigned long agn_transcript_clique_size(AgnTranscriptClique *clique);

/**
 * Get an array containing all the transcripts in this clique.
 *
 * @param[in] clique    the clique
 * @returns             an array of all the transcripts
 */
GtArray* agn_transcript_clique_to_array(AgnTranscriptClique *clique);

/**
 * Print the transcript clique to the given outstream in GFF3 format.
 *
 * @param[in] clique       the clique
 * @param[in] outstream    file to which output will be written
 * @param[in] prefix       optional prefix to print at the beginning of each
 *                         line of output; if NULL is provided, not prefix is
 *                         printed
 */
void agn_transcript_clique_to_gff3(AgnTranscriptClique *clique, FILE *outstream,
                                   const char *prefix);

/**
 * Apply the given function to each transcript in the clique.
 *
 * @param[in] clique      transcript clique
 * @param[in] func        function to be applied to each transcript
 * @param[in] funcdata    optional pointer to data to be used in the function
 */
void agn_transcript_clique_traverse(AgnTranscriptClique *clique,
                                    AgnCliqueVisitFunc func, void *funcdata);

/**
 * Run unit tests for this class.
 *
 * @param[out] test    object for storing test results
 * @returns            true for success, false for failure
 */
bool agn_transcript_clique_unit_test(AgnUnitTest *test);

#endif
