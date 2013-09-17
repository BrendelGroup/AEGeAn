#ifndef AEGEAN_TRANSCRIPT_CLIQUE
#define AEGEAN_TRANSCRIPT_CLIQUE

#include "genometools.h"
#include "AgnUnitTest.h"

/**
 * @class AgnTranscriptClique
 *
 * The purpose of the AgnTranscriptClique class is to store data pertaining to
 * an individual maximal transcript clique. This clique may only contain a
 * single transcript, or it may contain many. The only stipulation is that the
 * transcripts do not overlap.
 */
typedef struct AgnTranscriptClique AgnTranscriptClique;

/**
 * @functype
 * The signature that functions must match to be applied to each transcript in
 * the given clique. The function will be called once for each transcript in the
 * clique. The transcript will be passed as the first argument, and a second
 * argument is available for an optional pointer to supplementary data (if
 * needed). See :c:func:`agn_transcript_clique_traverse`.
 */
typedef void (*AgnCliqueVisitFunc)(GtFeatureNode*, void*);

/**
 * @function Add a transcript to this clique.
 */
void agn_transcript_clique_add(AgnTranscriptClique *clique,
                               GtFeatureNode *transcript);

/**
 * @function Get the CDS length (in amino acids) for this transcript clique.
 */
unsigned long agn_transcript_clique_cds_length(AgnTranscriptClique *clique);

/**
 * @function Make a shallow copy of this transcript clique.
 */
AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique);

/**
 * @function Class destructor.
 */
void agn_transcript_clique_delete(AgnTranscriptClique *clique);

/**
 * @function Determine whether any of the transcript IDs associated with this
 * clique are keys in the given hash map.
 */
bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique,
                                          GtHashmap *map);

/**
 * @function Retrieve the ID attribute of the transcript associated with this
 * clique. Will cause an assertion error if there is more than one trancript
 * associated with the clique.
 */
const char *agn_transcript_clique_id(AgnTranscriptClique *clique);

/**
 * @function Class constructor.
 */
AgnTranscriptClique* agn_transcript_clique_new();

/**
 * @function Get the number of exons in this clique.
 */
unsigned long agn_transcript_clique_num_exons(AgnTranscriptClique *clique);

/**
 * @function Get the number of UTR segments in this clique.
 */
unsigned long agn_transcript_clique_num_utrs(AgnTranscriptClique *clique);

/**
 * @function Print the IDs of this clique's transcripts to the given output
 * stream.
 */
void agn_transcript_clique_print_ids(AgnTranscriptClique *clique,
                                     FILE *outstream);

/**
 * @function Add all of the IDs associated with this clique to the given hash
 * map.
 */
void agn_transcript_clique_put_ids_in_hash(AgnTranscriptClique *clique,
                                           GtHashmap *map);

/**
 * @function Get the number of transcripts in this clique.
 */
unsigned long agn_transcript_clique_size(AgnTranscriptClique *clique);

/**
 * @function Get an array containing all the transcripts in this clique. User is
 * responsible for deleting the array.
 */
GtArray* agn_transcript_clique_to_array(AgnTranscriptClique *clique);

/**
 * @function Print the transcript clique to the given outstream in GFF3 format,
 * optionally with a prefix.
 */
void agn_transcript_clique_to_gff3(AgnTranscriptClique *clique, FILE *outstream,
                                   const char *prefix);

/**
 * @function Apply ``func`` to each transcript in the clique. See
 * :c:type:`AgnCliqueVisitFunc`.
 */
void agn_transcript_clique_traverse(AgnTranscriptClique *clique,
                                    AgnCliqueVisitFunc func, void *funcdata);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_transcript_clique_unit_test(AgnUnitTest *test);

#endif
