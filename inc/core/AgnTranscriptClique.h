#ifndef AEGEAN_TRANSCRIPT_CLIQUE
#define AEGEAN_TRANSCRIPT_CLIQUE

#include "extended/feature_node_api.h"
#include "core/hashmap_api.h"
#include "AgnUnitTest.h"
#include "AgnUtils.h"

/**
 * @class AgnTranscriptClique
 *
 * The purpose of the AgnTranscriptClique class is to store data pertaining to
 * an individual maximal transcript clique. This clique may only contain a
 * single transcript, or it may contain many. The only stipulation is that the
 * transcripts do not overlap. Under the hood, each ``AgnTranscriptClique``
 * instance is a pseudo node (a GtFeatureNode object) with one or more
 * transcript features as direct children.
 */
typedef GtGenomeNode AgnTranscriptClique;

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
 * @function Get the combined CDS length (in base pairs) for all transcripts in
 * this clique.
 */
GtUword agn_transcript_clique_cds_length(AgnTranscriptClique *clique);

/**
 * @function Make a shallow copy of this transcript clique.
 */
AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique);

/**
 * @function Class destructor.
 */
void agn_transcript_clique_delete(AgnTranscriptClique *clique);

/**
 * @function Get a pointer to the string representing this clique's transcript
 * structure.
 */
const char *agn_transcript_clique_get_model_vector(AgnTranscriptClique *clique);

/**
 * @function Determine whether any of the transcript IDs associated with this
 * clique are keys in the given hash map.
 */
bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique,
                                          GtHashmap *map);

/**
 * @function Retrieve the ID attribute of the transcript associated with this
 * clique. User is responsible to free the string.
 */
char *agn_transcript_clique_id(AgnTranscriptClique *clique);

/**
 * @function Retrieve the ID attributes of all transcripts associated with this
 * clique.
 */
GtArray *agn_transcript_clique_ids(AgnTranscriptClique *clique);

/**
 * @function Class constructor. ``locusrange`` should be a pointer to the
 * genomic coordinates of the locus to which this transcript clique belongs.
 */
AgnTranscriptClique *agn_transcript_clique_new(AgnSequenceRegion *region);

/**
 * @function Get the number of exons in this clique.
 */
GtUword agn_transcript_clique_num_exons(AgnTranscriptClique *clique);

/**
 * @function Get the number of UTR segments in this clique.
 */
GtUword agn_transcript_clique_num_utrs(AgnTranscriptClique *clique);

/**
 * @function Add all of the IDs associated with this clique to the given hash
 * map.
 */
void agn_transcript_clique_put_ids_in_hash(AgnTranscriptClique *clique,
                                           GtHashmap *map);

/**
 * @function Get the number of transcripts in this clique.
 */
GtUword agn_transcript_clique_size(AgnTranscriptClique *clique);

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
