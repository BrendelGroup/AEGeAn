#ifndef AEGEAN_LOCUS_STREAM
#define AEGEAN_LOCUS_STREAM

#include "core/hashmap_api.h"
#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnLocusStream
 *
 * Implements the ``GtNodeStream`` interface. The only feature nodes delivered
 * by this stream have type ``locus``, and the only direct children of these
 * features are transcript features (of types mRNA, rRNA, or tRNA) present in
 * the input stream. Any overlapping transcripts are children of the same locus
 * feature.
 */
typedef struct AgnLocusStream AgnLocusStream;


/**
 * @function This constructor searches the complete feature graph of each
 * feature node in the input stream for transcript features.
 */
GtNodeStream *agn_locus_stream_new(GtNodeStream *in_stream, GtLogger *logger);

/**
 * @function This constructor accepts two :c:type:`AgnTranscriptStream` objects
 * as input. Locus features are created as per the class description, with
 * additional data stored to track the source (reference vs prediction) of each
 * transcript in each locus.
 */
GtNodeStream *agn_locus_stream_new_pairwise(GtNodeStream *refr_stream,
                                            GtNodeStream *pred_stream,
                                            GtLogger *logger);

/**
 * @function Set the source value to be used for all iLoci created by this
 * stream. Default value is 'AEGeAn::AgnLocusStream'.
 */
void agn_locus_stream_set_source(AgnLocusStream *stream, GtStr *source);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_locus_stream_unit_test(AgnUnitTest *test);

#endif
