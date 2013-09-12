#ifndef AEGEAN_FILTER_STREAM
#define AEGEAN_FILTER_STREAM

#include "extended/node_stream_api.h"
#include "core/hashmap_api.h"

/**
 * @class AgnFilterStream
 *
 * A node stream used to filter out unwanted nodes. The user can either specify
 * the types of nodes to keep (all others will be deleted and not passed), or
 * the types of nodes not to keep (these will be deleted, all others will be
 * passed).
 */
typedef struct AgnFilterStream AgnFilterStream;

/**
 * Constructor for the AgnCanonGeneStream class.
 *
 * @param[in] in_stream        the node stream feeding nodes to this stream
 * @param[in] typestokeep      hashmap containing node types (as char *s) to
 *                             keep; if this argument is provided, typestofilter
 *                             must be NULL
 * @param[in] typestofilter    hashmap containing node types (as char *s) to
 *                             discard; if this argument is provided,
 *                             typestokeep must be NULL
 * @returns                    a new node stream that will feed only feature
 *                             nodes that meet the given type criteria
 */
GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream,
                                    GtHashmap *typestokeep,
                                    GtHashmap *typestofilter);

#endif
