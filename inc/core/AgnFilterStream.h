#ifndef AEGEAN_FILTER_STREAM
#define AEGEAN_FILTER_STREAM

#include "extended/node_stream_api.h"
#include "core/hashmap_api.h"

/**
 * @class AgnFilterStream
 *
 * Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream
 * used to select features of a certain type from a node stream.
 */
typedef struct AgnFilterStream AgnFilterStream;

/**
 * @function Class constructor. The keys of the ``typestokeep`` hashmap should
 * be the type(s) to be kept from the node stream. Any non-NULL value can be
 * associated with those keys.
 */
GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream,
                                    GtHashmap *typestokeep);

#endif
