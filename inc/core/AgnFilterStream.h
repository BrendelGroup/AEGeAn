#ifndef AEGEAN_FILTER_STREAM
#define AEGEAN_FILTER_STREAM

#include "extended/node_stream_api.h"
#include "core/hashmap_api.h"

/**
 * @class AgnFilterStream
 *
 * Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream
 * used to filter out unwanted nodes. The user can either specify the types of
 * nodes to keep (all others will be deleted and not passed), or the types of
 * nodes not to keep (these will be deleted, all others will be passed).
 */
typedef struct AgnFilterStream AgnFilterStream;

/**
 * @function Class constructor. Provide ``typestokeep`` or ``typestofilter``,
 * not both.
 */
GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream,
                                    GtHashmap *typestokeep,
                                    GtHashmap *typestofilter);

#endif
