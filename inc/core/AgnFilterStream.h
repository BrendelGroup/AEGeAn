#ifndef AEGEAN_FILTER_STREAM
#define AEGEAN_FILTER_STREAM

#include "extended/node_stream_api.h"
#include "core/hashmap_api.h"

/**
 * FIXME
 */
typedef struct AgnFilterStream AgnFilterStream;

/**
 * FIXME
 */
const GtNodeStreamClass* agn_filter_stream_class(void);

/**
 * FIXME
 */
GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream,
                                    GtHashmap *typestokeep,
                                    GtHashmap *typestofilter);

#endif