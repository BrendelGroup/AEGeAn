/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_ATTRIBUTE_FILTER_STREAM
#define AEGEAN_ATTRIBUTE_FILTER_STREAM

#include "extended/node_stream_api.h"
#include "core/hashmap_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnFilterStream
 *
 * Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream
 * used to remove features with certain attributes from a node stream.
 */
typedef struct AgnAttributeFilterStream AgnAttributeFilterStream;

/**
 * @function Class constructor. The keys of the `filters` hashmap should be the
 * attribute keys/value pairs (such as `partial=true` or `pseudo=true`) to test
 * each feature node against. The values associated with each key in the hashmap
 * can be any non-NULL value. Any feature node having an attribute key/value
 * pair matching an entry in the hashmap will be discarded.
 */
GtNodeStream* agn_attribute_filter_stream_new(GtNodeStream *in_stream,
                                              GtHashmap *filters);

/**
 * @function Run unit tests for this class.
 */
bool agn_attribute_filter_stream_unit_test(AgnUnitTest *test);

#endif
