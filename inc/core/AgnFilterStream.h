/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_FILTER_STREAM
#define AEGEAN_FILTER_STREAM

#include "extended/node_stream_api.h"
#include "core/hashmap_api.h"
#include "AgnUnitTest.h"

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

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_filter_stream_unit_test(AgnUnitTest *test);

#endif
