/**

Copyright (c) 2017, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_ID_FILTER_STREAM
#define AEGEAN_ID_FILTER_STREAM

#include "extended/node_stream_api.h"
#include "core/hashmap_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnIdFilterStream
 *
 * Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream
 * used to select features from a node stream using a pre-specified list of IDs.
 */
typedef struct AgnIdFilterStream AgnIdFilterStream;

/**
 * @function Class constructor. The keys of the ``ids2keep`` hashmap should
 * be strings of the IDs of features to be kept from the node stream. Any
 * non-NULL value can be associated with those keys.
 */
GtNodeStream* agn_id_filter_stream_new(GtNodeStream *in_stream,
                                       GtHashmap *ids2keep);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_id_filter_stream_unit_test(AgnUnitTest *test);

#endif
