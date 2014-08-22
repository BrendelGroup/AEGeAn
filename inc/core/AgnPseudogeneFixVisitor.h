/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_PSEUDOGENE_FIX_VISITOR
#define AEGEAN_PSEUDOGENE_FIX_VISITOR

#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnPseudogeneFixVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for correcting the ``type`` value for pseudogene features
 * erroneously using the ``gene`` type instead of the more appropriate
 * ``pseudogene`` type.
 */
typedef struct AgnPseudogeneFixVisitor AgnPseudogeneFixVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_pseudogene_fix_stream_new(GtNodeStream *in);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_pseudogene_fix_visitor_new();

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_pseudogene_fix_visitor_unit_test(AgnUnitTest *test);

#endif
