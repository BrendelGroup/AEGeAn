/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_REMOVE_CHILDREN_VISITOR
#define AEGEAN_REMOVE_CHILDREN_VISITOR

#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnRemoveChildrenVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for correcting removing all children of each top-level feature.
 * Psuedo-features are not modified.
 */
typedef struct AgnRemoveChildrenVisitor AgnRemoveChildrenVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_remove_children_stream_new(GtNodeStream *in);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_remove_children_visitor_new();

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_remove_children_visitor_unit_test(AgnUnitTest *test);

#endif
