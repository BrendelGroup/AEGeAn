/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_NODE_DELETE_VISITOR
#define AEGEAN_NODE_DELETE_VISITOR

#include "extended/node_stream_api.h"

/**
 * @class AgnNodeDeleteVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used to decrement the reference count to all feature nodes passing
 * through the node stream.
 */
typedef struct AgnNodeDeleteVisitor AgnNodeDeleteVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_node_delete_stream_new(GtNodeStream *in);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_node_delete_visitor_new();

#endif
