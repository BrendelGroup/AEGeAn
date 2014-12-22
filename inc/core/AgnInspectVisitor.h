/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_INSPECT_VISITOR
#define AEGEAN_INSPECT_VISITOR

#include "extended/node_stream_api.h"

/**
 * @class AgnInspectVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. Prints out
 * debugging information for each node it encounters in the node stream.
 */
typedef struct AgnInspectVisitor AgnInspectVisitor;

/**
 * @function Constructor for a node stream based on this node visitor. See
 * :c:func:`agn_inspect_visitor_new` for a description of the function
 * arguments.
 */
GtNodeStream*
agn_inspect_stream_new(GtNodeStream *in, FILE *outstream);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_inspect_visitor_new(FILE *outstream);

#endif

