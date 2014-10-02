/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef AEGEAN_AS_INSPECT_IR_VISITOR
#define AEGEAN_AS_INSPECT_IR_VISITOR

#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnASInspectIRVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for summarizing the extent of intron retention in an annotation.
 */
typedef struct AgnASInspectIRVisitor AgnASInspectIRVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_as_inspect_ir_stream_new(GtNodeStream *in, FILE *report);

/**
 * @function Constructor. If ``report`` == NULL, output will be written to
 * terminal (stdout).
 */
GtNodeVisitor *agn_as_inspect_ir_visitor_new(FILE *report);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_as_inspect_ir_visitor_unit_test(AgnUnitTest *test);

#endif
