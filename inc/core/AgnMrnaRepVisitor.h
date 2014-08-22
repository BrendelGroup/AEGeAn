/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_MRNA_REP_VISITOR
#define AEGEAN_MRNA_REP_VISITOR

#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnMrnaRepVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for filtering out all but the longest mRNA (as measured by CDS
 * length) from alternatively spliced genes.
 */
typedef struct AgnMrnaRepVisitor AgnMrnaRepVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_mrna_rep_stream_new(GtNodeStream *in);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_mrna_rep_visitor_new();

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_mrna_rep_visitor_unit_test(AgnUnitTest *test);

#endif
