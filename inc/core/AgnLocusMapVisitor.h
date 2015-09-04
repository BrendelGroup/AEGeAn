/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_LOCUS_MAP_VISITOR
#define AEGEAN_LOCUS_MAP_VISITOR

#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnLocusMapVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for printing out gene --> locus and mRNA --> locus relationships
 * as part of a locus/iLocus processing stream.
 */
typedef struct AgnLocusMapVisitor AgnLocusMapVisitor;

/**
 * @function Constructor for a node stream based on this node visitor. See
 * :c:func:`agn_locus_map_visitor_new` for a description of the function
 * arguments.
 */
GtNodeStream*
agn_locus_map_stream_new(GtNodeStream *in, FILE *genefh, FILE *mrnafh,
                         bool usename);

/**
 * @function Constructor for the node visitor. Gene-to-locus relationships are
 * printed to the ``genefh`` file handle, while mRNA-to-locus relationships are
 * printed to the ``mrnafh`` file handle. Setting either file handle to NULL
 * will disable printing the corresponding output.
 */
GtNodeVisitor *agn_locus_map_visitor_new(FILE *genefh, FILE *mrnafh);

/**
 * @function Report gene/mRNA `Name` attributes rather than ID attributes.
 */
void agn_locus_map_visitor_use_name(AgnLocusMapVisitor *mv);

#endif
