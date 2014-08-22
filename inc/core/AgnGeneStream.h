/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_GENE_STREAM
#define AEGEAN_GENE_STREAM

#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnGeneStream
 *
 * Implements the ``GtNodeStream`` interface. Searches the complete feature
 * graph of each feature node in the input for canonical protein-coding gene
 * features. Some basic sanity checks are performed on the mRNA(s) associated
 * with each gene, and genes are only delivered to the output stream if they
 * include one or more valid mRNA subfeatures.
 */
typedef struct AgnGeneStream AgnGeneStream;

/**
 * @function Class constructor.
 */
GtNodeStream* agn_gene_stream_new(GtNodeStream *in_stream, GtLogger *logger);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_gene_stream_unit_test(AgnUnitTest *test);

#endif
