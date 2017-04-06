/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef AEGEAN_GAEVAL_VISITOR
#define AEGEAN_GAEVAL_VISITOR

#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnGaevalVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for calculating transcript coverage and integrity scores for
 * gene models using alignment data.
 */
typedef struct AgnGaevalVisitor AgnGaevalVisitor;

/**
 * @type Parameters used in calculating GAEVAL integrity.
 * See http://www.plantgdb.org/GAEVAL/docs/integrity.html
 */
struct AgnGaevalParams
{
  double alpha;
  double beta;
  double gamma;
  double epsilon;
  GtUword exp_cds_len;
  GtUword exp_5putr_len;
  GtUword exp_3putr_len;
};
typedef struct AgnGaevalParams AgnGaevalParams;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_gaeval_stream_new(GtNodeStream *in, GtNodeStream *astream,
                                    AgnGaevalParams gparams);

/**
 * @function Class constructor for the node visitor.
 */
GtNodeVisitor*
agn_gaeval_visitor_new(GtNodeStream *astream, AgnGaevalParams gparams);

/**
* @function Indicate a file to be used for printing TSV output.
*/
void agn_gaeval_visitor_tsv_out(AgnGaevalVisitor *v, GtStr *tsvfilename);

/**
 * @function Run unit tests for this class.
 */
bool agn_gaeval_visitor_unit_test(AgnUnitTest *test);

#endif
