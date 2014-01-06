#ifndef AEGEAN_INFER_EXONS_VISITOR
#define AEGEAN_INFER_EXONS_VISITOR

#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnInferExonsVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for inferring exon features when only CDS and UTR features are
 * provided explicitly. 
 */
typedef struct AgnInferExonsVisitor AgnInferExonsVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_infer_exons_stream_new(GtNodeStream *in, GtLogger *logger);

/**
 * @function Class constructor for the node visitor.
 */
GtNodeVisitor* agn_infer_exons_visitor_new(GtLogger *logger);

/**
 * @function Run unit tests for this class.
 */
bool agn_infer_exons_visitor_unit_test(AgnUnitTest *test);

#endif
