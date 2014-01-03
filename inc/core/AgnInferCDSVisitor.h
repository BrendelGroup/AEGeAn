#ifndef AEGEAN_INFER_CDS_STREAM
#define AEGEAN_INFER_CDS_STREAM

#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnInferCDSVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for inferring an mRNA's CDS from explicitly defined exon and
 * start/stop codon features.
 */
typedef struct AgnInferCDSVisitor AgnInferCDSVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream* agn_infer_cds_stream_new(GtNodeStream *in, GtLogger *logger);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_infer_cds_visitor_new(GtLogger *logger);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test);

#endif
