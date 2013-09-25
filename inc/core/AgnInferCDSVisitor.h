#ifndef AEGEAN_INFER_CDS_VISITOR
#define AEGEAN_INFER_CDS_VISITOR

#include "genometools.h"
#include "AgnLogger.h"
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
 * @function Class constructor.
 */
GtNodeVisitor* agn_infer_cds_visitor_new(AgnLogger *logger);

/**
 * @function Run unit tests for this class.
 */
bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test);

#endif
