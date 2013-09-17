#ifndef AEGEAN_INFER_EXONS_VISITOR
#define AEGEAN_INFER_EXONS_VISITOR

#include "genometools.h"
#include "AgnLogger.h"

/**
 * @class AgnInferExonsVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node
 * visitor used for inferring exon features when only CDS and UTR features are
 * provided explicitly. 
 */
typedef struct AgnInferExonsVisitor AgnInferExonsVisitor;

/**
 * @function Class constructor.
 */
GtNodeVisitor* agn_infer_exons_visitor_new(AgnLogger *logger);

#endif
