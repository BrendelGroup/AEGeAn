#ifndef AEGEAN_INFER_EXONS_VISITOR
#define AEGEAN_INFER_EXONS_VISITOR

#include "genometools.h"
#include "AgnLogger.h"

/**
 * @class AgnInferExonsVisitor
 *
 * FIXME
 */
typedef struct AgnInferExonsVisitor AgnInferExonsVisitor;

/**
 * A node visitor used for inferring an mRNA's CDS from explicitly defined exon
 * and start/stop codon features.
 *
 * @param[in]  logger    object for storing error/warning messages
 * @returns              a node visitor object
 */
GtNodeVisitor* agn_infer_exons_visitor_new(AgnLogger *logger);

#endif
