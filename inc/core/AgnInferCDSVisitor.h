#ifndef AEGEAN_INFER_CDS_VISITOR
#define AEGEAN_INFER_CDS_VISITOR

#include "genometools.h"
#include "AgnLogger.h"

typedef struct AgnInferCDSVisitor AgnInferCDSVisitor;

/**
 * A node visitor used for inferring an mRNA's CDS from explicitly defined exon
 * and start/stop codon features.
 *
 * @param[in]  logger    object for storing error/warning messages
 * @returns              a node visitor object
 */
GtNodeVisitor* agn_infer_cds_visitor_new(AgnLogger *logger);

#endif
