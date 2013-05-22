#ifndef AEGEAN_CANONICAL_NODE_VISITOR
#define AEGEAN_CANONICAL_NODE_VISITOR

#include "genometools.h"
#include "AgnGeneValidator.h"

typedef struct AgnCanonNodeVisitor AgnCanonNodeVisitor;

/**
 * Allocate memory for a genome node visitor used to validate and load data into
 * memory.
 *
 * @param[out] index        the feature index into which data will be loaded
 * @param[in]  validator    a gene validator object
 * @param[in]  logger       object for storing error/warning messages
 * @returns                 a node visitor object
 */
GtNodeVisitor* agn_canon_node_visitor_new(GtFeatureIndex *index,
                                          AgnGeneValidator *validator,
                                          AgnLogger *logger);

#endif
