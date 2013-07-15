#ifndef AEGEAN_SIMPLE_NODE_VISITOR
#define AEGEAN_SIMPLE_NODE_VISITOR

#include "genometools.h"
#include "AgnGeneValidator.h"

typedef struct AgnSimpleNodeVisitor AgnSimpleNodeVisitor;

/**
 * Allocate memory for a genome node visitor used to validate and load data into
 * memory.
 *
 * @param[out] index     the feature index into which data will be loaded
 * @param[in]  type      the feature type to load into memory
 * @param[in]  logger    object for storing error/warning messages
 * @returns              a node visitor object
 */
GtNodeVisitor* agn_simple_node_visitor_new(GtFeatureIndex *index,
                                           const char *type,
                                           AgnLogger *logger);

#endif
