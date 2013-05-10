#include "genometools.h"
#include "AgnGeneValidator.h"

typedef struct PeNodeVisitor PeNodeVisitor;

/**
 * Allocate memory for a genome node visitor used to validate and load data into
 * memory.
 *
 * @param[out] index        the feature index into which data will be loaded
 * @param[in]  validator    a gene validator object
 * @param[in]  settings     feature inference settings
 * @param[in]  options      values for ParsEval's command-line options
 * @returns                 a node visitor object
 */
GtNodeVisitor* pe_node_visitor_new( GtFeatureIndex *index,
                                    AgnGeneValidator *validator );
