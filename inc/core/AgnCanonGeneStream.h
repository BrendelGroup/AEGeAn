#ifndef AEGEAN_CANON_GENE_STREAM
#define AEGEAN_CANON_GENE_STREAM

#include "extended/node_stream_api.h"
#include "AgnLogger.h"

/**
 * A node stream that returns only canonical protein-coding genes that pass
 * stringent validation.
 */
typedef struct AgnCanonGeneStream AgnCanonGeneStream;

/**
 * Constructor for the AgnCanonGeneStream class.
 *
 * @param[in] in_stream    the node stream feeding nodes to this stream
 * @param[in] logger       stores error and warning messages
 * @returns                a new node stream that will feed only gene feature
 *                         nodes
 */
GtNodeStream* agn_canon_gene_stream_new(GtNodeStream *in_stream,
                                        AgnLogger *logger);

#endif