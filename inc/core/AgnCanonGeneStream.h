#ifndef AEGEAN_CANON_GENE_STREAM
#define AEGEAN_CANON_GENE_STREAM

#include "extended/node_stream_api.h"
#include "AgnLogger.h"

/**
 * @class AgnCanonGeneStream
 * 
 * Implements the GenomeTools ``GtNodeStream`` interfact. This is a node stream
 * that returns only canonical protein-coding genes that pass stringent
 * validation.
 */
typedef struct AgnCanonGeneStream AgnCanonGeneStream;

/**
 * @function Class constructor.
 */
GtNodeStream* agn_canon_gene_stream_new(GtNodeStream *in_stream,
                                        AgnLogger *logger);

#endif
