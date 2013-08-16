#ifndef AEGEAN_CANON_GENE_STREAM
#define AEGEAN_CANON_GENE_STREAM

#include "extended/node_stream_api.h"

/**
 * FIXME
 */
typedef struct AgnCanonGeneStream AgnCanonGeneStream;

/**
 * FIXME
 */
const GtNodeStreamClass* agn_canon_gene_stream_class(void);

/**
 * FIXME
 */
GtNodeStream* agn_canon_gene_stream_new(GtNodeStream *in_stream);

#endif