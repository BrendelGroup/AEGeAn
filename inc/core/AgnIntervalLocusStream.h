#ifndef AEGEAN_INTERVAL_LOCUS_STREAM
#define AEGEAN_INTERVAL_LOCUS_STREAM

#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnIntervalLocusStream
 *
 * Implements the ``GtNodeStream`` interface. Input is a stream of
 * gene/transcript loci and output is a stream of interval loci. See online docs
 * for more information about interval loci (iLoci).
 */
typedef struct AgnIntervalLocusStream AgnIntervalLocusStream;

/**
 * @function Class constructor. The delta parameter specifies how far beyond
 * each transcript the iLocus boundaries should extend, and the minimum length
 * of an iLocus containing no transcripts. See the online docs for a complete
 * description of iLoci.
 */
GtNodeStream *agn_interval_locus_stream_new(GtNodeStream *locus_stream,
                                            GtUword delta,
                                            bool skipterminal,
                                            GtLogger *logger);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_interval_locus_stream_unit_test(AgnUnitTest *test);

#endif