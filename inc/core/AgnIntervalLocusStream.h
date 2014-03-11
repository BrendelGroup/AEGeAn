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
 * of an iLocus containing no transcripts. If ``endmode == 0``, all iLoci will
 * be included in the output; if ``endmode < 0``, terminal iLoci will not be
 * included in the output; and if ``endmode > 0``, then _only_ terminal iLoci
 * will be included in the output. See the online docs for a complete
 * description of iLoci.
 */
GtNodeStream *agn_interval_locus_stream_new(GtNodeStream *locus_stream,
                                            GtUword delta, int endmode,
                                            GtLogger *logger);

/**
 * @function Set the source value to be used for all iLoci created by this
 * stream. Default value is 'AEGeAn::AgnIntervalLocusStream'.
 */
void agn_interval_locus_stream_set_source(AgnIntervalLocusStream *stream,
                                          GtStr *source);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_interval_locus_stream_unit_test(AgnUnitTest *test);

#endif