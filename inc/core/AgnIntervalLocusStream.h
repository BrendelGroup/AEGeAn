/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
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
 * included in the output; and if ``endmode > 0``, then *only* terminal iLoci
 * will be included in the output. See the online docs for a complete
 * description of iLoci.
 */
GtNodeStream *agn_interval_locus_stream_new(GtNodeStream *locus_stream,
                                            GtUword delta, int endmode,
                                            GtLogger *logger);

/**
 * @function iLoci created by this stream are assigned an ID with an arbitrary
 * number. The default format is 'iLocus%lu' (that is, iLocus1, iLocus2, etc).
 * Use this function to override the default ID format.
 */
void agn_interval_locus_stream_set_idformat(AgnIntervalLocusStream *stream,
                                            const char *format);

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
