#ifndef AEGEAN_TRANSCRIPT_STREAM
#define AEGEAN_TRANSCRIPT_STREAM

#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnTranscriptStream
 *
 * Implements the ``GtNodeStream`` interface. Searches the complete
 * feature graph of each feature node in the input for transcript features
 * (mRNA, rRNA, or tRNA). Only transcripts that pass validation are delivered--
 * warning messages for all other transcripts are printed to the console.
 */
typedef struct AgnTranscriptStream AgnTranscriptStream;

/**
 * @function Class constructor.
 */
GtNodeStream* agn_transcript_stream_new(GtNodeStream *in_stream,
                                        GtLogger *logger);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_transcript_stream_unit_test(AgnUnitTest *test);

#endif
