/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef AEGEAN_LOCUS_STREAM
#define AEGEAN_LOCUS_STREAM

#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnUnitTest.h"

/**
 * @class AgnLocusStream
 *
 * Implements the ``GtNodeStream`` interface. The only feature nodes delivered
 * by this stream have type ``locus``, and the only direct children of these
 * features are gene features present in the input stream. Any overlapping genes
 * are children of the same locus feature.
 */
typedef struct AgnLocusStream AgnLocusStream;


/**
 * @function Use the given filenames to label the direct children of each iLocus
 * as a 'reference' feature or a 'prediction' feature, to facilitate pairwise
 * comparison. Note that these labels carry no connotation as to the relative
 * quality of the respective annotation sources.
 */
void agn_locus_stream_label_pairwise(AgnLocusStream *stream,
                                     const char *refrfile,const char *predfile);

/**
 * @function Calculate iLoci from a node stream which may or may not include
 * data from multiple sources. Extend each iLocus boundary as far as possible
 * without overlapping a gene from another iLocus, or by `delta` nucleotides,
 * whichever is shorter.
 */
GtNodeStream *agn_locus_stream_new(GtNodeStream *in_stream, GtUword delta);

/**
 * @function Terminal iLoci or 'end loci' are empty iLoci at either end of a
 * sequence. To exclude terminal iLoci from the output, set `endmode` < 0. To
 * output only terminal iLoci, set `endmode` > 0. By default (`endmode == 0`),
 * terminal iLoci are reported along with all other iLoci.
 */
void agn_locus_stream_set_endmode(AgnLocusStream *stream, int endmode);

/**
 * @function Loci created by this stream are assigned an ID with a serial
 * number. The default format is 'locus%lu' (that is, locus1, locus2, etc).
 * Use this function to override the default ID format.
 */
void agn_locus_stream_set_idformat(AgnLocusStream *stream, const char *format);

/**
 * @function By default, the locus stream will produce loci containing features
 * and loci containing no features. This function disables reporting of the
 * latter.
 */
void agn_locus_stream_skip_empty_loci(AgnLocusStream *stream);

/**
 * @function Set the source value to be used for all iLoci created by this
 * stream. Default value is 'AEGeAn::AgnLocusStream'.
 */
void agn_locus_stream_set_source(AgnLocusStream *stream, const char *source);

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_locus_stream_unit_test(AgnUnitTest *test);

#endif
