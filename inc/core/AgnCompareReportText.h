#ifndef AEGEAN_COMPARE_REPORT_TEXT
#define AEGEAN_COMPARE_REPORT_TEXT

#include "AgnCompareReport.h"

/**
 * @class AgnCompareReportText
 *
 * The ``AgnCompareReportText`` class is a special case of the ``GtNodeStream``
 * class that uses the ``AgnCompareReport`` node visitor to process a stream of
 * ``AgnLocus`` objects (containing two alternative sources of annotation to be
 * compared). The ``AgnCompareReportText`` class then produces textual reports
 * of the comparison statistics.
 */
typedef struct AgnCompareReportText AgnCompareReportText;


/**
 * @function After the node stream has been processed, call this function to
 * write a summary of all locus comparisons to ``outstream``.
 */
void agn_compare_report_text_create_summary(AgnCompareReportText *rpt,
                                            FILE *outstream);

/**
 * @function Class constructor. Input is a stream of ``AgnLocus`` objects
 * containing two sources of annotation to be compared. Reports will be written
 * to ``outstream`` and status messages will be written to the logger.
 */
GtNodeStream *agn_compare_report_text_new(GtNodeStream *locus_stream,
                                          FILE *outstream, GtLogger *logger);

#endif
