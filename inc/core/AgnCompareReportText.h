/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_COMPARE_REPORT_TEXT
#define AEGEAN_COMPARE_REPORT_TEXT

#include "core/logger_api.h"
#include "extended/node_visitor_api.h"

/**
 * @class AgnCompareReportText
 *
 * The ``AgnCompareReportText`` class is an extension of the
 * ``AgnCompareReport`` class. This node visitor relies on its parent class to
 * process a stream of ``AgnLocus`` objects (containing two alternative sources
 * of annotation to be compared) and then produces textual reports of the
 * comparison statistics.
 */
typedef struct AgnCompareReportText AgnCompareReportText;


/**
 * @function After the node stream has been processed, call this function to
 * write a summary of all locus comparisons to ``outstream``.
 */
void agn_compare_report_text_create_summary(AgnCompareReportText *rpt,
                                            FILE *outstream);

/**
 * @function Class constructor. Creates a node visitor used to process a stream
 * of ``AgnLocus`` objects containing two sources of annotation to be compared.
 * Reports will be written to ``outstream`` and status messages will be written
 * to the logger.
 */
GtNodeVisitor *agn_compare_report_text_new(FILE *outstream, GtLogger *logger);

#endif
