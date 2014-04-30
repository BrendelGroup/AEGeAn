#ifndef AEGEAN_COMPARE_REPORT_HTML
#define AEGEAN_COMPARE_REPORT_HTML

#include "AgnCompareReport.h"

/**
 * @class AgnCompareReportHTML
 *
 * The ``AgnCompareReportHTML`` class is an extension of the
 * ``AgnCompareReport`` class. This node visitor relies on its parent class to
 * process a stream of ``AgnLocus`` objects (containing two alternative sources
 * of annotation to be compared) and then produces hyperlinked HTML reports of
 * the comparison statistics.
 */
typedef AgnCompareReport AgnCompareReportHTML;


/**
 * @function After the node stream has been processed, call this function to
 * write a summary of all locus comparisons to ``outstream``.
 */
void agn_compare_report_html_create_summary(AgnCompareReportHTML *rpt,
                                            const char *outdir);

/**
 * @function Class constructor. Creates a node visitor used to process a stream
 * of ``AgnLocus`` objects containing two sources of annotation to be compared.
 * Reports will be written to ``outdir`` and status messages will be written
 * to the logger. An assertion 
 */
GtNodeVisitor *agn_compare_report_text_new(const char *outdir,GtLogger *logger);

#endif
