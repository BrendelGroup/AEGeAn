#ifndef AEGEAN_COMPARE_REPORT_HTML
#define AEGEAN_COMPARE_REPORT_HTML

#include "core/logger_api.h"
#include "extended/node_visitor_api.h"

/**
 * @class AgnCompareReportHTML
 *
 * The ``AgnCompareReportHTML`` class is an extension of the
 * ``AgnCompareReport`` class. This node visitor relies on its parent class to
 * process a stream of ``AgnLocus`` objects (containing two alternative sources
 * of annotation to be compared) and then produces textual reports of the
 * comparison statistics.
 */
typedef struct AgnCompareReportHTML AgnCompareReportHTML;


/**
 * @function After the node stream has been processed, call this function to
 * write a summary of all locus comparisons to the output directory.
 */
void agn_compare_report_html_create_summary(AgnCompareReportHTML *rpt,
                                            int argc, char **argv,
                                            const char *refrlabel,
                                            const char *predlabel,
                                            const char *start_time);

/**
 * @function Class constructor. Creates a node visitor used to process a stream
 * of ``AgnLocus`` objects containing two sources of annotation to be compared.
 * Reports will be written in ``outdir`` and status messages will be written
 * to the logger.
 */
GtNodeVisitor *agn_compare_report_html_new(const char *outdir,GtLogger *logger);

#endif
