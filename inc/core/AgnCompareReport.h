#ifndef AEGEAN_COMPARE_REPORT
#define AEGEAN_COMPARE_REPORT

#include "core/array_api.h"
#include "core/logger_api.h"
#include "extended/node_stream_api.h"
#include "AgnComparison.h"
#include "AgnLocus.h"

/**
 * @class AgnCompareReport
 *
 * The ``AgnCompareReport`` class implements the ``GtNodeVisitor`` class. It is
 * used to process a stream of ``AgnLocus`` objects (containing two alternative
 * sources of annotation to be compared) and then produce reports of the
 * comparison statistics.
 */
typedef struct AgnCompareReport AgnCompareReport;

/**
 * @functype Signature of callback functions used to create locus-level reports,
 * if any.
 */
typedef void (*AgnCompareReportLocusFunc)(AgnLocus *locus, void *data);

/**
 * @functype Signature of callback functions used to create sequence-level
 * summary reports, if any.
 */
typedef void (*AgnCompareReportSequenceFunc)(const AgnComparisonData *cd,
                                             const char *seqid, void *data);

/**
 * @function Return a pointer to the statistics collected by the node visitor.
 */
AgnComparisonData *agn_compare_report_data(AgnCompareReport *rpt);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_compare_report_new(GtArray *locusfilters, GtLogger *logger);

/**
 * @function Return a string array containing IDs of all sequences for which the
 * node stream contained annotations.
 */
GtStrArray *agn_compare_report_seqids(AgnCompareReport *rpt);

/**
 * @function Set the callback function to be used for generating locus reports.
 * The ``data`` pointer is optional but will be available to the callback
 * function. Set ``func`` to NULL to disable locus report generation (default).
 */
void agn_compare_report_set_locus_callback(AgnCompareReport *rpt,
                                           AgnCompareReportLocusFunc func,
                                           void *data);

/**
 * @function Set the callback function to be used for generating sequence
 * reports. The ``data`` pointer is optional but will be available to the
 * callback function. Set ``func`` to NULL to disable sequence report generation
 * (default).
 */
void agn_compare_report_set_sequence_callback(AgnCompareReport *rpt,
                                              AgnCompareReportSequenceFunc func,
                                              void *data);

#endif