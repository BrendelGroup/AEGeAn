#ifndef AEGEAN_COMPARE_TEXT_REPORT_VISITOR
#define AEGEAN_COMPARE_TEXT_REPORT_VISITOR

#include "core/logger_api.h"
#include "extended/node_stream_api.h"

/**
 * @class AgnCompareTextReportVisitor
 *
 * Implements the GenomeTools ``GtNodeVisitor`` interface. This node visitor
 * will perform comparative analysis on each locus it encounters and write the
 * comparison results to the given text file.
 */
typedef struct AgnCompareTextReportVisitor AgnCompareTextReportVisitor;

/**
 * @function Constructor for a node stream based on this node visitor.
 */
GtNodeStream *agn_compare_text_report_stream_new(GtNodeStream *instream,
                                                 FILE *reports,
                                                 FILE *summary,
                                                 GtLogger *logger);

/**
 * @function Constructor for the node visitor.
 */
GtNodeVisitor *agn_compare_text_report_visitor_new(FILE *reports,
                                                   FILE *summary,
                                                   GtLogger *logger);

/**
 * @function If ``max_locus_transcripts`` > 0, loci containing more than
 * ``max_locus_transcripts`` reference transcripts or prediction transcripts
 * will not be analyzed or reported.
 */
void agn_compare_text_report_visitor_trans_max(AgnCompareTextReportVisitor *v,
                                               GtUword max_locus_transcripts);

#endif
