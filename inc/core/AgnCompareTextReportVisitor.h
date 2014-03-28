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
 * @function As part of a comparative analysis, all reference transcripts (or
 * transcript cliques) associated with a locus will be compared to all
 * prediction transcripts from the locus. If ``max_comparisons`` > 0, loci
 * containing more than ``max_comparisons`` transcript (clique) pairs will be
 * ignored. This cutoff is motivated by the observation that some individual
 * loci require an inordinate amount of runtime and memory for comparative
 * analysis, and due to the complexity of these loci the comparison statistics
 * don't provide much insight.
 */
void agn_compare_text_report_visitor_compare_max(AgnCompareTextReportVisitor *v,
                                                 GtUword max_comparisons);

/**
 * @function Include GFF3 corresponding to each clique pair in that clique
 * pair's comparison report.
 */
void
agn_compare_text_report_visitor_enable_gff3(AgnCompareTextReportVisitor *v);

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
