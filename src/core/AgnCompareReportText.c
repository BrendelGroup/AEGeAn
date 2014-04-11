#include "extended/visitor_stream_api.h"
#include "AgnCompareReportText.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnCompareReportText
{
  GtNodeStream *ns;
  AgnCompareReport *nv;
};

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Create a report for each locus.
 */
static void compare_report_text_locus_handler(AgnLocus *locus, void *data);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_compare_report_text_create_summary(AgnCompareReportText *rpt,
                                            FILE *outstream)
{
  GT_UNUSED AgnComparisonData *data = agn_compare_report_data(rpt->nv); 
  /*
       implement method here
   */
}

void agn_compare_report_text_delete(AgnCompareReportText *rpt)
{
  gt_assert(rpt);
  gt_node_stream_delete(rpt->ns);
  gt_free(rpt);
}

GtNodeStream *agn_compare_report_text_new(GtNodeStream *locus_stream,
                                          FILE *outstream, GtLogger *logger)
{
  AgnCompareReportText *rpt;
  GtArray *filters;

  filters = gt_array_new( sizeof(AgnLocusFilter) );
  rpt = gt_malloc( sizeof(AgnCompareReportText) );
  rpt->nv = (AgnCompareReport *)agn_compare_report_new(filters, logger);
  agn_compare_report_set_locus_callback(rpt->nv,
                                        compare_report_text_locus_handler,
                                        outstream);
  rpt->ns = gt_visitor_stream_new(locus_stream, (GtNodeVisitor *)rpt->nv);
  return (GtNodeStream *)rpt;
}

static void compare_report_text_locus_handler(AgnLocus *locus, void *data)
{
  GT_UNUSED FILE *outstream = data;
  /*
       implement method here
   */
}
