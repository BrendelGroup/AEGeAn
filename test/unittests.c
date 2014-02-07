#include <string.h>
#include "AgnCliquePair.h"
#include "AgnFilterStream.h"
#include "AgnInferCDSVisitor.h"
#include "AgnInferExonsVisitor.h"
#include "AgnInferParentStream.h"
#include "AgnIntervalLocusStream.h"
#include "AgnLocus.h"
#include "AgnLocusStream.h"
#include "AgnPseudogeneFixVisitor.h"
#include "AgnTranscriptClique.h"
#include "AgnTranscriptStream.h"

int main(int argc, char **argv)
{
  puts("AEGeAn Unit Tests");
  gt_lib_init();

  GtQueue *tests = gt_queue_new();
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnPseudogeneFixVisitor",
                                        agn_pseudogene_fix_visitor_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnInferParentStream",
                                        agn_infer_parent_stream_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnTranscriptClique",
                                        agn_transcript_clique_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnCliquePair",
                                        agn_clique_pair_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnLocus",
                                        agn_locus_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnFilterStream",
                                        agn_filter_stream_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnInferCDSVisitor",
                                        agn_infer_cds_visitor_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnInferExonsVisitor",
                                        agn_infer_exons_visitor_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnTranscriptStream",
                                        agn_transcript_stream_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnLocusStream",
                                        agn_locus_stream_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnIntervalLocusStream",
                                        agn_interval_locus_stream_unit_test));

  while(gt_queue_size(tests) > 0)
  {
    AgnUnitTest *test = gt_queue_get(tests);
    agn_unit_test_run(test);
    agn_unit_test_print(test, stdout);
    agn_unit_test_delete(test);
  }

  gt_queue_delete(tests);
  gt_lib_clean();
  return 0;
}
