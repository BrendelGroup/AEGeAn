/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include <string.h>
#include "AgnCliquePair.h"
#include "AgnFilterStream.h"
#include "AgnGeneStream.h"
#include "AgnInferCDSVisitor.h"
#include "AgnInferExonsVisitor.h"
#include "AgnInferParentStream.h"
#include "AgnIntervalLocusStream.h"
#include "AgnLocus.h"
#include "AgnLocusStream.h"
#include "AgnMrnaRepVisitor.h"
#include "AgnPseudogeneFixVisitor.h"
#include "AgnRemoveChildrenVisitor.h"
#include "AgnTranscriptClique.h"

int main(int argc, char **argv)
{
  puts("AEGeAn Unit Tests");
  gt_lib_init();

  GtQueue *tests = gt_queue_new();
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnPseudogeneFixVisitor",
                                        agn_pseudogene_fix_visitor_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnInferParentStream",
                                        agn_infer_parent_stream_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnMrnaRepVisitor",
                                        agn_mrna_rep_visitor_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnRemoveChildrenVisitor",
                                        agn_remove_children_visitor_unit_test));
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
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnGeneStream",
                                        agn_gene_stream_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnLocusStream",
                                        agn_locus_stream_unit_test));
  gt_queue_add(tests, agn_unit_test_new("AEGeAn::AgnIntervalLocusStream",
                                        agn_interval_locus_stream_unit_test));

  unsigned passes   = 0;
  unsigned failures = 0;
  while(gt_queue_size(tests) > 0)
  {
    AgnUnitTest *test = gt_queue_get(tests);
    agn_unit_test_run(test);
    agn_unit_test_print(test, stdout);
    if(agn_unit_test_success(test))
      passes++;
    else
      failures++;
    agn_unit_test_delete(test);
  }

  bool returnval = 0;
  if(failures == 0)
  {
    printf("\n===== Unit tests passed for all %u classes! =====\n\n", passes);
  }
  else
  {
    printf("\n===== Unit tests failed for %u/%u classes! ===== \n\n", failures,
           passes+failures);
    returnval = 1;
  }

  gt_queue_delete(tests);
  gt_lib_clean();
  return returnval;
}
