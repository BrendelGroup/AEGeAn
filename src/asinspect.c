/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include "genometools.h"
#include "AgnASInspectVisitor.h"

int main(int argc, const char **argv)
{
  gt_lib_init();
  GtNodeStream *gtf = gt_gtf_in_stream_new(NULL);
  GtNodeStream *as  = agn_as_inspect_stream_new(gtf);
  GtError *error = gt_error_new();
  int result = gt_node_stream_pull(as, error);
  if(result == -1)
    fprintf(stderr, "error processing node stream: %s\n", gt_error_get(error));

  gt_error_delete(error);
  gt_node_stream_delete(as);
  gt_node_stream_delete(gtf);
  gt_lib_clean();
  return result;
}
