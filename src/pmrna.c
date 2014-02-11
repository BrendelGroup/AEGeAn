#include <getopt.h>
#include "genometools.h"
#include "AgnMrnaRepVisitor.h"
#include "AgnPseudogeneFixVisitor.h"

int main()
{
  gt_lib_init();

  GtNodeStream *in = gt_gff3_in_stream_new_unsorted(0, NULL);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)in);
  GtNodeStream *pfs = agn_pseudogene_fix_stream_new(in);
  GtNodeStream *ais = gt_add_introns_stream_new(pfs);
  GtNodeStream *rep_stream = agn_mrna_rep_stream_new(ais);
  GtFile *fh = gt_file_new_from_fileptr(stdout);
  GtNodeStream *out = gt_gff3_out_stream_new(rep_stream, fh);
  gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream *)out);

  GtError *error = gt_error_new();
  int result = gt_node_stream_pull(out, error);
  if(result == -1)
  {
    fprintf(stderr, "Error processing node stream: %s\n", gt_error_get(error));
  }
  gt_error_delete(error);

  gt_node_stream_delete(out);
  gt_node_stream_delete(rep_stream);
  gt_node_stream_delete(ais);
  gt_node_stream_delete(pfs);
  gt_node_stream_delete(in);
  gt_file_delete_without_handle(fh);
  gt_lib_clean();
  return 0;
}
