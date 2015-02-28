/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include <getopt.h>
#include "genometools.h"
#include "AgnMrnaRepVisitor.h"
#include "AgnPseudogeneFixVisitor.h"

typedef struct
{
  bool infer_introns;
  FILE *outstream;
  bool fix_pseudogenes;
} PmrnaOptions;

static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\npmrna: filter out all but the primary isoform from each gene of the input\n"
"Usage: pmrna [options] < annot.gff3 > new.gff3\n"
"  Options:\n"
"    -h|--help           print this help message and exit\n"
"    -i|--introns        flag indicating that introns are declared explicitly\n"
"                        and do not need to be inferred from exon features;\n"
"                        default is to infer introns\n"
"    -p|--pseudogenes    disable pseudogene detection and correction\n\n");
}

static void parse_options(int argc, char **argv, PmrnaOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hip";
  const struct option pmrna_options[] =
  {
    { "help",        no_argument,       NULL, 'h' },
    { "introns",     no_argument,       NULL, 'i' },
    { "pseudogenes", no_argument,       NULL, 'o' },
  };
  for(opt  = getopt_long(argc, argv + 0, optstr, pmrna_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv + 0, optstr, pmrna_options, &optindex))
  {
    switch(opt)
    {
      case 'h':
        print_usage(stdout);
        exit(0);
        break;
      case 'i':
        options->infer_introns = false;
        break;
      case 'p':
        options->fix_pseudogenes = false;
        break;
    }
  }
}

int main(int argc, char **argv)
{
  GtError *error;
  GtNodeStream *stream, *last_stream;
  GtQueue *streams;
  PmrnaOptions options = { true, NULL, true };
  parse_options(argc, argv, &options);

  //----------
  // Set up the processing stream
  //----------
  gt_lib_init();
  streams = gt_queue_new();

  stream = gt_gff3_in_stream_new_unsorted(0, NULL);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  if(options.fix_pseudogenes)
  {
    stream = agn_pseudogene_fix_stream_new(last_stream);
    gt_queue_add(streams, stream);
    last_stream = stream;
  }

  if(options.infer_introns)
  {
    stream = gt_add_introns_stream_new(last_stream);
    gt_queue_add(streams, stream);
    last_stream = stream;
  }

  stream = agn_mrna_rep_stream_new(last_stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  stream = gt_gff3_out_stream_new(last_stream, NULL);
  gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream *)stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  //----------
  // Execute the processing stream
  //----------
  error = gt_error_new();
  int had_err = gt_node_stream_pull(last_stream, error);
  if(had_err)
    fprintf(stderr, "Error processing node stream: %s\n", gt_error_get(error));
  gt_error_delete(error);

  //----------
  // Free memory
  //----------
  while(gt_queue_size(streams) > 0)
  {
    stream = gt_queue_get(streams);
    gt_node_stream_delete(stream);
  }
  gt_queue_delete(streams);
  gt_lib_clean();
  return had_err;
}
