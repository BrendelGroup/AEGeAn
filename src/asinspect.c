/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include "genometools.h"
#include "AgnASInspectVisitor.h"

typedef struct
{
  bool gtf;
} ASInspectOptions;

static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\nasinspect: summarize alternative splicing annotated in the given "
"GFF3/GTF file\n"
"Usage: asinspect [options] annot.gff3\n"
"  Options:\n"
"    -g|--gtf:     specifies that input is in GTF format; defaut is GFF3\n"
"    -h|--help:    print this help message and exit\n\n");
}

static void parse_options(int argc, char **argv, ASInspectOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "gh";
  const struct option pmrna_options[] =
  {
    { "gtf",  no_argument, NULL, 'g' },
    { "help", no_argument, NULL, 'h' },
  };
  for(opt  = getopt_long(argc, argv + 0, optstr, pmrna_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv + 0, optstr, pmrna_options, &optindex))
  {
    switch(opt)
    {
      case 'g':
        options->gtf = true;
        break;
      case 'h':
        print_usage(stdout);
        exit(0);
        break;
    }
  }
}

int main(int argc, char **argv)
{
  GtNodeStream *annot, *as;
  GtError *error;
  ASInspectOptions options = { false };
  parse_options(argc, argv, &options);
  int numargs = argc - optind;

  gt_lib_init();
  if(options.gtf)
  {
    const char *infile = NULL;
    if(numargs > 0)
    {
      infile = argv[optind];
      if(numargs > 1)
      {
        fprintf(stderr, "[asinspect] warning: only one input file allowed in "
                "GTF mode; ignoring %d file(s)\n", numargs - 1);
      }
    }
    annot = gt_gtf_in_stream_new(infile);
  }
  else
  {
    const char **infiles = NULL;
    if(numargs > 0)
      infiles = (const char **)argv + optind;
    annot = gt_gff3_in_stream_new_unsorted(numargs, infiles);
    gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)annot);
  }

  as  = agn_as_inspect_stream_new(annot);
  error = gt_error_new();
  int result = gt_node_stream_pull(as, error);
  if(result == -1)
    fprintf(stderr, "error processing node stream: %s\n", gt_error_get(error));

  gt_error_delete(error);
  gt_node_stream_delete(as);
  gt_node_stream_delete(annot);
  gt_lib_clean();
  return result;
}
