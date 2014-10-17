/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include "genometools.h"
#include "aegean.h"

typedef struct
{
  bool gtf;
  FILE *outstream;
  const char *refrfile;
} ASInspectOptions;

static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\nasinspect: summarize alternative splicing annotated in the given "
"GFF3/GTF file\n"
"Usage: asinspect [options] annot.gff3\n"
"  Options:\n"
"    -g|--gtf:          specifies that input is in GTF format; defaut is GFF3\n"
"    -h|--help:         print this help message and exit\n"
"    -o|--out: FILE     file to which output will be written; default is\n"
"                       terminal (standard output)\n"
"    -r|--refr: FILE    do not search for alternative splicing events within\n"
"                       the input file, but by comparison against the\n"
"                       provided reference\n\n");
}

static void parse_options(int argc, char **argv, ASInspectOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "gho:r:";
  const struct option pmrna_options[] =
  {
    { "gtf",  no_argument,       NULL, 'g' },
    { "help", no_argument,       NULL, 'h' },
    { "out",  required_argument, NULL, 'o' },
    { "refr", required_argument, NULL, 'r' },
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
      case 'o':
        options->outstream = fopen(optarg, "w");
        if(options->outstream == NULL)
        {
          fprintf(stderr, "error: could not open output file '%s'\n", optarg);
          exit(1);
        }
        break;
      case 'r':
        options->refrfile = optarg;
        break;
    }
  }
}

int main(int argc, char **argv)
{
  GtQueue *streams;
  GtNodeStream *current_stream, *last_stream;
  GtError *error;
  GtLogger *logger;
  ASInspectOptions options = { false, stdout, NULL };
  parse_options(argc, argv, &options);
  int numargs = argc - optind;

  gt_lib_init();
  streams = gt_queue_new();
  logger = gt_logger_new(true, "", stderr);

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
    current_stream = gt_gtf_in_stream_new(infile);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }
  else
  {
    const char **infiles = NULL;
    if(numargs > 0)
      infiles = (const char **)argv + optind;
    current_stream = gt_gff3_in_stream_new_unsorted(numargs, infiles);
    gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  if(options.refrfile != NULL)
  {
    GtNodeStream *rstream = gt_gff3_in_stream_new_unsorted(1,&options.refrfile);
    gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)rstream);
    current_stream = agn_locus_stream_new_pairwise(rstream,last_stream,logger);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  current_stream = agn_as_inspect_ce_stream_new(last_stream, options.outstream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_as_inspect_ir_stream_new(last_stream,
                                                options.outstream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  error = gt_error_new();
  fprintf(options.outstream, "Mechanism\tGeneID\tRefrMrna\tTestMrna\tSeqID\t"
          "Description\n");
  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
    fprintf(stderr, "error processing node stream: %s\n", gt_error_get(error));

  gt_error_delete(error);
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_queue_delete(streams);
  gt_logger_delete(logger);
  gt_lib_clean();
  return result;
}
