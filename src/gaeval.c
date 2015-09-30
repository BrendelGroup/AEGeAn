/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include <math.h>
#include "genometools.h"
#include "AgnGaevalVisitor.h"
#include "AgnInferCDSVisitor.h"
#include "AgnInferExonsVisitor.h"
#include "AgnUtils.h"

typedef struct
{
  const char *alignfile;
  const char **genefiles;
  int numgenefiles;
  AgnGaevalParams params;
} GaevalOptions;

static void default_params(AgnGaevalParams *params)
{
  params->alpha = 0.6;
  params->beta = 0.3;
  params->gamma = 0.05;
  params->epsilon = 0.05;
  params->exp_cds_len = 400;
  params->exp_5putr_len = 200;
  params->exp_3putr_len = 100;
}

static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\ngaeval: calculate coverage and intergrity scores for gene models based on "
"transcript alignments\n"
"Usage: gaeval [options] alignments.gff3 genes.gff3 [moregenes.gff3 ...]\n"
"  Basic options:\n"
"    -h|--help               print this help message and exit\n"
"    -v|--version            print version number and exit\n\n"
"  Weights for calculating integrity score (must add up to 1.0):\n"
"    -a|--alpha: DOUBLE      introns confirmed, or %% expected CDS length for\n"
"                            single-exon genes; default is 0.6\n"
"    -b|--beta: DOUBLE       exon coverage; default is 0.3\n"
"    -g|--gamma: DOUBLE      %% expected 5' UTR length; default is 0.05\n"
"    -e|--epsilon: DOUBLE    %% expected 3' UTR length; default is 0.05\n\n"
"  Expected feature lengths for calculating integrity score:\n"
"    -c|--exp-cds: INT       expected CDS length (in bp); default is 400\n"
"    -5|--exp-5putr: INT     expected 5' UTR length; default is 200\n"
"    -3|--exp-3putr: INT     expected 3' UTR length; default is 100\n\n");
}

static void parse_options(int argc, char **argv, GaevalOptions *options)
{
  default_params(&options->params);
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hva:b:g:e:c:5:3:";
  const struct option gaeval_options[] =
  {
    { "help",      no_argument,       NULL, 'h' },
    { "version",   no_argument,       NULL, 'v' },
    { "alpha",     required_argument, NULL, 'a' },
    { "beta",      required_argument, NULL, 'b' },
    { "gamma",     required_argument, NULL, 'g' },
    { "epsilon",   required_argument, NULL, 'e' },
    { "exp-cds",   required_argument, NULL, 'c' },
    { "exp-5putr", required_argument, NULL, '5' },
    { "exp-3putr", required_argument, NULL, '3' },
    { NULL,        no_argument,       NULL,  0  },
  };
  for(opt  = getopt_long(argc, argv + 0, optstr, gaeval_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv + 0, optstr, gaeval_options, &optindex))
  {
    if(opt == '3')
      options->params.exp_3putr_len = atoi(optarg);
    else if(opt == '5')
      options->params.exp_5putr_len = atoi(optarg);
    else if(opt == 'a')
      options->params.alpha = atof(optarg);
    else if(opt == 'b')
      options->params.beta = atof(optarg);
    else if(opt == 'c')
      options->params.exp_cds_len = atoi(optarg);
    else if(opt == 'e')
      options->params.epsilon = atof(optarg);
    else if(opt == 'h')
    {
      print_usage(stdout);
      exit(0);
    }
    else if(opt == 'g')
      options->params.gamma = atof(optarg);
    else if(opt == 'v')
    {
      agn_print_version("GAEVAL", stdout);
      exit(0);
    }
  }
  int numargs = argc - optind;
  if(numargs < 2)
  {
    print_usage(stderr);
    fprintf(stderr, "error: must provide at least 2 input files, %d provided\n",
            numargs);
    exit(1);
  }

  double weight_total = options->params.alpha + options->params.beta +
                        options->params.gamma + options->params.epsilon;
  if(fabs(weight_total - 1.0) > 0.00001)
  {
    print_usage(stderr);
    fprintf(stderr, "error: integrity score weights must add up to 1.0; "
            "specified weights total %.2lf\n", weight_total);
    exit(1);
  }

  options->alignfile = argv[optind + 0];
  options->genefiles = (const char **)argv + optind + 1;
  options->numgenefiles = numargs - 1;
}

int main(int argc, char **argv)
{
  GtError *error;
  GtNodeStream *stream, *last_stream, *align_stream;
  GtQueue *streams;
  GaevalOptions options;
  parse_options(argc, argv, &options);

  //----------
  // Set up the processing stream
  //----------
  gt_lib_init();
  streams = gt_queue_new();

  stream = gt_gff3_in_stream_new_unsorted(1, &options.alignfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)stream);
  gt_queue_add(streams, stream);
  align_stream = stream;

  stream = gt_gff3_in_stream_new_unsorted(options.numgenefiles,
                                          options.genefiles);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  GtStr *source = gt_str_new_cstr("AEGeAn::GAEVAL");
  GtLogger *logger = gt_logger_new(true, "", stderr);
  stream = agn_infer_cds_stream_new(last_stream, source, logger);
  gt_queue_add(streams, stream);
  last_stream = stream;

  stream = agn_infer_exons_stream_new(last_stream, source, logger);
  gt_queue_add(streams, stream);
  last_stream = stream;
  gt_str_delete(source);

  stream = agn_gaeval_stream_new(last_stream, align_stream, options.params);
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
  gt_logger_delete(logger);
  gt_lib_clean();
  return had_err;
}
