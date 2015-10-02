/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include "ga_commands.h"
#include "genometools.h"
#include "aegean.h"

typedef struct
{
  const char *source;
  const char *repo;
  const char **filenames;
  int numfiles;
} GaCleanOptions;

static void ga_clean_print_usage(FILE *outstream)
{
  fputs("\n"
"[GeneAnnoLogy::clean] remove existing (committed) annotations from a\n"
"                      repository and start over from scratch with a new\n"
"                      annotation\n\n"
"Usage: geneannology clean [options] repo annot.gff3\n"
"  Options:\n"
"    -h|--help              print this help message and exit\n"
"    -s|--source: STRING    replace the 'source' label (GFF3 2nd column) of\n"
"                           the input with the given value\n"
"    -v|--version           print version number and exit\n\n", outstream);
}

static void ga_clean_parse_options(int argc, char * const *argv,
                                   GaCleanOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hs:v";
  const struct option clean_options[] =
  {
    { "help",    no_argument,       NULL, 'h' },
    { "source",  required_argument, NULL, 's' },
    { "version", no_argument,       NULL, 'v' },
  };

  options->source = NULL;
  options->repo = NULL;
  options->filenames = NULL;
  options->numfiles = 0;

  for(opt  = getopt_long(argc, argv, optstr, clean_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv, optstr, clean_options, &optindex))
  {
    if(opt == 'h')
    {
      ga_clean_print_usage(stdout);
      exit(0);
    }
    else if(opt == 's')
    {
      options->source = optarg;
    }
    else if(opt == 'v')
    {
      agn_print_version("GeneAnnoLogy::clean", stdout);
      exit(0);
    }
    else
    {
      ga_clean_print_usage(stderr);
      fprintf(stderr, "error: unknown option '%c'", opt);
      exit(1);
    }
  }

  int numargs = argc - optind;
  if(numargs == 0)
    return;

  options->repo = argv[optind];
  options->numfiles = numargs - 1;
  if(numargs == 1)
    return;

  options->filenames = (const char **)(argv + optind + 1);
}

int ga_clean(int argc, char * const *argv)
{
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams = gt_queue_new();
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtError *error = gt_error_new();

  agn_assert(argc > 0);
  GaCleanOptions options;
  ga_clean_parse_options(argc, argv, &options);
  if(options.repo == NULL)
  {
    ga_clean_print_usage(stderr);
    return 0;
  }
  else if(options.numfiles == 0)
  {
    fprintf(stderr, "[GeneAnnoLogy::clean] warning: repo=%s, reading input from"
            " stdin\n", options.repo);
    current_stream = gt_gff3_in_stream_new_unsorted(0, NULL);
  }
  else
  {
    current_stream = gt_gff3_in_stream_new_unsorted(options.numfiles,
                                                    options.filenames);
  }
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, 0);
  agn_locus_stream_set_source((AgnLocusStream *)current_stream, "GeneAnnoLogy");
  agn_locus_stream_skip_iiLoci((AgnLocusStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(options.source)
  {
    GtStr *src = gt_str_new_cstr(options.source);
    GtNodeVisitor *nv = gt_set_source_visitor_new(src);
    current_stream = gt_visitor_stream_new(last_stream, nv);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
    gt_str_delete(src);
  }

  current_stream = agn_repo_stream_open_clean(last_stream, options.repo, error);
  if(current_stream == NULL)
  {
    fprintf(stderr, "[GeneAnnoLogy::clean] error setting up repo: %s\n",
            gt_error_get(error));
    return -1;
  }
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[GeneAnnoLogy::clean] error processing node stream: %s\n",
            gt_error_get(error));
  }

  gt_error_delete(error);
  gt_logger_delete(logger);
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_queue_delete(streams);
  return 0;
}
