/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include "ga_commands.h"
#include "genometools.h"
#include "aegean.h"

static void ga_clean_print_usage(FILE *outstream)
{
  fputs("\n"
"[GeneAnnoLogy::clean] remove existing (committed) annotations from a\n"
"                      repository and start over from scratch with a new\n"
"                      annotation\n\n"
"Usage: geneannology clean [options] repo annot.gff3\n"
"  Options:\n"
"    -h|--help       print this help message and exit\n"
"    -v|--version    print version number and exit\n\n", outstream);
}

static void ga_clean_parse_options(int argc, char * const *argv)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hv";
  const struct option clean_options[] =
  {
    { "help",    no_argument, NULL, 'h' },
    { "version", no_argument, NULL, 'v' },
  };

  for( opt = getopt_long(argc, argv, optstr, clean_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv, optstr, clean_options, &optindex) )
  {
    if(opt == 'h')
    {
      ga_clean_print_usage(stdout);
      exit(0);
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
    }
  }
}

int ga_clean(int argc, char * const *argv)
{
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams = gt_queue_new();
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtError *error = gt_error_new();

  agn_assert(argc > 0);
  ga_clean_parse_options(argc, argv);
  if(argc == 1)
  {
    ga_clean_print_usage(stderr);
    return -1;
  }
  else if(argc == 2)
  {
    fprintf(stderr, "[GeneAnnoLogy] warning: repo=%s, reading input from "
            "stdin\n", argv[1]);
    current_stream = gt_gff3_in_stream_new_unsorted(0, NULL);
  }
  else
  {
    const char **filenames = (const char **)argv + 2;
    current_stream = gt_gff3_in_stream_new_unsorted(argc - 2, filenames);
  }
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, 0);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_repo_stream_open_clean(last_stream, argv[1], error);
  if(current_stream == NULL)
  {
    fprintf(stderr, "[GeneAnnoLogy] error setting up repo: %s\n",
            gt_error_get(error));
    return -1;
  }
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[GeneAnnoLogy] error processing node stream: %s\n",
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
