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
  GtFile *outstream;
  GtStr *source;
  bool infer;
} CanonGFF3Options;

static void print_usage(FILE *outstream)
{
  fputs("\nUsage: canon-gff3 [options] gff3file1 [gff3file2 ...]\n"
"  Options:\n"
"     -h|--help               print this help message and exit\n"
"     -i|--infer              for transcript features lacking an explicitly\n"
"                             declared gene feature as a parent, create this\n"
"                             feature on-they-fly\n"
"     -o|--outfile: STRING    name of file to which GFF3 data will be\n"
"                             written; default is terminal (stdout)\n"
"     -s|--source: STRING     reset the source of each feature to the given\n"
"                             value\n"
"     -v|--version            print version number and exit\n\n",
        outstream);
}

static void canon_gff3_parse_options(int argc, char * const *argv,
                                     CanonGFF3Options *options, GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hio:s:v";
  const struct option init_options[] =
  {
    { "help",    no_argument,       NULL, 'h' },
    { "infer",   no_argument,       NULL, 'i' },
    { "outfile", required_argument, NULL, 'o' },
    { "source",  required_argument, NULL, 's' },
    { "version", no_argument,       NULL, 'v' },
    { NULL,      no_argument,       NULL, 0 },
  };

  for(opt = getopt_long(argc, argv, optstr, init_options, &optindex);
      opt != -1;
      opt = getopt_long(argc, argv, optstr, init_options, &optindex))
  {
    switch(opt)
    {
      case 'h':
        print_usage(stdout);
        exit(0);
        break;

      case 'i':
        options->infer = true;
        break;

      case 'o':
        if(options->outstream != NULL)
          gt_file_delete(options->outstream);
        options->outstream = gt_file_new(optarg, "w", error);
        break;

      case 's':
        if(options->source != NULL)
          gt_str_delete(options->source);
        options->source = gt_str_new_cstr(optarg);
        break;

      case 'v':
        agn_print_version("CanonGFF3", stdout);
        exit(0);
        break;

      default:
        break;
    }
  }
}

// Main method
int main(int argc, char * const *argv)
{
  GtError *error;
  GtLogger *logger;
  GtQueue *streams;
  GtNodeStream *stream, *last_stream;
  CanonGFF3Options options = { NULL, NULL, false };

  gt_lib_init();
  error = gt_error_new();
  canon_gff3_parse_options(argc, argv + 0, &options, error);

  streams = gt_queue_new();
  logger = gt_logger_new(true, "", stderr);

  stream = gt_gff3_in_stream_new_unsorted(argc - optind, (const char **)
                                                          argv+optind);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  if(options.infer)
  {
    GtHashmap *type_parents = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                             gt_free_func);
    gt_hashmap_add(type_parents, gt_cstr_dup("mRNA"), gt_cstr_dup("gene"));
    gt_hashmap_add(type_parents, gt_cstr_dup("tRNA"), gt_cstr_dup("gene"));
    stream = agn_infer_parent_stream_new(last_stream,
                                                 type_parents);
    gt_hashmap_delete(type_parents);
    gt_queue_add(streams, stream);
    last_stream = stream;
  }

  stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, stream);
  last_stream = stream;

  if(options.source != NULL)
  {
    GtNodeVisitor *ssv = gt_set_source_visitor_new(options.source);
    stream = gt_visitor_stream_new(last_stream, ssv);
    gt_queue_add(streams, stream);
    last_stream = stream;
  }

  stream = gt_gff3_out_stream_new(last_stream, options.outstream);
  if(!options.infer)
    gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream *)stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  if(gt_node_stream_pull(last_stream, error) == -1)
  {
    fprintf(stderr, "[CanonGFF3] error processing node stream: %s",
            gt_error_get(error));
  }

  while(gt_queue_size(streams) > 0)
  {
    stream = gt_queue_get(streams);
    gt_node_stream_delete(stream);
  }
  gt_queue_delete(streams);
  if(options.source != NULL)
    gt_str_delete(options.source);
  if(options.outstream != NULL)
    gt_file_delete(options.outstream);
  gt_error_delete(error);
  gt_logger_delete(logger);
  gt_lib_clean();

  return 0;
}
