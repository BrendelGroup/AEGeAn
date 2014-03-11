#include <getopt.h>
#include "genometools.h"
#include "aegean.h"

// Simple data structure for program options
typedef struct
{
  FILE *idfile;
  FILE *outfile;
  bool typeoverride;
  GtHashmap *typestoextract;
  unsigned width;
} XtractoreOptions;

static void set_option_defaults(XtractoreOptions *options)
{
  options->idfile = NULL;
  options->outfile = stdout;
  options->typeoverride = false;
  options->typestoextract = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
  char *defaulttype = gt_cstr_dup("gene");
  gt_hashmap_add(options->typestoextract, defaulttype, defaulttype);
  options->width = 100;
}

static void free_option_memory(XtractoreOptions *options)
{
  if(options->idfile != NULL)
    fclose(options->idfile);
  fclose(options->outfile);
  gt_hashmap_delete(options->typestoextract);
}

static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\nxtractore: extract sequences corresponding to annotated features from the\n"
"           given sequence file\n\n"
"Usage: xtractore [options] features.gff3 sequences.fasta"
"  Options:\n"
"    -h|--help             print this help message and exit\n"
"    -i|--idfile: FILE     file containing a list of feature IDs (1 per line\n"
"                          with no spaces); if provided, only features with\n"
"                          IDs in this file will be extracted\n"
"    -o|--outfile: FILE    file to which output sequences will be written;\n"
"                          default is terminal (stdout)\n"
"    -t|--type: FILE       feature type to extract; can be used multiple\n"
"                          times to extract features of multiple types\n"
"    -w|--width: INT       width of each line of sequence in the Fasta\n"
"                          output; default is 100\n\n");
}

static void
parse_options(int argc, char **argv, XtractoreOptions *options, GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hi:o:t:w:";
  char *type;
  const struct option xtractore_options[] =
  {
    { "help",     no_argument,      NULL, 'h' },
    { "idfile",  required_argument, NULL, 'i' },
    { "outfile", required_argument, NULL, 'o' },
    { "type",    required_argument, NULL, 't' },
    { "width",   required_argument, NULL, 'w' },
  };
  for( opt = getopt_long(argc, argv + 0, optstr, xtractore_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, xtractore_options, &optindex))
  {
    switch(opt)
    {
      case 'h':
        print_usage(stdout);
        exit(0);
        break;
      case 'i':
        options->idfile = fopen(optarg, "r");
        if(options->idfile == NULL)
          gt_error_set(error, "could not open ID file '%s'", optarg);
        break;
      case 'o':
        options->outfile = fopen(optarg, "w");
        if(options->outfile == NULL)
          gt_error_set(error, "could not open output file '%s'", optarg);
        break;
      case 't':
        if(options->typeoverride == false)
        {
          gt_hashmap_delete(options->typestoextract);
          gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
          options->typeoverride = true;
        }
        type = gt_cstr_dup(optarg);
        gt_hashmap_add(options->typestoextract, type, type);
        break;
      case 'w':
        if(sscanf(optarg, "%u", &options->width) == EOF)
        {
          gt_error_set(error, "could not convert width '%s' to an integer",
                       optarg);
        }
        break;
    }
  }
}

int main(int argc, char **argv)
{
  const char *featfile, *seqfile;
  GtFeatureIndex *features;
  GtError *error;
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams;
  GtSeqIterator *seqiter;
  GtStrArray *seqfastas;
  gt_lib_init();

  XtractoreOptions options;
  set_option_defaults(&options);
  error = gt_error_new();
  parse_options(argc, argv, &options, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[xtractore] error: %s", gt_error_get(error));
    return 1;
  }
  int numfiles = argc - optind;
  if(numfiles < 2)
  {
    fprintf(stderr, "[xtractore] error: must provide a feature annotation file "
            "(GFF3 format) and a sequence file (Fasta format)\n");
    print_usage(stderr);
    return 1;
  }
  featfile = argv[optind + 0];
  seqfile  = argv[optind + 1];

  streams = gt_queue_new();

  current_stream = gt_gff3_in_stream_new_unsorted(1, &featfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_filter_stream_new(last_stream, options.typestoextract);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  features = gt_feature_index_memory_new();
  current_stream = gt_feature_out_stream_new(last_stream, features);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[xtractore] error processing GFF3: %s",
            gt_error_get(error));
  }

  seqfastas = gt_str_array_new();
  gt_str_array_add_cstr(seqfastas, seqfile);
  seqiter = gt_seq_iterator_sequence_buffer_new(seqfastas, error);

  gt_seq_iterator_delete(seqiter);
  gt_str_array_delete(seqfastas);
  while(gt_queue_size(streams) > 0)
  {
    GtNodeStream *stream = gt_queue_get(streams);
    gt_node_stream_delete(stream);
  }
  gt_queue_delete(streams);
  gt_feature_index_delete(features);
  gt_error_delete(error);
  free_option_memory(&options);
  gt_lib_clean();
  return 0;
}
