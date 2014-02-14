#include <getopt.h>
#include "genometools.h"
#include "aegean.h"

// Simple data structure for program options
typedef struct
{
  bool debug;
  FILE *genestream;
  bool intloci;
  unsigned long delta;
  GtFile *outstream;
  void (*filefreefunc)(GtFile *);
  bool retainids;
  bool skipends;
  FILE *transstream;
  bool pseudofix;
  bool verbose;
} LocusPocusOptions;

// Set default values for program
static void set_option_defaults(LocusPocusOptions *options)
{
  options->debug = false;
  options->genestream = NULL;
  options->intloci = false;
  options->delta = 500;
  options->outstream = gt_file_new_from_fileptr(stdout);
  options->filefreefunc = gt_file_delete_without_handle;
  options->retainids = false;
  options->skipends = false;
  options->transstream = NULL;
  options->pseudofix = false;
  options->verbose = false;
}

// Usage statement
static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\nLocusPocus: calculate locus coordinates for the given gene annotation\n"
"Usage: locuspocus [options] gff3file1 [gff3file2 gff3file3 ...]\n"
"  Options:\n"
"    -d|--debug             print detailed debugging messages to terminal\n"
"                           (standard error)\n"
"    -g|--genemap: FILE     print a mapping from each gene annotation to its\n"
"                           corresponding locus to the given file\n"
"    -h|--help              print this help message and exit\n"
"    -i|--intloci           parse and report interval loci rather than gene\n"
"                           loci\n"
"    -l|--delta: INT        when parsing interval loci, use the following\n"
"                           delta to extend gene loci and include potential\n"
"                           regulatory regions; default is 500\n"
"    -o|--outfile: FILE     name of file to which results will be written;\n"
"                           default is terminal (standard output)\n"
"    -r|--retainids         retain IDs from input data; selecting 'genemap'\n"
"                           and/or 'transmap' will activate this setting\n"
"    -s|--skipends          when enumerating interval loci, exclude gene-less\n"
"                           iloci at either end of the sequence\n"
"    -t|--transmap: FILE    print a mapping from each transcript annotation\n"
"                           to its corresponding locus to the given file\n"
"    -u|--pseudo            correct erroneously labeled pseudogenes\n"
"    -v|--verbose           include all locus subfeatures (genes, RNAs, etc)\n"
"                           in the GFF3 output; default includes only locus\n"
"                           features\n\n");
}

// Adjust program settings from command-line arguments/options
static void
parse_options(int argc, char **argv, LocusPocusOptions *options, GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "dg:hil:o:rst:uv";
  const struct option locuspocus_options[] =
  {
    { "debug",     no_argument,       NULL, 'd' },
    { "genemap",   required_argument, NULL, 'g' },
    { "help",      no_argument,       NULL, 'h' },
    { "intloci",   no_argument,       NULL, 'i' },
    { "delta",     required_argument, NULL, 'l' },
    { "outfile",   required_argument, NULL, 'o' },
    { "retainids", no_argument,       NULL, 'r' },
    { "skipends",  no_argument,       NULL, 's' },
    { "transmap",  required_argument, NULL, 't' },
    { "pseudo",    no_argument,       NULL, 'u' },
    { "verbose",   no_argument,       NULL, 'v' },
  };
  for( opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex))
  {
    switch(opt)
    {
      case 'd':
        options->debug = 1;
        break;
      case 'g':
        options->retainids = 1;
        options->genestream = fopen(optarg, "w");
        if(options->genestream == NULL)
          gt_error_set(error, "could not open genemap file '%s'", optarg);
        break;
      case 'h':
        print_usage(stdout);
        exit(0);
        break;
      case 'i':
        options->intloci = 1;
        break;
      case 'l':
        if(sscanf(optarg, "%lu", &options->delta) == EOF)
        {
          gt_error_set(error, "could not convert delta '%s' to an integer",
                       optarg);
        }
        break;
      case 'o':
        options->filefreefunc(options->outstream);
        options->outstream = gt_file_new(optarg, "w", error);
        options->filefreefunc = gt_file_delete;
        break;
      case 'r':
        options->retainids = 1;
        break;
      case 's':
        options->skipends = 1;
        break;
      case 't':
        options->retainids = 1;
        options->transstream = fopen(optarg, "w");
        if(options->transstream == NULL)
          gt_error_set(error, "could not open genemap file '%s'", optarg);
        break;
      case 'u':
        options->pseudofix = 1;
        break;
      case 'v':
        options->verbose = 1;
        break;
    }
  }
}

// Main program
int main(int argc, char **argv)
{
  GtError *error;
  GtLogger *logger;
  GtQueue *streams;
  GtNodeStream *current_stream, *last_stream;
  gt_lib_init();

  // Parse command-line options
  LocusPocusOptions options;
  set_option_defaults(&options);
  error = gt_error_new();
  parse_options(argc, argv, &options, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[LocusPocus] error: %s", gt_error_get(error));
    return 1;
  }
  int numfiles = argc - optind;
  if(numfiles < 1)
  {
    fprintf(stderr, "[LocusPocus] error: must provide at least one GFF3 file "
            "as input; use '-' if you want to provide data via standard "
            "input\n");
    print_usage(stderr);
    return 1;
  }

  logger = gt_logger_new(true, "", stderr);
  streams = gt_queue_new();


  //----- Set up the node processing stream -----//
  //---------------------------------------------//

  current_stream = gt_gff3_in_stream_new_unsorted(numfiles,
                                                  (const char **)argv + optind);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(options.pseudofix)
  {
    current_stream = agn_pseudogene_fix_stream_new(last_stream);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  current_stream = agn_locus_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(options.intloci)
  {
    current_stream = agn_interval_locus_stream_new(last_stream, options.delta,
                                                   options.skipends, logger);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  /* FIXME I don't understand why this is needed, but without it memory is
   * leaked; if it's not in this precise location, no memory is leaked but
   * segfaults result */
  current_stream = agn_node_delete_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(options.genestream != NULL || options.transstream != NULL)
  {
    current_stream = agn_locus_map_stream_new(last_stream, options.genestream,
                                              options.transstream);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  if(options.verbose == 0)
  {
    current_stream = agn_remove_children_stream_new(last_stream);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

  current_stream = gt_gff3_out_stream_new(last_stream, options.outstream);
  if(options.retainids)
    gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;


  //----- Execute the node processing stream -----//
  //----------------------------------------------//

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
    fprintf(stderr, "[LocusPocus] error: %s", gt_error_get(error));


  // Free memory and terminate
  options.filefreefunc(options.outstream);
  if(options.genestream != NULL)  fclose(options.genestream);
  if(options.transstream != NULL) fclose(options.transstream);
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_queue_delete(streams);
  gt_logger_delete(logger);
  gt_error_delete(error);
  gt_lib_clean();
  return 0;
}
