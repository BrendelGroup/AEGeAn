/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include <string.h>
#include "genometools.h"
#include "aegean.h"

// Simple data structure for program options
typedef struct
{
  bool debug;
  GtHashmap *filter;
  FILE *genestream;
  char *idformat;
  unsigned long delta;
  GtFile *outstream;
  void (*filefreefunc)(GtFile *);
  GtHashmap *type_parents;
  int endmode;
  FILE *transstream;
  bool pseudofix;
  bool verbose;
  bool skipempty;
  bool by_cds;
  GtUword minoverlap;
} LocusPocusOptions;

// Set default values for program
static void set_option_defaults(LocusPocusOptions *options)
{
  options->debug = false;
  options->filter = gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
  gt_hashmap_add(options->filter, gt_cstr_dup("gene"), gt_cstr_dup("gene"));
  options->genestream = NULL;
  options->idformat = NULL;
  options->delta = 500;
  options->outstream = gt_file_new_from_fileptr(stdout);
  options->filefreefunc = gt_file_delete_without_handle;
  options->type_parents = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                         gt_free_func);
  options->endmode = 0;
  options->transstream = NULL;
  options->pseudofix = false;
  options->verbose = false;
  options->skipempty = false;
  options->by_cds = false;
  options->minoverlap = 1;
}

static void free_option_memory(LocusPocusOptions *options)
{
  options->filefreefunc(options->outstream);
  gt_hashmap_delete(options->type_parents);
  gt_hashmap_delete(options->filter);
  if(options->genestream != NULL)
    fclose(options->genestream);
  if(options->transstream != NULL)
    fclose(options->transstream);
  if(options->idformat != NULL)
    gt_free(options->idformat);
}

// Usage statement
static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\nLocusPocus: calculate locus coordinates for the given gene annotation\n"
"Usage: locuspocus [options] gff3file1 [gff3file2 gff3file3 ...]\n"
"  Basic options:\n"
"    -d|--debug             print detailed debugging messages to terminal\n"
"                           (standard error)\n"
"    -h|--help              print this help message and exit\n"
"    -v|--version           print version number and exit\n\n"
"  iLocus parsing:\n"
"    -l|--delta: INT        when parsing interval loci, use the following\n"
"                           delta to extend gene loci and include potential\n"
"                           regulatory regions; default is 500\n"
"    -r|--overlap: INT      specify the minimum overlap (in bp) required to\n"
"                           group two genes together in the same locus;\n"
"                           default is 1\n"
"    -c|--cds               when applicable, use coding sequence rather than\n"
"                           untranslated regions when determining gene\n"
"                           overlap\n"
"    -s|--skipends          when enumerating interval loci, exclude gene-less\n"
"                           iLoci at either end of the sequence\n"
"    -e|--endsonly          report only empty iLoci at the ends of sequences\n"
"                           (complement of --skipends)\n"
"    -y|--skipempty         do not report empty iLoci\n\n"
"  Output options:\n"
"    -I|--idformat: STR     provide a printf-style format string to override\n"
"                           the default ID format for newly created loci;\n"
"                           default is 'locus%%lu' (locus1, locus2, etc) for\n"
"                           loci and 'iLocus%%lu' (iLocus1, iLocus2, etc) for\n"
"                           interval loci; note the format string should\n"
"                           include a single %%lu specifier to be filled in\n"
"                           with a long unsigned integer value\n"
"    -g|--genemap: FILE     print a mapping from each gene annotation to its\n"
"                           corresponding locus to the given file\n"
"    -o|--outfile: FILE     name of file to which results will be written;\n"
"                           default is terminal (standard output)\n"
"    -t|--transmap: FILE    print a mapping from each transcript annotation\n"
"                           to its corresponding locus to the given file\n"
"    -V|--verbose           include all locus subfeatures (genes, RNAs, etc)\n"
"                           in the GFF3 output; default includes only locus\n"
"                           features\n\n"
"  Input options:\n"
"    -f|--filter: TYPE      comma-separated list of feature types to use in\n"
"                           constructing loci/iLoci; default is 'gene'\n"
"    -p|--parent: CT:PT     if a feature of type $CT exists without a parent,\n"
"                           create a parent for this feature with type $PT;\n"
"                           for example, mRNA:gene will create a gene feature\n"
"                           as a parent for any top-level mRNA feature;\n"
"                           this option can be specified multiple times\n"
"    -u|--pseudo            correct erroneously labeled pseudogenes\n\n");
}

// Adjust program settings from command-line arguments/options
static void
parse_options(int argc, char **argv, LocusPocusOptions *options, GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "cdef:g:hI:l:o:p:r:st:uVvy";
  const char *key, *value, *oldvalue;
  const struct option locuspocus_options[] =
  {
    { "cds",       no_argument,       NULL, 'c' },
    { "debug",     no_argument,       NULL, 'd' },
    { "endsonly",  no_argument,       NULL, 'e' },
    { "filter",    required_argument, NULL, 'f' },
    { "genemap",   required_argument, NULL, 'g' },
    { "help",      no_argument,       NULL, 'h' },
    { "idformat",  required_argument, NULL, 'I' },
    { "delta",     required_argument, NULL, 'l' },
    { "outfile",   required_argument, NULL, 'o' },
    { "parent",    required_argument, NULL, 'p' },
    { "overlap",   required_argument, NULL, 'r' },
    { "skipends",  no_argument,       NULL, 's' },
    { "transmap",  required_argument, NULL, 't' },
    { "pseudo",    no_argument,       NULL, 'u' },
    { "version",   no_argument,       NULL, 'v' },
    { "verbose",   no_argument,       NULL, 'V' },
    { "skipempty", no_argument,       NULL, 'y' },
  };
  for( opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex))
  {
    switch(opt)
    {
      case 'c':
        options->by_cds = 1;
        break;
      case 'd':
        options->debug = 1;
        break;
      case 'e':
        if(options->endmode < 0)
        {
          gt_error_set(error, "cannot set 'skipends' and 'endsonly' options"
                       "simultaneously; stubbornly refusing to proceed");
        }
        options->endmode = 1;
        break;
      case 'f':
        gt_hashmap_delete(options->filter);
        options->filter = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                         gt_free_func);

        value = strtok(optarg, ",");
        gt_hashmap_add(options->filter, gt_cstr_dup(value), gt_cstr_dup(value));
        while((value = strtok(NULL, ",")) != NULL)
        {
          oldvalue = gt_hashmap_get(options->filter, value);
          if(!oldvalue)
          {
            gt_hashmap_add(options->filter, gt_cstr_dup(value),
                           gt_cstr_dup(value));
          }
        }
        break;
      case 'g':
        options->genestream = fopen(optarg, "w");
        if(options->genestream == NULL)
          gt_error_set(error, "could not open genemap file '%s'", optarg);
        break;
      case 'h':
        print_usage(stdout);
        exit(0);
        break;
      case 'I':
        if(options->idformat != NULL)
          gt_free(options->idformat);
        options->idformat = gt_cstr_dup(optarg);
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
      case 'p':
        key   = strtok(optarg, ":");
        value = strtok(NULL, ":");
        oldvalue = gt_hashmap_get(options->type_parents, "key");
        if(oldvalue)
        {
          gt_error_set(error, "cannot set parent type for '%s' to '%s'; "
                       "already set to '%s'", key, value, oldvalue);
        }
        else
        {
          gt_hashmap_add(options->type_parents, gt_cstr_dup(key),
                         gt_cstr_dup(value));
        }
        break;
      case 'r':
        if(sscanf(optarg, "%lu", &options->minoverlap) == EOF)
        {
          gt_error_set(error, "could not convert overlap '%s' to an integer",
                       optarg);
        }
        break;
      case 's':
        if(options->endmode > 0)
        {
          gt_error_set(error, "cannot set 'skipends' and 'endsonly' options"
                       "simultaneously; stubbornly refusing to proceed");
        }
        options->endmode = -1;
        break;
      case 't':
        options->transstream = fopen(optarg, "w");
        if(options->transstream == NULL)
          gt_error_set(error, "could not open transmap file '%s'", optarg);
        break;
      case 'u':
        options->pseudofix = 1;
        break;
      case 'v':
        agn_print_version("LocusPocus", stdout);
        exit(0);
        break;
      case 'V':
        options->verbose = 1;
        break;
      case 'y':
        options->skipempty = true;
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

  current_stream = agn_infer_parent_stream_new(last_stream,
                                               options.type_parents);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_filter_stream_new(last_stream, options.filter);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = gt_sort_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, options.delta);
  agn_locus_stream_set_source((AgnLocusStream *)current_stream,
                              "AEGeAn::LocusPocus");
  agn_locus_stream_set_endmode((AgnLocusStream*)current_stream,options.endmode);
  agn_locus_stream_set_overlap((AgnLocusStream*)current_stream,
                               options.minoverlap);
  if(options.idformat != NULL)
  {
    agn_locus_stream_set_idformat((AgnLocusStream *)current_stream,
                                  options.idformat);
  }
  if(options.skipempty)
  {
    agn_locus_stream_skip_empty_loci((AgnLocusStream *)current_stream);
  }
  if(options.by_cds)
  {
    agn_locus_stream_by_cds((AgnLocusStream *)current_stream);
  }
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
  gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;


  //----- Execute the node processing stream -----//
  //----------------------------------------------//

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
    fprintf(stderr, "[LocusPocus] error: %s", gt_error_get(error));


  // Free memory and terminate
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_queue_delete(streams);
  gt_logger_delete(logger);
  gt_error_delete(error);
  free_option_memory(&options);
  gt_lib_clean();
  return 0;
}
