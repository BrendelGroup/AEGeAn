/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include <string.h>
#include "genometools.h"
#include "aegean.h"

// Data structure for program options
typedef struct
{
  bool debug;
  GtHashmap *filter;
  FILE *genestream;
  char *nameformat;
  unsigned long delta;
  GtFile *outstream;
  void (*filefreefunc)(GtFile *);
  GtHashmap *type_parents;
  int endmode;
  FILE *transstream;
  bool pseudofix;
  bool verbose;
  bool skipiiLoci;
  bool refine;
  bool by_cds;
  GtUword minoverlap;
  FILE *ilenfile;
  bool retain;
} LocusPocusOptions;

// Set default values for program
static void set_option_defaults(LocusPocusOptions *options)
{
  options->debug = false;
  options->filter = gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
  gt_hashmap_add(options->filter, gt_cstr_dup("gene"), gt_cstr_dup("gene"));
  options->genestream = NULL;
  options->nameformat = NULL;
  options->delta = 500;
  options->outstream = gt_file_new_from_fileptr(stdout);
  options->filefreefunc = gt_file_delete_without_handle;
  options->type_parents = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                         gt_free_func);
  options->endmode = 0;
  options->transstream = NULL;
  options->pseudofix = false;
  options->verbose = false;
  options->skipiiLoci = false;
  options->refine = false;
  options->by_cds = false;
  options->minoverlap = 1;
  options->ilenfile = NULL;
  options->retain = false;
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
  if(options->nameformat != NULL)
    gt_free(options->nameformat);
  if(options->ilenfile != NULL)
    fclose(options->ilenfile);
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
"    -s|--skipends          when enumerating interval loci, exclude\n"
"                           unannotated (and presumably incomplete) iLoci at\n"
"                           either end of the sequence\n"
"    -e|--endsonly          report only incomplete iLocus fragments at the\n"
"                           unannotated ends of sequences (complement of\n"
"                           --skipends)\n"
"    -y|--skipiiloci        do not report intergenic iLoci\n\n"
"  Refinement options:\n"
"    -r|--refine            by default genes are grouped in the same iLocus\n"
"                           if they have any overlap; 'refine' mode allows\n"
"                           for a more nuanced handling of overlapping genes\n"
"    -c|--cds               use CDS rather than UTRs for determining gene\n"
"                           overlap; implies 'refine' mode\n"
"    -m|--minoverlap: INT   the minimum number of nucleotides two genes must\n"
"                           overlap to be grouped in the same iLocus; default\n"
"                           is 1\n\n"
"  Output options:\n"
"    -n|--namefmt: STR     provide a printf-style format string to override\n"
"                           the default ID format for newly created loci;\n"
"                           default is 'locus%%lu' (locus1, locus2, etc) for\n"
"                           loci and 'iLocus%%lu' (iLocus1, iLocus2, etc) for\n"
"                           interval loci; note the format string should\n"
"                           include a single %%lu specifier to be filled in\n"
"                           with a long unsigned integer value\n"
"    -i|--ilens: FILE       create a file with the lengths of each intergenic\n"
"                           iLocus\n"
"    -g|--genemap: FILE     print a mapping from each gene annotation to its\n"
"                           corresponding locus to the given file\n"
"    -o|--outfile: FILE     name of file to which results will be written;\n"
"                           default is terminal (standard output)\n"
"    -T|--retainids         retain original feature IDs from input files;\n"
"                           conflicts will arise if input contains duplicated\n"
"                           ID values\n"
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
  const char *optstr = "cdef:g:hi:l:m:n:o:p:rsTt:uVvy";
  const char *key, *value, *oldvalue;
  const struct option locuspocus_options[] =
  {
    { "cds",        no_argument,       NULL, 'c' },
    { "debug",      no_argument,       NULL, 'd' },
    { "endsonly",   no_argument,       NULL, 'e' },
    { "filter",     required_argument, NULL, 'f' },
    { "genemap",    required_argument, NULL, 'g' },
    { "help",       no_argument,       NULL, 'h' },
    { "ilens",      required_argument, NULL, 'i' },
    { "delta",      required_argument, NULL, 'l' },
    { "minoverlap", required_argument, NULL, 'm' },
    { "namefmt",    required_argument, NULL, 'n' },
    { "outfile",    required_argument, NULL, 'o' },
    { "parent",     required_argument, NULL, 'p' },
    { "refine",     no_argument,       NULL, 'r' },
    { "skipends",   no_argument,       NULL, 's' },
    { "retainids",  no_argument,       NULL, 'T' },
    { "transmap",   required_argument, NULL, 't' },
    { "pseudo",     no_argument,       NULL, 'u' },
    { "version",    no_argument,       NULL, 'v' },
    { "verbose",    no_argument,       NULL, 'V' },
    { "skipiiloci", no_argument,       NULL, 'y' },
    { NULL,         no_argument,       NULL,  0  },
  };
  for( opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex))
  {
    if(opt == 'c')
    {
      options->by_cds = 1;
      options->refine = 1;
    }
    else if(opt == 'd')
      options->debug = 1;
    else if(opt == 'e')
    {
      if(options->endmode < 0)
      {
        gt_error_set(error, "cannot set 'skipends' and 'endsonly' options"
                     "simultaneously; stubbornly refusing to proceed");
      }
      options->endmode = 1;
    }
    else if(opt == 'f')
    {
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
    }
    else if(opt == 'g')
    {
      options->genestream = fopen(optarg, "w");
      if(options->genestream == NULL)
        gt_error_set(error, "could not open genemap file '%s'", optarg);
    }
    else if(opt == 'h')
    {
      print_usage(stdout);
      exit(0);
    }
    else if(opt == 'i')
    {
      options->ilenfile = fopen(optarg, "w");
      if(options->ilenfile == NULL)
        gt_error_set(error, "could not open ilenfile file '%s'", optarg);
    }
    else if(opt == 'l')
    {
      if(sscanf(optarg, "%lu", &options->delta) == EOF)
      {
        gt_error_set(error, "could not convert delta '%s' to an integer",
                     optarg);
      }
    }
    else if(opt == 'm')
    {
      if(sscanf(optarg, "%lu", &options->minoverlap) == EOF)
      {
        gt_error_set(error, "could not convert overlap '%s' to an integer",
                     optarg);
      }
    }
    else if(opt == 'n')
    {
      if(options->nameformat != NULL)
        gt_free(options->nameformat);
      options->nameformat = gt_cstr_dup(optarg);
    }
    else if(opt == 'o')
    {
      options->filefreefunc(options->outstream);
      options->outstream = gt_file_new(optarg, "w", error);
      options->filefreefunc = gt_file_delete;
    }
    else if(opt == 'p')
    {
      key = strtok(optarg, ":");
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
    }
    else if(opt == 'r')
      options->refine = 1;
    else if(opt == 's')
    {
      if(options->endmode > 0)
      {
        gt_error_set(error, "cannot set 'skipends' and 'endsonly' options"
                     "simultaneously; stubbornly refusing to proceed");
      }
      options->endmode = -1;
    }
    else if(opt == 'T')
      options->retain = true;
    else if(opt == 't')
    {
      options->transstream = fopen(optarg, "w");
      if(options->transstream == NULL)
        gt_error_set(error, "could not open transmap file '%s'", optarg);
    }
    else if(opt == 'u')
      options->pseudofix = 1;
    else if(opt == 'v')
    {
      agn_print_version("LocusPocus", stdout);
      exit(0);
    }
    else if(opt == 'V')
      options->verbose = 1;
    else if(opt == 'y')
      options->skipiiLoci = true;
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
  AgnLocusStream *ls = (AgnLocusStream*)current_stream;
  agn_locus_stream_set_source(ls, "AEGeAn::LocusPocus");
  agn_locus_stream_set_endmode(ls, options.endmode);
  agn_locus_stream_track_ilens(ls, options.ilenfile);
  if(options.nameformat != NULL)
    agn_locus_stream_set_name_format(ls, options.nameformat);
  if(options.skipiiLoci)
    agn_locus_stream_skip_iiLoci(ls);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(options.refine)
  {
    current_stream = agn_locus_refine_stream_new(last_stream, options.delta,
                                                 options.minoverlap,
                                                 options.by_cds);
    AgnLocusRefineStream *lrs = (AgnLocusRefineStream *)current_stream;
    agn_locus_refine_stream_set_source(lrs, "AEGeAn::LocusPocus");
    agn_locus_refine_stream_track_ilens(lrs, options.ilenfile);
    if(options.nameformat != NULL)
      agn_locus_refine_stream_set_name_format(lrs, options.nameformat);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
  }

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
  if(options.retain)
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
