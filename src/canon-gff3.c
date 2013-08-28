#include <getopt.h>
#include <string.h>
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

typedef struct
{
  FILE *outstream;
  GtStr *source;
  bool read_stdin;
  const char **gff3files;
  int numfiles;
} CanonGFF3Options;

/**
 * Print the usage statement for CanonGFF3
 *
 * @param[out] outstream    stream to which the statement will be written
 */
void print_usage(FILE *outstream)
{
  fputs("Usage: canon-gff3 [options] gff3file1 [gff3file2 ...]\n"
"  Options:\n"
"     -h|--help               print this help message and exit\n"
"     -o|--outfile: STRING    name of file to which GFF3 data will be\n"
"                             written; default is terminal (stdout)\n"
"     -s|--source: STRING     reset the source of each feature to the given\n"
"                             value\n"
"     -t|--stdin              read input from terminal (stdin)\n",
        outstream);
}

int canon_gff3_parse_options(int argc, char * const *argv,
                             CanonGFF3Options *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "ho:s:t";
  const struct option init_options[] =
  {
    { "help",    no_argument,       NULL, 'h' },
    { "outfile", required_argument, NULL, 'o' },
    { "source",  required_argument, NULL, 's' },
    { "stdin",   no_argument,       NULL, 't' },
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
        return -1;
        break;

      case 'o':
        options->outstream = agn_fopen(optarg, "w", stderr);
        break;

      case 's':
        options->source = gt_str_new_cstr(optarg);
        break;

      case 't':
        options->read_stdin = true;
        break;

      default:
        break;
    }
  }

  // Validate options and arguments
  options->numfiles = argc - optind;
  if(options->read_stdin)
  {
    if(options->numfiles > 0)
    {
      fprintf(stderr, "[CanonGFF3] warning: reading from standard input, "
              "ignoring %d files\n", options->numfiles);
      options->numfiles = 0;
    }
  }
  else
  {
    if(options->numfiles < 1)
    {
      fprintf(stderr, "[CanonGFF3] error: must provide 1 or more GFF3 files to "
              "process\n\n");
      print_usage(stderr);
      return 1;
    }
  }
  
  // Create a char ** of the GFF3 filenames
  int x;
  if(options->numfiles > 0)
  {
    options->gff3files = gt_malloc( options->numfiles * sizeof(char *) );
    for(x = 0; x < options->numfiles; x++)
    {
      options->gff3files[x] = gt_cstr_dup(argv[x+optind]);
    }
  }

  return 0;
}

// Main method
int main(int argc, char * const *argv)
{
  // Options
  gt_lib_init();
  CanonGFF3Options options = { stdout, NULL, false, NULL, 0 };
  int result = canon_gff3_parse_options(argc, argv, &options);
  if(result)
  {
    if(result < 1)
      return 0;
    else
      return result;
  }

  // Load GFF3 data into memory feature index
  AgnLogger *logger = agn_logger_new();
  GtFeatureIndex *features = agn_import_canonical(options.numfiles,
                                                  options.gff3files, logger);
  unsigned long i;
  for(i = 0; i < options.numfiles; i++)
  {
    gt_free((char *)options.gff3files[i]);
  }
  gt_free(options.gff3files);
  bool haderror = agn_logger_print_all(logger, stderr, "[CanonGFF3] load "
                                       "annotation data into memory");
  if(haderror) return EXIT_FAILURE;

  // Write header for new GFF3 file: version pragma, sequence regions
  GtError *error = gt_error_new();
  fputs("##gff-version   3\n", options.outstream);
  GtStrArray *seqids = gt_feature_index_get_seqids(features, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[CanonGFF3] error fetching sequence IDs: %s",
            gt_error_get(error));
    return 1;
  }
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtRange seqrange;
    int code = gt_feature_index_get_orig_range_for_seqid(features, &seqrange,
                                                         seqid, error);
    if(gt_error_is_set(error))
    {
      fprintf(stderr, "[CanonGFF3] error fetching range for sequence '%s': %s",
              seqid, gt_error_get(error));
      return code;
    }
    fprintf(options.outstream, "##sequence-region   %s %lu %lu\n", seqid,
            seqrange.start, seqrange.end);
  }

  // Write features to new GFF3 file
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *genes = gt_feature_index_get_features_for_seqid(features, seqid,
                                                             error);
    if(gt_error_is_set(error))
    {
      fprintf(stderr, "[CanonGFF3] error fetching features for sequence '%s': "
              "%s", seqid, gt_error_get(error));
      return 1;
    }
    gt_array_sort(genes, (GtCompare)agn_gt_genome_node_compare);

    unsigned long j;
    for(j = 0; j < gt_array_size(genes); j++)
    {
      GtFeatureNode *gene = *(GtFeatureNode **)gt_array_get(genes, j);
      if(options.source != NULL)
        agn_gt_feature_node_set_source_recursive(gene, options.source);
      agn_gt_feature_node_to_gff3(gene, options.outstream, true, NULL, NULL);
      fputs("###\n", options.outstream);
    }
  }

  // Clean up
  fclose(options.outstream);
  gt_str_delete(options.source);
  gt_error_delete(error);
  if(gt_lib_clean() != 0)
  {
    fputs("error: issue cleaning GenomeTools library\n", stderr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
