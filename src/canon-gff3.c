#include <getopt.h>
#include <string.h>
#include "AgnCanonGeneStream.h"
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
  int code = canon_gff3_parse_options(argc, argv, &options);
  if(code)
  {
    if(code < 1)
      return 0;
    else
      return code;
  }

  // Create input and output streams to load, process, and write the data
  AgnLogger *logger = agn_logger_new();
  GtFile *outfile = gt_file_new_from_fileptr(options.outstream);
  GtNodeVisitor *ssv = gt_set_source_visitor_new(options.source);

  GtNodeStream *gff3in;
  gff3in = gt_gff3_in_stream_new_unsorted(options.numfiles, options.gff3files);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtNodeStream *cgstream = agn_canon_gene_stream_new(gff3in, logger);
  GtNodeStream *ssstream = gt_visitor_stream_new(cgstream, ssv);
  GtNodeStream *gff3out  = gt_gff3_out_stream_new(ssstream, outfile);
  gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream *)gff3out);

  GtError *error = gt_error_new();
  int result = gt_node_stream_pull(gff3out, error);
  if(result == -1)
  {
    fprintf(stderr, "error processing node stream: %s", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(cgstream);
  gt_node_stream_delete(ssstream);
  gt_node_stream_delete(gff3out);
  gt_file_delete_without_handle(outfile);
  gt_error_delete(error);

  bool haderror = agn_logger_print_all(logger, stderr, "[CanonGFF3] processing "
                                       "annotation data");
  if(haderror) return EXIT_FAILURE;
  agn_logger_delete(logger);

  // Clean up
  fclose(options.outstream);
  gt_str_delete(options.source);
  if(gt_lib_clean() != 0)
  {
    fputs("error: issue cleaning GenomeTools library\n", stderr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
