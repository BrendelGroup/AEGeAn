#include <getopt.h>
#include <string.h>
#include "genometools.h"
#include "aegean.h"

// Simple data structure for program options
typedef struct
{
  bool debug;
  FILE *outfile;
  const char *outfilename;
  bool gff3;
  bool summary_only;
  bool graphics;
  const char *refrfile;
  const char *predfile;
  const char *refrlabel;
  const char *predlabel;
  const char *outfmt;
  bool overwrite;
  const char *data_path;
  bool makefilter;
  bool usefilter;
  const char *filterfile;
  // AgnCompareFilters filters;
  bool verbose;
} ParsEvalOptions;

// Set default values for program
static void set_option_defaults(ParsEvalOptions *options)
{
  options->debug = false;
  options->outfile = stdout;
  options->gff3 = true;
  options->summary_only = false;
  options->graphics = false;
  options->refrfile = NULL;
  options->predfile = NULL;
  options->refrlabel = NULL;
  options->predlabel = NULL;
  options->outfmt = "text";
  options->overwrite = false;
  options->data_path = AGN_DATA_PATH;
  options->makefilter = false;
  options->usefilter = false;
  options->filterfile = NULL;
  // AgnCompareFilters filters;
  options->verbose = false;
}

static void free_option_memory(ParsEvalOptions *options)
{
  fclose(options->outfile);
}

// Usage statement
static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\nParsEval: comparative analysis of two alternative sources of annotation\n"
"Usage: parseval [options] reference.gff3 prediction.gff3\n"
"  Basic options:\n"
"    -d|--debug:                 Print debugging messages\n"
"    -h|--help:                  Print help message and exit\n"
"    -v|--verbose:               Print verbose warning messages\n\n"
"  Output options:\n"
"    -a|--datashare: STRING      Location from which to copy shared data for\n"
"                                HTML output (if `make install' has not yet\n"
"                                been run)\n"
"    -f|--outformat: STRING      Indicate desired output format; possible\n"
"                                options: 'csv', 'text', or 'html'\n"
"                                (default='text'); in 'text' or 'csv' mode,\n"
"                                will create a single file; in 'html' mode,\n"
"                                will create a directory\n"
"    -g|--nogff3:                Do no print GFF3 output corresponding to each\n"
"                                comparison\n"
"    -o|--outfile: FILENAME      File/directory to which output will be\n"
"                                written; default is the terminal (STDOUT)\n"
"    -p|--png:                   Generate individual PNG graphics for each\n"
"                                gene locus (HTML mode only)\n"
"    -s|--summary:               Only print summary statistics, do not print\n"
"                                individual comparisons\n"
"    -w|--overwrite:             Force overwrite of any existing output files\n"
"    -x|--refrlabel: STRING      Optional label for reference annotations\n"
"    -y|--predlabel: STRING      Optional label for prediction annotations\n\n"
"  Filtering options:\n"
"    -k|--makefilter             Create a default configuration file for\n"
"                                filtering reported results and quit,\n"
"                                performing no comparisons\n"
"    -r|--filterfile: STRING     Use the indicated configuration file to\n"
"                                filter reported results;\n\n");
}

// Adjust program settings from command-line arguments/options
static int
parse_options(int argc, char **argv, ParsEvalOptions *options, GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "a:df:ghko:pr:svwx:y:";
  const struct option parseval_options[] =
  {
    { "datashare",  required_argument, NULL, 'a' },
    { "debug",      no_argument,       NULL, 'd' },
    { "outformat",  required_argument, NULL, 'f' },
    { "printgff3",  no_argument,       NULL, 'g' },
    { "help",       no_argument,       NULL, 'h' },
    { "makefilter", no_argument,       NULL, 'k' },
    { "outfile",    required_argument, NULL, 'o' },
    { "png",        no_argument,       NULL, 'p' },
    { "filterfile", required_argument, NULL, 'r' },
    { "summary",    no_argument,       NULL, 's' },
    { "verbose",    no_argument,       NULL, 'v' },
    { "overwrite",  no_argument,       NULL, 'w' },
    { "refrlabel",  required_argument, NULL, 'x' },
    { "predlabel",  required_argument, NULL, 'y' },
    { NULL,         no_argument,       NULL,  0  },
  };

  for(opt  = getopt_long(argc, argv, optstr, parseval_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv, optstr, parseval_options, &optindex))
  {
    switch(opt)
    {
      case 'a':
        options->data_path = optarg;
        break;

      case 'd':
        options->debug = true;
        break;

      case 'f':
        if(strcmp(optarg,  "csv") != 0 &&
           strcmp(optarg, "text") != 0 &&
           strcmp(optarg, "html") != 0)
        {
          fprintf(stderr, "error: unknown value '%s' for '-f|--outformat' "
                  "option\n\n", optarg);
          print_usage(stderr);
          exit(1);
        }
        options->outfmt = optarg;
        break;

      case 'g':
        options->gff3 = false;
        break;

      case 'h':
        print_usage(stdout);
        exit(0);
        break;

      case 'k':
        options->makefilter = true;
        break;

      case 'o':
        options->outfilename = optarg;
        break;

      case 'p':
        options->graphics = true;
#ifdef WITHOUT_CAIRO
        fputs("error: AEGeAn was compiled without graphics support. Please "
              "recompile to enable this feature.\n", stderr);
        exit(1);
#endif
        break;

      case 'r':
        options->usefilter = true;
        options->filterfile = optarg;
        /*if(options->usefilter)
        {
          if(options->debug)
            fprintf(stderr, "debug: opening filter file '%s'\n",
                    options->filterfile);

          FILE *filterfile = agn_fopen(options->filterfile, "r", stderr);
          AgnLogger *logger = agn_logger_new();
          agn_compare_filters_parse(&options->filters, filterfile, logger);
          bool haderrors = agn_logger_print_all(logger, stderr,
                                                "[ParsEval] parsing filters");
          if(haderrors)
            exit(1);
          agn_logger_delete(logger);
          if(options->debug)
            fprintf(stderr, "debug: closing filter file\n");
          fclose(filterfile);
        }*/
        break;

      case 's':
        options->summary_only = true;
        break;

      case 'v':
        options->verbose = true;
        break;

      case 'w':
        options->overwrite = true;
        break;

      case 'x':
        options->refrlabel = optarg;
        break;

      case 'y':
        options->predlabel = optarg;
        break;

      default:
        break;
    }
  }

  // For debugging
  // pe_option_print(options, stderr);

  if(options->makefilter)
  {
    char cmd[512];
    sprintf(cmd, "cp %s/pe.filter pe.filter", options->data_path);
    if(options->debug)
      fprintf(stderr, "debug: creating filter file '%s'\n", cmd);
    fputs("Created new filter file 'pe.filter'\n", stderr);
    if(system(cmd) != 0)
    {
      fprintf(stderr, "error: could not create filter file 'pe.filter'\n");
      exit(1);
    }
    exit(0);
  }

  if(argc - optind != 2)
  {
    fprintf(stderr, "error: must provide 2 (and only 2) input files, you "
            "provided %d\n\n", argc - optind);
    print_usage(stderr);
    exit(1);
  }

  if(options->outfilename)
  {
    if(strcmp(options->outfmt, "html") == 0)
    {
      char dircmd[1024];
      sprintf(dircmd, "test -d %s", options->outfilename);
      if(system(dircmd) == 0)
      {
        if(options->overwrite)
        {
          char rmcmd[1024];
          sprintf(rmcmd, "rm -r %s", options->outfilename);
          if(system(rmcmd) != 0)
          {
            fprintf(stderr, "error: could not overwrite output directory '%s'",
                    options->outfilename);
            exit(1);
          }
        }
        else
        {
          fprintf(stderr, "error: outfile '%s' exists; use '-w' to force "
                  "overwrite\n", options->outfilename);
          exit(1);
        }
      }
    }
    else
    {
      char filecmd[1024];
      sprintf(filecmd, "test -f %s", options->outfilename);
      if(system(filecmd) == 0 && !options->overwrite)
      {
          fprintf(stderr, "error: outfile '%s' exists; use '-w' to force "
                  "overwrite\n", options->outfilename);
          exit(1);
      }
    }
  }

  if(strcmp(options->outfmt, "html") == 0)
  {
    if(options->outfilename == NULL)
    {
      fputs("error: will not print results to terminal in HTML mode; must "
            "provide outfile\n\n", stderr);
      print_usage(stderr);
      exit(1);
    }
    char dircmd[1024];
    sprintf(dircmd, "mkdir %s", options->outfilename);
    if(system(dircmd) != 0)
    {
      fprintf(stderr, "error: cannot open output directory '%s'\n",
              options->outfilename);
      exit(1);
    }
    char outname[1024];
    sprintf(outname, "%s/index.html", options->outfilename);
    options->outfile = fopen(outname, "w");
    if(!options->outfile)
    {
      fprintf(stderr, "error: could not open output file '%s'\n", outname);
      exit(1);
    }

    char copy_cmd[1024];
    sprintf(copy_cmd,"cp -r %s/* %s", options->data_path, options->outfilename);
    if(options->debug)
      fprintf(stderr, "debug: copying shared data: '%s'\n", copy_cmd);
    if(system(copy_cmd) != 0)
    {
      fprintf(stderr, "error: could not copy data files '%s'\n", copy_cmd);
      exit(1);
    }

    if(options->summary_only && options->graphics)
    {
      fprintf(stderr, "warning: cannot print PNG graphics in summary only "
              "mode; ignoring\n");
      options->graphics = false;
    }
  }
  else
  {
    if(options->graphics)
    {
      fputs("warning: will only generate PNG graphics when outformat='html'; "
            "ignoring\n\n", stderr);
      options->graphics = false;
    }
    if(options->outfilename)
    {
      options->outfile = fopen(options->outfilename, "w");
      if(options->outfile == NULL)
      {
        fprintf(stderr, "error: cannot open output file '%s'\n",
                options->outfilename);
        exit(1);
      }
    }
  }

  options->refrfile = argv[optind];
  options->predfile = argv[optind + 1];
  return optind;
}

// Main program
int main(int argc, char **argv)
{
  GtError *error;
  GtLogger *logger;
  GtQueue *streams;
  GtNodeStream *current_stream, *last_stream, *refrgff3, *predgff3;
  gt_lib_init();

  // Parse command-line options
  ParsEvalOptions options;
  set_option_defaults(&options);
  error = gt_error_new();
  parse_options(argc, argv, &options, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[ParsEval] error: %s", gt_error_get(error));
    return 1;
  }
  int numfiles = argc - optind;
  if(numfiles != 2)
  {
    fprintf(stderr, "[ParsEval] error: must provide two GFF3 files as input");
    print_usage(stderr);
    return 1;
  }

  logger = gt_logger_new(true, "", stderr);
  streams = gt_queue_new();


  //----- Set up the node processing stream -----//
  //---------------------------------------------//

  refrgff3 = gt_gff3_in_stream_new_unsorted(numfiles, (const char **)argv+optind);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)refrgff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)refrgff3);
  gt_queue_add(streams, refrgff3);

  predgff3 = gt_gff3_in_stream_new_unsorted(numfiles, (const char **)argv+optind+1);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)predgff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)predgff3);
  gt_queue_add(streams, predgff3);

  current_stream = agn_locus_stream_new_pairwise(refrgff3, predgff3, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  /*if(options.pseudofix)
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
  last_stream = current_stream;*/

  /*GtStr *source = gt_str_new_cstr("AEGeAn::ParsEval");
  current_stream = agn_locus_stream_new(last_stream, logger);
  agn_locus_stream_set_source((AgnLocusStream *)current_stream, source);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;
  gt_str_delete(source);*/

  /* FIXME I don't understand why this is needed, but without it memory is
   * leaked; if it's not in this precise location, no memory is leaked but
   * segfaults result */
  current_stream = agn_node_delete_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_compare_report_text_new(last_stream, options.outfile,
                                               logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;



  //----- Execute the node processing stream -----//
  //----------------------------------------------//

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
    fprintf(stderr, "[ParsEval] error: %s", gt_error_get(error));

  agn_compare_report_text_create_summary((AgnCompareReportText *)last_stream,
                                         options.outfile);

  // Free memory and terminate
  free_option_memory(&options);
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
