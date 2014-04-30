#include "pe_options.h"

void pe_free_option_memory(ParsEvalOptions *options)
{
  fclose(options->outfile);
  gt_array_delete(options->filters);
}

char *pe_get_start_time()
{
  time_t start_time;
  struct tm *start_time_info;
  time(&start_time);
  start_time_info = localtime(&start_time);

  char timestr[128];
  strftime(timestr, 128, "%d %b %Y, %I:%M%p", start_time_info);
  return gt_cstr_dup(timestr);
}

int pe_parse_options(int argc, char **argv, ParsEvalOptions *options,
                     GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "a:df:ghko:pr:st:vwx:y:";
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
    { "maxtrans",   required_argument, NULL, 't' },
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
          pe_print_usage(stderr);
          exit(1);
        }
        options->outfmt = optarg;
        break;

      case 'g':
        options->gff3 = false;
        break;

      case 'h':
        pe_print_usage(stdout);
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
        if(true)
        {
          FILE *filterfile = fopen(optarg, "r");
          if(filterfile == NULL)
          {
            gt_error_set(error, "unable to open filter file '%s'", optarg);
            return -1;
          }
          agn_locus_filter_parse(filterfile, options->filters);
          fclose(filterfile);
        }
        break;

      case 's':
        options->summary_only = true;
        break;

      case 't':
        if(sscanf(optarg, "%d", &options->max_transcripts) == EOF)
        {
          fprintf(stderr, "error: could not convert maxtrans '%s' to an "
                          "integer\n", optarg);
          exit(1);
        }
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

  if(options->max_transcripts > 0)
  {
    AgnLocusFilter filter;
    filter.testvalue = options->max_transcripts;
    filter.operator = AGN_LOCUS_FILTER_LE;
    filter.src = DEFAULTSOURCE;
    filter.function = agn_locus_mrna_num;
    gt_array_add(options->filters, filter);
  }

  if(argc - optind != 2)
  {
    fprintf(stderr, "error: must provide 2 (and only 2) input files, you "
            "provided %d\n\n", argc - optind);
    pe_print_usage(stderr);
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
      pe_print_usage(stderr);
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

void pe_print_usage(FILE *outstream)
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
"                                filter reported results;\n"
"    -t|--maxtrans: INT          Maximum transcripts allowed per locus; use 0\n"
"                                to disable limit; default is 32\n\n");
}

void pe_set_option_defaults(ParsEvalOptions *options)
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
  options->filters = gt_array_new( sizeof(AgnLocusFilter) );
  options->verbose = false;
  options->max_transcripts = 32;
}

void pe_summary_header(ParsEvalOptions *options, FILE *outstream,
                       char *start_time, int argc, char **argv)
{
  fprintf(outstream,
          "============================================================\n"
          "========== ParsEval Summary\n"
          "============================================================\n\n");

  fprintf(outstream, "Started:                %s\n", start_time);

  if(options->refrlabel != NULL)
    fprintf(outstream, "Reference annotations:  %s\n", options->refrlabel);
  else
    fprintf(outstream, "Reference annotations:  %s\n", options->refrfile);
  if(options->predlabel != NULL)
    fprintf(outstream, "Prediction annotations: %s\n", options->predlabel);
  else
    fprintf(outstream, "Prediction annotations: %s\n", options->predfile);
  fprintf(outstream, "Executing command:      ");

  int x;
  for(x = 0; x < argc; x++)
  {
    fprintf(outstream, "%s ", argv[x]);
  }
  fprintf(outstream, "\n\n");
}
