#include <getopt.h>
#include <omp.h>
#include <string.h>
#include "PeOptions.h"

int pe_parse_options(int argc, char * const argv[], PeOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "a:c:df:ghkl:mn:o:pr:t:svwx:y:";
  const struct option parseval_options[] =
  {
    { "datashare", required_argument, NULL, 'a' },
    { "complimit", required_argument, NULL, 'c' },
    { "debug", no_argument, NULL, 'd' },
    { "outformat", required_argument, NULL, 'f' },
    { "printgff3", required_argument, NULL, 'g' },
    { "help", no_argument, NULL, 'h' },
    { "makefilter", no_argument, NULL, 'k' },
    { "locusgff3", required_argument, NULL, 'l' },
    { "vectors", no_argument, NULL, 'm' },
    { "numprocs", required_argument, NULL, 'n' },
    { "outfile", required_argument, NULL, 'o' },
    { "png", no_argument, NULL, 'p' },
    { "filterfile", required_argument, NULL, 'r' },
    { "summary", no_argument, NULL, 's' },
    { "maxtrans", no_argument, NULL, 't' },
    { "verbose", no_argument, NULL, 'v' },
    { "overwrite", no_argument, NULL, 'w' },
    { "refrlabel", required_argument, NULL, 'x' },
    { "predlabel", required_argument, NULL, 'y' },
    { NULL, no_argument, NULL, 0 },
  };

  for( opt = getopt_long(argc, argv + 0, optstr, parseval_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, parseval_options, &optindex) )
  {
    switch(opt)
    {
      case 'a':
        options->data_path = optarg;
        break;

      case 'c':
        if( sscanf(optarg, "%d", &options->complimit) == EOF )
        {
          fprintf(stderr, "error: could not convert comparison limit '%s' to an integer", optarg);
          exit(1);
        }
        break;

      case 'd':
        options->debug = true;
        break;

      case 'f':
        if( strcmp(optarg,  "csv") != 0 &&
            strcmp(optarg, "text") != 0 &&
            strcmp(optarg, "html") != 0 )
        {
            fprintf(stderr, "error: unknown value '%s' for '-f|--outformat' option\n\n", optarg);
            pe_print_usage();
            exit(1);
        }
        options->outfmt = optarg;
        if(strcmp(options->outfmt, "html") == 0)
          options->html = true;
        break;

      case 'g':
        options->gff3 = true;
        break;

      case 'h':
        pe_print_usage();
        exit(0);
        break;

      case 'k':
        if(true)
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
        break;

      case 'l':
        options->locusfilename = optarg;
        options->locusgff3 = 1;
        break;

      case 'm':
        options->vectors = true;
        break;

      case 'n':
        if( sscanf(optarg, "%d", &options->numprocs) == EOF )
        {
          fprintf(stderr, "error: could not convert number of processors '%s' to an integer", optarg);
          exit(1);
        }
        break;

      case 'o':
        options->outfilename = optarg;
        break;

      case 'p':
        options->locus_graphics = true;
#ifdef WITHOUT_CAIRO
        fputs( "error: ParsEval was compiled without graphics support. Please recompile to enable "
               "this feature.\n", stderr);
        exit(1);
#endif
        break;

      case 'r':
        options->usefilter = true;
        options->filterfile = optarg;
        if(options->usefilter)
        {
          if(options->debug)
            fprintf(stderr, "debug: opening filter file '%s'\n", options->filterfile);
          FILE *filterfile = fopen(options->filterfile, "r");
          if(!filterfile)
          {
            fprintf(stderr, "error: could not open filter file '%s'\n", options->filterfile);
            exit(1);
          }
          agn_compare_filters_parse(&options->filters, filterfile);
          if(options->debug)
            fprintf(stderr, "debug: closing filter file\n");
          fclose(filterfile);
        }
        break;

      case 's':
        options->summary_only = true;
        break;

      case 't':
        if( sscanf(optarg, "%d", &options->trans_per_locus) == EOF )
        {
          fprintf(stderr, "error: could not convert transcript limit '%s' to an integer", optarg);
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

  if(argc - optind != 2)
  {
    fprintf( stderr, "error: must provide 2 (and only 2) input files, you provided %d\n\n",
             argc - optind );
    pe_print_usage();
    exit(1);
  }

  if(strcmp(options->outfilename, "STDOUT") != 0)
  {
    if(options->html)
    {
      char dircmd[256];
      sprintf(dircmd, "test -d %s", options->outfilename);
      if(system(dircmd) == 0)
      {
        if(options->overwrite)
        {
          char rmcmd[256];
          sprintf(rmcmd, "rm -r %s", options->outfilename);
          if(system(rmcmd) != 0)
          {
            fprintf( stderr, "error: could not overwrite output directory '%s'",
                     options->outfilename );
            exit(1);
          }
        }
        else
        {
          fprintf( stderr, "error: outfile '%s' exists; use '-w' to force overwrite\n",
                   options->outfilename );
          exit(1);
        }
      }
    }
    else
    {
      char filecmd[256];
      sprintf(filecmd, "test -f %s", options->outfilename);
      if(system(filecmd) == 0 && !options->overwrite)
      {
          fprintf( stderr, "error: outfile '%s' exists; use '-w' to force overwrite\n",
                   options->outfilename );
          exit(1);
      }
    }
  }

  if(options->html)
  {
    if(strcmp(options->outfilename, "STDOUT") == 0)
    {
      fputs( "error: will not print results to terminal in HTML mode; must provide outfile\n\n",
             stderr );
      pe_print_usage();
      exit(1);
    }
    char dircmd[256];
    sprintf(dircmd, "mkdir %s", options->outfilename);
    if(system(dircmd) != 0)
    {
      fprintf(stderr, "error: cannot open output directory '%s'\n", options->outfilename);
      exit(1);
    }
    char outname[256];
    sprintf(outname, "%s/index.html", options->outfilename);
    options->outfile = fopen(outname, "w");
    if(!options->outfile)
    {
      fprintf(stderr, "error: could not open output file '%s'\n", outname);
      exit(1);
    }

    char copy_cmd[1024];
    sprintf(copy_cmd, "cp -r %s/* %s", options->data_path, options->outfilename);
    if(options->debug)
      fprintf(stderr, "debug: copying shared data: '%s'\n", copy_cmd);
    if(system(copy_cmd) != 0)
    {
      fprintf(stderr, "error: could not copy data files '%s'\n", copy_cmd);
      exit(1);
    }

    if(options->summary_only && options->locus_graphics)
    {
      fprintf(stderr, "warning: cannot print PNG graphics in summary only mode; ignoring\n");
      options->locus_graphics = false;
    }
  }
  else
  {
    if(options->locus_graphics)
    {
      fputs("warning: will only generate PNG graphics when outformat='html'; ignoring\n\n", stderr);
      options->locus_graphics = false;
    }
    if(strcmp(options->outfilename, "STDOUT") != 0)
    {
      options->outfile = fopen(options->outfilename, "w");
      if(options->outfile == NULL)
      {
        fprintf(stderr, "error: cannot open output file '%s'\n", options->outfilename);
        exit(1);
      }
    }
  }

  return optind;
}

void pe_print_usage()
{
  fprintf( stderr,
           "Usage: parseval [options] reference prediction\n"
           "  Options:\n"
           "    -a|--datashare: STRING      Location from which to copy shared data for\n"
           "                                HTML output (if `make install' has not yet\n"
           "                                been run)\n"
           "    -c|--complimit: INT         Maximum number of comparisons per locus; set\n"
           "                                to 0 for no limit (default=512)\n"
           "    -d|--debug:                 Print debugging messages\n"
           "    -f|--outformat: STRING      Indicate desired output format; possible\n"
           "                                options: 'csv', 'text', or 'html'\n"
           "                                (default='text'); in 'text' or 'csv' mode,\n"
           "                                will create a single file; in 'html' mode,\n"
           "                                will create a directory\n"
           "    -g|--printgff3:             Include GFF3 output corresponding to each\n"
           "                                comparison\n"
           "    -h|--help:                  Print help message and exit\n"
           "    -k|--makefilter             Create a default configuration file for\n"
           "                                filtering reported results\n"
           "    -l|--locusgff3: STRING      Print out locus coordinates to the given GFF3\n"
           "                                file (ideal for genome browser integration)\n"
           "    -m|--vectors:               Print model vectors in output file\n"
           "    -n|--numprocs: INT          Number of processors to use (default=1)\n"
           "    -o|--outfile: FILENAME      File/directory to which output will be written;\n"
           "                                default is the terminal (STDOUT)\n"
           "    -p|--png:                   Generate individual PNG graphics for each gene\n"
           "                                locus\n"
           "    -r|--filterfile: STRING     Use the indicated configuration file to filter\n"
           "                                reported results;\n"
           "    -s|--summary:               Only print summary statistics, do not print\n"
           "                                individual comparisons\n"
           "    -t|--maxtrans: INT          The maximum number of transcripts that can be\n"
           "                                annotated at a given gene locus; set to 0 for\n"
           "                                no limit (default=32)\n"
           "    -v|--verbose:               Print verbose warning messages\n"
           "    -w|--overwrite:             Force overwrite of any existing output files\n"
           "    -x|--refrlabel: STRING      Optional label for reference annotations\n"
           "    -y|--predlabel: STRING      Optional label for prediction annotations\n" );
}

void pe_set_option_defaults(PeOptions *options)
{
  options->debug = false;
  options->outfile = stdout;
  options->outfilename = "STDOUT";
  options->gff3 = false;
  options->verbose = false;
  options->complimit = 512;
  options->summary_only = false;
  options->vectors = false;
  options->locus_graphics = false;
  options->outfmt = "text";
  options->overwrite = false;
  options->html = false;
  options->data_path = AGN_DATA_PATH;
  options->usefilter = false;
  options->filterfile = "";
  agn_compare_filters_init(&options->filters);
  options->numprocs = 1;
  options->trans_per_locus = 32;
  options->locusgff3 = 0;
  options->locusfilename = "";
  options->refrlabel = "";
  options->predlabel = "";
}
