/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include "pe_options.h"

void pe_free_option_memory(ParsEvalOptions *options)
{
  fclose(options->outfile);
  gt_array_delete(options->filters);
}

int pe_parse_options(int argc, char **argv, ParsEvalOptions *options,
                     GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "a:df:ghkl:o:pr:st:Vvwx:y:";
  const struct option parseval_options[] =
  {
    { "datashare",  required_argument, NULL, 'a' },
    { "debug",      no_argument,       NULL, 'd' },
    { "outformat",  required_argument, NULL, 'f' },
    { "printgff3",  no_argument,       NULL, 'g' },
    { "help",       no_argument,       NULL, 'h' },
    { "makefilter", no_argument,       NULL, 'k' },
    { "delta",      required_argument, NULL, 'l' },
    { "outfile",    required_argument, NULL, 'o' },
    { "nopng",      no_argument,       NULL, 'p' },
    { "filterfile", required_argument, NULL, 'r' },
    { "summary",    no_argument,       NULL, 's' },
    { "maxtrans",   required_argument, NULL, 't' },
    { "verbose",    no_argument,       NULL, 'V' },
    { "version",    no_argument,       NULL, 'v' },
    { "overwrite",  no_argument,       NULL, 'w' },
    { "refrlabel",  required_argument, NULL, 'x' },
    { "predlabel",  required_argument, NULL, 'y' },
    { NULL,         no_argument,       NULL,  0  },
  };

  for(opt  = getopt_long(argc, argv, optstr, parseval_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv, optstr, parseval_options, &optindex))
  {
    if(opt == 'a')
    {
      options->data_path = optarg;
    }
    else if(opt == 'd')
    {
      options->debug = true;
    }
    else if(opt == 'f')
    {
      if      (strcmp(optarg, "csv")  == 0) options->outfmt = CSVMODE;
      else if (strcmp(optarg, "text") == 0) options->outfmt = TEXTMODE;
      else if (strcmp(optarg, "html") == 0) options->outfmt = HTMLMODE;
      else
      {
        fprintf(stderr, "error: unknown value '%s' for '-f|--outformat' "
                "option\n\n", optarg);
        pe_print_usage(stderr);
        exit(1);
      }
    }
    else if(opt == 'g')
    {
      options->gff3 = false;
    }
    else if(opt == 'h')
    {
      pe_print_usage(stdout);
      exit(0);
    }
    else if(opt == 'k')
    {
      options->makefilter = true;
    }
    else if(opt == 'l')
    {
      if(sscanf(optarg, "%ld", &options->delta) == EOF)
      {
        fprintf(stderr, "error: could not convert delta '%s' to an integer\n",
                optarg);
        exit(1);
      }
    }
    else if(opt == 'o')
    {
      options->outfilename = optarg;
    }
    else if(opt == 'p')
    {
      options->graphics = false;
    }
    else if(opt == 'r')
    {
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
    }
    else if(opt == 's')
    {
      options->summary_only = true;
    }
    else if(opt == 't')
    {
      if(sscanf(optarg, "%d", &options->max_transcripts) == EOF)
      {
        fprintf(stderr, "error: could not convert maxtrans '%s' to an "
                        "integer\n", optarg);
        exit(1);
      }
    }
    else if(opt == 'V')
    {
      options->verbose = true;
    }
    else if(opt == 'v')
    {
      agn_print_version("ParsEval", stdout);
      exit(0);
    }
    else if(opt == 'w')
    {
      options->overwrite = true;
    }
    else if(opt == 'x')
    {
      options->refrlabel = optarg;
    }
    else if(opt == 'y')
    {
      options->predlabel = optarg;
    }
  }
  
#ifdef WITHOUT_CAIRO
  if(options->graphics)
  {
    fputs("error: AEGeAn was compiled without graphics support. Please "
          "recompile to enable this feature.\n", stderr);
    exit(1);
  }
#endif

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
    pe_print_usage(stderr);
    fprintf(stderr, "error: must provide 2 (and only 2) input files, you "
            "provided %d\n\n", argc - optind);
    exit(1);
  }

  if(options->outfmt == HTMLMODE && options->summary_only)
  {
    fprintf(stderr, "warning: summary-only mode requires text output format; "
            "ignoring\n");
    options->summary_only = false;
  }

  if(options->outfilename)
  {
    if(options->outfmt == HTMLMODE)
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

  if(options->outfmt == HTMLMODE)
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

    if(options->graphics && options->summary_only)
    {
      fprintf(stderr, "warning: cannot print PNG graphics in summary only "
              "mode; ignoring\n");
      options->graphics = false;
    }
  }
  else
  {
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
  if(options->outfmt != HTMLMODE && options->graphics)
    options->graphics = false;

  if(options->graphics)
  {
    sprintf(options->pngdata.filename_template, "%s/%%s/%%s_%%lu-%%lu.png",
            options->outfilename);
    sprintf(options->pngdata.stylefile, "%s/pe.style", options->data_path);
    options->pngdata.refrfile  = options->refrfile;
    options->pngdata.predfile  = options->predfile;
    options->pngdata.refrlabel = options->refrlabel;
    options->pngdata.predlabel = options->predlabel;
  }
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
"    -l|--delta: INT             Extend gene loci by this many nucleotides;\n"
"                                default is 0\n"
"    -V|--verbose:               Print verbose warning messages\n"
"    -v|--version:               Print version number and exit\n\n"
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
"    -p|--nopng:                 In HTML output mode, skip generation of PNG\n"
"                                graphics for each gene locus\n"
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
  options->outfilename = NULL;
  options->gff3 = true;
  options->summary_only = false;
  options->graphics = true;
#ifdef WITHOUT_CAIRO
  options->graphics = false;
#endif
  options->refrfile = NULL;
  options->predfile = NULL;
  options->refrlabel = NULL;
  options->predlabel = NULL;
  options->outfmt = TEXTMODE;
  options->overwrite = false;
  options->data_path = AGN_DATA_PATH;
  options->makefilter = false;
  options->filters = gt_array_new( sizeof(AgnLocusFilter) );
  options->verbose = false;
  options->max_transcripts = 32;
  options->delta = 0;
}
