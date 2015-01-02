/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include <unistd.h>
#include "ga_commands.h"
#include "genometools.h"
#include "aegean.h"

typedef struct
{
  const char *author;
  bool dryrun;
  bool fix;
  const char *message;
  const char *repo;
} GaCommitOptions;

static void ga_commit_print_usage(FILE *outstream)
{
  fputs("\n"
"[GeneAnnoLogy::commit] save a snapshot of the annotation repository\n\n"
"Usage: geneannology commit [options] repo\n"
"  Options:\n"
"    -a|--author: STRING     override author for commit\n"
"    -d|--dry-run            do not commit, but show what would be committed\n"
"    -f|--fix                overwrite previous commit (a la 'git --amend')\n"
"    -h|--help               print this help message and exit\n"
"    -m|--message: STRING    commit message\n"
"    -v|--version            print version number and exit\n\n", outstream);
}

static void ga_commit_parse_options(int argc, char * const *argv,
                                  GaCommitOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hs:v";
  const struct option commit_options[] =
  {
    { "author",  required_argument, NULL, 'a' },
    { "dry-run", no_argument,       NULL, 'd' },
    { "fix",     no_argument,       NULL, 'f' },
    { "help",    no_argument,       NULL, 'h' },
    { "message", required_argument, NULL, 'm' },
  };

  options->author = NULL;
  options->dryrun = false;
  options->fix = false;
  options->message = NULL;
  options->repo = NULL;

  for(opt  = getopt_long(argc, argv, optstr, commit_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv, optstr, commit_options, &optindex))
  {
    if(opt == 'a')
    {
      options->author = optarg;
    }
    else if(opt == 'd')
    {
      options->dryrun = true;
    }
    else if(opt == 'f')
    {
      options->fix = true;
    }
    else if(opt == 'h')
    {
      ga_commit_print_usage(stdout);
      exit(0);
    }
    else if(opt == 'm')
    {
      options->message = optarg;
    }
    else if(opt == 'v')
    {
      agn_print_version("GeneAnnoLogy::commit", stdout);
      exit(0);
    }
    else
    {
      ga_commit_print_usage(stderr);
      fprintf(stderr, "error: unknown option '%c'", opt);
      exit(1);
    }
  }

  int numargs = argc - optind;
  if(numargs == 0)
    return;

  options->repo = argv[optind];
}

int ga_commit(int argc, char * const *argv)
{
  agn_assert(argc > 0);
  GaCommitOptions options;
  ga_commit_parse_options(argc, argv, &options);
  if(options.repo == NULL)
  {
    ga_commit_print_usage(stderr);
    return 0;
  }
  
  char cwd[PATH_MAX];
  if(getcwd(cwd, PATH_MAX) == NULL)
  {
    fprintf(stderr, "[GeneAnnoLogy::commit] error: could not retrieve current "
            "working directory\n");
    exit(1);
  }
  chdir(options.repo);

  GtStr *command = gt_str_new_cstr("git commit");
  if(options.author)
  {
    gt_str_append_cstr(command, " --author ");
    gt_str_append_cstr(command, options.author);
  }
  if(options.dryrun)
    gt_str_append_cstr(command, " --dry-run");
  if(options.fix)
    gt_str_append_cstr(command, " --amend");
  if(options.message)
  {
    gt_str_append_cstr(command, " -m \"");
    gt_str_append_cstr(command, options.message);
    gt_str_append_char(command, '"');
  }
  if(system(gt_str_get(command)) != 0)
  {
    fprintf(stderr, "[GeneAnnoLogy::commit] error: could not commit "
            "repository\n");
    exit(1);
  }

  chdir(cwd);
  return 0;
}
