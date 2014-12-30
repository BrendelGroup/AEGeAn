/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "genometools.h"
#include "aegean.h"
#include "ga_commands.h"

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Remind user that certain commands can be delegated entirely to git.
 */
static int ga_delegate(int argc, const char **argv);

/**
 * @function Indicate that the command is not yet implemented.
 */
static int ga_notyetimplemented(int argc, const char **argv);

/**
 * @function Print a usage statement.
 */
static void ga_print_usage(FILE *outstream);

/**
 * @function Indicate that the git command is not supported in GeneAnnoLogy.
 */
static int ga_unsupported(int argc, const char **argv);


//------------------------------------------------------------------------------
// Function implementations
//------------------------------------------------------------------------------

static int ga_delegate(int argc, const char **argv)
{
  printf("command '%s' not supported, should be delegated to git\n\n", argv[0]);
  printf(
"A GeneAnnoLogy repository is in every sense a git repository. While this\n"
"does not necessarily give the user free license to directly modify data in\n"
"the repository, it *does* mean that you can use certain git commands (such\n"
"as '%s') to manage your GeneAnnoLogy repository. Enter the command\n\n"
"  git %s --help\n\n"
"for more details.\n",
    argv[0], argv[0]);
  return 0;
}

static int ga_notyetimplemented(int argc, const char **argv)
{
  fprintf(stderr, "[GeneAnnoLogy] error: the command '%s' has not yet been "
          "implemented.\n", argv[0]);
  ga_print_usage(stderr);
  return 1;
}

static void ga_print_usage(FILE *outstream)
{
  fputs("Usage: geneannology <command> [options] repo [<args>]\n"
    "  Commands:\n"
    "    cat\n"
    "    commit\n"
    "    delete\n"
    "    diff\n"
    "    init\n"
    "    status\n"
    "    update\n",
    outstream);
}

static int ga_unsupported(int argc, const char **argv)
{
  fprintf(stderr,
"[GeneAnnoLogy] error: the git command '%s' is unsupported in GeneAnnoLogy\n\n"
"The git commands 'add', 'mv', and 'rm' have been replaced with the\n"
"GeneAnnoLogy commands 'update' and 'delete'. While some git commands can be\n"
"used to manage entire repositories (clone, pull, etc), please only use\n"
"GeneAnnoLogy commands for modifying repository data.\n", argv[0]);
  ga_print_usage(stderr);
  return 1;
}


//------------------------------------------------------------------------------
// Main method
//------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
  // Initialize the GenomeTools library
  gt_lib_init();

  // Commands associated with git and/or GeneAnnoLogy
  GtHashmap *commands = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);

  //   - commands defined in git but not supported in GeneAnnoLogy
  gt_hashmap_add(commands, "add",  ga_unsupported);
  gt_hashmap_add(commands, "mv",   ga_unsupported);
  gt_hashmap_add(commands, "rm",   ga_unsupported);

  //   - commands that should be delegated entirely to git
  gt_hashmap_add(commands, "bisect",   ga_delegate);
  gt_hashmap_add(commands, "branch",   ga_delegate);
  gt_hashmap_add(commands, "checkout", ga_delegate);
  gt_hashmap_add(commands, "clone",    ga_delegate);
  gt_hashmap_add(commands, "fetch",    ga_delegate);
  gt_hashmap_add(commands, "grep",     ga_delegate);
  gt_hashmap_add(commands, "log",      ga_delegate);
  gt_hashmap_add(commands, "merge",    ga_delegate);
  gt_hashmap_add(commands, "pull",     ga_delegate);
  gt_hashmap_add(commands, "push",     ga_delegate);
  gt_hashmap_add(commands, "rebase",   ga_delegate);
  gt_hashmap_add(commands, "reset",    ga_delegate);
  gt_hashmap_add(commands, "tag",      ga_delegate);

  //   - commands defined by git but overridden by GeneAnnoLogy
  gt_hashmap_add(commands, "commit", ga_notyetimplemented);
  gt_hashmap_add(commands, "diff",   ga_notyetimplemented);
  gt_hashmap_add(commands, "init",   ga_init);
  gt_hashmap_add(commands, "show",   ga_notyetimplemented);
  gt_hashmap_add(commands, "status", ga_notyetimplemented);

  //   - commands unique to GeneAnnoLogy
  gt_hashmap_add(commands, "clean",  ga_clean);
  gt_hashmap_add(commands, "delete", ga_notyetimplemented);
  gt_hashmap_add(commands, "update", ga_notyetimplemented);


  // Parse command name, run the associated command
  if(argc < 2)
  {
    fprintf(stderr, "[GeneAnnoLogy] error: no command specified\n");
    ga_print_usage(stderr);
    return 1;
  }

  if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0 ||
     strcmp(argv[1], "--help") == 0)
  {
    ga_print_usage(stdout);
    return 0;
  }

  const char *command = argv[1];
  int (*action)(int, const char **) = gt_hashmap_get(commands, command);
  if(!action)
  {
    fprintf(stderr, "[GeneAnnoLogy] error: unknown command '%s'\n", command);
    ga_print_usage(stderr);
    return 1;
  }

  int code = action(argc - 2, argv + 2);

  // Free memory and exit
  gt_hashmap_delete(commands);
  gt_lib_clean();
  return code;
}
