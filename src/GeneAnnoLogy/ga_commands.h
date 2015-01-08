/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef GENEANNOLOGY_COMMANDS
#define GENEANNOLOGY_COMMANDS

#include "genometools.h"
#include "aegean.h"

// Remove all existing annotations from the repository and start over fresh
int ga_clean(int argc, char * const *argv);

// Save a snapshot of the repository
int ga_commit(int argc, char * const *argv);

// Initialize a new repository
int ga_init(int argc, char * const *argv);

// Merge two sets of annotations into one non-redundant set
int ga_merge(int argc, char * const *argv);

// Add new annotations to a repository. Keep all old and new annotations: do
// not attempt to resolve redundant or overlapping annotations.
int ga_union(int argc, char * const *argv);


//------------------------------------------------------------------------------
// Utility commands
//------------------------------------------------------------------------------

/**
 * @function Create a string array containing filenames for any loci associated
 * with the sequences specified by `seqids`. Silently ignore any sequences in
 * `seqids` not present in the repository.
 */
GtStrArray *ga_repo_get_filenames_for_seqids(const char *repo,
                                             GtStrArray *seqids,GtError *error);

#endif
