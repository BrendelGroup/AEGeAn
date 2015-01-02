/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef GENEANNOLOGY_COMMANDS
#define GENEANNOLOGY_COMMANDS

// Remove all existing annotations from the repository and start over fresh
int ga_clean(int argc, char * const *argv);

// Save a snapshot of the repository
int ga_commit(int argc, char * const *argv);

// Initialize a new repository
int ga_init(int argc, char * const *argv);

#endif
