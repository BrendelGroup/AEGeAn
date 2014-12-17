/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef AEGEAN_REPO_VISITOR
#define AEGEAN_REPO_VISITOR

#include "extended/node_stream_api.h"

/**
 * @class AgnRepoStream
 *
 * Node visitor that handles populating GFF3 files in a git annotation
 * repository.
 */
typedef struct AgnRepoVisitor AgnRepoVisitor;

/**
 * @class Create a repo node stream.
 */
GtNodeStream *agn_repo_stream_new(GtNodeStream *instream, const char *path,
                                  GtError *error);

/**
 * @class Constructor: allocate memory for the node visitor.
 */
GtNodeVisitor *agn_repo_visitor_new(const char *path, GtError *error);

#endif
