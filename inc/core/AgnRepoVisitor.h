/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

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
 * @class Create a stream that will initialize a new repo.
 */
GtNodeStream *agn_repo_stream_new(GtNodeStream *instream, const char *path,
                                  GtError *error);

/**
 * @class Open a stream to write to an existing repo.
 */
GtNodeStream *agn_repo_stream_open(GtNodeStream *instream, const char *path,
                                   GtError *error);

/**
 * @class Create a stream that will remove existing annotations from the repo
 * before filling it in with new annotations.
 */
GtNodeStream *agn_repo_stream_open_clean(GtNodeStream *instream,
                                         const char *path, GtError *error);

/**
 * @class Constructor: create a visitor that will initialize a new repo.
 */
GtNodeVisitor *agn_repo_visitor_new(const char *path, GtError *error);

/**
 * @class Constructor: create a visitor to write to an existing repo.
 */
GtNodeVisitor *agn_repo_visitor_open(const char *path, GtError *error);

/**
 * @class Constructor: create a visitor that will remove existing annotations
 * from the repo before filling it in with new annotations.
 */
GtNodeVisitor *agn_repo_visitor_open_clean(const char *path, GtError *error);

#endif
