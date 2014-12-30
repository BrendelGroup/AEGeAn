/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>
#include "AgnRepoVisitor.h"
#include "AgnGeneLocusMapping.h"
#include "AgnSeqGeneMapping.h"

#define repo_visitor_cast(GV)\
        gt_node_visitor_cast(repo_visitor_class(), GV)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnRepoVisitor
{
  const GtNodeVisitor parent_instance;
  char *repopath;
  char cwd[PATH_MAX];
  AgnGeneLocusMapping *glmap;
  AgnSeqGeneMapping *sgmap;
};

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the GtNodeVisitor interface.
 */
static const GtNodeVisitorClass *repo_visitor_class();

/**
 * @function Destructor: free instance memory.
 */
static void repo_visitor_free(GtNodeVisitor *nv);

/**
 * @function Process feature nodes in the node stream: verify that each feature
 * node is a locus feature, and write it to disk.
 */
static int repo_visitor_fn_handler(GtNodeVisitor *nv, GtFeatureNode *fn,
                                   GtError *error);

/**
 * @function Create a new repo and initialize it.
 */
static int repo_visitor_init(AgnRepoVisitor *rv, GtError *error);

/**
 * @function Determine the filename of the locus using the following format:
 * `repopath/seqid/seqid_start-end.gff3`.
 */
static GtStr *repo_visitor_locus_filename(AgnRepoVisitor *rv, AgnLocus *locus);

/**
 * @function Process region nodes in the node stream: delete features from the
 * given sequence.
 */
static int repo_visitor_rn_handler(GtNodeVisitor *nv, GtRegionNode *rn,
                                   GtError *error);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream *agn_repo_stream_new(GtNodeStream *instream, const char *path,
                                  GtError *error)
{
  agn_assert(instream && path && error);
  GtNodeVisitor *nv = agn_repo_visitor_new(path, error);
  if(nv == NULL)
    return NULL;
  return gt_visitor_stream_new(instream, nv);
}

GtNodeStream *agn_repo_stream_open_clean(GtNodeStream *instream,
                                         const char *path, GtError *error)
{
  agn_assert(instream && path && error);
  GtNodeVisitor *nv = agn_repo_visitor_open_clean(path, error);
  if(nv == NULL)
    return NULL;
  return gt_visitor_stream_new(instream, nv);
}

GtNodeVisitor *agn_repo_visitor_new(const char *path, GtError *error)
{
  struct stat buffer;
  agn_assert(path && error);
  if(stat(path, &buffer) == 0)
  {
    gt_error_set(error, "path '%s' exists, will not create repository", path);
    return NULL;
  }

  GtNodeVisitor *nv = gt_node_visitor_create(repo_visitor_class());
  AgnRepoVisitor *rv = repo_visitor_cast(nv);
  rv->repopath = gt_cstr_dup(path);
  if(getcwd(rv->cwd, PATH_MAX) == NULL)
  {
    gt_error_set(error, "error getting cwd");
    return NULL;
  }
  if(repo_visitor_init(rv, error) != 0)
    return NULL;

  GtStr *gl_filename = gt_str_new_cstr(path);
  gt_str_append_cstr(gl_filename, "/gene-locus.map");
  GtStr *sg_filename = gt_str_new_cstr(path);
  gt_str_append_cstr(sg_filename, "/sequence-gene.map");
  rv->glmap = agn_gene_locus_mapping_new(gt_str_get(gl_filename));
  rv->sgmap = agn_seq_gene_mapping_new(gt_str_get(sg_filename));
  gt_str_delete(gl_filename);
  gt_str_delete(sg_filename);

  return nv;
}

GtNodeVisitor *agn_repo_visitor_open_clean(const char *path, GtError *error)
{
  struct stat buffer;
  agn_assert(path && error);
  if(stat(path, &buffer))
  {
    gt_error_set(error, "repo '%s' does not exist", path);
    return NULL;
  }

  GtNodeVisitor *nv = gt_node_visitor_create(repo_visitor_class());
  AgnRepoVisitor *rv = repo_visitor_cast(nv);
  rv->repopath = gt_cstr_dup(path);
  if(getcwd(rv->cwd, PATH_MAX) == NULL)
  {
    gt_error_set(error, "error getting cwd");
    return NULL;
  }

  chdir(rv->repopath);
  int result = system("git rm -r *");
  if(result != 0)
  {
    gt_error_set(error, "error removing existing annotations");
    return NULL;
  }
  chdir(rv->cwd);

  GtStr *gl_filename = gt_str_new_cstr(path);
  gt_str_append_cstr(gl_filename, "/gene-locus.map");
  GtStr *sg_filename = gt_str_new_cstr(path);
  gt_str_append_cstr(sg_filename, "/sequence-gene.map");
  rv->glmap = agn_gene_locus_mapping_new(gt_str_get(gl_filename));
  rv->sgmap = agn_seq_gene_mapping_new(gt_str_get(sg_filename));
  gt_str_delete(gl_filename);
  gt_str_delete(sg_filename);

  return nv;
}

static const GtNodeVisitorClass *repo_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnRepoVisitor), repo_visitor_free,
                                    NULL, repo_visitor_fn_handler,
                                    repo_visitor_rn_handler, NULL, NULL);
  }
  return nvc;
}

static int repo_visitor_fn_handler(GtNodeVisitor *nv, GtFeatureNode *fn,
                                   GtError *error)
{
  agn_assert(nv);
  AgnRepoVisitor *rv = repo_visitor_cast(nv);
  agn_assert(gt_feature_node_has_type(fn, "locus"));
  AgnLocus *locus = (AgnLocus *)fn;
  GtStr *filename = repo_visitor_locus_filename(rv, locus);
  GtFile *file = gt_file_new(gt_str_get(filename), "w", error);
  if(file == NULL)
    return -1;
  GtNodeVisitor *gff3 = gt_gff3_visitor_new(file);
  gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor *)gff3);
  if(gt_genome_node_accept(locus, gff3, error) == -1)
    return -1;

  agn_gene_locus_mapping_add(rv->glmap, locus);
  agn_seq_gene_mapping_add(rv->sgmap, locus);
  gt_str_delete(filename);
  gt_file_delete(file);
  gt_node_visitor_delete(gff3);
  return 0;
}

static void repo_visitor_free(GtNodeVisitor *nv)
{
  agn_assert(nv);
  AgnRepoVisitor *rv = repo_visitor_cast(nv);
  agn_gene_locus_mapping_close(rv->glmap);
  agn_seq_gene_mapping_close(rv->sgmap);
  chdir(rv->repopath);
  if(system("git add .") != 0)
    fprintf(stderr, "warning: file staging failed\n");
  chdir(rv->cwd);
  gt_free(rv->repopath);
}

static int repo_visitor_init(AgnRepoVisitor *rv, GtError *error)
{
  int status = mkdir(rv->repopath, 0755);
  if(status != 0)
  {
    gt_error_set(error, "could not create directory '%s', error code %d",
                 rv->repopath, status);
    return status;
  }

#ifdef DEBUG
  fprintf(stderr, "[AgnRepoVisitor::repo_visitor_init] chdir: %s -> %s\n",
          rv->cwd, rv->repopath);
#endif
  chdir(rv->repopath);

  status = system("git init");
  if(status != 0)
  {
    gt_error_set(error, "error initializing git repository");
    return -1;
  }
#ifdef DEBUG
  fprintf(stderr, "[AgnRepoVisitor::repo_visitor_init] chdir: %s -> %s\n",
          rv->repopath, rv->cwd);
#endif
  chdir(rv->cwd);
  return 0;
}

static GtStr *repo_visitor_locus_filename(AgnRepoVisitor *rv, AgnLocus *locus)
{
  agn_assert(rv && locus);
  char filename[PATH_MAX];
  const char *seqid = gt_str_get(gt_genome_node_get_seqid(locus));
  GtRange range = gt_genome_node_get_range(locus);

  sprintf(filename, "%s/%s/%s_%lu-%lu.gff3", rv->repopath, seqid, seqid,
          range.start, range.end);
  return gt_str_new_cstr(filename);
}

static int repo_visitor_rn_handler(GtNodeVisitor *nv, GtRegionNode *rn,
                                   GtError *error)
{
  agn_assert(nv && rn && error);
  AgnRepoVisitor *rv = repo_visitor_cast(nv);
  const char *seqid = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode *)rn));

  chdir(rv->repopath);
  struct stat buffer;
  if(stat(seqid, &buffer) == 0)
  {
    GtStr *command = gt_str_new_cstr("git rm -r ");
    gt_str_append_cstr(command, seqid);
    if(system(gt_str_get(command)) != 0)
    {
      gt_error_set(error, "error removing seq directory '%s'", seqid);
      return -1;
    }
    gt_str_delete(command);

    GtStrArray *geneids = agn_seq_gene_mapping_unmap_seqid(rv->sgmap, seqid);
    agn_assert(geneids);
    GtUword i;
    for(i = 0; i < gt_str_array_size(geneids); i++)
    {
      const char *geneid = gt_str_array_get(geneids, i);
      GtStr *locuspos = agn_gene_locus_mapping_unmap_gene(rv->glmap, geneid);
      gt_str_delete(locuspos);
    }
    gt_str_array_delete(geneids);
  }

  if(mkdir(seqid, 0755) != 0)
  {
    gt_error_set(error, "error creating seq directory '%s'", seqid);
    return -1;
  }
  chdir(rv->cwd);

  return 0;
}
