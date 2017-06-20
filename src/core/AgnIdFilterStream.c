/**

Copyright (c) 2017, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "core/queue_api.h"
#include "extended/array_in_stream_api.h"
#include "extended/array_out_stream_api.h"
#include "extended/feature_node_iterator_api.h"
#include "AgnIdFilterStream.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnIdFilterStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtHashmap *ids2keep;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

#define id_filter_stream_cast(GS)\
        gt_node_stream_cast(id_filter_stream_class(), GS)

/**
 * @function Implements the GtNodeStream interface for this class.
 */
static const GtNodeStreamClass* id_filter_stream_class(void);

/**
 * @function Class destructor.
 */
static void id_filter_stream_free(GtNodeStream *ns);

/**
 * @function Pulls nodes from the input stream and feeds them to the output
 * stream if they are in the list of IDs to keep.
 */
static int id_filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                GtError *error);

/**
 * @function Generate data for unit testing.
 */
static void id_filter_stream_test_data(GtQueue *queue);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_id_filter_stream_new(GtNodeStream *in_stream,
                                       GtHashmap *ids2keep)
{
  GtNodeStream *ns;
  AgnIdFilterStream *stream;
  agn_assert(in_stream && ids2keep);
  ns = gt_node_stream_create(id_filter_stream_class(), false);
  stream = id_filter_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);
  stream->ids2keep = gt_hashmap_ref(ids2keep);
  return ns;
}

static const GtNodeStreamClass *id_filter_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnIdFilterStream),
                                   id_filter_stream_free,
                                   id_filter_stream_next);
  }
  return nsc;
}

static void id_filter_stream_free(GtNodeStream *ns)
{
  AgnIdFilterStream *stream = id_filter_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_hashmap_delete(stream->ids2keep);
}

static int id_filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                 GtError *error)
{
  AgnIdFilterStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = id_filter_stream_cast(ns);

  while(1)
  {
    had_err = gt_node_stream_next(stream->in_stream, gn, error);
    if(had_err)
      return had_err;
    if(!*gn)
      return 0;

    fn = gt_feature_node_try_cast(*gn);
    if(!fn)
      return 0;

    const char *featureid = gt_feature_node_get_attribute(fn, "ID");
    if(gt_hashmap_get(stream->ids2keep, featureid) != NULL)
    {
      return 0;
    }
    else
    {
      // gt_genome_node_delete(*gn);
      continue;
    }
  }

  // *gn = NULL;
  return 0;
}

bool agn_id_filter_stream_unit_test(AgnUnitTest *test)
{
    GtArray *source, *sink;
    GtHashmap *ids;
    GtNodeStream *aos, *ais, *ifs;
    GtUword progress;

    GtError *error = gt_error_new();
    GtQueue *queue = gt_queue_new();
    id_filter_stream_test_data(queue);
    agn_assert(gt_queue_size(queue) == 1);

    source = gt_queue_get(queue);
    sink = gt_array_new( sizeof(GtFeatureNode *) );
    ids = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
    gt_hashmap_add(ids, "gene2", "gene2");
    ais = gt_array_in_stream_new(source, &progress, error);
    ifs = agn_id_filter_stream_new(ais, ids);
    aos = gt_array_out_stream_new(ifs, sink, error);
    gt_node_stream_pull(aos, error);
    bool test1 = gt_array_size(sink) == 1;
    if(test1)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(sink);
      const char *featureid = gt_feature_node_get_attribute(fn, "ID");
      test1 = featureid != NULL && strcmp(featureid, "gene2") == 0;
      gt_genome_node_delete((GtGenomeNode *)fn);
    }
    agn_unit_test_result(test, "gene2", test1);
    gt_array_delete(source);
    gt_array_delete(sink);
    gt_hashmap_delete(ids);
    gt_node_stream_delete(ais);
    gt_node_stream_delete(ifs);
    gt_node_stream_delete(aos);

    return agn_unit_test_success(test);
}

static void id_filter_stream_test_data(GtQueue *queue)
{
  GtError *error = gt_error_new();
  const char *refrfile = "data/gff3/bogus-three-genes.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &refrfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(gff3in, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnFilterStream::id_filter_stream_test_data] error "
            "processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(arraystream);
  gt_queue_add(queue, feats);

  gt_error_delete(error);
}
