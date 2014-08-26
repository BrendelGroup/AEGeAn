/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include "core/queue_api.h"
#include "extended/array_in_stream_api.h"
#include "extended/array_out_stream_api.h"
#include "extended/feature_node_iterator_api.h"
#include "AgnFilterStream.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnFilterStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtQueue *cache;
  GtHashmap *typestokeep;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

#define filter_stream_cast(GS)\
        gt_node_stream_cast(filter_stream_class(), GS)

/**
 * @function Implements the GtNodeStream interface for this class.
 */
static const GtNodeStreamClass* filter_stream_class(void);

/**
 * @function Class destructor.
 */
static void filter_stream_free(GtNodeStream *ns);

/**
 * @function Pulls nodes from the input stream and feeds them to the output
 * stream if they pass the provided filtering criteria.
 */
static int filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GtError *error);

/**
 * @function Generate data for unit testing.
 */
static void filter_stream_test_data(GtQueue *queue);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream,
                                    GtHashmap *typestokeep)
{
  GtNodeStream *ns;
  AgnFilterStream *stream;
  agn_assert(in_stream && typestokeep);
  ns = gt_node_stream_create(filter_stream_class(), false);
  stream = filter_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);
  stream->cache = gt_queue_new();
  stream->typestokeep = gt_hashmap_ref(typestokeep);
  return ns;
}

static const GtNodeStreamClass *filter_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnFilterStream),
                                   filter_stream_free,
                                   filter_stream_next);
  }
  return nsc;
}

static void filter_stream_free(GtNodeStream *ns)
{
  AgnFilterStream *stream = filter_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_queue_delete(stream->cache);
  gt_hashmap_delete(stream->typestokeep);
}

static int filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GtError *error)
{
  AgnFilterStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = filter_stream_cast(ns);

  if(gt_queue_size(stream->cache) > 0)
  {
    *gn = gt_queue_get(stream->cache);
    return 0;
  }

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

    GtFeatureNode *current;
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
    GtHashmap *multi_parents = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    for(current  = gt_feature_node_iterator_next(iter);
        current != NULL;
        current  = gt_feature_node_iterator_next(iter))
    {
      const char *type = gt_feature_node_get_type(current);
      if(gt_hashmap_get(stream->typestokeep, type) == NULL)
        continue;

      gt_genome_node_ref((GtGenomeNode *)current);
      if(gt_feature_node_is_multi(current))
      {
        GtFeatureNode *rep = gt_feature_node_get_multi_representative(current);
        GtFeatureNode *parent = gt_hashmap_get(multi_parents, rep);
        GtRange currentrange = gt_genome_node_get_range((GtGenomeNode*)current);
        if(parent == NULL)
        {
          GtGenomeNode *parentgn = gt_feature_node_new_pseudo_template(rep);
          const char *id = gt_feature_node_get_attribute(rep, "ID");
          const char *parentid = gt_feature_node_get_attribute(rep, "Parent");
          parent = gt_feature_node_cast(parentgn);
          if(id)
            gt_feature_node_set_attribute(parent, "ID", id);
          if(parentid)
            gt_feature_node_set_attribute(parent, "Parent", parentid);
          gt_hashmap_add(multi_parents, rep, parent);
          gt_queue_add(stream->cache, parent);
        }
        GtRange parentrange = gt_genome_node_get_range((GtGenomeNode*)parent);
        GtRange newrange = gt_range_join(&parentrange, &currentrange);
        gt_genome_node_set_range((GtGenomeNode *)parent, &newrange);
        gt_feature_node_add_child(parent, current);
      }
      else
      {
        gt_queue_add(stream->cache, current);
      }
    }
    gt_feature_node_iterator_delete(iter);
    gt_genome_node_delete((GtGenomeNode *)fn);
    gt_hashmap_delete(multi_parents);
    if(gt_queue_size(stream->cache) > 0)
    {
      *gn = gt_queue_get(stream->cache);
      return 0;
    }
  }

  return 0;
}

bool agn_filter_stream_unit_test(AgnUnitTest *test)
{
  GtArray *source, *sink;
  GtHashmap *types;
  GtNodeStream *aos, *ais, *fs;
  GtUword progress;

  GtError *error = gt_error_new();
  GtQueue *queue = gt_queue_new();
  filter_stream_test_data(queue);
  agn_assert(gt_queue_size(queue) == 6);

  source = gt_queue_get(queue);
  sink = gt_array_new( sizeof(GtFeatureNode *) );
  types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(types, "mRNA", "mRNA");
  ais = gt_array_in_stream_new(source, &progress, error);
  fs = agn_filter_stream_new(ais, types);
  aos = gt_array_out_stream_new(fs, sink, error);
  gt_node_stream_pull(aos, error);
  bool test1 = gt_array_size(sink) == 1;
  if(test1)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(sink);
    test1 = agn_typecheck_mrna(fn);
    gt_genome_node_delete((GtGenomeNode *)fn);
  }
  agn_unit_test_result(test, "mRNA", test1);
  gt_array_delete(source);
  gt_array_delete(sink);
  gt_hashmap_delete(types);
  gt_node_stream_delete(ais);
  gt_node_stream_delete(fs);
  gt_node_stream_delete(aos);

  source = gt_queue_get(queue);
  sink = gt_array_new( sizeof(GtFeatureNode *) );
  types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(types, "exon", "exon");
  ais = gt_array_in_stream_new(source, &progress, error);
  fs = agn_filter_stream_new(ais, types);
  aos = gt_array_out_stream_new(fs, sink, error);
  gt_node_stream_pull(aos, error);
  bool test2 = gt_array_size(sink) == 3;
  if(test2)
  {
    while(gt_array_size(sink) > 0)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(sink);
      test2 = test2 && agn_typecheck_exon(fn);
      gt_genome_node_delete((GtGenomeNode *)fn);
    }
  }
  agn_unit_test_result(test, "exon", test2);
  gt_array_delete(source);
  gt_array_delete(sink);
  gt_hashmap_delete(types);
  gt_node_stream_delete(ais);
  gt_node_stream_delete(fs);
  gt_node_stream_delete(aos);

  source = gt_queue_get(queue);
  sink = gt_array_new( sizeof(GtFeatureNode *) );
  types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(types, "CDS", "CDS");
  ais = gt_array_in_stream_new(source, &progress, error);
  fs = agn_filter_stream_new(ais, types);
  aos = gt_array_out_stream_new(fs, sink, error);
  gt_node_stream_pull(aos, error);
  bool test3 = gt_array_size(sink) == 12;
  if(test3)
  {
    while(gt_array_size(sink) > 0)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(sink);
      test3 = test3 && agn_typecheck_cds(fn);
      gt_genome_node_delete((GtGenomeNode *)fn);
    }
  }
  agn_unit_test_result(test, "CDS", test3);
  gt_array_delete(source);
  gt_array_delete(sink);
  gt_hashmap_delete(types);
  gt_node_stream_delete(ais);
  gt_node_stream_delete(fs);
  gt_node_stream_delete(aos);

  source = gt_queue_get(queue);
  sink = gt_array_new( sizeof(GtFeatureNode *) );
  types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(types, "intron", "intron");
  ais = gt_array_in_stream_new(source, &progress, error);
  fs = agn_filter_stream_new(ais, types);
  aos = gt_array_out_stream_new(fs, sink, error);
  gt_node_stream_pull(aos, error);
  bool test4 = gt_array_size(sink) == 6;
  if(test4)
  {
    while(gt_array_size(sink) > 0)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(sink);
      test4 = test4 && agn_typecheck_intron(fn);
      gt_genome_node_delete((GtGenomeNode *)fn);
    }
  }
  agn_unit_test_result(test, "intron", test4);
  gt_array_delete(source);
  gt_array_delete(sink);
  gt_hashmap_delete(types);
  gt_node_stream_delete(ais);
  gt_node_stream_delete(fs);
  gt_node_stream_delete(aos);

  source = gt_queue_get(queue);
  sink = gt_array_new( sizeof(GtFeatureNode *) );
  types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(types, "exon", "exon");
  gt_hashmap_add(types, "intron", "intron");
  ais = gt_array_in_stream_new(source, &progress, error);
  fs = agn_filter_stream_new(ais, types);
  aos = gt_array_out_stream_new(fs, sink, error);
  gt_node_stream_pull(aos, error);
  bool test5 = gt_array_size(sink) == 13;
  if(test5)
  {
    unsigned ecount = 0, icount = 0;
    while(gt_array_size(sink) > 0)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(sink);
      if(agn_typecheck_exon(fn))
        ecount++;
      else if(agn_typecheck_intron(fn))
        icount++;
      gt_genome_node_delete((GtGenomeNode *)fn);
    }
    test5 = ecount == 7 && icount == 6;
  }
  agn_unit_test_result(test, "exon+intron", test5);
  gt_array_delete(source);
  gt_array_delete(sink);
  gt_hashmap_delete(types);
  gt_node_stream_delete(ais);
  gt_node_stream_delete(fs);
  gt_node_stream_delete(aos);

  source = gt_queue_get(queue);
  sink = gt_array_new( sizeof(GtFeatureNode *) );
  types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(types, "CDS", "CDS");
  gt_hashmap_add(types, "intron", "intron");
  ais = gt_array_in_stream_new(source, &progress, error);
  fs = agn_filter_stream_new(ais, types);
  aos = gt_array_out_stream_new(fs, sink, error);
  gt_node_stream_pull(aos, error);
  bool test6 = gt_array_size(sink) == 5;
  if(test6)
  {
    unsigned ccount = 0, icount = 0;
    while(gt_array_size(sink) > 0)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(sink);
      if(agn_typecheck_cds(fn))
        ccount++;
      else if(agn_typecheck_intron(fn))
        icount++;
      gt_genome_node_delete((GtGenomeNode *)fn);
    }
    test6 = ccount == 3 && icount == 2;
  }
  agn_unit_test_result(test, "CDS+intron", test6);
  gt_array_delete(source);
  gt_array_delete(sink);
  gt_hashmap_delete(types);
  gt_node_stream_delete(ais);
  gt_node_stream_delete(fs);
  gt_node_stream_delete(aos);

  gt_error_delete(error);
  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static void filter_stream_test_data(GtQueue *queue)
{
  GtError *error = gt_error_new();
  const char *refrfile = "data/gff3/grape-refr-mrnas.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &refrfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(gff3in, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnFilterStream::filter_stream_test_data] error "
            "processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(arraystream);
  gt_array_sort(feats, (GtCompare)agn_genome_node_compare);

  GtUword i;
  for(i = 0; i < gt_array_size(feats); i++)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_get(feats, i);
    if(i >= 6)
    {
      gt_genome_node_delete((GtGenomeNode *)fn);
      continue;
    }
    GtArray *array = gt_array_new( sizeof(GtFeatureNode *) );
    gt_array_add(array, fn);
    gt_queue_add(queue, array);
  }

  gt_array_delete(feats);
  gt_error_delete(error);
}
