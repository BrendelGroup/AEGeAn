/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "extended/array_in_stream_api.h"
#include "extended/array_out_stream_api.h"
#include "extended/feature_node_iterator_api.h"
#include "AgnAttributeFilterStream.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnAttributeFilterStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtHashmap *filters;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

#define attribute_filter_stream_cast(GS)\
        gt_node_stream_cast(attribute_filter_stream_class(), GS)

/**
 * @function Implements the GtNodeStream interface for this class.
 */
static const GtNodeStreamClass* attribute_filter_stream_class(void);

/**
 * @function Class destructor.
 */
static void attribute_filter_stream_free(GtNodeStream *ns);

/**
 * @function Pulls nodes from the input stream and feeds them to the output
 * stream if they pass the provided filtering criteria.
 */
static int attribute_filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                        GtError *error);

/**
 * @function Generate data for unit testing.
 */
static void attribute_filter_stream_test_data(GtQueue *queue);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_attribute_filter_stream_new(GtNodeStream *in_stream,
                                              GtHashmap *filters)
{
  GtNodeStream *ns;
  AgnAttributeFilterStream *stream;
  agn_assert(in_stream && filters);
  ns = gt_node_stream_create(attribute_filter_stream_class(), false);
  stream = attribute_filter_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);
  stream->filters = gt_hashmap_ref(filters);
  return ns;
}

static const GtNodeStreamClass *attribute_filter_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnAttributeFilterStream),
                                   attribute_filter_stream_free,
                                   attribute_filter_stream_next);
  }
  return nsc;
}

bool agn_attribute_filter_stream_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  attribute_filter_stream_test_data(queue);
  bool test1 = gt_queue_size(queue) == 2;
  if(test1)
  {
    GtGenomeNode *gene = gt_queue_get(queue);
    test1 = test1 && gt_genome_node_get_start(gene) == 646485;
    gt_genome_node_delete(gene);

    gene = gt_queue_get(queue);
    test1 = test1 && gt_genome_node_get_start(gene) == 654911;
    gt_genome_node_delete(gene);
  }

  agn_unit_test_result(test, "Pogonomyrmex barbatus (partial)", test1);
  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static void attribute_filter_stream_free(GtNodeStream *ns)
{
  AgnAttributeFilterStream *stream = attribute_filter_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_hashmap_delete(stream->filters);
}

static int attribute_filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                        GtError *error)
{
  AgnAttributeFilterStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = attribute_filter_stream_cast(ns);

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

    bool dodelete = false;
    GtFeatureNode *current;
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
    for(current  = gt_feature_node_iterator_next(iter);
        current != NULL;
        current  = gt_feature_node_iterator_next(iter))
    {
      GtStrArray *attrkeys = gt_feature_node_get_attribute_list(current);
      GtUword i;
      for(i = 0; i < gt_str_array_size(attrkeys); i++)
      {
        const char *key = gt_str_array_get(attrkeys, i);
        const char *value = gt_feature_node_get_attribute(current, key);

        GtStr *keyvalue = gt_str_new_cstr(key);
        gt_str_append_char(keyvalue, '=');
        gt_str_append_cstr(keyvalue, value);

        void *test = gt_hashmap_get(stream->filters, gt_str_get(keyvalue));
        gt_str_delete(keyvalue);
        if(test != NULL)
        {
          dodelete = true;
          break;
        }
      }
      gt_str_array_delete(attrkeys);
    }
    gt_feature_node_iterator_delete(iter);
    if(dodelete)
      gt_genome_node_delete((GtGenomeNode *)fn);
    else
      return 0;
  }

  return 0;
}

static void attribute_filter_stream_test_data(GtQueue *queue)
{
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams = gt_queue_new();

  const char *infile = "data/gff3/pbar-partial.gff3";
  current_stream = gt_gff3_in_stream_new_unsorted(1, &infile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  GtHashmap *attrs = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
  char *keyvalue = gt_cstr_dup("partial=true");
  gt_hashmap_add(attrs, keyvalue, keyvalue);
  current_stream = agn_attribute_filter_stream_new(last_stream, attrs);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;
  gt_hashmap_delete(attrs);

  GtError *error = gt_error_new();
  GtArray *genes = gt_array_new( sizeof(GtGenomeNode *) );
  current_stream = gt_array_out_stream_new(last_stream, genes, error);
  agn_assert(!gt_error_is_set(error));
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "error loading unit test data: %s\n", gt_error_get(error));
    exit(1);
  }

  gt_array_reverse(genes);
  while(gt_array_size(genes) > 0)
  {
    GtGenomeNode **gene = gt_array_pop(genes);
    gt_queue_add(queue, *gene);
  }

  while(gt_queue_size(streams) > 0)
  {
    GtNodeStream *ns = gt_queue_get(streams);
    gt_node_stream_delete(ns);
  }
  gt_queue_delete(streams);
  gt_error_delete(error);
  gt_array_delete(genes);
}
