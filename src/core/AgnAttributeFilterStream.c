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

bool agn_attribute_filter_stream_unit_test(AgnUnitTest *test)
{
  attribute_filter_stream_test_data(NULL);
  return false;
}

static void attribute_filter_stream_test_data(GtQueue *queue)
{
  
}
