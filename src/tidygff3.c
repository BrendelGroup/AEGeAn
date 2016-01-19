/**

Copyright (c) 2010-2016, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include "genometools.h"
#include "aegean.h"

int main(int argc, char **argv)
{
    // Set up the processing stream
    //----------
    gt_lib_init();
    GtQueue *streams = gt_queue_new();

    GtNodeStream *stream = gt_gff3_in_stream_new_unsorted(0, NULL);
    gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)stream);
    gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)stream);
    gt_queue_add(streams, stream);
    GtNodeStream *last_stream = stream;

    stream = agn_pseudogene_fix_stream_new(last_stream);
    gt_queue_add(streams, stream);
    last_stream = stream;

    GtHashmap *filters =
        gt_hashmap_new(GT_HASH_STRING, (GtFree)gt_free_func, NULL);
    char *filter;
    filter = gt_cstr_dup("exception=unclassified transcription discrepancy");
    gt_hashmap_add(filters, filter, filter);
    filter = gt_cstr_dup("exception=unclassified translation discrepancy");
    gt_hashmap_add(filters, filter, filter);
    filter = gt_cstr_dup("gene_biotype=other");
    gt_hashmap_add(filters, filter, filter);
    stream = agn_attribute_filter_stream_new(last_stream, filters);
    gt_queue_add(streams, stream);
    last_stream = stream;
    gt_hashmap_delete(filters);

    GtHashmap *types =
        gt_hashmap_new(GT_HASH_STRING, gt_free_func, gt_free_func);
    GtStr *source = gt_str_new_cstr("AEGeAn::tidygff3");
    gt_hashmap_add(types, gt_cstr_dup("mRNA"), gt_cstr_dup("gene"));
    gt_hashmap_add(types, gt_cstr_dup("rRNA"), gt_cstr_dup("gene"));
    gt_hashmap_add(types, gt_cstr_dup("tRNA"), gt_cstr_dup("gene"));
    stream = agn_infer_parent_stream_new(last_stream, types);
    agn_infer_parent_stream_set_source((AgnInferParentStream *)stream, source);
    gt_queue_add(streams, stream);
    last_stream = stream;
    gt_hashmap_delete(types);
    gt_str_delete(source);

    stream = gt_gff3_out_stream_new(last_stream, NULL);
    gt_queue_add(streams, stream);
    last_stream = stream;

    //----------
    // Execute the processing stream
    //----------
    GtError *error = gt_error_new();
    int had_err = gt_node_stream_pull(last_stream, error);
    if (had_err)
    {
        fprintf(stderr, "Error processing node stream: %s\n",
                gt_error_get(error));
    }
    gt_error_delete(error);

    //----------
    // Free memory
    //----------
    while (gt_queue_size(streams) > 0)
    {
        stream = gt_queue_get(streams);
        gt_node_stream_delete(stream);
    }
    gt_queue_delete(streams);
    gt_lib_clean();
    return had_err;
}
