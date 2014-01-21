#ifndef AGN_UTILS
#define AGN_UTILS

#include "core/array_api.h"
#include "core/str_api.h"
#include "extended/feature_index_api.h"
#include "extended/genome_node_api.h"

/**
 * @module AgnUtils
 *
 * Collection of assorted functions that are otherwise unrelated.
 */ //;


/**
 * @type This data structure combines sequence coordinates with a sequence ID to
 * facilitate their usage together.
 */
struct AgnSequenceRegion
{
  GtStr *seqid;
  GtRange range;
};
typedef struct AgnSequenceRegion AgnSequenceRegion;


/**
 * @function Similar to ``gt_array_copy``, except that array elements are
 * treated as pointers and dereferenced before being added to the new array.
 */
GtArray* agn_array_copy(GtArray *source, size_t size);

/**
 * @function Determine the splice complexity of the given set of transcripts.
 */
double agn_calc_splice_complexity(GtArray *transcripts);

/**
 * @function Copy the sequence regions from ``src`` to ``dest``. If ``use_orig``
 * is true, regions specified by input region nodes (such as those parsed from
 * ``##sequence-region`` pragmas in GFF3) are used. Otherwise, regions inferred
 * directly from the feature nodes are used.
 */
GtUword
agn_feature_index_copy_regions(GtFeatureIndex *dest, GtFeatureIndex *src,
                               bool use_orig, GtError *error);

/**
 * @function Copy the sequence regions from ``refrsrc`` and ``predsrc`` to
 * ``dest``. If ``use_orig`` is true, regions specified by input region nodes
 * (such as those parsed from ``##sequence-region`` pragmas in GFF3) are used.
 * Otherwise, regions inferred directly from the feature nodes are used.
 */
GtUword
agn_feature_index_copy_regions_pairwise(GtFeatureIndex *dest,
                                        GtFeatureIndex *refrsrc,
                                        GtFeatureIndex *predsrc,
                                        bool use_orig, GtError *error);

/**
 * @function Compare function for data type ``GtGenomeNode **``, needed for
 * sorting ``GtGenomeNode *`` stored in ``GtArray`` objects.
 */
int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b);

/**
 * @function Dereference the given pointers and compare the resulting strings
 * (*a la* ``strcmp``).
 */
int agn_string_compare(const void *p1, const void *p2);

/**
 * @function Find the strings that are present in either (or both) of the string
 * arrays.
 */
GtStrArray* agn_str_array_union(GtStrArray *a1, GtStrArray *a2);

#endif
