/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AGN_UTILS
#define AGN_UTILS

#include "core/array_api.h"
#include "core/str_api.h"
#include "extended/feature_index_api.h"
#include "extended/genome_node_api.h"

#define AGN_MAX_FILENAME_SIZE 4096

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

#ifndef NDEBUG
/* Stolen shamelessley from gt_assert() */
#define agn_assert(expression)                                               \
        do {                                                                 \
          if (!(expression)) {                                               \
            fprintf(stderr, "Assertion failed: (%s), function %s, file %s, " \
                    "line %d.\nThis is a bug, please report it at "          \
                    "https://github.com/standage/AEGeAn/issues.\n",          \
                    #expression, __func__, __FILE__, __LINE__);              \
            /*@ignore@*/                                                     \
            exit(1);                                                         \
            /*@end@*/                                                        \
          }                                                                  \
        } while (false)
#else
#define agn_assert(expression) ((void) 0)
#endif


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
 * @function Traverse the given feature and its subfeatures and find the
 * range occupied by coding sequence, or {0,0} if there is no coding sequence.
 */
GtRange agn_feature_node_get_cds_range(GtFeatureNode *fn);

/**
 * @function Return an array of the given feature node's direct children.
 */
GtArray *agn_feature_node_get_children(GtFeatureNode *fn);

/**
 * @function Remove feature ``fn`` and all its subfeatures from ``root``.
 * Analogous to ``gt_feature_node_remove_leaf`` with the difference that ``fn``
 * need not be a leaf feature.
 */
void agn_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn);

/**
 * @function Returns true if any of the features in ``feats`` overlaps, false
 * otherwise.
 */
bool agn_feature_overlap_check(GtArray *feats);

/**
 * @function Compare function for data type ``GtGenomeNode **``, needed for
 * sorting ``GtGenomeNode *`` stored in ``GtArray`` objects.
 */
int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b);

/**
 * @function Determine the length of an mRNA's 3' UTR.
 */
GtUword agn_mrna_3putr_length(GtFeatureNode *mrna);

/**
 * @function Determine the length of an mRNA's 5'UTR.
 */
GtUword agn_mrna_5putr_length(GtFeatureNode *mrna);

/**
 * @function Determine the length of an mRNA's coding sequence.
 */
GtUword agn_mrna_cds_length(GtFeatureNode *mrna);

/**
 * @function If a top-level feature ``top`` contains a multifeature child (with
 * multi representative ``rep``), use this function to get the complete range of
 * the multifeature.
 */
GtRange agn_multi_child_range(GtFeatureNode *top, GtFeatureNode *rep);

/**
 * @function Determine if two features overlap such that they should be
 * assigned to the same iLocus. Specify the minimum overlap (in bp) required
 * and whether the location of the feature's coding sequence (CDS) should be
 * used or not.
 */
bool agn_overlap_ilocus(GtGenomeNode *f1, GtGenomeNode *f2,
                        GtUword minoverlap, bool by_cds);

/**
 * @function CLI function: provide the name of the program, and this function
 * prints out the AEGeAn version number to the specified outstream.
 */
void agn_print_version(const char *progname, FILE *outstream);

/**
 * @function Format the given non-negative number with commas as the thousands
 * separator. The resulting string will be written to ``buffer``.
 */
int agn_sprintf_comma(GtUword n, char *buffer);

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
