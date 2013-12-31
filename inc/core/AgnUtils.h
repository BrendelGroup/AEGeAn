#ifndef AGN_UTILS
#define AGN_UTILS

#include "core/array_api.h"
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
typedef struct
{
  char *seqid;
  GtRange range;
} AgnSequenceRegion;


/**
 * @function Similar to ``gt_array_copy``, except that array elements are
 * treated as pointers and dereferenced before being added to the new array.
 */
GtArray* agn_array_copy(GtArray *source, size_t size);

/**
 * @function FIXME
 */
double agn_calc_splice_complexity(GtArray *transcripts);

/**
 * @function Compare function for data type ``GtGenomeNode **``, needed for
 * sorting ``GtGenomeNode *`` stored in ``GtArray`` objects.
 */
int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b);

#endif
