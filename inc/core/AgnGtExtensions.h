#ifndef AEGEAN_GT_EXTENSIONS
#define AEGEAN_GT_EXTENSIONS

#include "genometools.h"

/**
 * @module AgnGtExtensions
 *
 * A collection of extensions to core GenomeTools classes.
 */

/**
 * @function This function makes a copy of an array.
 */
GtArray* agn_gt_array_copy(GtArray *source, size_t size);

/**
 * @function Write the given feature index to GFF3 format.
 */
void agn_gt_feature_index_to_gff3(GtFeatureIndex *index, FILE *outstream);

/**
 * @function Calculate the length of the given transcript's coding sequence in
 * amino acids.
 */
unsigned long agn_gt_feature_node_cds_length(GtFeatureNode *transcript);

/**
 * @function Gather the children of a given feature that have a certain type.
 * Type is tested by ``typetestfunc``, which accepts a single ``GtFeatureNode``
 * object.
 */
GtArray *agn_gt_feature_node_children_of_type(GtFeatureNode *fn,
                                         bool (*typetestfunc)(GtFeatureNode *));

/**
 * @function When a feature has multiple parents but only one of them is valid,
 * this function will fix the ``Parent`` attribute so that it only points to the
 * valid parent.
 */
bool agn_gt_feature_node_fix_parent_attribute(GtFeatureNode *feature,
                                              GtFeatureNode *parent);

/**
 * @function Determine whether the given feature belongs to a CDS.
 */
bool agn_gt_feature_node_is_cds_feature(GtFeatureNode *fn);

/**
 * @function Determine whether the given feature is an exon.
 */
bool agn_gt_feature_node_is_exon_feature(GtFeatureNode *fn);

/**
 * @function Determine whether the given feature is a gene.
 */
bool agn_gt_feature_node_is_gene_feature(GtFeatureNode *fn);

/**
 * @function Determine whether the given feature is an intron.
 */
bool agn_gt_feature_node_is_intron_feature(GtFeatureNode *fn);

/**
 * @function Determine whether the given feature is an mRNA.
 */
bool agn_gt_feature_node_is_mrna_feature(GtFeatureNode *fn);

/**
 * @function Determine whether the given feature is a start codon.
 */
bool agn_gt_feature_node_is_start_codon_feature(GtFeatureNode *fn);

/**
 * @function Determine whether the given feature is a stop codon.
 */
bool agn_gt_feature_node_is_stop_codon_feature(GtFeatureNode *fn);

/**
 * @function Determine whether the given feature is part of a UTR.
 */
bool agn_gt_feature_node_is_utr_feature(GtFeatureNode *fn);

/**
 * @function Determine the number of transcripts for the given gene feature.
 */
unsigned long agn_gt_feature_node_num_transcripts(GtFeatureNode *gene);

/**
 * @function Determine whether the given features overlap.
 */
bool agn_gt_feature_node_overlap(GtFeatureNode *first, GtFeatureNode *second);

/**
 * @function Determine whether the range of n2 falls within the range of n1.
 */
bool agn_gt_feature_node_range_contains(GtFeatureNode *n1, GtFeatureNode *n2);

/**
 * @function The ``gt_feature_node_remove_leaf`` function only allows removal of
 * a leaf node. This function will remove all of a node's children so that it is
 * a leaf node, which can then be removed.
 */
void agn_gt_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn);

/**
 * @function Reset the source of the feature and all its children to the given
 * value.
 */
void agn_gt_feature_node_set_source_recursive(GtFeatureNode *feature,
                                              GtStr *source);

/**
 * @function Print the given feature in GFF3 format. Use ``filtered_types`` to
 * provide a list of feature types to exclude when printing. If
 * ``filtered_types`` is NULL, a default exclusion list will be used.
 */
void agn_gt_feature_node_to_gff3(GtFeatureNode *feature, FILE *outstream,
                                 bool printchildren, char *prefix,
                                 GtHashmap *filtered_types);

/**
 * @function Comparison function to be used for sorting GtGenomeNode objects
 * stored in a GtArray (for GtDlist, use gt_genome_node_cmp).
 */
int agn_gt_genome_node_compare(const void *n1, const void *n2);

/**
 * @function Convert a GtPhase object into its corresponding character
 * representation.
 */
char agn_gt_phase_to_char(GtPhase phase);

/**
 * @function Convert a GtStrand object into its corresponding character
 * representation.
 */
char agn_gt_strand_to_char(GtStrand strand);

/**
 * @function Find the strings that are common to both string arrays.
 */
GtStrArray* agn_gt_str_array_intersection(GtStrArray *a1, GtStrArray *a2);

/**
 * @function Find the strings that are present in either (or both) of the string
 * arrays.
 */
GtStrArray* agn_gt_str_array_union(GtStrArray *a1, GtStrArray *a2);

#endif
