#ifndef AEGEAN_GT_EXTENSIONS_MODULE
#define AEGEAN_GT_EXTENSIONS_MODULE

/**
 * @module AgnGtExtensions
 *
 * A collection of extensions to core GenomeTools classes.
 */

#include "genometools.h"

/**
 * This function makes a copy of an array.
 *
 * @param[in] source    array to be copied
 * @param[in] size      size of elements in the array
 * @returns             an array with all the same elements as the source
 */
GtArray* agn_gt_array_copy(GtArray *source, size_t size);

/**
 * Write the given feature index to GFF3 format.
 *
 * @param[in]  index        feature index containing the data
 * @param[out] outstream    file to which output will be written
 */
void agn_gt_feature_index_to_gff3(GtFeatureIndex *index, FILE *outstream);

/**
 * Calculate the length of the given transcript's coding sequence
 *
 * @param[in] transcript    the transcript
 * @returns                 the length of its CDS (in codons/amino acids)
 */
unsigned long agn_gt_feature_node_cds_length(GtFeatureNode *transcript);

/**
 * Gather the children of a given feature that have a certain type.
 *
 * @param[in] fn              a feature
 * @param[in] typetestfunc    function used to select children of interest
 * @returns                   an array containing the children that match the
 *                            selection criteria specified by typetestfunc
 */
GtArray *agn_gt_feature_node_children_of_type(GtFeatureNode *fn,
                                         bool (*typetestfunc)(GtFeatureNode *));

/**
 * When a feature has multiple parents but only one of them is valid, this
 * function will fix the Parent attribute so that it only points to the valid
 * parent.
 */
bool agn_gt_feature_node_fix_parent_attribute(GtFeatureNode *feature,
                                              GtFeatureNode *parent);

/**
 * Determine whether the given feature belongs to a CDS.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the feature belongs to a CDS
 */
bool agn_gt_feature_node_is_cds_feature(GtFeatureNode *fn);

/**
 * Determine whether the given feature is an exon.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the is an exon
 */
bool agn_gt_feature_node_is_exon_feature(GtFeatureNode *fn);

/**
 * Determine whether the given feature is a gene.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the is an exon
 */
bool agn_gt_feature_node_is_gene_feature(GtFeatureNode *fn);

/**
 * Determine whether the given feature is an intron.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the is an intron
 */
bool agn_gt_feature_node_is_intron_feature(GtFeatureNode *fn);

/**
 * Determine whether the given feature is a transcript.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the is a transcript
 */
bool agn_gt_feature_node_is_mrna_feature(GtFeatureNode *fn);

/**
 * Determine whether the given feature is a start codon.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the is a start codon
 */
bool agn_gt_feature_node_is_start_codon_feature(GtFeatureNode *fn);

/**
 * Determine whether the given feature is a stop codon.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the is a stop codon
 */
bool agn_gt_feature_node_is_stop_codon_feature(GtFeatureNode *fn);

/**
 * Determine whether the given feature is a UTR.
 *
 * @param[in] fn    the feature to be tested
 * @returns         a boolean indicating whether the is a UTR
 */
bool agn_gt_feature_node_is_utr_feature(GtFeatureNode *fn);

/**
 * Determine the number of transcripts for the given gene feature.
 *
 * @param[in] gene    the gene
 * @returns           the number of transcripts associated with this gene
 */
unsigned long agn_gt_feature_node_num_transcripts(GtFeatureNode *gene);

/**
 * Determine whether the following features overlap.
 *
 * @param[in] first     a gene feature
 * @param[in] second    another gene feature
 * @returns             boolean indicating whether the features overlap
 */
bool agn_gt_feature_node_overlap(GtFeatureNode *first, GtFeatureNode *second);

/**
 * Determine whether the range of n2 falls within the range of n1.
 *
 * @param[in] n1    a feature node
 * @param[in] n2    another feature node
 * @retuns          true if n2 falls within n1, false otherwise
 */
bool agn_gt_feature_node_range_contains(GtFeatureNode *n1, GtFeatureNode *n2);

/**
 * The 'gt_feature_node_remove_leaf' function only allows removal of a leaf
 * node. This function will remove all of a node's children so that it is a
 * leaf node, which can then be removed.
 *
 * @param[in] root    a root or parent node to the node to be deleted
 * @param[in] fn      the feature node to be removed
 */
void agn_gt_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn);

/**
 * Pseudo nodes are used when connected features have more than one top-level
 * element. This function takes any feature node as input; if it is a pseudo
 * node, it adds all the node's direct children to the provided array. If it is
 * not a pseudo node, the node itself is added to the array. This could be a
 * recursive function, but for now it will simply fail if any of the root's
 * direct children are pseudo nodes.
 *
 * @param[in]  root     the node, possibly a pseudo node
 * @param[out] nodes    array to which top-level features will be added
 */
void agn_gt_feature_node_resolve_pseudo_node(GtFeatureNode *root, GtArray *nodes);

/**
 * Reset the source of the feature and all its children to the given value
 *
 * @param[out] feature    the feature whose source will be updated
 * @param[in]  source     the new source value
 */
void agn_gt_feature_node_set_source_recursive( GtFeatureNode *feature,
                                               GtStr *source );

/**
 * Print the given feature in GFF3 format
 *
 * @param[in]  feature           the feature to be printed
 * @param[out] outstream         file to which GFF3 data will be written
 * @param[in]  printchildren     flag indicating whether to print the feature's
 *                               children as well
 * @param[in]  filtered_types    collection of feature types to exclude from
 *                               printing; if NULL is provided, will use a
 *                               default exclusion list
 * @param[in]  prefix            optional prefix
 */
void agn_gt_feature_node_to_gff3(GtFeatureNode *feature, FILE *outstream,
                                 bool printchildren, char *prefix,
                                 GtHashmap *filtered_types);

/**
 * Comparison function to be used for sorting GtGenomeNode objects stored in a
 * GtArray (for GtDlist, use gt_genome_node_cmp).
 *
 * @param[in] n1    a node
 * @param[in] n2    another node
 * @returns         1 if n1 comes before n2, -1 if n2 comes before n1, and 0 if
 *                  n1 and n2 are sort-equal
 */
int agn_gt_genome_node_compare(const void *n1, const void *n2);

/**
 * Convert a GtPhase object into its corresponding character representation
 *
 * @param[in] phase    an object representing the phase of a CDS segment
 * @returns            the phase's character representation
 */
char agn_gt_phase_to_char(GtPhase phase);

/**
 * Convert a GtStrand object into its corresponding character representation
 *
 * @param[in] phase    an object representing the strand of a sequence feature
 * @returns            the strand's character representation
 */
char agn_gt_strand_to_char(GtStrand strand);

/**
 * Find the strings that are common to both string arrays.
 *
 * @param[in] a1    an array of strings
 * @param[in] a2    another array of strings
 * @returns         array of shared strings
 */
GtStrArray* agn_gt_str_array_intersection(GtStrArray *a1, GtStrArray *a2);

/**
 * Find the strings that are present in either (or both) of the string arrays.
 *
 * @param[in] a1    an array of strings
 * @param[in] a2    another array of strings
 * @returns         array of each string contained in either array
 */
GtStrArray* agn_gt_str_array_union(GtStrArray *a1, GtStrArray *a2);

#endif
