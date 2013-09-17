#ifndef AEGEAN_UTILS
#define AEGEAN_UTILS

/**
 * @module AgnUtils
 * A collection of assorted core utility functions.
 */

#include "genometools.h"
#include "AgnLogger.h"

/**
 * @type AgnSequenceRegion
 * Simple data structure for referencing genomic locations.
 *
 * @member [char *] seqid identifier for the sequence
 * @member [GtRange] range start and stop coordinates for the region
 */
typedef struct
{
  char *seqid;
  GtRange range;
} AgnSequenceRegion;

/**
 * @function The Bron-Kerbosch algorithm is an algorithm for enumerating all
 * maximal cliques in an undirected graph. See the `algorithm's Wikipedia entry
 * <http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm>`_
 * for a description of ``R``, ``P``, and ``X``. All maximal cliques will be
 * stored in ``cliques``. If ``skipsimplecliques`` is true, cliques containing a
 * single item will not be stored.
 */
void agn_bron_kerbosch( GtArray *R, GtArray *P, GtArray *X, GtArray *cliques,
                        bool skipsimplecliques );

/**
 * @function A paper by Eilbeck `et al`
 * (http://dx.doi.org/10.1186/1471-2105-10-67) described the annotation edit
 * distance as a measure for comparative evaluation of annotations. This
 * function calculates the AED between the two given annotations.
 */
double agn_calc_edit_distance(GtFeatureNode *t1, GtFeatureNode *t2);

/**
 * @function Determine the splice complexity of the given set of transcripts.
 */
double agn_calc_splice_complexity(GtArray *transcripts);

/**
 * @function Create the canonical gene structure in memory.
 */
GtFeatureNode *agn_eden();

/**
 * @function If reference transcripts belonging to the same locus overlap, they
 * must be separated before comparison with prediction transcript models (and
 * vice versa). This is an instance of the maximal clique enumeration problem
 * (NP-complete), for which the Bron-Kerbosch algorithm provides a solution.
 */
GtArray* agn_enumerate_feature_cliques(GtArray *feature_set);

/**
 * @function For a set of features, we can construct a graph where each node
 * represents a feature and where two nodes are connected if the corresponding
 * features do not overlap. This function returns the intersection of
 * feature_set with the neighbors of feature (where a "neighbor" refers to an
 * adjacent node).
 */
GtArray* agn_feature_neighbors(GtGenomeNode *feature, GtArray *feature_set);

/**
 * @function Wrapper around the stdio.h function that will exit in case of an IO
 * error.
 */
FILE *agn_fopen(const char *filename, const char *mode, FILE *errstream);

/**
 * @function Load canonical protein-coding genes from the given GFF3 files into
 * memory.
 */
GtFeatureIndex *agn_import_canonical(int numfiles, const char **filenames,
                                     AgnLogger *logger);

/**
 * @function Load features whose type is equal to ``type`` into memory from the
 * given GFF3 files.
 */
GtFeatureIndex *agn_import_simple(int numfiles, const char **filenames,
                                  char *type, AgnLogger *logger);

/**
 * @function Given an exon and the start/stop codons associated with its
 * corresponding mRNA, determine which parts of the exon (if any) correspond to
 * coding sequence. If the exon contains coding sequence, the range of that
 * coding sequence will be stored in ``cds_range`` and the function will return
 * true. Otherwise, the function will return false. If the mRNA is on the
 * forward strand, ``left_codon_range`` should contain the coordinates for the
 * start codon and ``right_codon_range`` should contain coordinates for the stop
 * codon. If the mRNA is on the reverse strand, these should be swapped.
 */
bool agn_infer_cds_range_from_exon_and_codons(GtRange *exon_range,
                                              GtRange *leftcodon_range,
                                              GtRange *rightcodon_range,
                                              GtRange *cds_range);

/**
 * @function Given two feature indices, determine which sequences are common
 * between them and return those sequences' IDs as a string array.
 */
GtStrArray* agn_seq_intersection(GtFeatureIndex *refrfeats,
                                 GtFeatureIndex *predfeats, AgnLogger *logger);

/**
 * @function Given two feature indices, determine all of the sequences that are
 * annotated by either of them and return those sequences' IDs as a string
 * array.
 */
GtStrArray* agn_seq_union(GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats,
                          AgnLogger *logger);

/**
 * @function Format the given non-negative number with commas as the thousands
 * separator. The resulting string will be written to ``buffer``.
 */
int agn_sprintf_comma(unsigned long n, char *buffer);

/**
 * @function Dereference the given pointers and compare the resulting strings
 * (*a la* ``strcmp``).
 */
int agn_string_compare(const void *p1, const void *p2);

/**
 * @function Determine the start and end coordinates of the given transcript's
 * CDS.
 */
GtRange agn_transcript_cds_range(GtFeatureNode *transcript);

/**
 * @function Write the structure of a gene transcript in GenBank format to
 * ``outstream``.
 */
void agn_transcript_structure_gbk(GtFeatureNode *transcript, FILE *outstream);

#endif
