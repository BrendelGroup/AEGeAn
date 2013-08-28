#ifndef AEGEAN_UTILS_MODULE
#define AEGEAN_UTILS_MODULE

/**
 * A collection of assorted utility functions for ParsEval.
 */

#include "genometools.h"
#include "AgnLogger.h"

/**
 * Simple data structure for referencing genomic locations: sequence identifier
 * plus start and end coordinates.
 */
typedef struct
{
  char *seqid;
  GtRange range;
} AgnSequenceRegion;

/**
 * The Bron-Kerbosch algorithm is an algorithm for enumerating all maximal
 * cliques in an undirected graph. See
 * http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm.
 *
 * @param[in]  R                    see algorithm
 * @param[in]  P                    see algorithm
 * @param[in]  X                    see algorithm
 * @param[out] cliques              each maximal clique is stored in this array
 * @param[in]  skipsimplecliques    boolean indicating whether or not to store
 *                                  transcript cliques containing a single
 *                                  transcript
 */
void agn_bron_kerbosch( GtArray *R, GtArray *P, GtArray *X, GtArray *cliques,
                        bool skipsimplecliques );

/**
 * A paper by Eilbeck et al (doi:10.1186/1471-2105-10-67) described the
 * annotation edit distance as a measure for comparative evaluation of
 * annotations. This function calculates the AED between the two given
 * annotations.
 */
double agn_calc_edit_distance(GtFeatureNode *t1, GtFeatureNode *t2);

/**
 * Determine the splice complexity of the given transcripts.
 *
 * @param[in] transcripts    a set of transcripts
 * @returns                  the calculated splice complexity
 */
double agn_calc_splice_complexity(GtArray *transcripts);

/**
 * Create the canonical gene structure in memory.
 *
 * @returns    a pointer to the top-level gene feature
 */
GtFeatureNode *agn_eden();

/**
 * If reference transcripts belonging to the same locus overlap, they must be
 * separated before comparison with prediction transcript models (and vice
 * versa). This is an instance of the maximal clique enumeration problem
 * (NP-complete), for which the Bron-Kerbosch algorithm provides a solution.
 *
 * @param[in] feature_set    a set of features (typically transcripts)
 * @returns                  a set of all maximal cliques from feature_set
 */
GtArray* agn_enumerate_feature_cliques(GtArray *feature_set);

/**
 * For a set of features, we can construct a graph where each node represents a
 * feature and where two nodes are connected if the corresponding features do
 * not overlap. This function returns the intersection of feature_set with the
 * neighbors of feature (where a "neighbor" refers to an adjacent node).
 *
 * @param[in] feature        a feature
 * @param[in] feature_set    a set of features
 * @returns                  the features in feature_set that are adjacent to
 *                           feature in the graph
 */
GtArray* agn_feature_neighbors(GtGenomeNode *feature, GtArray *feature_set);

/**
 * Wrapper around the stdio.h function that will exit in case of an IO error
 *
 * @param[in] filename    name of the file
 * @param[in] mode        access mode for this file handle
 * @param[in] errstream   output stream to which error messages will be written
 * @returns               a file handle to the opened file; will ungracefully
 *                        exit if file handle cannot be created
 */
FILE *agn_fopen(const char *filename, const char *mode, FILE *errstream);

/**
 * Load canonical protein-coding genes from the given GFF3 files into memory.
 *
 * @param[in] numfiles     number of files to be parsed
 * @param[in] filenames    list of filenames
 * @param[in] logger       object for logging error/warning messages
 * @returns                a feature index containing annotations for all files
 */
GtFeatureIndex *agn_import_canonical(int numfiles, const char **filenames,
                                     AgnLogger *logger);

/**
 * Load simple features from the given GFF3 files into memory.
 *
 * @param[in] numfiles     number of files to be parsed
 * @param[in] filenames    list of filenames
 * @param[in] type         feature type to load into memory
 * @param[in] logger       object for logging error/warning messages
 * @returns                a feature index containing annotations for all files
 */
GtFeatureIndex *agn_import_simple(int numfiles, const char **filenames,
                                  char *type, AgnLogger *logger);

/**
 * Given an exon and the start/stop codons associated with its corresponding
 * mRNA, determine which parts of the exon (if any) correspond to coding
 * sequence.
 *
 * @param[in]  exon_range          the coordinates of an exon
 * @param[in]  leftcodon_range     the coordinates of the corresponding mRNA's
 *                                 start codon (or stop codon if on the reverse
 *                                 strand)
 * @param[in]  rightcodon_range    the coordinates of the corresponding mRNA's
 *                                 stop codon (or start codon if on the reverse
 *                                 strand)
 * @param[out] cds_range           if this exon includes any coding sequence,
 *                                 the coordinates of this coding sequence will
 *                                 be written to this; if not, {0,0} will be
 *                                 written
 * @returns                        true if the exon contains coding sequence,
 *                                 false otherwise
 */
bool agn_infer_cds_range_from_exon_and_codons( GtRange *exon_range,
                                               GtRange *leftcodon_range,
                                               GtRange *rightcodon_range,
                                               GtRange *cds_range );

/**
 * Given two feature indices, determine which sequences are common and return
 * them as an array.
 *
 * @param[in]  refrfeats    index of reference features
 * @param[in]  predfeats    index of prediction feature
 * @param[out] error        object to which warning/error messages will be
 *                          written if necessary
 * @returns                 a string array containing IDs of shared sequences
 */
GtStrArray* agn_seq_intersection(GtFeatureIndex *refrfeats,
                                 GtFeatureIndex *predfeats, AgnLogger *logger);

/**
 * Given two feature indices, determine all of the sequence that are annotated
 * between them and return the sequence IDs as an array.
 *
 * @param[in]  refrfeats    index of reference features
 * @param[in]  predfeats    index of prediction feature
 * @param[out] error        object to which warning/error messages will be
 *                          written if necessary
 * @returns                 a string array containing IDs of all sequences
 */
GtStrArray* agn_seq_union(GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats,
                          AgnLogger *logger);

/**
 * Format the given non-negative number with commas as the thousands separator.
 *
 * @param[in]  n         the number
 * @param[out] buffer    the string to which the formatted number will be
 *                       written
 */
int agn_sprintf_comma(unsigned long n, char *buffer);

/**
 * Dereference the given pointers and compare the resulting strings.
 *
 * @param[in] p1    pointer to a char array
 * @param[in] p2    pointer to another char array
 * @returns         1 if p1 comes lexicographically before p2, -1 if p2 comes
 *                  lexicographically before p1, and 0 if p1 and p2 are
 *                  lexicographically equal
 */
int agn_string_compare(const void *p1, const void *p2);

/**
 * Determine the start and end coordinates of the given transcript's CDS
 *
 * @param[in] transcript    the transcript feature
 * @returns                 a range object with the start and end positions of
 *                          the transcript's CDS
 */
GtRange agn_transcript_cds_range(GtFeatureNode *transcript);

/**
 * Write the structure of a gene transcript in GenBank format
 *
 * @param[in]  transcript    the transcript feature
 * @param[out] outstream     file to which structure will be written
 */
void agn_transcript_structure_gbk(GtFeatureNode *transcript, FILE *outstream);

#endif
