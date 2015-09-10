/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_LOCUS
#define AEGEAN_LOCUS

#include "annotationsketch/block_api.h"
#include "core/array_api.h"
#include "extended/feature_node_api.h"
#include "AgnCliquePair.h"

/**
 * @class AgnLocus
 *
 * The AgnLocus class represents gene loci and interval loci in memory and can
 * be used to facilitate comparison of two different sources of annotation.
 * Under the hood, each ``AgnLocus`` object is a feature node with one or more
 * gene features as direct children.
 */
typedef GtGenomeNode AgnLocus;

/**
 * @type When tracking the source of an annotation for comparison purposes, use
 * this enumerated type to refer to reference (``REFERENCESOURCE``) vs
 * prediction (``PREDICTIONSOURCE``) annotations. ``DEFAULTSOURCE`` is for when
 * the source is not a concern.
 */
enum AgnComparisonSource
{
  REFERENCESOURCE,
  PREDICTIONSOURCE,
  DEFAULTSOURCE
};
typedef enum AgnComparisonSource AgnComparisonSource;

/**
 * @type This data structure provides a convenient container for metadata needed
 * to produce a PNG graphic for pairwise comparison
 * loci.
 */
struct AgnLocusPngMetadata
{
  char filename_template[512];
  char stylefile[512];
  const char *refrfile;
  const char *predfile;
  const char *refrlabel;
  const char *predlabel;
};
typedef struct AgnLocusPngMetadata AgnLocusPngMetadata;

/**
 * @type Comparison operators to use when filtering loci.
 */
enum AgnLocusFilterOp
{
  AGN_LOCUS_FILTER_EQ,
  AGN_LOCUS_FILTER_NE,
  AGN_LOCUS_FILTER_GT,
  AGN_LOCUS_FILTER_GE,
  AGN_LOCUS_FILTER_LT,
  AGN_LOCUS_FILTER_LE,
  AGN_LOCUS_FILTER_NO, // no-op
};
typedef enum AgnLocusFilterOp AgnLocusFilterOp;

/**
 * @type Data by which to filter a locus. If the value returned by ``function``
 * satisfies the criterion specified by ``testvalue`` and ``operator``, then
 * the locus is to be kept.
 */
struct AgnLocusFilter
{
  GtUword (*function) (AgnLocus *, AgnComparisonSource);
  GtUword testvalue;
  AgnLocusFilterOp operator;
  AgnComparisonSource src;
};
typedef struct AgnLocusFilter AgnLocusFilter;

/**
 * @function Associate the given annotation with this locus. Rather
 * than calling this function directly, users are recommended to use one of the
 * following macros: ``agn_locus_add_pred_feature(locus, gene)`` and
 * ``agn_locus_add_refr_feature(locus, gene)``, to be used when keeping
 * track of an annotation's source is important (i.e. for pairwise comparison);
 * and ``agn_locus_add_feature(locus, gene)`` otherwise.
 */
void agn_locus_add(AgnLocus *locus, GtFeatureNode *feature,
                   AgnComparisonSource source);
#define agn_locus_add_pred_feature(LC, GN)\
        agn_locus_add(LC, GN, PREDICTIONSOURCE)
#define agn_locus_add_refr_feature(LC, GN)\
        agn_locus_add(LC, GN, REFERENCESOURCE)
#define agn_locus_add_feature(LC, GN)\
        agn_locus_add(LC, GN, DEFAULTSOURCE)

/**
 * @function Do a semi-shallow copy of this data structure--for members whose
 * data types support reference counting, the same pointer is used and the
 * reference is incremented. For the other members a new object is created and
 * populated with the same content.
 */
AgnLocus *agn_locus_clone(AgnLocus *locus);

/**
 * @function The combined length of all coding sequences associated with this
 * locus. Rather than calling this function directly, users are encouraged to
 * use one of the following macros: ``agn_locus_refr_cds_length(locus)``
 * for the combined length of all reference CDSs,
 * ``agn_locus_pred_cds_length(locus)`` for the combined length of all
 * prediction CDSs, and ``agn_locus_get_cds_length(locus)`` for the
 * combined length of all CDSs.
 */
GtUword agn_locus_cds_length(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_cds_length(LC)\
        agn_locus_cds_length(LC, PREDICTIONSOURCE)
#define agn_locus_refr_cds_length(LC)\
        agn_locus_cds_length(LC, REFERENCESOURCE)
#define agn_locus_get_cds_length(LC)\
        agn_locus_cds_length(LC, DEFAULTSOURCE)

/**
 * @function Compare every reference transcript clique with every prediction
 * transcript clique. For gene loci with multiple transcript cliques, each
 * comparison is not necessarily reported. Instead, we report the set of clique
 * pairs that provides the optimal pairing of reference and prediction
 * transcripts. If there are more reference transcript cliques than prediction
 * cliques (or vice versa), these unmatched cliques are reported separately.
 */
void agn_locus_comparative_analysis(AgnLocus *locus, GtLogger *logger);

/**
 * @function Analog of ``strcmp`` for sorting AgnLocus objects. Loci are first
 * sorted lexicographically by sequence ID, and then spatially by genomic
 * coordinates.
 */
int agn_locus_array_compare(const void *p1, const void *p2);

/**
 * @function Add this locus' internal comparison stats to a larger set of
 * aggregate stats.
 */
void agn_locus_comparison_aggregate(AgnLocus *locus, AgnComparison *comp);

/**
 * @function Add this locus' internal comparison stats to a larger set of
 * aggregate stats.
 */
void agn_locus_data_aggregate(AgnLocus *locus, AgnComparisonData *data);

/**
 * @function Class destructor.
 */
void agn_locus_delete(AgnLocus *locus);

/**
 * @function Get the number of exons for the locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * ``agn_locus_num_pred_exons(locus)`` for the number of prediction exons,
 * ``agn_locus_num_refr_exons(locus)`` for the number of reference exons,
 * or ``agn_locus_num_exons(locus)`` if the source of annotation is
 * undesignated or irrelevant.
 */
GtUword agn_locus_exon_num(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_num_pred_exons(LC)\
        agn_locus_exon_num(LC, PREDICTIONSOURCE)
#define agn_locus_num_refr_exons(LC)\
        agn_locus_exon_num(LC, REFERENCESOURCE)
#define agn_locus_num_exons(LC)\
        agn_locus_exon_num(LC, DEFAULTSOURCE)

/**
 * @function Parse filters from ``filterfile`` and place ``AgnLocusFilter``
 * objects in ``filters``.
 */
void agn_locus_filter_parse(FILE *filterfile, GtArray *filters);

/**
 * @function Return true if ``locus`` satisfies the given filtering
 * criterion.
 */
bool agn_locus_filter_test(AgnLocus *locus, AgnLocusFilter *filter);


/**
 * @function Return an array of the locus' top-level children, regardless of
 * their type.
 */
GtArray *agn_locus_get(AgnLocus *locus);

/**
 * @function Get a list of all the prediction transcript cliques that have no
 * corresponding reference transcript clique.
 */
GtArray *agn_locus_get_unique_pred_cliques(AgnLocus *locus);

/**
 * @function Get a list of all the reference transcript cliques that have no
 * corresponding prediction transcript clique.
 */
GtArray *agn_locus_get_unique_refr_cliques(AgnLocus *locus);

/**
 * @function Get the genes associated with this locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_locus_pred_genes(locus)`` to retrieve prediction genes,
 * ``agn_locus_refr_genes(locus)`` to retrieve reference genes, or
 * ``agn_locus_get_genes(locus)`` if the source of annotation is undesignated or
 * irrelevant.
 */
GtArray *agn_locus_genes(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_genes(LC)\
        agn_locus_genes(LC, PREDICTIONSOURCE)
#define agn_locus_refr_genes(LC)\
        agn_locus_genes(LC, REFERENCESOURCE)
#define agn_locus_get_genes(LC)\
        agn_locus_genes(LC, DEFAULTSOURCE)

/**
 * @function Get the gene IDs associated with this locus. Rather than
 * calling this function directly, users are encouraged to use one of the
 * following macros: ``agn_locus_pred_gene_ids(locus)`` to retrieve
 * prediction IDs, ``agn_locus_refr_gene_ids(locus)`` to retrieve
 * reference IDs, or ``agn_locus_get_gene_ids(locus)`` if the source of
 * annotation is undesignated or irrelevant.
 */

GtArray *agn_locus_gene_ids(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_gene_ids(LC)\
        agn_locus_gene_ids(LC, PREDICTIONSOURCE)
#define agn_locus_refr_gene_ids(LC)\
        agn_locus_gene_ids(LC, REFERENCESOURCE)
#define agn_locus_get_gene_ids(LC)\
        agn_locus_gene_ids(LC, DEFAULTSOURCE)

/**
 * @function Get the number of genes for the locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_locus_num_pred_genes(locus)`` for the number of prediction
 * genes, ``agn_locus_num_refr_genes(locus)`` for the number of reference
 * genes, or ``agn_locus_num_genes(locus)`` if the source of annotation is
 * undesignated or irrelevant.
 */
GtUword agn_locus_gene_num(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_num_pred_genes(LC)\
        agn_locus_gene_num(LC, PREDICTIONSOURCE)
#define agn_locus_num_refr_genes(LC)\
        agn_locus_gene_num(LC, REFERENCESOURCE)
#define agn_locus_num_genes(LC)\
        agn_locus_gene_num(LC, DEFAULTSOURCE)

/**
 * @function Given two adjacent gene-containing iLoci, determine their
 * orientation: 0 for both forward ('>>'), 1 for inner ('><'), 2 for outer
 * ('<>'), and 3 for reverse ('<<').
 */
int agn_locus_inner_orientation(AgnLocus *left, AgnLocus *right);

/**
 * @function Get the mRNAs associated with this locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_locus_pred_mrnas(locus)`` to retrieve prediction mRNAs,
 * ``agn_locus_refr_mrnas(locus)`` to retrieve reference mRNAs, or
 * ``agn_locus_get_mrnas(locus)`` if the source of annotation is undesignated or
 * irrelevant.
 */
GtArray *agn_locus_mrnas(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_mrnas(LC)\
        agn_locus_mrnas(LC, PREDICTIONSOURCE)
#define agn_locus_refr_mrnas(LC)\
        agn_locus_mrnas(LC, REFERENCESOURCE)
#define agn_locus_get_mrnas(LC)\
        agn_locus_mrnas(LC, DEFAULTSOURCE)

/**
 * @function Get the mRNA IDs associated with this locus. Rather than
 * calling this function directly, users are encouraged to use one of the
 * following macros: ``agn_locus_pred_mrna_ids(locus)`` to retrieve
 * prediction IDs, ``agn_locus_refr_mrna_ids(locus)`` to retrieve
 * reference IDs, or ``agn_locus_get_mrna_ids(locus)`` if the source of
 * annotation is undesignated or irrelevant.
 */
GtArray *agn_locus_mrna_ids(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_mrna_ids(LC)\
        agn_locus_mrna_ids(LC, PREDICTIONSOURCE)
#define agn_locus_refr_mrna_ids(LC)\
        agn_locus_mrna_ids(LC, REFERENCESOURCE)
#define agn_locus_get_mrna_ids(LC)\
        agn_locus_mrna_ids(LC, DEFAULTSOURCE)

/**
 * @function Get the number of mRNAs for the locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_locus_num_pred_mrnas(locus)`` for the number of prediction
 * mRNAs, ``agn_locus_num_refr_mrnas(locus)`` for the number of reference
 * mRNAs, or ``agn_locus_num_mrnas(locus)`` if the source of annotation is
 * undesignated or irrelevant.
 */
GtUword agn_locus_mrna_num(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_num_pred_mrnas(LC)\
        agn_locus_mrna_num(LC, PREDICTIONSOURCE)
#define agn_locus_num_refr_mrnas(LC)\
        agn_locus_mrna_num(LC, REFERENCESOURCE)
#define agn_locus_num_mrnas(LC)\
        agn_locus_mrna_num(LC, DEFAULTSOURCE)

/**
 * @function Class constructor.
 */
AgnLocus* agn_locus_new(GtStr *seqid);

/**
 * @function Return the clique pairs to be reported for this locus.
 */
GtArray *agn_locus_pairs_to_report(AgnLocus *locus);

#ifndef WITHOUT_CAIRO
/**
 * @function Track selector function for generating PNG graphics of pairwise
 * comparison loci. The track name to will be written to ``track``.
 */
void agn_locus_png_track_selector(GtBlock *block, GtStr *track,void *data);
#endif

#ifndef WITHOUT_CAIRO
/**
 * @function Print a PNG graphic for this locus.
 */
void agn_locus_print_png(AgnLocus *locus, AgnLocusPngMetadata *metadata);
#endif

/**
 * @function Print a mapping of the transcript(s) associated with this locus in
 * a two-column tab-delimited format: ``transcriptId<tab>locusId``.
 */
void agn_locus_print_transcript_mapping(AgnLocus *locus, FILE *outstream);

/**
 * @function Set the start and end coordinates for this locus.
 */
void agn_locus_set_range(AgnLocus *locus, GtUword start, GtUword end);

/**
 * @function Calculate the splice complexity of this gene locus. Rather than
 * calling this method directly, users are recommended to use one of the
 * following macros: ``agn_locus_prep_splice_complexity(locus)`` to
 * calculate the splice complexity of just the prediction transcripts,
 * ``agn_locus_refr_splice_complexity(locus)`` to calculate the splice
 * complexity of just the reference transcripts, and
 * ``agn_locus_calc_splice_complexity(locus)`` to calculate the splice
 * complexity taking into account all transcripts.
 */
double agn_locus_splice_complexity(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_splice_complexity(LC)\
        agn_locus_splice_complexity(LC, PREDICTIONSOURCE)
#define agn_locus_refr_splice_complexity(LC)\
        agn_locus_splice_complexity(LC, REFERENCESOURCE)
#define agn_locus_calc_splice_complexity(LC)\
        agn_locus_splice_complexity(LC, DEFAULTSOURCE)

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_locus_unit_test(AgnUnitTest *test);

#endif
