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
 * transcript features as direct children.
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
  char filename[512];
  char stylefile[512];
  const char *refrfile;
  const char *predfile;
  const char *refrlabel;
  const char *predlabel;
  GtUword graphic_width;
  int (*track_order_func)(const char *s1, const char *s2, void *data);
};
typedef struct AgnLocusPngMetadata AgnLocusPngMetadata;

/**
 * @function Associate the given transcript annotation with this locus. Rather
 * than calling this function directly, users are recommended to use one of the
 * following macros: ``agn_locus_add_pred_transcript(locus, trans)`` and
 * ``agn_locus_add_refr_transcript(locus, trans)``, to be used when keeping
 * track of an annotation's source is important (i.e. for pairwise comparison);
 * and ``agn_locus_add_transcript(locus, trans)`` otherwise.
 */
void agn_locus_add(AgnLocus *locus, GtFeatureNode *transcript,
                   AgnComparisonSource source);
#define agn_locus_add_pred_transcript(LC, TR)\
        agn_locus_add(LC, TR, PREDICTIONSOURCE)
#define agn_locus_add_refr_transcript(LC, TR)\
        agn_locus_add(LC, TR, REFERENCESOURCE)
#define agn_locus_add_transcript(LC, TR)\
        agn_locus_add(LC, TR, DEFAULTSOURCE)

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
void agn_locus_comparative_analysis(AgnLocus *locus);

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
 * @function Class constructor.
 */
AgnLocus* agn_locus_new(const char *seqid);

/**
 * @function Return the number of clique pairs to be reported for this locus.
 */
GtUword agn_locus_num_clique_pairs(AgnLocus *locus);

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
 * @function Get the transcripts associated with this locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_locus_pred_transcripts(locus)`` to retrieve prediction
 * transcripts, ``agn_locus_refr_transcripts(locus)`` to retrieve reference
 * transcripts, or ``agn_locus_get_genes(locus)`` if the source of
 * annotation is undesignated or irrelevant.
 */
GtArray *agn_locus_transcripts(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_transcripts(LC)\
        agn_locus_transcripts(LC, PREDICTIONSOURCE)
#define agn_locus_refr_transcripts(LC)\
        agn_locus_transcripts(LC, REFERENCESOURCE)
#define agn_locus_get_transcripts(LC)\
        agn_locus_transcripts(LC, DEFAULTSOURCE)

/**
 * @function Get the transcript IDs associated with this locus. Rather than
 * calling this function directly, users are encouraged to use one of the
 * following macros: ``agn_locus_pred_transcripts(locus)`` to retrieve
 * prediction IDs, ``agn_locus_refr_transcripts(locus)`` to retrieve
 * reference IDs, or ``agn_locus_get_genes(locus)`` if the source of
 * annotation is undesignated or irrelevant.
 */
GtArray *agn_locus_transcript_ids(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_pred_transcript_ids(LC)\
        agn_locus_transcript_ids(LC, PREDICTIONSOURCE)
#define agn_locus_refr_transcript_ids(LC)\
        agn_locus_transcript_ids(LC, REFERENCESOURCE)
#define agn_locus_get_transcript_ids(LC)\
        agn_locus_transcript_ids(LC, DEFAULTSOURCE)

/**
 * @function Get the number of transcripts for the locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_transcript_locus_num_pred_transcripts(locus)`` for the number
 * of prediction transcripts,
 * ``agn_transcript_locus_num_refr_transcripts(locus)`` for the number of
 * reference transcripts, or ``agn_transcript_locus_num_transcripts(locus)``
 * if the source of annotation is undesignated or irrelevant.
 */
GtUword agn_locus_transcript_num(AgnLocus *locus, AgnComparisonSource src);
#define agn_locus_num_pred_transcripts(LC)\
        agn_locus_transcript_num(LC, PREDICTIONSOURCE)
#define agn_locus_num_refr_transcripts(LC)\
        agn_locus_transcript_num(LC, REFERENCESOURCE)
#define agn_locus_num_transcripts(LC)\
        agn_locus_transcript_num(LC, DEFAULTSOURCE)

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_locus_unit_test(AgnUnitTest *test);

#endif
