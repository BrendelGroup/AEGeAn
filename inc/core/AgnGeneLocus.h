#ifndef AEGEAN_GENE_LOCUS
#define AEGEAN_GENE_LOCUS

#include "genometools.h"
#include "AgnCliquePair.h"

/**
 * @class AgnGeneLocus
 *
 * The purpose of the AgnGeneLocus class is to store all of the data
 * associated with a distinct locus, in many cases to facilitate the comparison
 * of two sets of gene structure annotations for that locus.
 */
typedef struct AgnGeneLocus AgnGeneLocus;

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
struct AgnGeneLocusPngMetadata
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
typedef struct AgnGeneLocusPngMetadata AgnGeneLocusPngMetadata;

/**
* @type This data structure provides a summary of the data and comparisons
* associated with a given locus.
*/
struct AgnGeneLocusSummary
{
  GtUword start;
  GtUword end;
  GtUword length;
  GtUword refrtrans;
  GtUword predtrans;
  GtUword reported;
  GtUword total;
  AgnCompSummary counts;
};
typedef struct AgnGeneLocusSummary AgnGeneLocusSummary;

/**
 * @function Associate the given gene annotation with this gene locus. Rather
 * than calling this function directly, users are recommended to use one of the
 * following macros: ``agn_gene_locus_add_pred_gene(locus, gene)`` and
 * ``agn_gene_locus_add_refr_gene(locus, gene)``, to be used when keeping track
 * of an annotation's source is important (i.e. for pairwise comparison); and
 * ``agn_gene_locus_add_gene(locus, gene)`` otherwise.
 */
void agn_gene_locus_add(AgnGeneLocus *locus, GtFeatureNode *gene,
                        AgnComparisonSource source);
#define agn_gene_locus_add_pred_gene(LC, GN)\
        agn_gene_locus_add(LC, GN, PREDICTIONSOURCE)
#define agn_gene_locus_add_refr_gene(LC, GN)\
        agn_gene_locus_add(LC, GN, REFERENCESOURCE)
#define agn_gene_locus_add_gene(LC, GN)\
        agn_gene_locus_add(LC, GN, DEFAULTSOURCE)

/**
 * @function Add the locus' comparison statistics to a set of aggregate
 * statistics.
 */
void agn_gene_locus_aggregate_results(AgnGeneLocus *locus,
                                      AgnCompEvaluation *eval);

/**
 * @function Do a semi-shallow copy of this data structure--for members whose
 * data types support reference counting, the same pointer is used and the
 * reference is incremented. For the other members a new object is created and
 * populated with the same content.
 */
AgnGeneLocus *agn_gene_locus_clone(AgnGeneLocus *locus);

/**
 * @function Analog of ``strcmp`` for comparing AgnGeneLocus objects, used for
 * sorting GtArray objects containing AgnGeneLocus objects.
 */
int agn_gene_locus_array_compare(const void *p1, const void *p2);

/**
 * @function The combined length of all coding sequences associated with this
 * locus. Rather than calling this function directly, users are encouraged to
 * use one of the following macros: ``agn_gene_locus_refr_cds_length(locus)``
 * for the combined length of all reference CDSs,
 * ``agn_gene_locus_pred_cds_length(locus)`` for the combined length of all
 * prediction CDSs, and ``agn_gene_locus_get_cds_length(locus)`` for the
 * combined length of all CDSs.
 */
GtUword agn_gene_locus_cds_length(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_pred_cds_length(LC)\
        agn_gene_locus_cds_length(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_cds_length(LC)\
        agn_gene_locus_cds_length(LC, REFERENCESOURCE)
#define agn_gene_locus_get_cds_length(LC)\
        agn_gene_locus_cds_length(LC, DEFAULTSOURCE)

/**
 * @function Compare every reference transcript clique with every prediction
 * transcript clique. For gene loci with multiple transcript cliques, each
 * comparison is not necessarily reported. Instead, we report the set of clique
 * pairs that provides the optimal pairing of reference and prediction
 * transcripts. If there are more reference transcript cliques than prediction
 * cliques (or vice versa), these unmatched cliques are reported separately.
 */
GtArray *agn_gene_locus_comparative_analysis(AgnGeneLocus *locus);
#define agn_gene_locus_pairs_to_report(LOC)\
        agn_gene_locus_comparative_analysis(LOC)

/**
 * @function Class destructor.
 */
void agn_gene_locus_delete(AgnGeneLocus *locus);

/**
 * @function Get the number of exons for the locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * ``agn_gene_locus_num_pred_exons(locus)`` for the number of prediction exons,
 * ``agn_gene_locus_num_refr_exons(locus)`` for the number of reference exons,
 * or ``agn_gene_locus_num_exons(locus)`` if the source of annotation is
 * undesignated or irrelevant.
 */
GtUword agn_gene_locus_exon_num(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_num_pred_exons(LC)\
        agn_gene_locus_exon_num(LC, PREDICTIONSOURCE)
#define agn_gene_locus_num_refr_exons(LC)\
        agn_gene_locus_exon_num(LC, REFERENCESOURCE)
#define agn_gene_locus_num_exons(LC)\
        agn_gene_locus_exon_num(LC, DEFAULTSOURCE)

/**
 * @function Given a set of filtering criteria, determine whether a locus meets
 * those criteria. Returns true if the locus should be filtered (if it does not
 * meet the criteria), false otherwise.
 */
bool agn_gene_locus_filter(AgnGeneLocus *locus, AgnCompareFilters *filters);

/**
 * @function Get the genes associated with this locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * ``agn_gene_locus_pred_genes(locus)`` to retrieve prediction genes,
 * ``agn_gene_locus_refr_genes(locus)`` to retrieve reference genes, or
 * ``agn_gene_locus_get_genes(locus)`` if the source of annotation is
 * undesignated or irrelevant.
 */
GtArray *agn_gene_locus_genes(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_pred_genes(LC)\
        agn_gene_locus_genes(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_genes(LC)\
        agn_gene_locus_genes(LC, REFERENCESOURCE)
#define agn_gene_locus_get_genes(LC)\
        agn_gene_locus_genes(LC, DEFAULTSOURCE)

/**
 * @function Get IDs of the genes associated with this locus. Rather than
 * calling this function directly, users are encouraged to use one of the
 * following macros: ``agn_gene_locus_pred_gene_ids(locus)`` to retrieve
 * prediction genes IDs, ``agn_gene_locus_refr_gene_ids(locus)`` to retrieve
 * reference genes IDs, or ``agn_gene_locus_get_gene_ids(locus)`` if the source
 * of annotation is undesignated or irrelevant.
 */
GtArray *agn_gene_locus_gene_ids(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_pred_gene_ids(LC)\
        agn_gene_locus_gene_ids(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_gene_ids(LC)\
        agn_gene_locus_gene_ids(LC, REFERENCESOURCE)
#define agn_gene_locus_get_gene_ids(LC)\
        agn_gene_locus_gene_ids(LC, DEFAULTSOURCE)

/**
 * @function Get the number of genes for the locus. Rather than calling this
 * function directly, users are encouraged to use one of the following macros:
 * ``agn_gene_locus_num_pred_genes(locus)`` for the number of prediction genes,
 * ``agn_gene_locus_num_refr_genes(locus)`` for the number of reference genes,
 * or ``agn_gene_locus_num_genes(locus)`` if the source of annotation is
 * undesignated or irrelevant.
 */
GtUword agn_gene_locus_gene_num(AgnGeneLocus *locus, AgnComparisonSource src);
#define agn_gene_locus_num_pred_genes(LC)\
        agn_gene_locus_gene_num(LC, PREDICTIONSOURCE)
#define agn_gene_locus_num_refr_genes(LC)\
        agn_gene_locus_gene_num(LC, REFERENCESOURCE)
#define agn_gene_locus_num_genes(LC)\
        agn_gene_locus_gene_num(LC, DEFAULTSOURCE)

/**
 * @function Get this locus' end coordinate.
 */
GtUword agn_gene_locus_get_end(AgnGeneLocus *locus);

/**
 * @function Get this locus' length.
 */
GtUword agn_gene_locus_get_length(AgnGeneLocus *locus);

/**
 * @function Get this locus' sequence ID.
 */
const char* agn_gene_locus_get_seqid(AgnGeneLocus *locus);

/**
 * @function Get this locus' start coordinate.
 */
GtUword agn_gene_locus_get_start(AgnGeneLocus *locus);

/**
 * @function Get a list of all the prediction transcript cliques that have no
 * corresponding reference transcript clique.
 */
GtArray *agn_gene_locus_get_unique_pred_cliques(AgnGeneLocus *locus);

/**
 * @function Get a list of all the reference transcript cliques that have no
 * corresponding prediction transcript clique.
 */
GtArray *agn_gene_locus_get_unique_refr_cliques(AgnGeneLocus *locus);

/**
 * @function Class constructor.
 */
AgnGeneLocus* agn_gene_locus_new(const char *seqid);

/**
 * @function Report the number of clique pairs to be reported for this locus.
 */
GtUword agn_gene_locus_num_clique_pairs(AgnGeneLocus *locus);

#ifndef WITHOUT_CAIRO
/**
 * @function Track selector function for generating PNG graphics of pairwise
 * comparison loci. The track name to will be written to ``track``.
 */
void agn_gene_locus_png_track_selector(GtBlock *block, GtStr *track,void *data);
#endif

/**
 * @function Print a mapping of the gene(s) associated with this locus in a two-
 * column tab-delimited format: ``geneId<tab>locusId``.
 */
void agn_gene_locus_print_gene_mapping(AgnGeneLocus *locus, FILE *outstream);

#ifndef WITHOUT_CAIRO
/**
 * @function Print a PNG graphic for this locus.
 */
void agn_gene_locus_print_png(AgnGeneLocus *locus,
                              AgnGeneLocusPngMetadata *metadata);
#endif

/**
 * @function Print a mapping of the transcript(s) associated with this locus in
 * a two-column tab-delimited format: ``transcriptId<tab>locusId``.
 */
void agn_gene_locus_print_transcript_mapping(AgnGeneLocus *locus,
                                             FILE *outstream);

/**
 * @function Return the coordinates of this locus.
 */
GtRange agn_gene_locus_range(AgnGeneLocus *locus);

/**
 * @function Set the range of this locus, no questions asked.
 */
void agn_gene_locus_set_range(AgnGeneLocus *locus, GtUword start, GtUword end);

/**
 * @function Calculate the splice complexity of this gene locus. Rather than
 * calling this method directly, users are recommended to use one of the
 * following macros: ``agn_gene_locus_prep_splice_complexity(locus)`` to
 * calculate the splice complexity of just the prediction transcripts,
 * ``agn_gene_locus_refr_splice_complexity(locus)`` to calculate the splice
 * complexity of just the reference transcripts, and
 * ``agn_gene_locus_calc_splice_complexity(locus)`` to calculate the splice
 * complexity taking into account all transcripts.
 */
double agn_gene_locus_splice_complexity(AgnGeneLocus *locus,
                                        AgnComparisonSource src);
#define agn_gene_locus_pred_splice_complexity(LC)\
        agn_gene_locus_splice_complexity(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_splice_complexity(LC)\
        agn_gene_locus_splice_complexity(LC, REFERENCESOURCE)
#define agn_gene_locus_calc_splice_complexity(LC)\
        agn_gene_locus_splice_complexity(LC, DEFAULTSOURCE)

/**
 * @function Class constructor.
 */
void agn_gene_locus_summary_init(AgnGeneLocusSummary *summary);


/**
 * @function Print the locus in GFF3 format. If ``source`` is NULL, the string
 * "AEGeAn" will be used.
 */
void agn_gene_locus_to_gff3(AgnGeneLocus *locus, FILE *outstream,
                            const char *source);

/**
 * @function Get the transcripts associated with this locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_gene_locus_pred_transcripts(locus)`` to retrieve prediction
 * transcripts, ``agn_gene_locus_refr_transcripts(locus)`` to retrieve reference
 * transcripts, or ``agn_gene_locus_get_genes(locus)`` if the source of
 * annotation is undesignated or irrelevant.
 */
GtArray *agn_gene_locus_transcripts(AgnGeneLocus *locus,
                                    AgnComparisonSource src);
#define agn_gene_locus_pred_transcripts(LC)\
        agn_gene_locus_transcripts(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_transcripts(LC)\
        agn_gene_locus_transcripts(LC, REFERENCESOURCE)
#define agn_gene_locus_get_transcripts(LC)\
        agn_gene_locus_transcripts(LC, DEFAULTSOURCE)

/**
 * @function Get the transcript IDs associated with this locus. Rather than
 * calling this function directly, users are encouraged to use one of the
 * following macros: ``agn_gene_locus_pred_transcripts(locus)`` to retrieve
 * prediction IDs, ``agn_gene_locus_refr_transcripts(locus)`` to retrieve
 * reference IDs, or ``agn_gene_locus_get_genes(locus)`` if the source of
 * annotation is undesignated or irrelevant.
 */
GtArray *agn_gene_locus_transcript_ids(AgnGeneLocus *locus,
                                       AgnComparisonSource src);
#define agn_gene_locus_pred_transcript_ids(LC)\
        agn_gene_locus_transcript_ids(LC, PREDICTIONSOURCE)
#define agn_gene_locus_refr_transcript_ids(LC)\
        agn_gene_locus_transcript_ids(LC, REFERENCESOURCE)
#define agn_gene_locus_get_transcript_ids(LC)\
        agn_gene_locus_transcript_ids(LC, DEFAULTSOURCE)

/**
 * @function Get the number of transcripts for the locus. Rather than calling
 * this function directly, users are encouraged to use one of the following
 * macros: ``agn_transcript_locus_num_pred_transcripts(locus)`` for the number
 * of prediction transcripts,
 * ``agn_transcript_locus_num_refr_transcripts(locus)`` for the number of
 * reference transcripts, or ``agn_transcript_locus_num_transcripts(locus)``
 * if the source of annotation is undesignated or irrelevant.
 */
GtUword agn_gene_locus_transcript_num(AgnGeneLocus *locus,
                                      AgnComparisonSource src);
#define agn_gene_locus_num_pred_transcripts(LC)\
        agn_gene_locus_transcript_num(LC, PREDICTIONSOURCE)
#define agn_gene_locus_num_refr_transcripts(LC)\
        agn_gene_locus_transcript_num(LC, REFERENCESOURCE)
#define agn_gene_locus_num_transcripts(LC)\
        agn_gene_locus_transcript_num(LC, DEFAULTSOURCE)

/**
 * @function Run unit tests for this class. Returns true if all tests passed.
 */
bool agn_gene_locus_unit_test(AgnUnitTest *test);

#endif
