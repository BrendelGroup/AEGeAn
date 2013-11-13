#ifndef AEGEAN_LOCUS_INDEX
#define AEGEAN_LOCUS_INDEX

#include "genometools.h"
#include "AgnComparEval.h"
#include "AgnGeneLocus.h"
#include "AgnLogger.h"

/**
 * @class AgnLocusIndex
 *
 * FIXME
 */
typedef struct AgnLocusIndex AgnLocusIndex;

/**
 * @functype Signature functions must match to be applied to each locus in the
 * index. The function will be called once for each locus, which will be passed
 * as the first argument to the function. a second argument is available for an
 * optional pointer to supplementary data (if needed). See
 * :c:func:`agn_locus_index_comparative_analysis`.
 */
typedef void (*AgnLocusIndexVisitFunc)(AgnGeneLocus *, void *);

/**
 * @function Perform a comparative analysis of each locus associated with
 * ``seqid`` in this index. If ``preanalyfunc`` is not NULL, it will be applied
 * to each locus immediately before comparative analysis. If ``postanalyfunc``
 * is not NULL, it will be applied to each locus immediately following
 * comparative analysis. ``analyfuncdata`` will be passed as supplementary data
 * to both functions.
 */
void agn_locus_index_comparative_analysis(AgnLocusIndex *idx, const char *seqid,
                                          AgnLocusIndexVisitFunc preanalyfunc,
                                          AgnLocusIndexVisitFunc postanalyfunc,
                                          void *analyfuncdata,
                                          AgnLogger *logger);

/**
 * @function Class destructor.
 */
void agn_locus_index_delete(AgnLocusIndex *idx);

/**
 * @functype Find all overlapping features in the given range stored in this
 * locus index and store them in ``loci``.
 */
void agn_locus_index_find(AgnLocusIndex *idx, const char *seqid, GtRange *range,
                          GtArray *loci);

/**
 * @function Retrieve all loci corresponding to the specified sequence ID.
 */
GtArray *agn_locus_index_get(AgnLocusIndex *idx, const char *seqid);

/**
 * @functype Compute interval loci with the given ``delta``. If running on
 * incomplete (contig/scaffold) genomic sequences, consider setting
 * ``skipterminal`` to true to ignore the ends of the sequence.
 */
GtArray *agn_locus_index_interval_loci(AgnLocusIndex *idx, const char *seqid,
                                       GtUword delta, bool skipterminal);

/**
 * @function Class constructor
 */
AgnLocusIndex *agn_locus_index_new(bool freeondelete);

/**
 * @function Given a pair of annotation feature sets in memory, identify loci
 * while keeping the two sources of annotation separate (to enable comparison).
 */
GtUword agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx,
                                              GtFeatureIndex *refrfeats,
                                              GtFeatureIndex *predfeats,
                                              int numprocs,
                                              AgnCompareFilters *filters,
                                              AgnLogger *logger);

/**
 * @function Given a pair of annotation feature sets in memory, identify loci
 * while keeping the two sources of annotation separate (to enable comparison).
 */
GtUword agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx,
                                            const char *refrfile,
                                            const char *predfile, int numprocs,
                                            AgnCompareFilters *filters,
                                            AgnLogger *logger);

/**
 * @function Identify loci given an index of annotation features.
 */
GtUword agn_locus_index_parse_memory(AgnLocusIndex *idx,
                                     GtFeatureIndex *features, int numprocs,
                                     AgnLogger *logger);

/**
 * @function Identify loci from the given set of annotation files.
 */
GtUword agn_locus_index_parse_disk(AgnLocusIndex *idx, int numfiles,
                                   const char **filenames, int numprocs,
                                   AgnLogger *logger);

/**
 * @function Get a list of the seqids stored in this locus index.
 */
GtStrArray *agn_locus_index_seqids(AgnLocusIndex *idx);

#endif
