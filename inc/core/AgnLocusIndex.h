#ifndef AEGEAN_LOCUS_INDEX
#define AEGEAN_LOCUS_INDEX

#include "genometools.h"
#include "AgnComparEval.h"
#include "AgnGeneLocus.h"
#include "AgnLogger.h"

typedef struct AgnLocusIndex AgnLocusIndex;
typedef void (*AgnLocusIndexVisitFunc)(AgnGeneLocus *, void *);

/**
 *
 */
void agn_locus_index_comparative_analysis(AgnLocusIndex *idx, const char *seqid,
                                          int numprocs,
                                          AgnLocusIndexVisitFunc preanalyfunc,
                                          AgnLocusIndexVisitFunc postanalyfunc,
                                          void *analyfuncdata,
                                          AgnLogger *logger);

/**
 * Free memory allocated to the given locus index.
 *
 * @param[in] idx    the locus index object
 */
void agn_locus_index_delete(AgnLocusIndex *idx);

/**
 * Find all overlapping features in the given range stored in this locus index.
 *
 * @param[in]  idx      locus index object
 * @param[in]  seqid    ID of the sequence of interest
 * @param[in]  range    range of interest
 * @param[out] loci     array in which all loci overlapping the specified range
 *                      will be stored
 */
void agn_locus_index_find(AgnLocusIndex *idx, const char *seqid, GtRange *range,
                          GtArray *loci);

/**
 * Retrieve all loci corresponding to the specified sequence ID.
 *
 * @param[in] idx      locus index object
 * @param[in] seqid    ID of the sequence of interest
 * @returns            array containing all loci corresponding to given seqid
 */
GtArray *agn_locus_index_get(AgnLocusIndex *idx, const char *seqid);

/**
 * Compute interval loci (more detailed description forthcoming).
 *
 * @param[in] idx      locus index object
 * @param[in] seqid    ID of the sequence of interest
 * @param[in] delta    number of nucleotides to extend gene loci up- and down-
 *                     stream assuming this does not cause overlap with a gene
 *                     belonging to another locus
 * @returns            a list of interval loci corresponding to this sequence
 */
GtArray *agn_locus_index_interval_loci(AgnLocusIndex *idx, const char *seqid,
                                       unsigned long delta);

/**
 * Allocate memory for a new locus index object.
 *
 * @param[in] freeondelete    indicate whether the locus index object should
 *                            deallocate locus memory on delete
 * @returns                   a pointer to the new object
 */
AgnLocusIndex *agn_locus_index_new(bool freeondelete);

/**
 * Given a pair of annotation feature sets in memory, identify loci while
 * keeping the two sources of annotation separate (to enable comparison).
 *
 * @param[out] idx          locus index object to be populated
 * @param[in]  refrfeats    reference annotations
 * @param[in]  predfeats    prediction annotations
 * @param[in]  numprocs     number of processors to use while identifying loci
 * @param[in]  filters      filtering criteria
 * @param[out] logger       object to which warning/error messages will be
 *                          written if necessary
 * @returns                 the number of loci identified
 */
unsigned long agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx,
                  GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats,
                  int numprocs, AgnCompareFilters *filters, AgnLogger *logger);

/**
 * Given a pair of annotation feature sets in memory, identify loci while
 * keeping the two sources of annotation separate (to enable comparison).
 *
 * @param[out] idx          locus index object to be populated
 * @param[in]  refrfile     path to GFF3 file containing reference annotations
 * @param[in]  predfile     path to GFF3 file containing prediction annotations
 * @param[in]  numprocs     number of processors to use while identifying loci
 * @param[in]  filters      filtering criteria
 * @param[out] logger       object to which warning/error messages will be
 *                          written if necessary
 * @returns                 the number of loci identified
 */
unsigned long agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx,
                  const char *refrfile, const char *predfile, int numprocs,
                  AgnCompareFilters *filters, AgnLogger *logger);

/**
 * Identify loci given an index of annotation features.
 *
 * @param[out] idx          locus index object to be populated
 * @param[in]  features     index containing gene features
 * @param[in]  numprocs     number of processors to use while identifying loci
 * @param[out] logger       object to which warning/error messages will be
 *                          written if necessary
 * @returns                 the number of loci identified
 */
unsigned long agn_locus_index_parse_memory(AgnLocusIndex *idx,
                  GtFeatureIndex *features, int numprocs, AgnLogger *logger);

/**
 * Identify loci given an index of annotation features.
 *
 * @param[out] idx          locus index object to be populated
 * @param[in]  filenames    list of filenames corresponding to GFF3 files
 * @param[in]  numprocs     number of processors to use while identifying loci
 * @param[out] logger       object to which warning/error messages will be
 *                          written if necessary
 * @returns                 the number of loci identified
 */
unsigned long agn_locus_index_parse_disk(AgnLocusIndex *idx, int numfiles,
                  const char **filenames, int numprocs, AgnLogger *logger);

/**
 * Get a list of the seqids stored in this locus index.
 *
 * @param[in] idx    the locus index
 * @returns          pointer to a string array containing the seqids
 */
GtStrArray *agn_locus_index_seqids(AgnLocusIndex *idx);

#endif
