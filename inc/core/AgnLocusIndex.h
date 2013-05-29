#ifndef AEGEAN_LOCUS_INDEX
#define AEGEAN_LOCUS_INDEX

#include "genometools.h"
#include "AgnLogger.h"

typedef struct AgnLocusIndex AgnLocusIndex;

/**
 *
 */
void agn_locus_index_delete(AgnLocusIndex *idx);

/**
 *
 */
AgnLocusIndex *agn_locus_index_new();

/**
 * Given a pair of annotation feature sets in memory, identify loci while
 * keeping the two sources of annotation separate (to enable comparison).
 *
 * @param[out] idx          locus index object to be populated
 * @param[in]  refrfeats    reference annotations
 * @param[in]  predfeats    prediction annotations
 * @param[in]  numprocs     number of processors to use while identifying loci
 * @param[out] logger       object to which warning/error messages will be
 *                          written if necessary
 * @returns                 the number of loci identified
 */
unsigned long agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx,
                  GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats,
                  int numprocs, AgnLogger *logger);

/**
 * Given a pair of annotation feature sets in memory, identify loci while
 * keeping the two sources of annotation separate (to enable comparison).
 *
 * @param[out] idx          locus index object to be populated
 * @param[in]  refrfile     path to GFF3 file containing reference annotations
 * @param[in]  predfile     path to GFF3 file containing prediction annotations
 * @param[in]  numprocs     number of processors to use while identifying loci
 * @param[out] logger       object to which warning/error messages will be
 *                          written if necessary
 * @returns                 the number of loci identified
 */
unsigned long agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx,
                  const char *refrfile, const char *predfile, int numprocs,
                  AgnLogger *logger);

/**
 * Identify loci given an index of annotation features.
 *
 * @param[out] idx          locus index object to be populated
 * @param[in]  features     index containing valid annotation features (FIXME)
 * @param[in]  numprocs     number of processors to use while identifying loci
 * @param[out] logger       object to which warning/error messages will be
 *                          written if necessary
 * @returns                 the number of loci identified
 */
unsigned long agn_locus_index_parse_memory(AgnLocusIndex * idx,
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
unsigned long agn_locus_index_parse_disk(AgnLocusIndex * idx, int numfiles,
                  const char **filenames, int numprocs, AgnLogger *logger);

#endif
