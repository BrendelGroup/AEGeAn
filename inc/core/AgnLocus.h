#ifndef AEGEAN_LOCUS
#define AEGEAN_LOCUS

#include "genometools.h"

/**
 * A simple data structure for a gene locus
 */
typedef struct
{
  GtDlist *genes;
  GtRange range;
  const char *seqid;
} AgnLocus;

/**
 * Associate the given gene with this locus
 *
 * @param[in] locus    the locus
 * @param[in] gene     the new gene to be associated
 */
void agn_locus_add(AgnLocus *locus, GtFeatureNode *gene);

/**
 * Free the memory occupied by this locus object
 *
 * @param[in] locus    the locus
 */
void agn_locus_delete(AgnLocus *locus);

/**
 * Allocate some memory for a new locus object
 *
 * @param[in] seqid    the sequence to which the locus belongs
 * @returns            a new locus object
 */
AgnLocus *agn_locus_new(const char *seqid);

/**
 * Print the locus in GFF3 format
 *
 * @param[in] locus        the locus
 * @param[in] outstream    the file to which the locus will be printed
 */
void agn_locus_print(AgnLocus *locus, FILE *outstream);

/**
 * Write a string representation of the locus (seq_start-end)
 *
 * @param[in]  locus     the locus
 * @param[out] string    the string to which the label will be written
 */
void agn_locus_stringify(AgnLocus *locus, char *string);

#endif
