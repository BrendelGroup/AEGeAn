#ifndef AEGEAN_GENE_VALIDATOR
#define AEGEAN_GENE_VALIDATOR

#include "genometools.h"
#include "AgnError.h"

/**
 * The purpose of the AgnGeneValidator class is to validate gene structure
 * annotations. This structure maintains lists of features needed for reference
 * when performing the validation.
 */
typedef struct AgnGeneValidator AgnGeneValidator;

/**
 * Free the memory previously used by the validator object.
 *
 * @param[in] v    the validator object to be deleted
 */
void agn_gene_validator_delete(AgnGeneValidator *v);

/**
 * Allocate some memory for the gene validator object.
 *
 * @returns               pointer to a validator object
 */
AgnGeneValidator* agn_gene_validator_new();

/**
 * Validate the given gene structure annotation.
 *
 * @param[in] v           validator object
 * @param[in] gene        gene feature to be validated
 * @param[in] settings    list of feature types to infer
 * @param[in] error       stores error messages in case of an error
 * @returns               a boolean indicating whether the gene is valid
 */
bool agn_gene_validator_validate_gene( AgnGeneValidator *v, GtFeatureNode *gene,
                                       AgnError *error );

#endif
