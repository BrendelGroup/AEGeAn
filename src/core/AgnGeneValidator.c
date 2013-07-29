#include <string.h>
#include "AgnGeneValidator.h"
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnGeneValidator
{
  GtArray *cdss;
  GtArray *exons;
  GtArray *introns;
  GtArray *utrs;
  GtArray *segments;
  GtRange start_codon_range;
  GtRange stop_codon_range;
  GtRange *left_codon_range;
  GtRange *right_codon_range;
};


//----------------------------------------------------------------------------//
// Prototypes for private method(s).
//----------------------------------------------------------------------------//

/**
 * Verify that the CDS features have correctly-defined phases.
 *
 * @param[in] v    validator object
 */
void agn_gene_validator_check_cds_phase(AgnGeneValidator *v);

/**
 * Determine whether any of the given features overlap with each other.
 *
 * @param[in]  v           validator object
 * @param[in]  features    list of features
 * @param[out] logger      object to which error/warning messages will be
 *                         written if necessary
 * @returns                false if any overlaps are detected, true otherwise
 */
bool agn_gene_validator_check_feature_overlap(AgnGeneValidator *v,
                                              GtArray *features,
                                              AgnLogger *logger);

/**
 * If exons and start/stop codons are provided explicitly, the (possibly)
 * disjoint segments(s) of the CDS can be inferred using this function.
 *
 * @param[in]  v        validator object
 * @param[in]  mrna     mRNA feature
 * @param[out] logger   object to which error/warning messages will be
 *                      written if necessary
 * @returns             true on success, false otherwise
 */
bool agn_gene_validator_infer_cds(AgnGeneValidator *v, GtFeatureNode *mrna,
                                  AgnLogger *logger);

/**
 * Infer exon features from explicitly defined CDS and UTR features.
 *
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] logger   object to which message will be written, if necessary
 * @returns            true if exons are successfully infered, false otherwise
 */
bool agn_gene_validator_infer_exons(AgnGeneValidator *v, GtFeatureNode *mrna,
                                    AgnLogger *logger);

/**
 * Infer intron features from explicitly defined exon features
 *
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] logger   object to which message will be written, if necessary
 * @returns            true if introns are successfully inferred, false
 *                     otherwise
 */
bool agn_gene_validator_infer_introns(AgnGeneValidator *v, GtFeatureNode *mrna,
                                      AgnLogger *logger);

/**
 * Determine which features have not been provided explicitly and infer these
 * from the features that were provided.
 *
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] logger   object to which message will be written, if necessary
 * @returns            true upon success, false otherwise
 */
bool agn_gene_validator_infer_mrna_features(AgnGeneValidator *v,
                                            GtFeatureNode *mrna,
                                            AgnLogger *logger);

/**
 * Infer UTR features from explicitly defined exon features
 *
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] logger   object to which messages will be written, if necessary
 * @returns            true if UTRs are successfully infered, false otherwise
 */
bool agn_gene_validator_infer_utrs(AgnGeneValidator *v, GtFeatureNode *mrna,
                                   AgnLogger *logger);

/**
 * Classify and store the given mRNA's subfeature, or delete it if it does not
 * have a valid feature type.
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] child    a subfeature of the mRNA
 * @param[in] logger   object to which messages will be written, if necessary
 */
void agn_gene_validator_mrna_typecheck(AgnGeneValidator *v, GtFeatureNode *mrna,
                                       GtFeatureNode *child, AgnLogger *logger);

/**
 * Reset the validator object to be ready to validate the next feature(s).
 *
 * @param[in] v        validator object
 */
void agn_gene_validator_reset(AgnGeneValidator *v);

/**
 * Make sure that each UTR is explicitly designated as 3' or 5'.
 *
 * @param[in] v    validator object
 */
void agn_gene_validator_set_utrs(AgnGeneValidator *v);

/**
 * Sort all of the currently stored mRNA subfeatures.
 *
 * @param[in] v    validator object
 */
void agn_gene_validator_sort_features(AgnGeneValidator *v);

/**
 * Validate the given feature as an mRNA
 *
 * @param[in] v           validator object
 * @param[in] gene        the gene from which the mRNA is derived
 * @param[in] mRNA        the feature to be validated
 * @param[in] logger      object to which messages will be written, ifnecessary
 * @returns               true if the feature is successfully validated as an
 *                        mRNA, false otherwise
 */
bool agn_gene_validator_validate_mrna(AgnGeneValidator *v,
                                      GtFeatureNode *gene,
                                      GtFeatureNode *mRNA,
                                      AgnLogger *logger);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

void agn_gene_validator_check_cds_phase(AgnGeneValidator *v)
{
  unsigned long num_cds_feats = gt_array_size(v->cdss);
  if(num_cds_feats == 0)
    return;

  GtFeatureNode *cdsf1 = *(GtFeatureNode **)gt_array_get(v->cdss, 0);
  gt_feature_node_set_phase(cdsf1, GT_PHASE_ZERO);
  if(num_cds_feats == 1)
    return;

  GtStrand strand = gt_feature_node_get_strand(cdsf1);
  unsigned long cds_length = gt_genome_node_get_length((GtGenomeNode *)cdsf1);
  unsigned long i;
  if(strand == GT_STRAND_REVERSE)
  {
    for(i = num_cds_feats - 1; i > 0; i--)
    {
      GtFeatureNode *cds = *(GtFeatureNode **)gt_array_get(v->cdss, i);
      int phasenum = cds_length % 3;
      GtPhase phase = GT_PHASE_ZERO;
      if(phasenum == 1)
        phase = GT_PHASE_TWO;
      else if(phasenum == 2)
        phase = GT_PHASE_ONE;
      gt_feature_node_set_phase(cds, phase);
      cds_length += gt_genome_node_get_length((GtGenomeNode *)cds);
    }
  }
  else
  {
    for(i = 1; i < num_cds_feats; i++)
    {
      GtFeatureNode *cds = *(GtFeatureNode **)gt_array_get(v->cdss, i);
      int phasenum = cds_length % 3;
      GtPhase phase = GT_PHASE_ZERO;
      if(phasenum == 1)
        phase = GT_PHASE_TWO;
      else if(phasenum == 2)
        phase = GT_PHASE_ONE;
      gt_feature_node_set_phase(cds, phase);
      cds_length += gt_genome_node_get_length((GtGenomeNode *)cds);
    }
  }
}

bool agn_gene_validator_check_feature_overlap(AgnGeneValidator *v,
                                              GtArray *features,
                                              AgnLogger *logger)
{
  gt_assert(features != NULL);
  if(gt_array_size(features) <= 1)
    return true;

  bool isvalid = true;
  unsigned long i;
  for(i = 0; i < gt_array_size(features) - 1; i++)
  {
    GtFeatureNode *fni = *(GtFeatureNode **)gt_array_get(features, i);
    GtFeatureNode *fnj = *(GtFeatureNode **)gt_array_get(features, i+1);
    GtRange ri = gt_genome_node_get_range((GtGenomeNode *)fni);
    GtRange rj = gt_genome_node_get_range((GtGenomeNode *)fnj);
    if(gt_range_overlap(&ri, &rj))
    {
      const char *type = gt_feature_node_get_type(fni);
      const char *tid = gt_feature_node_get_attribute(fni, "Parent");
      agn_logger_log_error(logger, "overlapping %ss for mRNA '%s'", type, tid);
      isvalid = false;
    }
  }
  return isvalid;
}

void agn_gene_validator_delete(AgnGeneValidator *v)
{
  gt_free(v);
  v = NULL;
}

bool agn_gene_validator_infer_cds(AgnGeneValidator *v, GtFeatureNode *mrna,
                                  AgnLogger *logger)
{
  bool left_codon_found = v->left_codon_range->start != 0;
  bool right_codon_found  = v->right_codon_range->start  != 0;
  unsigned long num_exons = gt_array_size(v->exons);
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");
  unsigned int tln = gt_genome_node_get_line_number((GtGenomeNode *)mrna);
  bool isvalid = true;

  gt_assert(left_codon_found && right_codon_found && num_exons > 0);
  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  unsigned long i;
  for(i = 0; i < num_exons; i++)
  {
    GtFeatureNode *exon = *(GtFeatureNode **)gt_array_get(v->exons, i);
    GtGenomeNode *exon_gn = (GtGenomeNode *)exon;
    GtRange exon_range = gt_genome_node_get_range(exon_gn);
    GtStrand exon_strand = gt_feature_node_get_strand(exon);

    GtRange cdsrange;
    bool exon_includes_cds = agn_infer_cds_range_from_exon_and_codons(
                                 &exon_range, v->left_codon_range,
                                 v->right_codon_range, &cdsrange);
    if(exon_includes_cds)
    {
      GtGenomeNode *cds = gt_feature_node_new
      (
        gt_genome_node_get_seqid(exon_gn), "CDS", cdsrange.start,
        cdsrange.end, exon_strand
      );
      gt_feature_node_set_source((GtFeatureNode *)cds, aegean);
      gt_feature_node_add_child(mrna, (GtFeatureNode *)cds);
      gt_feature_node_add_attribute((GtFeatureNode *)cds, "Parent", tid);
      gt_array_add(v->cdss, cds);
    }
  }
  gt_str_delete(aegean);

  if(gt_array_size(v->cdss) == 0)
  {
    agn_logger_log_error(logger, "error inferring CDS from codons and exons "
                         "for mRNA '%s' (line %u)", tid, tln);
    isvalid = false;
  }

  return isvalid;
}

bool agn_gene_validator_infer_exons(AgnGeneValidator *v, GtFeatureNode *mrna,
                                    AgnLogger *logger)
{
  unsigned long i,j;
  unsigned long num_exons = gt_array_size(v->exons);
  unsigned long num_cds_segs = gt_array_size(v->cdss);
  gt_assert(num_exons == 0 && num_cds_segs > 0);
  GtHashmap *adjacent_utrs = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  GtArray *exons_to_add = gt_array_new( sizeof(GtRange) );

  GtGenomeNode *n = (GtGenomeNode *)mrna;
  unsigned int tln = gt_genome_node_get_line_number(n);
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");

  for(i = 0; i < num_cds_segs; i++)
  {
    GtGenomeNode *cds = *(GtGenomeNode **)gt_array_get(v->cdss, i);
    GtRange cds_range = gt_genome_node_get_range(cds);
    GtRange exon_range = cds_range;
    for(j = 0; j < gt_array_size(v->utrs); j++)
    {
      GtGenomeNode *utr = *(GtGenomeNode **)gt_array_get(v->utrs, j);
      GtRange utr_range = gt_genome_node_get_range(utr);

      // If the UTR is adjacent to the CDS, merge the ranges
      if( utr_range.end+1 == cds_range.start ||
          cds_range.end+1 == utr_range.start )
      {
        exon_range = gt_range_join(&exon_range, &utr_range);
        gt_hashmap_add(adjacent_utrs, utr, utr);
      }
    }

    gt_array_add(exons_to_add, exon_range);
  }

  // Now create UTR-only exons
  for(i = 0; i < gt_array_size(v->utrs); i++)
  {
    GtGenomeNode *utr = *(GtGenomeNode **)gt_array_get(v->utrs, i);
    GtRange utr_range = gt_genome_node_get_range(utr);
    if(gt_hashmap_get(adjacent_utrs, utr) == NULL)
    {
      gt_array_add(exons_to_add, utr_range);
    }
  }

  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  for(i = 0; i < gt_array_size(exons_to_add); i++)
  {
    GtRange *exon_range = (GtRange *)gt_array_get(exons_to_add, i);
    GtGenomeNode *firstcds = *(GtGenomeNode **)gt_array_get(v->cdss, 0);
    GtGenomeNode *exon = gt_feature_node_new
    (
      gt_genome_node_get_seqid(firstcds), "exon", exon_range->start,
      exon_range->end, gt_feature_node_get_strand((GtFeatureNode *)firstcds)
    );
    GtFeatureNode *fn_exon = (GtFeatureNode *)exon;
    gt_feature_node_set_source(fn_exon, aegean);
    gt_feature_node_add_child((GtFeatureNode *)mrna, fn_exon);
    gt_feature_node_add_attribute(fn_exon, "Parent", tid);
    gt_array_add(v->exons, exon);
  }
  gt_str_delete(aegean);
  gt_array_delete(exons_to_add);

  if(gt_array_size(v->exons) == 0)
  {
    agn_logger_log_error(logger,"unable to infer exons for mRNA '%s' (line %u)",
                         tid,tln);
    return false;
  }
  gt_hashmap_delete(adjacent_utrs);

  // Sort and return
  gt_array_sort(v->exons, (GtCompare)agn_gt_genome_node_compare);
  return true;
}

bool agn_gene_validator_infer_introns(AgnGeneValidator *v, GtFeatureNode *mrna,
                                      AgnLogger *logger)
{
  gt_assert(gt_array_size(v->introns) == 0 && gt_array_size(v->exons) > 0);

  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  unsigned long i;
  for(i = 1; i < gt_array_size(v->exons); i++)
  {
    GtGenomeNode *exon1 = *(GtGenomeNode**)gt_array_get(v->exons, i-1);
    GtGenomeNode *exon2 = *(GtGenomeNode**)gt_array_get(v->exons, i);
    GtRange first_range = gt_genome_node_get_range(exon1);
    GtRange second_range = gt_genome_node_get_range(exon2);

    if(first_range.end == second_range.start - 1)
    {
      agn_logger_log_error(logger, "mRNA '%s' has directly adjacent exons",
                           gt_feature_node_get_attribute(mrna, "ID"));
      gt_str_delete(aegean);
      return false;
    }
    else
    {
      GtStrand intronstrnd = gt_feature_node_get_strand((GtFeatureNode *)exon1);
      GtGenomeNode *intron = gt_feature_node_new
      (
        gt_genome_node_get_seqid(exon1), "intron", first_range.end + 1,
        second_range.start - 1, intronstrnd
      );
      gt_feature_node_set_source((GtFeatureNode *)intron, aegean);
      gt_feature_node_add_child((GtFeatureNode *)mrna, (GtFeatureNode *)intron);
      gt_feature_node_add_attribute((GtFeatureNode *)intron, "Parent",
                                    gt_feature_node_get_attribute(mrna, "ID"));
      gt_array_add(v->introns, intron);
    }
  }
  gt_str_delete(aegean);

  gt_array_sort(v->introns, (GtCompare)agn_gt_genome_node_compare);
  return true;
}

bool agn_gene_validator_infer_mrna_features(AgnGeneValidator *v,
                                            GtFeatureNode *mrna,
                                            AgnLogger *logger)
{
  bool isvalid = true;

  // If CDS segments are not explicitly provided, infer CDS from exons and
  // start/stop codons
  bool left_codon_found = v->left_codon_range->start != 0;
  bool right_codon_found  = v->right_codon_range->start  != 0;
  unsigned long num_exons = gt_array_size(v->exons);
  if(gt_array_size(v->cdss) == 0 && left_codon_found && right_codon_found &&
     num_exons > 0)
  {
    bool success = agn_gene_validator_infer_cds(v, mrna, logger);
    if(!success)
      isvalid = false;
  }
  agn_gene_validator_check_feature_overlap(v, v->cdss, logger);
  agn_gene_validator_check_cds_phase(v);

  // If exons are not explicitly provided, infer them from CDS and UTR segments
  if(isvalid && gt_array_size(v->exons) == 0 && gt_array_size(v->cdss) > 0)
  {
    bool success = agn_gene_validator_infer_exons(v, mrna, logger);
    if(!success)
      isvalid = false;
  }
  agn_gene_validator_check_feature_overlap(v, v->exons, logger);

  // If exons and CDS are provided explicitly, but UTR segments are
  // not, infer the UTRs
  if(isvalid && v->start_codon_range.start != 0 && gt_array_size(v->exons) > 0
     && gt_array_size(v->utrs) == 0)
  {
    bool success = agn_gene_validator_infer_utrs(v, mrna, logger);
    if(!success)
      isvalid = false;
  }
  agn_gene_validator_check_feature_overlap(v, v->utrs, logger);
  agn_gene_validator_set_utrs(v);

  // If there are multiple exons and introns are not explicitly provided, infer
  // introns from the exon features
  if(isvalid && gt_array_size(v->exons) > 1 && gt_array_size(v->introns) == 0)
  {
    bool success = agn_gene_validator_infer_introns(v, mrna, logger);
    if(!success)
      isvalid = false;
  }
  if(isvalid && gt_array_size(v->exons) != gt_array_size(v->introns) + 1)
  {
    agn_logger_log_error(logger, "mRNA '%s' (line %u) has %lu exons but %lu "
        "introns", gt_feature_node_get_attribute(mrna, "ID"),
        gt_genome_node_get_line_number((GtGenomeNode *)mrna),
        gt_array_size(v->exons), gt_array_size(v->introns)
    );
    isvalid = false;
  }

  if(gt_array_size(v->cdss) == 0)
  {
    agn_logger_log_warning(logger, "mRNA '%s' (line %u) has no CDS",
        gt_feature_node_get_attribute(mrna, "ID"),
        gt_genome_node_get_line_number((GtGenomeNode *)mrna)
    );
    isvalid = false;
  }

  return isvalid;
}

bool agn_gene_validator_infer_utrs(AgnGeneValidator *v, GtFeatureNode *mrna,
                                   AgnLogger *logger)
{
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");
  GtStrand strand = gt_feature_node_get_strand(mrna);
  const char *lefttype  = "five_prime_UTR";
  const char *righttype = "three_prime_UTR";
  if(gt_feature_node_get_strand(mrna) == GT_STRAND_REVERSE)
  {
    lefttype  = "three_prime_UTR";
    righttype = "five_prime_UTR";
  }

  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  unsigned long i;
  for(i = 0; i < gt_array_size(v->exons); i++)
  {
    GtFeatureNode *exon = *(GtFeatureNode **)gt_array_get(v->exons, i);
    GtRange exonrange = gt_genome_node_get_range((GtGenomeNode *)exon);

    if(exonrange.start < v->left_codon_range->start)
    {
      GtRange utrrange;
      if(gt_range_overlap(&exonrange, v->left_codon_range))
      {
        utrrange.start = exonrange.start;
        utrrange.end   = v->left_codon_range->start - 1;
      }
      else
      {
        utrrange = exonrange;
      }
      GtGenomeNode *utr = gt_feature_node_new(
                              gt_genome_node_get_seqid((GtGenomeNode *)exon),
                              lefttype, utrrange.start, utrrange.end, strand);
      gt_feature_node_set_source((GtFeatureNode *)utr, aegean);
      gt_feature_node_add_child(mrna, (GtFeatureNode *)utr);
      gt_feature_node_add_attribute((GtFeatureNode *)utr, "Parent", tid);
      gt_array_add(v->utrs, utr);
    }

    if(exonrange.end > v->right_codon_range->end)
    {
      GtRange utrrange;
      if(gt_range_overlap(&exonrange, v->right_codon_range))
      {
        utrrange.start = v->right_codon_range->end + 1;
        utrrange.end   = exonrange.end;
      }
      else
      {
        utrrange = exonrange;
      }
      GtGenomeNode *utr = gt_feature_node_new(
                              gt_genome_node_get_seqid((GtGenomeNode *)exon),
                              righttype, utrrange.start, utrrange.end, strand);
      gt_feature_node_set_source((GtFeatureNode *)utr, aegean);
      gt_feature_node_add_child(mrna, (GtFeatureNode *)utr);
      gt_feature_node_add_attribute((GtFeatureNode *)utr, "Parent", tid);
      gt_array_add(v->utrs, utr);
    }
  }
  gt_str_delete(aegean);

  return true;
}

AgnGeneValidator* agn_gene_validator_new()
{
  AgnGeneValidator *v = gt_malloc(sizeof(AgnGeneValidator));
  v->cdss    = NULL;
  v->exons   = NULL;
  v->introns = NULL;
  v->utrs    = NULL;

  v->start_codon_range.start = 0;
  v->start_codon_range.end   = 0;
  v->stop_codon_range.start  = 0;
  v->stop_codon_range.end    = 0;

  v->left_codon_range  = NULL;
  v->right_codon_range = NULL;

  return v;
}

void agn_gene_validator_mrna_typecheck(AgnGeneValidator *v, GtFeatureNode *mrna,
                                       GtFeatureNode *child, AgnLogger *logger)
{
  if(agn_gt_feature_node_is_cds_feature(child))
  {
    gt_array_add(v->cdss, child);
    GtRange cdsrange = gt_genome_node_get_range((GtGenomeNode *)child);
    if(v->left_codon_range->start == 0 ||
       cdsrange.start < v->left_codon_range->start)
    {
      v->left_codon_range->start = cdsrange.start;
      v->left_codon_range->end   = cdsrange.start + 2;
    }
    if(v->right_codon_range->end == 0 ||
       cdsrange.end > v->right_codon_range->end)
    {
      v->right_codon_range->end   = cdsrange.end;
      v->right_codon_range->start = cdsrange.end - 2;
    }
  }

  else if(agn_gt_feature_node_is_start_codon_feature(child))
  {
    v->start_codon_range = gt_genome_node_get_range((GtGenomeNode *)child);
    gt_feature_node_remove_leaf(mrna, child);
  }

  else if(agn_gt_feature_node_is_stop_codon_feature(child))
  {
    v->stop_codon_range = gt_genome_node_get_range((GtGenomeNode *)child);
    gt_feature_node_remove_leaf(mrna, child);
  }

  else if(agn_gt_feature_node_is_exon_feature(child))
    gt_array_add(v->exons, child);

  else if(agn_gt_feature_node_is_intron_feature(child))
    gt_array_add(v->introns, child);

  else if(agn_gt_feature_node_is_utr_feature(child))
    gt_array_add(v->utrs, child);

  else if(agn_gt_feature_node_is_mrna_feature(child))
  {
    if(child != mrna)
    {
      agn_logger_log_warning(logger, "%s feature (line %u) is not a valid child"
          " for mRNA '%s', ignoring", gt_feature_node_get_type(child),
          gt_genome_node_get_line_number((GtGenomeNode *)child),
          gt_feature_node_get_attribute(mrna, "ID")
      );
      gt_feature_node_remove_leaf(mrna, child);
      gt_genome_node_delete((GtGenomeNode *)child);
    }
  }
  else
  {
    agn_logger_log_warning(logger, "%s feature (line %u) is not a valid child "
        "for mRNA '%s', ignoring", gt_feature_node_get_type(child),
        gt_genome_node_get_line_number((GtGenomeNode *)child),
        gt_feature_node_get_attribute(mrna, "ID")
    );
    gt_feature_node_remove_leaf(mrna, child);
  }
}

void agn_gene_validator_reset(AgnGeneValidator *v)
{
  gt_array_delete(v->cdss);
  gt_array_delete(v->exons);
  gt_array_delete(v->introns);
  gt_array_delete(v->utrs);
  gt_array_delete(v->segments);

  v->start_codon_range.start = 0;
  v->start_codon_range.end   = 0;
  v->stop_codon_range.start  = 0;
  v->stop_codon_range.end    = 0;
  v->left_codon_range  = NULL;
  v->right_codon_range = NULL;
}

void agn_gene_validator_set_utrs(AgnGeneValidator *v)
{
  GtGenomeNode *gn;
  unsigned long i, cds_start, num_utrs;

  num_utrs = gt_array_size(v->utrs);
  if(num_utrs == 0)
    return;

  gn = *(GtGenomeNode **)gt_array_get(v->cdss, 0);
  cds_start = gt_genome_node_get_start(gn);

  for(i = 0; i < num_utrs; i++)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_get(v->utrs, i);
    gn = (GtGenomeNode *)fn;

    if( !gt_feature_node_has_type(fn, "five_prime_UTR") &&
        !gt_feature_node_has_type(fn, "three_prime_UTR") )
    {
      GtStrand strand = gt_feature_node_get_strand(fn);
      unsigned long utr_start = gt_genome_node_get_start(gn);

      if(strand == GT_STRAND_FORWARD)
      {
        if(utr_start < cds_start)
          gt_feature_node_set_type(fn, "five_prime_UTR");
        else
          gt_feature_node_set_type(fn, "three_prime_UTR");
      }
      else
      {
        if(utr_start < cds_start)
          gt_feature_node_set_type(fn, "three_prime_UTR");
        else
          gt_feature_node_set_type(fn, "five_prime_UTR");
      }
    }
  }
}

void agn_gene_validator_sort_features(AgnGeneValidator *v)
{
  gt_array_sort(v->cdss,     (GtCompare)agn_gt_genome_node_compare);
  gt_array_sort(v->exons,    (GtCompare)agn_gt_genome_node_compare);
  gt_array_sort(v->introns,  (GtCompare)agn_gt_genome_node_compare);
  gt_array_sort(v->utrs,     (GtCompare)agn_gt_genome_node_compare);
  gt_array_sort(v->segments, (GtCompare)agn_gt_genome_node_compare);
}

bool agn_gene_validator_validate_gene(AgnGeneValidator *v, GtFeatureNode *gene,
                                      AgnLogger *logger)
{
  if(!agn_gt_feature_node_is_gene_feature(gene))
  {
    agn_logger_log_warning(logger, "feature from line %u is a %s, expected a "
        "gene", gt_genome_node_get_line_number((GtGenomeNode *)gene),
        gt_feature_node_get_type(gene)
    );
    return false;
  }

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
  GtFeatureNode *childfeature;
  bool isvalid = true;
  unsigned int num_mrnas = 0;
  for(childfeature = gt_feature_node_iterator_next(iter);
      childfeature != NULL;
      childfeature = gt_feature_node_iterator_next(iter))
  {
    if(agn_gene_validator_validate_mrna(v, gene, childfeature, logger))
    {
      GtRange grange = gt_genome_node_get_range((GtGenomeNode *)gene);
      GtRange trange = gt_genome_node_get_range((GtGenomeNode *)childfeature);
      if(!gt_range_contains(&grange, &trange))
      {
        agn_logger_log_error(logger, "range of mRNA '%s' exceeds range of gene "
            "'%s'", gt_feature_node_get_attribute(childfeature, "ID"),
            gt_feature_node_get_attribute(gene, "ID")
        );
        isvalid = false;
      }

      num_mrnas++;
    }
    else
    {
      // gt_assert(agn_logger_has_error(logger));
      isvalid = false;
    }
  }
  gt_feature_node_iterator_delete(iter);

  if(num_mrnas == 0)
  {
    agn_logger_log_warning(logger, "found no valid mRNAs for gene '%s' "
        "(line %u)", gt_feature_node_get_attribute(gene, "ID"),
        gt_genome_node_get_line_number((GtGenomeNode *)gene)
    );
    return false;
  }

  return isvalid;
}

bool agn_gene_validator_validate_mrna(AgnGeneValidator *v, GtFeatureNode *gene,
                                      GtFeatureNode *mrna, AgnLogger *logger)
{
  if(!agn_gt_feature_node_is_mrna_feature(mrna))
  {
    agn_logger_log_warning(logger, "feature from line %u is a %s, expected an "
        "mRNA", gt_genome_node_get_line_number((GtGenomeNode *)mrna),
        gt_feature_node_get_type(mrna)
    );
    return false;
  }

  v->cdss     = gt_array_new(sizeof(GtFeatureNode *));
  v->exons    = gt_array_new(sizeof(GtFeatureNode *));
  v->introns  = gt_array_new(sizeof(GtFeatureNode *));
  v->utrs     = gt_array_new(sizeof(GtFeatureNode *));
  v->segments = gt_array_new(sizeof(GtFeatureNode *));

  v->left_codon_range  = &v->start_codon_range;
  v->right_codon_range = &v->stop_codon_range;
  if(gt_feature_node_get_strand(mrna) == GT_STRAND_REVERSE)
  {
    v->left_codon_range  = &v->stop_codon_range;
    v->right_codon_range = &v->start_codon_range;
  }

  bool isvalid = true;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(mrna);
  GtFeatureNode *childfeature;
  for(childfeature = gt_feature_node_iterator_next(iter);
      childfeature != NULL;
      childfeature = gt_feature_node_iterator_next(iter))
  {
    if(!agn_gt_feature_node_fix_parent_attribute(childfeature, mrna))
    {
      agn_logger_log_error(logger, "could not fix 'Parent' attribute for %s "
          "feature (child of %s '%s')", gt_feature_node_get_type(childfeature),
          gt_feature_node_get_type(mrna),
          gt_feature_node_get_attribute(mrna, "ID")
      );
      isvalid = false;
    }
    if(!agn_gt_feature_node_range_contains(mrna, childfeature))
    {
      GtRange crange = gt_genome_node_get_range((GtGenomeNode *)childfeature);
      agn_logger_log_error(logger, "%s at [%lu, %lu] exceeds range of its "
          "parent mRNA '%s'", gt_feature_node_get_type(childfeature),
          crange.start, crange.end,
          gt_feature_node_get_attribute(childfeature, "Parent")
      );
      isvalid = false;
    }
    agn_gene_validator_mrna_typecheck(v, mrna, childfeature, logger);
  }
  gt_feature_node_iterator_delete(iter);

  if(isvalid)
  {
    agn_gene_validator_sort_features(v);
    int success = agn_gene_validator_infer_mrna_features(v, mrna, logger);
    if(!success)
      isvalid = false;
  }

  agn_gene_validator_reset(v);
  return isvalid;
}
