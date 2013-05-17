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
  GtRange start_codon_range;
  GtRange stop_codon_range;
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
 * Infer exon features from explicitly defined CDS and UTR features.
 *
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] error    error object to which message will be written, if
 *                     necessary
 * @returns            true if exons are successfully infered, false otherwise
 */
bool agn_gene_validator_infer_exons( AgnGeneValidator *v, GtFeatureNode *mrna,
                                     AgnError *error );

/**
 * Infer intron features from explicitly defined exon features
 *
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] error    error object to which message will be written, if
 *                     necessary
 * @returns            true if introns are successfully infered, false otherwise
 */
bool agn_gene_validator_infer_introns( AgnGeneValidator *v, GtFeatureNode *mrna,
                                       AgnError *error );

/**
 * Infer UTR features from explicitly defined exon features
 *
 * @param[in] v        validator object
 * @param[in] mrna     an mRNA feature
 * @param[in] error    error object to which messages will be written, if
 *                     necessary
 * @returns            true if UTRs are successfully infered, false otherwise
 */
bool agn_gene_validator_infer_utrs( AgnGeneValidator *v, GtFeatureNode *mrna,
                                    AgnError *error );

/**
 * Make sure that each UTR is explicitly designated as 3' or 5'
 *
 * @param[in] v    validator object
 */
void agn_gene_validator_set_utrs(AgnGeneValidator *v);

/**
 * Validate the given feature as an mRNA
 *
 * @param[in] v           validator object
 * @param[in] gene        the gene from which the mRNA is derived
 * @param[in] mRNA        the feature to be validated
 * @param[in] error       error object to which messages will be written, if
 *                        necessary
 * @returns               true if the feature is successfully validated as an
 *                        mRNA, false otherwise
 */
bool agn_gene_validator_validate_mrna( AgnGeneValidator *v,
                                       GtFeatureNode *gene,
                                       GtFeatureNode *mRNA,
                                       AgnError *error );


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

void agn_gene_validator_check_cds_phase(AgnGeneValidator *v)
{
  unsigned long num_cds_feats = gt_array_size(v->cdss);
  gt_assert(num_cds_feats > 0);

  GtFeatureNode *first_cds_feature = *(GtFeatureNode **)gt_array_get(v->cdss, 0);
  gt_feature_node_set_phase(first_cds_feature, GT_PHASE_ZERO);
  if(num_cds_feats == 1)
    return;

  GtStrand strand = gt_feature_node_get_strand(first_cds_feature);
  unsigned long cds_length = gt_genome_node_get_length((GtGenomeNode *)first_cds_feature);
  unsigned long i;
  if(strand == GT_STRAND_REVERSE)
  {
    for(i = num_cds_feats - 1; i > 0; i--)
    {
      GtFeatureNode *current_feature = *(GtFeatureNode **)gt_array_get(v->cdss, i);
      int phasenum = cds_length % 3;
      GtPhase phase = GT_PHASE_ZERO;
      if(phasenum == 1)
        phase = GT_PHASE_TWO;
      else if(phasenum == 2)
        phase = GT_PHASE_ONE;
      gt_feature_node_set_phase(current_feature, phase);
      cds_length += gt_genome_node_get_length((GtGenomeNode *)current_feature);
    }
  }
  else
  {
    for(i = 1; i < num_cds_feats; i++)
    {
      GtFeatureNode *current_feature = *(GtFeatureNode **)gt_array_get(v->cdss, i);
      int phasenum = cds_length % 3;
      GtPhase phase = GT_PHASE_ZERO;
      if(phasenum == 1)
        phase = GT_PHASE_TWO;
      else if(phasenum == 2)
        phase = GT_PHASE_ONE;
      gt_feature_node_set_phase(current_feature, phase);
      cds_length += gt_genome_node_get_length((GtGenomeNode *)current_feature);
    }
  }
}

void agn_gene_validator_delete(AgnGeneValidator *v)
{
  gt_free(v);
  v = NULL;
}

AgnGeneValidator* agn_gene_validator_new()
{
  AgnGeneValidator *v = gt_malloc( sizeof(AgnGeneValidator) );
  v->cdss    = NULL;
  v->exons   = NULL;
  v->introns = NULL;
  v->utrs    = NULL;

  v->start_codon_range.start = 0;
  v->start_codon_range.end   = 0;
  v->stop_codon_range.start  = 0;
  v->stop_codon_range.end    = 0;

  return v;
}

bool agn_gene_validator_infer_exons( AgnGeneValidator *v, GtFeatureNode *mrna,
                                     AgnError *error )
{
  unsigned long i,j;
  unsigned long num_exons = gt_array_size(v->exons);
  unsigned long num_cds_segs = gt_array_size(v->cdss);
  gt_assert(num_exons == 0 && num_cds_segs > 0);

  GtGenomeNode *n = (GtGenomeNode *)mrna;
  unsigned int tln = gt_genome_node_get_line_number(n);
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");

  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  for(i = 0; i < num_cds_segs; i++)
  {
    GtGenomeNode *cds = *(GtGenomeNode **)gt_array_get(v->cdss, i);
    GtRange cds_range = gt_genome_node_get_range(cds);
    GtRange exon_range = cds_range;

    // For each CDS, iterate through all UTR segments
    for(j = 0; j < gt_array_size(v->utrs); j++)
    {
      GtGenomeNode *utr = *(GtGenomeNode **)gt_array_get(v->utrs, j);
      GtRange utr_range = gt_genome_node_get_range(utr);

      // If the UTR is adjacent to the CDS, merge the ranges
      if( utr_range.end+1 == cds_range.start ||
          cds_range.end+1 == utr_range.start )
      {
        exon_range = gt_range_join(&exon_range, &utr_range);
      }
    }

    // Store the exon feature
    GtStrand exonstrand = gt_feature_node_get_strand((GtFeatureNode *)cds);
    GtGenomeNode *exon = gt_feature_node_new
    (
      gt_genome_node_get_seqid(cds), "exon", exon_range.start,
      exon_range.end, exonstrand
    );
    GtFeatureNode *fn_exon = (GtFeatureNode *)exon;
    gt_feature_node_set_source(fn_exon, aegean);
    gt_feature_node_add_child((GtFeatureNode *)mrna, fn_exon);
    gt_feature_node_add_attribute(fn_exon, "Parent", tid);
    gt_array_add(v->exons, exon);
  }
  gt_str_delete(aegean);
  if(gt_array_size(v->exons) == 0)
  {
    agn_error_add( error,false, "unable to infer exons for mRNA '%s' (line %u)",
                   tid, tln );
    return false;
  }

  // Sort and return
  gt_array_sort(v->exons, (GtCompare)agn_gt_genome_node_compare);
  return true;
}

bool agn_gene_validator_infer_introns( AgnGeneValidator *v, GtFeatureNode *mrna,
                                       AgnError *error )
{
  unsigned long i;
  unsigned long num_introns = gt_array_size(v->introns);
  unsigned long num_exons = gt_array_size(v->exons);
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");
  gt_assert(num_introns == 0 && num_exons > 0);

  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  for(i = 1; i < gt_array_size(v->exons); i++)
  {
    GtGenomeNode *exon1 = *(GtGenomeNode**)gt_array_get(v->exons, i-1);
    GtGenomeNode *exon2 = *(GtGenomeNode**)gt_array_get(v->exons, i);
    GtRange first_range = gt_genome_node_get_range(exon1);
    GtRange second_range = gt_genome_node_get_range(exon2);

    if(first_range.end == second_range.start - 1)
    {
      agn_error_add(error, false, "mRNA '%s' has directly adjacent exons\n",
                    tid);
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
      GtFeatureNode *fn_intron = (GtFeatureNode *)intron;
      gt_feature_node_set_source(fn_intron, aegean);
      gt_feature_node_add_child((GtFeatureNode *)mrna, fn_intron);
      gt_feature_node_add_attribute(fn_intron, "Parent", tid);
      gt_array_add(v->introns, intron);
    }
  }
  gt_str_delete(aegean);

  gt_array_sort(v->introns, (GtCompare)agn_gt_genome_node_compare);
  return true;
}

// Function needs a prototype
void agn_gene_validator_infer_utr_from_partial_exon(AgnGeneValidator *v,
                                                    GtFeatureNode *mrna,
                                                    GtFeatureNode *exon,
                                                    GtFeatureNode *cds)
{
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");
  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  GtRange exon_range = gt_genome_node_get_range((GtGenomeNode *)exon);
  GtRange cds_range = gt_genome_node_get_range((GtGenomeNode *)cds);
  GtStrand exonstrand = gt_feature_node_get_strand(exon);
  
  if(exon_range.start < cds_range.start)
  {
    char *utr_type = "five_prime_UTR";
    if(exonstrand == GT_STRAND_REVERSE)
      utr_type = "three_prime_UTR";

    GtGenomeNode *utr = gt_feature_node_new(
                            gt_genome_node_get_seqid((GtGenomeNode *)exon),
                            utr_type, exon_range.start, cds_range.start - 1,
                            exonstrand);
    gt_feature_node_set_source((GtFeatureNode *)utr, aegean);
    gt_feature_node_add_child(mrna, (GtFeatureNode *)utr);
    gt_feature_node_add_attribute((GtFeatureNode *)utr, "Parent", tid);
    gt_array_add(v->utrs, utr);
  }
  if(exon_range.end > cds_range.end)
  {
    char *utr_type = "three_prime_UTR";
    if(exonstrand == GT_STRAND_REVERSE)
      utr_type = "five_prime_UTR";

    GtGenomeNode *utr = gt_feature_node_new(
                          gt_genome_node_get_seqid((GtGenomeNode *)exon),
                          utr_type, cds_range.end + 1, exon_range.end,
                          exonstrand);
    gt_feature_node_set_source((GtFeatureNode *)utr, aegean);
    gt_feature_node_add_child(mrna, (GtFeatureNode *)utr);
    gt_feature_node_add_attribute((GtFeatureNode *)utr, "Parent", tid);
    gt_array_add(v->utrs, utr);
  }
  gt_str_delete(aegean);
}

// Function needs a prototype
void agn_gene_validator_infer_utr_from_full_exon(AgnGeneValidator *v,
                                                 GtFeatureNode *mrna,
                                                 GtFeatureNode *exon,
                                                 GtFeatureNode *cds)
{
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");
  GtStr *aegean = gt_str_new_cstr("AEGeAn");
  GtRange exon_range = gt_genome_node_get_range((GtGenomeNode *)exon);
  GtRange cds_range = gt_genome_node_get_range((GtGenomeNode *)cds);
  char *utr_type = "five_prime_UTR";
  GtStrand exon_strand = gt_feature_node_get_strand(exon);
  bool is_3prime_rev = ( exon_strand == GT_STRAND_REVERSE &&
                         exon_range.start < cds_range.start );
  bool is_3prime_fwd = ( exon_strand == GT_STRAND_FORWARD &&
                         exon_range.end > cds_range.end );
  if(is_3prime_rev || is_3prime_fwd)
  {
    utr_type = "three_prime_UTR";
  }

  GtGenomeNode *utr = gt_feature_node_new(
                        gt_genome_node_get_seqid((GtGenomeNode *)exon),
                        utr_type, exon_range.start,exon_range.end, exon_strand);
  gt_feature_node_set_source((GtFeatureNode *)utr, aegean);
  gt_feature_node_add_child(mrna, (GtFeatureNode *)utr);
  gt_feature_node_add_attribute((GtFeatureNode *)utr, "Parent", tid);
  gt_array_add(v->utrs, utr);
  gt_str_delete(aegean);
}

bool agn_gene_validator_infer_utrs( AgnGeneValidator *v, GtFeatureNode *mrna,
                                    AgnError *error )
{
  unsigned long exon_index = 0;
  unsigned long cds_index = 0;
  unsigned long num_exons = gt_array_size(v->exons);
  unsigned long num_cds_segs = gt_array_size(v->cdss);
  gt_assert(num_exons > 0 && num_cds_segs > 0);

  while(exon_index < num_exons)
  {
    if(cds_index >= gt_array_size(v->cdss))
      cds_index = gt_array_size(v->cdss) - 1;
    GtFeatureNode *exon = *(GtFeatureNode **)gt_array_get(v->exons, exon_index);
    GtFeatureNode *cds  = *(GtFeatureNode **)gt_array_get(v->cdss, cds_index);
    GtRange exon_range = gt_genome_node_get_range((GtGenomeNode *)exon);
    GtRange cds_range = gt_genome_node_get_range((GtGenomeNode *)cds);

    if(gt_range_overlap(&exon_range, &cds_range))
    {
      agn_gene_validator_infer_utr_from_partial_exon(v, mrna, exon, cds);
      exon_index++;
      cds_index++;
    }
    else
    {
      agn_gene_validator_infer_utr_from_full_exon(v, mrna, exon, cds);
      exon_index++;
    }
  }

  gt_array_sort(v->utrs, (GtCompare)agn_gt_genome_node_compare);
  return true;
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

bool agn_gene_validator_validate_gene( AgnGeneValidator *v, GtFeatureNode *gene,
                                       AgnError *error )
{
  if(!agn_gt_feature_node_is_gene_feature(gene))
  {
    agn_error_add(error, false, "feature from line %u is a %s, expected a gene",
                  gt_genome_node_get_line_number((GtGenomeNode *)gene),
                  gt_feature_node_get_type(gene));
    return false;
  }

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
  GtFeatureNode *childfeature;
  bool isvalid = true;
  unsigned int num_mrnas = 0;
  for
  (
    childfeature = gt_feature_node_iterator_next(iter);
    childfeature != NULL;
    childfeature = gt_feature_node_iterator_next(iter)
  )
  {
    if(agn_gene_validator_validate_mrna(v, gene, childfeature, error))
    {
      GtRange grange = gt_genome_node_get_range((GtGenomeNode *)gene);
      GtRange trange = gt_genome_node_get_range((GtGenomeNode *)childfeature);
      if(!gt_range_contains(&grange, &trange))
      {
        agn_error_add(error, true, "range of mRNA '%s' exceeds range of gene "
                      "'%s'", gt_feature_node_get_attribute(childfeature, "ID"),
                      gt_feature_node_get_attribute(gene, "ID"));
        isvalid = false;
      }
      
      num_mrnas++;
    }
    else
    {
      gt_assert(agn_error_is_set(error));
      isvalid = false;
    }
  }
  gt_feature_node_iterator_delete(iter);

  if(num_mrnas == 0)
  {
    agn_error_add(error, false, "found no valid mRNAs for gene '%s' (line %u)",
                  gt_feature_node_get_attribute(gene, "ID"),
                  gt_genome_node_get_line_number((GtGenomeNode *)gene));
    return false;
  }

  return isvalid;
}

void agn_gene_validator_fix_parent(AgnGeneValidator *v, GtFeatureNode *p,
                                   GtFeatureNode *c)
{
  const char *parentattr = gt_feature_node_get_attribute(c,"Parent");
  const char *tid = gt_feature_node_get_attribute(p, "ID");
  if(strchr(parentattr, ','))
  {
    char *parentstr = gt_cstr_dup(parentattr);
    char *ptok = strtok(parentstr, ",");
    do
    {
      if(strcmp(tid, ptok) == 0)
      {
        gt_feature_node_set_attribute(c, "Parent", ptok);
        break;
      }
    } while ((ptok = strtok(NULL, ",")) != NULL);

    if(strchr(parentattr, ','))
    {
      fprintf(stderr, "could not fix Parent attribute for mRNA child feauture\n");
      exit(1);
    }

    gt_free(parentstr); 
  }
}

void agn_gene_validator_mrna_typecheck(AgnGeneValidator *v, GtFeatureNode *mrna,
                                       GtFeatureNode *child, AgnError *error)
{
  if(agn_gt_feature_node_is_cds_feature(child))
  {
    gt_array_add(v->cdss, child);
    GtRange cdsrange = gt_genome_node_get_range((GtGenomeNode *)child);
    if(v->start_codon_range.start == 0 || cdsrange.start < v->start_codon_range.start)
    {
      v->start_codon_range.start = cdsrange.start;
      v->start_codon_range.end   = cdsrange.start + 2;
    }
    if(v->stop_codon_range.end == 0 || cdsrange.end > v->stop_codon_range.end)
    {
      v->stop_codon_range.end   = cdsrange.end;
      v->stop_codon_range.start = cdsrange.end - 2;
    }
  }

  else if(agn_gt_feature_node_is_start_codon_feature(child))
  {
    v->start_codon_range = gt_genome_node_get_range((GtGenomeNode *)child);
    agn_gt_feature_node_remove_child(mrna, child);
  }

  else if(agn_gt_feature_node_is_stop_codon_feature(child))
  {
    v->stop_codon_range = gt_genome_node_get_range((GtGenomeNode *)child);
    agn_gt_feature_node_remove_child(mrna, child);
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
      agn_error_add(
        error, false, "%s feature (line %u) is not a valid child for mRNA '%s',"
        " ignoring", gt_feature_node_get_type(child),
        gt_genome_node_get_line_number((GtGenomeNode *)child),
        gt_feature_node_get_attribute(mrna, "ID")
      );
      agn_gt_feature_node_remove_child(mrna, child);
      gt_genome_node_delete((GtGenomeNode *)child);
    }
  }
  else
  {
    agn_error_add
    (
      error, false, "%s feature (line %u) is not a valid child for mRNA '%s', "
      "ignoring", gt_feature_node_get_type(child),
      gt_genome_node_get_line_number((GtGenomeNode *)child),
      gt_feature_node_get_attribute(mrna, "ID")
    );
    agn_gt_feature_node_remove_child(mrna, child);
    gt_genome_node_delete((GtGenomeNode *)child);
  }
}

bool agn_gene_validator_infer_cds(AgnGeneValidator *v, GtFeatureNode *mrna,
                                  AgnError *error)
{
  bool start_codon_found = v->start_codon_range.start != 0;
  bool stop_codon_found  = v->stop_codon_range.start  != 0;
  unsigned long num_exons = gt_array_size(v->exons);
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");
  unsigned int tln = gt_genome_node_get_line_number((GtGenomeNode *)mrna);
  bool isvalid = true;
  
  // Either infer CDS from start/stop codons and exons...
  if(start_codon_found && stop_codon_found && num_exons > 0)
  {
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
                                   &exon_range, exon_strand,
                                   &v->start_codon_range, &v->stop_codon_range,
                                   &cdsrange);

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
      agn_error_add(error, false, "error inferring CDS from codons and exons "
                    "for mRNA '%s' (line %u)", tid, tln);
      isvalid = false;
    }
  }
  // ...or return with error
  else
  {
    agn_error_add(error,false,"found no CDS for mRNA '%s' (line %u)", tid, tln);
    isvalid = false;
  }
  
  return isvalid;
}

bool agn_gene_validator_validate_mrna(AgnGeneValidator *v, GtFeatureNode *gene,
                                      GtFeatureNode *mrna, AgnError *error)
{
  const char *tid = gt_feature_node_get_attribute(mrna, "ID");
  if(!agn_gt_feature_node_is_mrna_feature(mrna))
  {
    agn_error_add(error, false,"feature from line %u is a %s, expected an mRNA",
                  gt_genome_node_get_line_number((GtGenomeNode *)mrna),
                  gt_feature_node_get_type(mrna) );
    return false;
  }

  v->cdss    = gt_array_new( sizeof(GtFeatureNode *) );
  v->exons   = gt_array_new( sizeof(GtFeatureNode *) );
  v->introns = gt_array_new( sizeof(GtFeatureNode *) );
  v->utrs    = gt_array_new( sizeof(GtFeatureNode *) );

  bool isvalid = true;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(mrna);
  GtFeatureNode *childfeature;
  for
  (
    childfeature = gt_feature_node_iterator_next(iter);
    childfeature != NULL;
    childfeature = gt_feature_node_iterator_next(iter)
  )
  {
    agn_gene_validator_fix_parent(v, mrna, childfeature);

    GtRange trange = gt_genome_node_get_range((GtGenomeNode *)mrna);
    GtRange crange = gt_genome_node_get_range((GtGenomeNode *)childfeature);
    if(!gt_range_contains(&trange, &crange))
    {
      agn_error_add(error, true, "%s at [%lu, %lu] exceeds range of its parent "
                    "mRNA '%s'", gt_feature_node_get_type(childfeature),
                    crange.start, crange.end, tid);
      isvalid = false;
    }

    agn_gene_validator_mrna_typecheck(v, mrna, childfeature, error);
  }
  gt_feature_node_iterator_delete(iter);

  if(isvalid)
  {
    gt_array_sort(v->cdss,    (GtCompare)agn_gt_genome_node_compare);
    gt_array_sort(v->exons,   (GtCompare)agn_gt_genome_node_compare);
    gt_array_sort(v->introns, (GtCompare)agn_gt_genome_node_compare);
    gt_array_sort(v->utrs,    (GtCompare)agn_gt_genome_node_compare);

    if(gt_array_size(v->cdss) == 0)
    {
      int success = agn_gene_validator_infer_cds(v, mrna, error);
      if(!success)
        isvalid = false;
    }
    agn_gene_validator_check_cds_phase(v);

    if(isvalid && gt_array_size(v->exons) == 0 && gt_array_size(v->cdss) > 0)
    {
      if(!agn_gene_validator_infer_exons(v, mrna, error))
        isvalid = false;
    }
    unsigned long i;
    for(i = 0; i < gt_array_size(v->exons) - 1; i++)
    {
      GtFeatureNode *fni = *(GtFeatureNode **)gt_array_get(v->exons, i);
      GtRange ri = gt_genome_node_get_range((GtGenomeNode *)fni);
      GtFeatureNode *fnj = *(GtFeatureNode **)gt_array_get(v->exons, i+1);
      GtRange rj = gt_genome_node_get_range((GtGenomeNode *)fnj);
      if(gt_range_overlap(&ri, &rj))
      {
        agn_error_add(error, true, "overlapping exons for mRNA '%s'", tid);
        isvalid = false;
      }
    }

    if(isvalid && gt_array_size(v->cdss) > 0 && gt_array_size(v->exons) > 0 &&
       gt_array_size(v->utrs) == 0)
    {
      if(!agn_gene_validator_infer_utrs(v, mrna, error))
        isvalid = false;
    }

    if(isvalid && gt_array_size(v->exons) > 1 && gt_array_size(v->introns) == 0)
    {
      if(!agn_gene_validator_infer_introns(v, mrna, error))
        isvalid = false;
    }

    if(isvalid && gt_array_size(v->exons) != gt_array_size(v->introns) + 1)
    {
      agn_error_add( error, false,
                     "mRNA '%s' (line %u) has %lu exons but %lu introns", tid,
                     gt_genome_node_get_line_number((GtGenomeNode *)mrna),
                     gt_array_size(v->exons), gt_array_size(v->introns) );
      isvalid = false;
    }

    agn_gene_validator_set_utrs(v);
  }

  gt_array_delete(v->cdss);
  gt_array_delete(v->exons);
  gt_array_delete(v->introns);
  gt_array_delete(v->utrs);

  return isvalid;
}
