#include <string.h>
#include "AgnGtExtensions.h"
#include "AgnUtils.h"

GtArray* agn_gt_array_copy(GtArray *source, size_t size)
{
  unsigned long i;
  GtArray *new = gt_array_new(size);
  for(i = 0; i < gt_array_size(source); i++)
  {
    void *data = *(void **)gt_array_get(source, i);
    gt_array_add(new, data);
  }
  return new;
}

void agn_gt_feature_index_to_gff3(GtFeatureIndex *index, FILE *outstream)
{
  fputs("##gff-version   3\n", outstream);
  GtError *error = gt_error_new();
  GtStrArray *seqids = gt_feature_index_get_seqids(index, error);
  unsigned long i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtRange seqrange;
    gt_feature_index_get_range_for_seqid(index, &seqrange, seqid, error);
    fprintf(outstream, "##sequence-region   %s %lu %lu\n", seqid,
            seqrange.start, seqrange.end);
  }

  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *features = gt_feature_index_get_features_for_seqid(index, seqid,
                                                                error);
    unsigned long j;
    for(j = 0; j < gt_array_size(features); j++)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_get(features, j);
      agn_gt_feature_node_to_gff3(fn, outstream, true, NULL, NULL);
      fputs("###\n", outstream);
    }
  }

  gt_error_delete(error);
  gt_str_array_delete(seqids);
}

unsigned long agn_gt_feature_node_cds_length(GtFeatureNode *transcript)
{
  unsigned long length = 0;
  GtFeatureNode *current;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    if(agn_gt_feature_node_is_cds_feature(current))
      length += gt_genome_node_get_length((GtGenomeNode *)current);
  }
  gt_feature_node_iterator_delete(iter);

  // FIXME how should this be reported?
  // if(length % 3 != 0)
  // {
  //   const char *tid = gt_feature_node_get_attribute(transcript, "ID");
  //   fprintf(stderr, "warning: CDS for mRNA '%s' has length of %lu, "
  //           "not a multiple of 3\n", tid, length);
  // }

  return length;
}

GtArray *agn_gt_feature_node_children_of_type(GtFeatureNode *fn,
                                         bool (*typetestfunc)(GtFeatureNode *))
{
  GtArray *children = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    if(typetestfunc(current))
      gt_array_add(children, current);
  }
  gt_feature_node_iterator_delete(iter);
  gt_array_sort(children, (GtCompare)agn_gt_genome_node_compare);
  return children;
}

bool agn_gt_feature_node_fix_parent_attribute(GtFeatureNode *feature,
                                              GtFeatureNode *parent)
{
  bool success = true;
  const char *pid = gt_feature_node_get_attribute(parent, "ID");
  const char *parentattr = gt_feature_node_get_attribute(feature, "Parent");
  if(strchr(parentattr, ','))
  {
    char *parentstr = gt_cstr_dup(parentattr);
    char *ptok = strtok(parentstr, ",");
    do
    {
      if(strcmp(pid, ptok) == 0)
      {
        gt_feature_node_set_attribute(feature, "Parent", ptok);
        break;
      }
    } while ((ptok = strtok(NULL, ",")) != NULL);

    if(strchr(parentattr, ','))
      success = false;
    gt_free(parentstr);
  }

  return success;
}

bool agn_gt_feature_node_is_cds_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "CDS") ||
         gt_feature_node_has_type(feature, "coding sequence") ||
         gt_feature_node_has_type(feature, "coding_sequence");
}

bool agn_gt_feature_node_is_exon_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "exon");
}

bool agn_gt_feature_node_is_gene_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "gene");
}

bool agn_gt_feature_node_is_intron_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "intron");
}

bool agn_gt_feature_node_is_mrna_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "mRNA") ||
         gt_feature_node_has_type(feature, "messenger RNA") ||
         gt_feature_node_has_type(feature, "messenger_RNA");
}

bool agn_gt_feature_node_is_start_codon_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "start_codon") ||
         gt_feature_node_has_type(feature, "start codon") ||
         gt_feature_node_has_type(feature, "initiation codon");
}

bool agn_gt_feature_node_is_stop_codon_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "stop_codon") ||
         gt_feature_node_has_type(feature, "stop codon");
         // What about 'termination codon'?
}

bool agn_gt_feature_node_is_utr_feature(GtFeatureNode *feature)
{
  return gt_feature_node_has_type(feature, "UTR") ||
         gt_feature_node_has_type(feature, "untranslated region") ||
         gt_feature_node_has_type(feature, "untranslated_region") ||
         gt_feature_node_has_type(feature, "5' UTR") ||
         gt_feature_node_has_type(feature, "5'UTR") ||
         gt_feature_node_has_type(feature, "five prime UTR") ||
         gt_feature_node_has_type(feature, "five_prime_UTR") ||
         gt_feature_node_has_type(feature, "five prime untranslated region") ||
         gt_feature_node_has_type(feature, "five_prime_untranslated_region") ||
         gt_feature_node_has_type(feature, "3' UTR") ||
         gt_feature_node_has_type(feature, "3'UTR") ||
         gt_feature_node_has_type(feature, "three prime UTR") ||
         gt_feature_node_has_type(feature, "three_prime_UTR") ||
         gt_feature_node_has_type(feature, "three prime untranslated region") ||
         gt_feature_node_has_type(feature, "three_prime_untranslated_region");
}

unsigned long agn_gt_feature_node_num_transcripts(GtFeatureNode *gene)
{
  unsigned long count = 0;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(gene);
  GtFeatureNode *feature;

  for
  (
    feature = gt_feature_node_iterator_next(iter);
    feature != NULL;
    feature = gt_feature_node_iterator_next(iter)
  )
  {
    if(agn_gt_feature_node_is_mrna_feature(feature))
      count++;
  }
  gt_feature_node_iterator_delete(iter);

  return count;
}

bool agn_gt_feature_node_overlap(GtFeatureNode *first, GtFeatureNode *second)
{
  GtGenomeNode *one = (GtGenomeNode *)first;
  GtGenomeNode *two = (GtGenomeNode *)second;
  GtRange r1 = gt_genome_node_get_range(one);
  GtRange r2 = gt_genome_node_get_range(two);
  return gt_range_overlap(&r1, &r2);
}

bool agn_gt_feature_node_range_contains(GtFeatureNode *n1, GtFeatureNode *n2)
{
  GtRange r1 = gt_genome_node_get_range((GtGenomeNode *)n1);
  GtRange r2 = gt_genome_node_get_range((GtGenomeNode *)n2);
  return gt_range_contains(&r1, &r2);
}

void agn_gt_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn)
{
  gt_assert(root && fn);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *child;
  for(child  = gt_feature_node_iterator_next(iter);
      child != NULL;
      child  = gt_feature_node_iterator_next(iter))
  {
    agn_gt_feature_node_remove_tree(root, child);
    gt_feature_node_remove_leaf(fn, child);
  }
  gt_feature_node_iterator_delete(iter);
  gt_feature_node_remove_leaf(root, fn);
}

void agn_gt_feature_node_set_source_recursive( GtFeatureNode *feature,
                                               GtStr *source )
{
  gt_assert(source != NULL);
  GtFeatureNode *current;
  GtFeatureNodeIterator *iter;
  iter = gt_feature_node_iterator_new(feature);
  for( current = gt_feature_node_iterator_next(iter);
       current != NULL;
       current = gt_feature_node_iterator_next(iter) )
  {
    gt_feature_node_set_source(current, source);
  }
}

void agn_gt_feature_node_to_gff3(GtFeatureNode *feature, FILE *outstream,
                                 bool printchildren, char *prefix,
                                 GtHashmap *filtered_types)
{
  bool default_filter = false;
  if(filtered_types == NULL)
  {
    default_filter = true;
    filtered_types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
    gt_hashmap_add(filtered_types, "intron", filtered_types);
  }

  GtFeatureNode *current;
  GtFeatureNodeIterator *iter;
  GtHashmap *printed = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  iter = gt_feature_node_iterator_new(feature);
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    const char *type = gt_feature_node_get_type(current);
    if(gt_hashmap_get(filtered_types, type) != NULL ||
       gt_hashmap_get(printed, current) != NULL ||
       (!printchildren && current != feature))
    {
      continue;
    }

    GtGenomeNode *gn = (GtGenomeNode *)current;
    GtStr *seqid = gt_genome_node_get_seqid(gn);
    const char *source = gt_feature_node_get_source(current);
    unsigned long start = gt_genome_node_get_start(gn);
    unsigned long end = gt_genome_node_get_end(gn);
    char score_string[16];
    if(gt_feature_node_score_is_defined(current))
    {
      float score = gt_feature_node_get_score(current);
      sprintf(score_string, "%.2f", score);
    }
    else
      sprintf(score_string, "%s", ".");
    GtStrand strand = gt_feature_node_get_strand(current);
    GtPhase phase = gt_feature_node_get_phase(current);

    if(prefix != NULL && strcmp(prefix, "") != 0)
      fputs(prefix, outstream);

    fprintf(outstream, "%s\t%s\t%s\t%lu\t%lu\t%s\t%c\t%c\t", gt_str_get(seqid),
            source, type, start, end, score_string,
            agn_gt_strand_to_char(strand), agn_gt_phase_to_char(phase));

    GtStrArray *attributes = gt_feature_node_get_attribute_list(current);
    unsigned long num_attrs = gt_str_array_size(attributes);
    unsigned long i;
    for(i = 0; i < num_attrs; i++)
    {
      const char *attr_name = gt_str_array_get(attributes, i);
      if(i > 0)
        fputs(";", outstream);
      fprintf(outstream, "%s=%s", attr_name,
              gt_feature_node_get_attribute(current, attr_name));
    }
    fputs("\n", outstream);
    gt_str_array_delete(attributes);
    gt_hashmap_add(printed, current, current);
  }
  gt_feature_node_iterator_delete(iter);
  if(default_filter)
  {
    gt_hashmap_delete(filtered_types);
  }
  gt_hashmap_delete(printed);
}

int agn_gt_genome_node_compare(const void *n1, const void *n2)
{
  GtGenomeNode *gn1 = *(GtGenomeNode **)n1;
  GtGenomeNode *gn2 = *(GtGenomeNode **)n2;
  return gt_genome_node_cmp(gn1, gn2);
}

char agn_gt_phase_to_char(GtPhase phase)
{
  switch (phase)
  {
    case GT_PHASE_ZERO: return '0';
    case GT_PHASE_ONE: return '1';
    case GT_PHASE_TWO: return '2';
    case GT_PHASE_UNDEFINED: return '.';
    default : gt_assert(0);
  }
}

char agn_gt_strand_to_char(GtStrand strand)
{
  switch(strand)
  {
    case GT_STRAND_FORWARD: return '+';
    case GT_STRAND_REVERSE: return '-';
    case GT_STRAND_BOTH:    return '.';
    case GT_STRAND_UNKNOWN: return '?';
    default:                return GT_NUM_OF_STRAND_TYPES;
  }
}

GtStrArray* agn_gt_str_array_intersection(GtStrArray *a1, GtStrArray *a2)
{
  GtStrArray *intersection = gt_str_array_new();
  GtHashmap *added = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  unsigned long i;
  for(i = 0; i < gt_str_array_size(a1); i++)
  {
    char *str = (char *)gt_str_array_get(a1, i);
    gt_hashmap_add(added, str, str);
  }
  for(i = 0; i < gt_str_array_size(a2); i++)
  {
    char *str = (char *)gt_str_array_get(a2, i);
    if(gt_hashmap_get(added, str) != NULL)
    {
      gt_str_array_add_cstr(intersection, str);
    }
  }
  gt_hashmap_delete(added);
  return intersection;
}

GtStrArray* agn_gt_str_array_union(GtStrArray *a1, GtStrArray *a2)
{
  GtArray *strings = gt_array_new( sizeof(char *) );
  GtHashmap *added = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  unsigned long i;
  for(i = 0; i < gt_str_array_size(a1); i++)
  {
    char *str = (char *)gt_str_array_get(a1, i);
    if(gt_hashmap_get(added, str) == NULL)
    {
      gt_hashmap_add(added, str, str);
      gt_array_add(strings, str);
    }
  }
  for(i = 0; i < gt_str_array_size(a2); i++)
  {
    char *str = (char *)gt_str_array_get(a2, i);
    if(gt_hashmap_get(added, str) == NULL)
    {
      gt_hashmap_add(added, str, str);
      gt_array_add(strings, str);
    }
  }
  gt_hashmap_delete(added);
  // The whole reason I'm going through this mess--GtStrArray class has no sort
  // function.
  gt_array_sort(strings, (GtCompare)agn_string_compare);

  GtStrArray *uniona = gt_str_array_new();
  for(i = 0; i < gt_array_size(strings); i++)
  {
    const char *str = *(const char **)gt_array_get(strings, i);
    gt_str_array_add_cstr(uniona, str);
  }
  return uniona;
}
