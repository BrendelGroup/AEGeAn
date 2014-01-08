#include <string.h>
#include "core/hashmap_api.h"
#include "AgnUtils.h"

GtArray* agn_array_copy(GtArray *source, size_t size)
{
  GtUword i;
  GtArray *new = gt_array_new(size);
  for(i = 0; i < gt_array_size(source); i++)
  {
    void *data = *(void **)gt_array_get(source, i);
    gt_array_add(new, data);
  }
  return new;
}

double agn_calc_splice_complexity(GtArray *transcripts)
{
  // FIXME
  return -1.0;
}

int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}

int agn_string_compare(const void *p1, const void *p2)
{
  const char *s1 = *(char **)p1;
  const char *s2 = *(char **)p2;
  return strcmp(s1, s2);
}

GtStrArray* agn_str_array_union(GtStrArray *a1, GtStrArray *a2)
{
  GtArray *strings = gt_array_new( sizeof(char *) );
  GtHashmap *added = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  GtUword i;
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
  gt_array_delete(strings);
  return uniona;
}