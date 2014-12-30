/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "core/error_api.h"
#include "core/hashmap_api.h"
#include "AgnSeqGeneMapping.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnSeqGeneMapping
{
  GtStr *filename;
  FILE *file;
  GtHashmap *mapping;
};

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * Callback function: print sequence-gene mapping.
 */
static int map_visit_print(void *key, void *value, void *data, GtError *error);

//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_seq_gene_mapping_add(AgnSeqGeneMapping *map, AgnLocus *locus)
{
  GtArray *genes;
  GtStr *seqid;
  GtStrArray *genelist;

  agn_assert(map && locus);
  seqid = gt_genome_node_get_seqid(locus);
  genelist = gt_hashmap_get(map->mapping, gt_str_get(seqid));
  if(genelist == NULL)
  {
    genelist = gt_str_array_new();
    char *key = gt_cstr_dup(gt_str_get(seqid));
    gt_hashmap_add(map->mapping, key, genelist);
  }

  genes = agn_locus_get_genes(locus);
  while(gt_array_size(genes) > 0)
  {
    GtFeatureNode **gene = gt_array_pop(genes);
    const char *geneid = gt_feature_node_get_attribute(*gene, "ID");
    gt_str_array_add_cstr(genelist, geneid);
  }
  gt_array_delete(genes);
}

void agn_seq_gene_mapping_close(AgnSeqGeneMapping *map)
{
  agn_assert(map);
  if(map->file == NULL)
  {
    map->file = fopen(gt_str_get(map->filename), "w");
    if(map->file == NULL)
    {
      fprintf(stderr, "unable to open sequence-gene map file '%s'\n",
              gt_str_get(map->filename));
      exit(1);
    }
  }

  GtError *error = gt_error_new();
  int result = gt_hashmap_foreach_in_key_order(map->mapping, map_visit_print,
                                               map->file, error);
  if(result)
  {
    fprintf(stderr, "error writing sequence-gene mapping to disk: %s\n",
            gt_error_get(error));
    exit(1);
  }
  agn_seq_gene_mapping_delete(map);
  gt_error_delete(error);
}

void agn_seq_gene_mapping_delete(AgnSeqGeneMapping *map)
{
  agn_assert(map);
  gt_str_delete(map->filename);
  if(map->file)
    fclose(map->file);
  gt_hashmap_delete(map->mapping);
}

GtStrArray *agn_seq_gene_mapping_get_genes(AgnSeqGeneMapping *map,
                                           const char *seqid)
{
  agn_assert(map && seqid);
  GtStrArray *genelist = gt_hashmap_get(map->mapping, seqid);
  agn_assert(genelist);
  return gt_str_array_ref(genelist);
}

bool agn_seq_gene_mapping_has_seqid(AgnSeqGeneMapping *map, const char *seqid)
{
  agn_assert(map && seqid);
  return gt_hashmap_get(map->mapping, seqid) != NULL;
}

AgnSeqGeneMapping *agn_seq_gene_mapping_open(const char *filepath)
{
  AgnSeqGeneMapping *map = agn_seq_gene_mapping_new(filepath);

  FILE *mapfile = fopen(gt_str_get(map->filename), "r");
  if(mapfile == NULL)
  {
    fprintf(stderr, "unable to open sequence-gene mapfile '%s'\n",
            gt_str_get(map->filename));
    exit(1);
  }

  char buffer[1048576];
  GtStrArray *geneids = gt_str_array_new();
  while(fgets(buffer, 1048576, mapfile) != NULL)
  {
    if(strlen(buffer) == 0)
      continue;
    char *seqid = strtok(buffer, "\t\n");
    agn_assert(seqid);

    char *geneid;
    while( (geneid = strtok(NULL, " \n")) != NULL)
      gt_str_array_add_cstr(geneids, geneid);

    gt_hashmap_add(map->mapping, gt_cstr_dup(seqid), geneids);
  }
  fclose(mapfile);

  return map;
}

AgnSeqGeneMapping *agn_seq_gene_mapping_new(const char *filepath)
{
  AgnSeqGeneMapping *map = gt_malloc( sizeof(AgnSeqGeneMapping) );
  map->filename = gt_str_new_cstr(filepath);
  map->file = NULL;
  map->mapping = gt_hashmap_new(GT_HASH_STRING, (GtFree)gt_free_func,
                                (GtFree)gt_str_array_delete);
  return map;
}

bool agn_seq_gene_mapping_unmap_gene(AgnSeqGeneMapping *map,
                                     const char *seqid,
                                     const char *geneid)
{
  agn_assert(map && seqid && geneid);
  unsigned long i;
  bool found;
  GtStrArray *genelist, *newgenelist;

  if(!agn_seq_gene_mapping_has_seqid(map, seqid))
    return false;

  genelist = gt_hashmap_get(map->mapping, seqid);
  newgenelist = gt_str_array_new();
  for(i = 0; i < gt_str_array_size(genelist); i++)
  {
    const char *test_geneid = gt_str_array_get(genelist, i);
    if(strcmp(test_geneid, geneid) == 0)
      found = true;
    else
      gt_str_array_add_cstr(newgenelist, test_geneid);
  }

  gt_hashmap_remove(map->mapping, (char *)seqid);
  if(gt_str_array_size(newgenelist) > 0)
    gt_hashmap_add(map->mapping, gt_cstr_dup(seqid), newgenelist);
  else
    gt_str_array_delete(newgenelist);

  return found;
}

GtStrArray *agn_seq_gene_mapping_unmap_seqid(AgnSeqGeneMapping *map,
                                             const char *seqid)
{
  agn_assert(map && seqid);
  if(agn_seq_gene_mapping_has_seqid(map, seqid))
  {
    GtStrArray *genelist = agn_seq_gene_mapping_get_genes(map, seqid);
    gt_hashmap_remove(map->mapping, seqid);
    return genelist;
  }
  return NULL;
}

static int map_visit_print(void *key, void *value, void *data, GtError *error)
{
  const char *seqid = key;
  GtStrArray *genelist = value;
  FILE *outstream = data;
  unsigned long i;
  fprintf(outstream, "%s\t", seqid);
  for(i = 0; i < gt_str_array_size(genelist); i++)
  {
    const char *geneid = gt_str_array_get(genelist, i);
    if(i > 0)
      fputs(" ", outstream);
    fputs(geneid, outstream);
  }
  fputs("\n", outstream);
  return 0;
}
