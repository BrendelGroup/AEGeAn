#include <string.h>
#include "genometools.h"
#include "VangRelation.h"
#include "VangSchemaEntry.h"

#define MAX_LINE_LENGTH 2048

//----------------------------------------------------------------------------//
// Private data structure definitions
//----------------------------------------------------------------------------//

struct VangSchemaEntry
{
  char *datatype;
  GtArray *relations;
  GtArray *relation_exclusions;
  GtArray *required_attributes;
};

struct VangRelationExclusion
{
  GtArray *exclusive_relations;
  char *note;
};


//----------------------------------------------------------------------------//
// Prototypes for private methods
//----------------------------------------------------------------------------//

VangRelation *vang_schema_parser_parse_relation(char *relation_string);
VangRelationExclusion *vang_schema_parser_parse_exclusion(char *exclusion_string);
GtArray *vang_schema_parser_parse_attributes(char *attribute_string);
VangSchemaEntry *vang_schema_entry_new(const char *datatype);
void vang_schema_entry_add_relation(VangSchemaEntry *entry, VangRelation *relation);
void vang_schema_entry_add_exclusion(VangSchemaEntry *entry, VangRelationExclusion *exclusion);
void vang_schema_entry_add_attributes(VangSchemaEntry *entry, GtArray *attributes);
VangRelationExclusion *vang_relation_exclusion_new();
void vang_relation_exclusion_add_relation(VangRelationExclusion *exclusion, char *relid);
void vang_relation_exclusion_set_note(VangRelationExclusion *exclusion, const char *note);


//----------------------------------------------------------------------------//
// Public method implementations
//----------------------------------------------------------------------------//

VangSchemaEntry *vang_schema_entry_next(FILE *schemafile)
{
  char buffer[MAX_LINE_LENGTH];
  VangSchemaEntry *entry = NULL;

  while( fgets(buffer, MAX_LINE_LENGTH, schemafile) != NULL )
  {
    if(buffer[0] == '#' || strcmp(buffer, "") == 0)
      continue;

    char *tok = strtok(buffer, "\t");
    entry = vang_schema_entry_new(tok);

    while((tok = strtok(NULL, "\t")) != NULL)
    {
      if(strncmp(tok, "Relation=", strlen("Relation=")) == 0)
      {
        VangRelation *relation = vang_schema_parser_parse_relation(tok);
        vang_schema_entry_add_relation(entry, relation);
      }
      else if(strncmp(tok, "Exclusive=", strlen("Exclusive=")) == 0)
      {
        VangRelationExclusion *exclusion = vang_schema_parser_parse_exclusion(tok);
        vang_schema_entry_add_exclusion(entry, exclusion);
      }
      else if(strncmp(tok, "Attribute=", strlen("Attribute=")) == 0)
      {
        GtArray *attributes = vang_schema_parser_parse_attributes(tok);
        vang_schema_entry_add_attributes(entry, attributes);
        gt_array_delete(attributes);
      }
      else
      {
        fprintf(stderr, "error: unsupported token '%s'", tok);
        exit(1);
      }
    }
    break;
  }

  return entry;
}

void vang_relation_exclusion_delete(VangRelationExclusion *exclusion)
{
  unsigned long i;

  for(i = 0; i < gt_array_size(exclusion->exclusive_relations); i++)
  {
    char *relid = *(char **)gt_array_get(exclusion->exclusive_relations, i);
    gt_free(relid);
  }
  gt_array_delete(exclusion->exclusive_relations);

  if(exclusion->note != NULL)
    gt_free(exclusion->note);

  gt_free(exclusion);
  exclusion = NULL;
}

void vang_schema_entry_delete(VangSchemaEntry *entry)
{
  unsigned int i;

  gt_free(entry->datatype);

  for(i = 0; i < gt_array_size(entry->relations); i++)
  {
    VangRelation *rel = *(VangRelation **)gt_array_get(entry->relations, i);
    vang_relation_delete(rel);
  }
  gt_array_delete(entry->relations);

  for(i = 0; i < gt_array_size(entry->relation_exclusions); i++)
  {
    VangRelationExclusion *exclusion = *(VangRelationExclusion **)gt_array_get(entry->relation_exclusions, i);
    vang_relation_exclusion_delete(exclusion);
  }
  gt_array_delete(entry->relation_exclusions);

  for(i = 0; i < gt_array_size(entry->required_attributes); i++)
  {
    char *attribute = *(char **)gt_array_get(entry->required_attributes, i);
    gt_free(attribute);
  }
  gt_array_delete(entry->required_attributes);

  entry = NULL;
}


//----------------------------------------------------------------------------//
// Private method implementations
//----------------------------------------------------------------------------//

VangRelation *vang_schema_parser_parse_relation(char *relation_string)
{
  unsigned long i;
  char *tok;
  char *relstrcpy = gt_cstr_dup(relation_string);
  GtArray *tokens = gt_array_new( sizeof(char *) );

  // Store all key/value token in an array
  tok = strtok(relstrcpy, ";");
  gt_array_add(tokens, tok);
  while((tok = strtok(NULL, ";")) != NULL)
  {
    gt_array_add(tokens, tok);
  }

  // Parse key/value pairs
  GtHashmap *attributes = gt_hashmap_new(GT_HASH_STRING, (GtFree)gt_free_func, (GtFree)gt_free_func);
  GtArray *mapkeys = gt_array_new( sizeof(char *) );
  for(i = 0; i < gt_array_size(tokens); i++)
  {
    char *key, *value;
    tok = *(char **)gt_array_get(tokens, i);
    key = strtok(tok, "=");
    key = gt_cstr_dup(key);
    value = strtok(NULL, "=");
    value = gt_cstr_dup(value);

    gt_hashmap_add(attributes, key, value);
    gt_array_add(mapkeys, key);
  }

  // Grab required attributes: relation ID and node type
  const char *id = gt_hashmap_get(attributes, "Relation");
  const char *nodetype = gt_hashmap_get(attributes, "Nodetype");
  if(id == NULL || nodetype == NULL)
  {
    fprintf(stderr, "error: relation ID and node type must be specified for each relation: '%s'\n", relation_string);
    exit(1);
  }

  // Set any other values
  VangRelation *rel = vang_relation_new(id, nodetype);
  for(i = 0; i < gt_array_size(mapkeys); i++)
  {
    const char *key = *(char **)gt_array_get(mapkeys, i);
    if(strcmp(key, "Relation") == 0 || strcmp(key, "Nodetype") == 0)
      continue;

    const char *value = gt_hashmap_get(attributes, key);
    if(strcmp(key, "Degree") == 0)
    {
      DegreeConstraint dc = vang_degree_constraint_parse(value);
      vang_relation_set_degree(rel, dc.context, dc.operator, dc.degree);
    }
    else if(strcmp(key, "Key") == 0)
    {
      vang_relation_set_key(rel, value);
    }
    if(strcmp(key, "Spatial") == 0)
    {
      SpatialConstraint constraint = vang_spatial_constraint_parse(value);
      vang_relation_set_spatial(rel, constraint);
    }
  }

  gt_free(relstrcpy);
  gt_array_delete(tokens);
  gt_array_delete(mapkeys);
  gt_hashmap_delete(attributes);

  return rel;
}

VangRelationExclusion *vang_schema_parser_parse_exclusion(char *exclusion_string)
{
  char *exclstrcpy = gt_cstr_dup(exclusion_string);
  char *excltok = strtok(exclstrcpy, ";");
  char *notetok = strtok(NULL, ";");

  VangRelationExclusion *exclusion = vang_relation_exclusion_new();

  char *tok = strtok(excltok, "=");
  tok = strtok(NULL, "=");
  char *rid = strtok(tok, ",");
  vang_relation_exclusion_add_relation(exclusion, rid);
  while( (rid = strtok(NULL, ",")) != NULL )
  {
    vang_relation_exclusion_add_relation(exclusion, rid);
  }

  if(notetok != NULL)
  {
    vang_relation_exclusion_set_note(exclusion, notetok + strlen("Note="));
  }

  gt_free(exclstrcpy);
  return exclusion;
}

GtArray *vang_schema_parser_parse_attributes(char *attribute_string)
{
  char *attrstrcpy = gt_cstr_dup(attribute_string);
  GtArray *attributes = gt_array_new( sizeof(char *) );

  char *tok = strtok(attrstrcpy + strlen("Attribute="), ",");
  gt_array_add_elem(attributes, gt_cstr_dup(tok), sizeof(char *));
  while( (tok = strtok(NULL, ",")) != NULL )
  {
    gt_array_add_elem(attributes, gt_cstr_dup(tok), sizeof(char *));
  }

  gt_free(attrstrcpy);

  return attributes;
}

VangSchemaEntry *vang_schema_entry_new(const char *datatype)
{
  VangSchemaEntry *entry = gt_malloc( sizeof(VangSchemaEntry) );

  entry->datatype = gt_cstr_dup(datatype);
  entry->relations = gt_array_new( sizeof(VangRelation *) );
  entry->relation_exclusions = gt_array_new( sizeof(VangRelationExclusion *) );
  entry->required_attributes = gt_array_new( sizeof(char *) );

  return entry;
}

void vang_schema_entry_add_relation(VangSchemaEntry *entry, VangRelation *relation)
{
  gt_array_add(entry->relations, relation);
}

void vang_schema_entry_add_exclusion(VangSchemaEntry *entry, VangRelationExclusion *exclusion)
{
  gt_array_add(entry->relation_exclusions, exclusion);
}

void vang_schema_entry_add_attributes(VangSchemaEntry *entry, GtArray *attributes)
{
  gt_array_add_array(entry->required_attributes, attributes);
}

VangRelationExclusion *vang_relation_exclusion_new()
{
  VangRelationExclusion *exclusion = gt_malloc( sizeof(VangRelationExclusion) );
  exclusion->exclusive_relations = gt_array_new( sizeof(char *) );
  exclusion->note = NULL;

  return exclusion;
}

void vang_relation_exclusion_add_relation(VangRelationExclusion *exclusion, char *relid)
{
  gt_array_add_elem(exclusion->exclusive_relations, gt_cstr_dup(relid), sizeof(char *));
}

void vang_relation_exclusion_set_note(VangRelationExclusion *exclusion, const char *note)
{
  if(exclusion->note != NULL)
    gt_free(exclusion->note);
  exclusion->note = gt_cstr_dup(note);
}

