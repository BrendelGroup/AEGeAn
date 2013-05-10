#include "VangSchemaEntry.h"

int main(int argc, char **argv)
{
  if(argc != 2)
  {
    fputs("usage: vang schema.file\n", stderr);
    return 1;
  }

  gt_lib_init();

  char *schemafile = argv[1];
  FILE *schema = fopen(schemafile, "r");
  if(schema == NULL)
  {
    fprintf(stderr, "error: unable to open schema file '%s'\n", schemafile);
    return 1;
  }

  GtArray *entry_datatypes = gt_array_new( sizeof(char *) );
  GtHashmap *entries = gt_hashmap_new( GT_HASH_STRING, NULL,
                                       (GtFree)vang_schema_entry_delete );
  VangSchemaEntry *entry;
  while( (entry = vang_schema_entry_next(schema)) != NULL )
  {
    char *datatype = (char *)vang_schema_entry_get_datatype(entry);
    VangSchemaEntry *testentry = gt_hashmap_get(entries, datatype);
    if(testentry != NULL)
    {
      fprintf( stderr, "warning: already have an entry for data type '%s'; "
               "replacing\n", datatype );
      vang_schema_entry_delete(testentry);
    }
    gt_hashmap_add(entries, datatype, entry);
    gt_array_add(entry_datatypes, datatype);
  }

  unsigned long i;
  for(i = 0; i < gt_array_size(entry_datatypes); i++)
  {
    const char *type = *(const char **)gt_array_get(entry_datatypes, i);
    VangSchemaEntry *entry = gt_hashmap_get(entries, type);
    vang_schema_entry_to_string(entry, stdout);
    puts("");
  }
  gt_array_delete(entry_datatypes);
  gt_hashmap_delete(entries);

  gt_lib_clean();

  return 0;
}

