#include "AgnError.h"

struct AgnError
{
  GtArray *errors;
  unsigned int numfatal;
  bool isset;
};

void agn_error_add(AgnError *error, bool isfatal, const char *format, ...)
{
  gt_assert(format);

  GtError *newerror = gt_error_new();
  va_list ap;
  va_start(ap, format);
  gt_error_vset(newerror, format, ap);
  va_end(ap);
  
  gt_array_add(error->errors, newerror);
  if(isfatal)
    error->numfatal++;
  error->isset = true;
}

void agn_error_delete(AgnError *error)
{
  unsigned long i;
  for(i = 0; i < gt_array_size(error->errors); i++)
  {
    GtError *suberror = *(GtError **)gt_array_get(error->errors, i);
    gt_error_delete(suberror);
  }
  gt_array_delete(error->errors);
  gt_free(error);
  error = NULL;
}

GtArray *agn_error_get_messages(AgnError *error)
{
  GtArray *messages = gt_array_new( sizeof(char *) );
  unsigned long i;
  for(i = 0; i < gt_array_size(error->errors); i++)
  {
    GtError *suberror = *(GtError **)gt_array_get(error->errors, i);
    const char *errmsg = gt_error_get(suberror);
    gt_array_add(messages, errmsg);
  }
  return error->errors;
}

bool agn_error_is_fatal(AgnError *error)
{
  return error->numfatal > 0;
}

bool agn_error_is_set(AgnError *error)
{
  return error->isset;
}

AgnError *agn_error_new()
{
  AgnError *error = gt_malloc( sizeof(AgnError) );
  error->errors = gt_array_new( sizeof(GtError *) );
  error->numfatal = 0;
  error->isset = false;
  return error;
}

void agn_error_print(AgnError *error, FILE *outstream, const char *format, ...)
{
  if(format)
  {
    va_list ap;
    va_start(ap, format);
    vfprintf(outstream, format, ap);
    va_end(ap);
    fputs("\n", outstream);
  }
  unsigned long i;
  for(i = 0; i < gt_array_size(error->errors); i++)
  {
    GtError *suberror = *(GtError **)gt_array_get(error->errors, i);
    fprintf(outstream, "    %s\n", gt_error_get(suberror));
  }
}

void agn_error_unset(AgnError *error)
{
  agn_error_delete(error);
  error = agn_error_new();
}
