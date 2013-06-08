#include "AgnLogger.h"

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//
struct AgnLogger
{
  GtArray *errors;
  GtArray *messages;
  GtArray *warnings;
};


//----------------------------------------------------------------------------//
// Prototypes of private methods
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
void agn_logger_delete(AgnLogger *logger)
{
  agn_logger_unset(logger);
  gt_array_delete(logger->errors);
  gt_array_delete(logger->messages);
  gt_array_delete(logger->warnings);
  gt_free(logger);
  logger = NULL;
}

GtArray *agn_logger_get_messages(AgnLogger *logger)
{
  GtArray *messages = gt_array_new( sizeof(char *) );
  unsigned long i;
  for(i = 0; i < gt_array_size(logger->errors); i++)
  {
    GtError *suberror = *(GtError **)gt_array_get(logger->errors, i);
    const char *errmsg = gt_error_get(suberror);
    gt_array_add(messages, errmsg);
  }
  return logger->errors;
}

bool agn_logger_has_error(AgnLogger *logger)
{
  return gt_array_size(logger->errors) > 0;
}

bool agn_logger_has_status(AgnLogger *logger)
{
  return gt_array_size(logger->messages) > 0;
}

bool agn_logger_has_warning(AgnLogger *logger)
{
  return gt_array_size(logger->warnings) > 0;
}

void agn_logger_log_error(AgnLogger *logger, const char *format, ...)
{
  gt_assert(format);

  GtError *message = gt_error_new();
  va_list ap;
  va_start(ap, format);
  gt_error_vset(message, format, ap);
  va_end(ap);

  gt_array_add(logger->errors, message);
}

void agn_logger_log_status(AgnLogger *logger, const char *format, ...)
{
  gt_assert(format);

  GtError *message = gt_error_new();
  va_list ap;
  va_start(ap, format);
  gt_error_vset(message, format, ap);
  va_end(ap);

  gt_array_add(logger->messages, message);
}

void agn_logger_log_warning(AgnLogger *logger, const char *format, ...)
{
  gt_assert(format);

  GtError *message = gt_error_new();
  va_list ap;
  va_start(ap, format);
  gt_error_vset(message, format, ap);
  va_end(ap);

  gt_array_add(logger->warnings, message);
}

AgnLogger *agn_logger_new()
{
  AgnLogger *logger = gt_malloc( sizeof(AgnLogger) );
  logger->errors   = gt_array_new( sizeof(GtError *) );
  logger->messages = gt_array_new( sizeof(GtError *) );
  logger->warnings = gt_array_new( sizeof(GtError *) );
  return logger;
}

bool agn_logger_print_all(AgnLogger *logger, FILE *outstream,
                          const char *format, ...)
{
  unsigned long nerrors   = gt_array_size(logger->errors);
  unsigned long nmessages = gt_array_size(logger->messages);
  unsigned long nwarnings = gt_array_size(logger->warnings);
  if(nerrors + nmessages + nwarnings == 0)
    return false;

  if(format)
  {
    va_list ap;
    va_start(ap, format);
    vfprintf(outstream, format, ap);
    va_end(ap);
    fputs("\n", outstream);
  }

  unsigned long i;
  if(nmessages > 0)
  {
    fputs("    status:\n", outstream);
    for(i = 0; i < nmessages; i++)
    {
      GtError *message = *(GtError **)gt_array_get(logger->messages, i);
      fprintf(outstream, "        %s\n", gt_error_get(message));
    }
  }

  if(nwarnings > 0)
  {
    fputs("    warning:\n", outstream);
    for(i = 0; i < nwarnings; i++)
    {
      GtError *message = *(GtError **)gt_array_get(logger->warnings, i);
      fprintf(outstream, "        %s\n", gt_error_get(message));
    }
  }

  if(nerrors > 0)
  {
    fputs("    error:\n", outstream);
    for(i = 0; i < nerrors; i++)
    {
      GtError *message = *(GtError **)gt_array_get(logger->errors, i);
      fprintf(outstream, "        %s\n", gt_error_get(message));
    }
  }

  return nerrors > 0;
}

bool agn_logger_print_error(AgnLogger *logger, FILE *outstream,
                            const char *format, ...)
{
  unsigned long nerrors = gt_array_size(logger->errors);
  if(nerrors == 0)
    return false;

  if(format)
  {
    va_list ap;
    va_start(ap, format);
    vfprintf(outstream, format, ap);
    va_end(ap);
    fputs("\n", outstream);
  }

  unsigned long i;
  for(i = 0; i < nerrors; i++)
  {
    GtError *message = *(GtError **)gt_array_get(logger->errors, i);
    fprintf(outstream, "    %s\n", gt_error_get(message));
  }

  return true;
}

bool agn_logger_print_status(AgnLogger *logger, FILE *outstream,
                             const char *format, ...)
{
  unsigned long nmessages = gt_array_size(logger->messages);
  if(nmessages == 0)
    return false;

  if(format)
  {
    va_list ap;
    va_start(ap, format);
    vfprintf(outstream, format, ap);
    va_end(ap);
    fputs("\n", outstream);
  }

  unsigned long i;
  for(i = 0; i < nmessages; i++)
  {
    GtError *message = *(GtError **)gt_array_get(logger->messages, i);
    fprintf(outstream, "    %s\n", gt_error_get(message));
  }

  return true;
}

bool agn_logger_print_warning(AgnLogger *logger, FILE *outstream,
                              const char *format, ...)
{
  unsigned long nwarnings = gt_array_size(logger->warnings);
  if(nwarnings == 0)
    return false;

  if(format)
  {
    va_list ap;
    va_start(ap, format);
    vfprintf(outstream, format, ap);
    va_end(ap);
    fputs("\n", outstream);
  }

  unsigned long i;
  for(i = 0; i < nwarnings; i++)
  {
    GtError *message = *(GtError **)gt_array_get(logger->warnings, i);
    fprintf(outstream, "    %s\n", gt_error_get(message));
  }

  return true;
}

void agn_logger_unset(AgnLogger *logger)
{
  while(gt_array_size(logger->errors) > 0)
  {
    GtError *message = *(GtError **)gt_array_pop(logger->errors);
    gt_error_delete(message);
  }
  while(gt_array_size(logger->messages) > 0)
  {
    GtError *message = *(GtError **)gt_array_pop(logger->messages);
    gt_error_delete(message);
  }
  while(gt_array_size(logger->warnings) > 0)
  {
    GtError *message = *(GtError **)gt_array_pop(logger->warnings);
    gt_error_delete(message);
  }
}
