#ifndef PARSEVAL_UTILS
#define PARSEVAL_UTILS

/**
 * @type FIXME
 */
struct PeHtmlOverviewData
{
  const char *start_time;
  const char *refrlabel;
  const char *predlabel;
  int         argc;
  char      **argv;
};
typedef struct PeHtmlOverviewData PeHtmlOverviewData;

char *pe_get_start_time();
void pe_summary_html_overview(FILE *outstream, void *data);
void pe_summary_header(ParsEvalOptions *options, FILE *outstream,
                       char *start_time, int argc, char **argv);

#endif
