/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef PARSEVAL_UTILS
#define PARSEVAL_UTILS

/**
 * @type Data required for ParsEval's default runtime summary.
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
