#ifndef AEGEAN_UNIT_TEST
#define AEGEAN_UNIT_TEST

#include "genometools.h"

/**
 * FIXME
 */
typedef struct AgnUnitTest AgnUnitTest;

/**
 * FIXME
 */
void agn_unit_test_delete(AgnUnitTest *test);

AgnUnitTest *agn_unit_test_new(const char *label,
                               bool (*testfunc)(AgnUnitTest *));

/**
 * FIXME
 */
void agn_unit_test_print(AgnUnitTest *test, FILE *outstream);

/**
 * FIXME
 */
void agn_unit_test_result(AgnUnitTest *test, const char *label, bool success);

/**
 * FIXME
 */
void agn_unit_test_run(AgnUnitTest *test);

#endif