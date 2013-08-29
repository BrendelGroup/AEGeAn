#ifndef AEGEAN_UNIT_TEST
#define AEGEAN_UNIT_TEST

#include "genometools.h"

/**
 * Class used for unit testing of classes and modules.
 */
typedef struct AgnUnitTest AgnUnitTest;

/**
 * Destructor.
 *
 * @param[in] test    object to be destroyed
 */
void agn_unit_test_delete(AgnUnitTest *test);

AgnUnitTest *agn_unit_test_new(const char *label,
                               bool (*testfunc)(AgnUnitTest *));

/**
 * Prints results of the unit test to the given output stream.
 *
 * @param[in]  test         the unit test
 * @param[out] outstream    the output file to which results should be written
 */
void agn_unit_test_print(AgnUnitTest *test, FILE *outstream);

/**
 * Add a result to this unit test.
 *
 * @param[in] test       the unit test
 * @param[in] label      label for the result
 * @param[in] success    true for a passing result, false for failure
 */
void agn_unit_test_result(AgnUnitTest *test, const char *label, bool success);

/**
 * Run the unit test.
 *
 * @param[in] test    the unit test
 */
void agn_unit_test_run(AgnUnitTest *test);

#endif
