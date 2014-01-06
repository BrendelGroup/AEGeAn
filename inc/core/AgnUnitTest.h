#ifndef AEGEAN_UNIT_TEST
#define AEGEAN_UNIT_TEST

#include "genometools.h"

/**
 * @class AgnUnitTest
 * Class used for unit testing of classes and modules.
 */
typedef struct AgnUnitTest AgnUnitTest;

/**
 * @function Destructor.
 */
void agn_unit_test_delete(AgnUnitTest *test);

/**
 * @function Class constructor, where ``label`` is a label for the test and
 * ``testfunc`` is a pointer to the function that will execute the test.
 */
AgnUnitTest *agn_unit_test_new(const char *label,
                               bool (*testfunc)(AgnUnitTest *));

/**
 * @function Prints results of the unit test to ``outstream``.
 */
void agn_unit_test_print(AgnUnitTest *test, FILE *outstream);

/**
 * @function Add a result to this unit test.
 */
void agn_unit_test_result(AgnUnitTest *test, const char *label, bool success);

/**
 * @function Returns true if all the results checked with this unit test passed,
 * false otherwise.
 */
bool agn_unit_test_success(AgnUnitTest *test);

/**
 * @function Run the unit test.
 */
void agn_unit_test_run(AgnUnitTest *test);

#endif
