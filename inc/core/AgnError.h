#ifndef AEGEAN_ERROR
#define AEGEAN_ERROR

#include "genometools.h"

/**
 * The AgnError class provides a few extensions to the GtError class from the
 * GenomeTools library.
 */
typedef struct AgnError AgnError;

/**
 * Add a message to this error.
 *
 * @param[in] error      the error object
 * @param[in] isfatal    indicates whether this is a fatal error
 * @param[in] format     format the error message a la printf
 * @param[in] ...        variable argument list for format
 */
void agn_error_add(AgnError *error, bool isfatal, const char *format, ...);

/**
 * Free memory previously allocated to the error.
 *
 * @param[in] error    the error object to be de-allocated
 */
void agn_error_delete(AgnError *error);

/**
 * Return an array containing all the messages associated with this error. The
 * user is responsible for freeing this array.
 *
 * @param[in] error    the error object
 * @returns            an array containing all associated error messages
 */
GtArray *agn_error_get_messages(AgnError *error);

/**
 * Is this error a fatal error?
 *
 * @param[in] error    the error object
 * @returns            true if any of the associated messages were marked as
 *                     fatal, false otherwise
 */
bool agn_error_is_fatal(AgnError *error);

/**
 * Is this error set?
 *
 * @param[in] error    the error object
 * @returns            true if the error is set, false otherwise
 */
bool agn_error_is_set(AgnError *error);

/**
 * Allocate memory for a new error object.
 *
 * @returns    a new error object
 */
AgnError *agn_error_new();

/**
 * Print the messages associated with this error to the give file stream.
 *
 * @param[in] error        the error object
 * @param[in] outstream    the file pointer to which the error message(s) will
 *                         be written
 * @param[in] format       optional format to print as a header message for the
 *                         error (useful when multiple messages are associated
 *                         with a single error object)
 * @param[in] ...          variable argument list for format
 */
void agn_error_print(AgnError *error, FILE *outstream, const char *format, ...);

/**
 * Reset this error object.
 *
 * @param[in] error    the error object
 */
void agn_error_unset(AgnError *error);

#endif
