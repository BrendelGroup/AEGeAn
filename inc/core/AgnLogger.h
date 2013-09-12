#ifndef AEGEAN_LOGGER
#define AEGEAN_LOGGER

#include "genometools.h"

/**
 * @class AgnLogger
 *
 * The AgnLogger class is desiged to store error, warning, and status messages.
 */
typedef struct AgnLogger AgnLogger;

/**
 * Free memory previously allocated to the logger.
 *
 * @param[in] logger   logger object to be freed
 */
void agn_logger_delete(AgnLogger *logger);

/**
 * Return an array containing all the error messages associated with this
 * logger. The user is responsible for freeing this array, but not its contents.
 *
 * @param[in] logger   logger object
 * @returns            an array containing all associated error messages
 */
GtArray *agn_logger_get_error_messages(AgnLogger *logger);

/**
 * Return an array containing all the status messages associated with this
 * logger. The user is responsible for freeing this array, but not its contents.
 *
 * @param[in] logger   logger object
 * @returns            an array containing all associated error messages
 */
GtArray *agn_logger_get_status_messages(AgnLogger *logger);

/**
 * Return an array containing all the warning messages associated with this
 * logger. The user is responsible for freeing this array, but not its contents.
 *
 * @param[in] logger   logger object
 * @returns            an array containing all associated error messages
 */
GtArray *agn_logger_get_warning_messages(AgnLogger *logger);

/**
 * Have any errors been logged?
 *
 * @param[in] logger   logger object
 * @returns            true if status message(s) have been logged, false
 *                     otherwise
 */
bool agn_logger_has_error(AgnLogger *logger);

/**
 * Have any status messages been logged?
 *
 * @param[in] logger   logger object
 * @returns            true if error(s) have been logged, false otherwise
 */
bool agn_logger_has_status(AgnLogger *logger);

/**
 * Have any warnings been logged?
 *
 * @param[in] logger   logger object
 * @returns            true if warning(s) have been logged, false otherwise
 */
bool agn_logger_has_warning(AgnLogger *logger);

/**
 * Add an error message to the logger.
 *
 * @param[in] logger     logger object
 * @param[in] format     format the error message a la printf
 * @param[in] ...        variable argument list for format
 */
void agn_logger_log_error(AgnLogger *logger, const char *format, ...);

/**
 * Add a status message/update to the logger.
 *
 * @param[in] logger     logger object
 * @param[in] format     format the status message a la printf
 * @param[in] ...        variable argument list for format
 */
void agn_logger_log_status(AgnLogger *logger, const char *format, ...);

/**
 * Add a warning message to the logger.
 *
 * @param[in] logger     logger object
 * @param[in] format     format the warning message a la printf
 * @param[in] ...        variable argument list for format
 */
void agn_logger_log_warning(AgnLogger *logger, const char *format, ...);

/**
 * Allocate memory for a new error object.
 *
 * @returns    a new logger object
 */
AgnLogger *agn_logger_new();

/**
 * Print the status messages, warnings, and errors that have been logged to the
 * given file stream.
 *
 * @param[in] logger       logger object
 * @param[in] outstream    the file pointer to which the error message(s) will
 *                         be written
 * @param[in] format       optional format to print as a header message for the
 *                         error (useful when multiple messages are associated
 *                         with a single error object)
 * @param[in] ...          variable argument list for format
 * @returns                true if errors were printed, false otherwise
 */
bool agn_logger_print_all(AgnLogger *logger, FILE *outstream,
                          const char *format, ...);

/**
 * Print the error messages associated with this logger to the given file stream.
 *
 * @param[in] logger       logger object
 * @param[in] outstream    the file pointer to which the error message(s) will
 *                         be written
 * @param[in] format       optional format to print as a header message for the
 *                         error (useful when multiple messages are associated
 *                         with a single error object)
 * @param[in] ...          variable argument list for format
 * @returns                true if errors were printed, false otherwise
 */
bool agn_logger_print_error(AgnLogger *logger, FILE *outstream,
                            const char *format, ...);

/**
 * Print the status messages associated with this logger to the given file
 * stream.
 *
 * @param[in] logger       logger object
 * @param[in] outstream    the file pointer to which the error message(s) will
 *                         be written
 * @param[in] format       optional format to print as a header message for the
 *                         error (useful when multiple messages are associated
 *                         with a single error object)
 * @param[in] ...          variable argument list for format
 * @returns                true if status messages were printed, false otherwise
 */
bool agn_logger_print_status(AgnLogger *logger, FILE *outstream,
                             const char *format, ...);

/**
 * Print the warning messages associated with this logger to the given file
 * stream.
 *
 * @param[in] logger       logger object
 * @param[in] outstream    the file pointer to which the error message(s) will
 *                         be written
 * @param[in] format       optional format to print as a header message for the
 *                         error (useful when multiple messages are associated
 *                         with a single error object)
 * @param[in] ...          variable argument list for format
 * @returns                true if warnings were printed, false otherwise
 */
bool agn_logger_print_warning(AgnLogger *logger, FILE *outstream,
                              const char *format, ...);

/**
 * Reset this logger object.
 *
 * @param[in] logger    logger object
 */
void agn_logger_unset(AgnLogger *logger);

#endif
