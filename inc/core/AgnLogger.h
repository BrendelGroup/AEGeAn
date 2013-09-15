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
 * @function Class destructor.
 */
void agn_logger_delete(AgnLogger *logger);

/**
 * @function Return an array containing all the error messages associated with
 * this logger. The user is responsible for freeing this array, but not its
 * contents.
 */
GtArray *agn_logger_get_error_messages(AgnLogger *logger);

/**
 * @function Return an array containing all the status messages associated with
 * this logger. The user is responsible for freeing this array, but not its
 * contents.
 */
GtArray *agn_logger_get_status_messages(AgnLogger *logger);

/**
 * @function Return an array containing all the warning messages associated with
 * this logger. The user is responsible for freeing this array, but not its
 * contents.
 */
GtArray *agn_logger_get_warning_messages(AgnLogger *logger);

/**
 * @function Have any errors been logged?
 */
bool agn_logger_has_error(AgnLogger *logger);

/**
 * @function Have any status messages been logged?
 */
bool agn_logger_has_status(AgnLogger *logger);

/**
 * @function Have any warnings been logged?
 */
bool agn_logger_has_warning(AgnLogger *logger);

/**
 * @function Add an error message to the logger using ``printf``-style string
 * formatting.
 */
void agn_logger_log_error(AgnLogger *logger, const char *format, ...);

/**
 * @function Add a status message/update to the logger using ``printf``-style
 * string formatting.
 */
void agn_logger_log_status(AgnLogger *logger, const char *format, ...);

/**
 * @function Add a warning message to the logger using ``print``-style string
 * formatting.
 */
void agn_logger_log_warning(AgnLogger *logger, const char *format, ...);

/**
 * @function Class constructor.
 */
AgnLogger *agn_logger_new();

/**
 * @function Print the status messages, warnings, and errors that have been
 * logged to the given file stream, ``printf``-style. Returns true if any errors
 * were printed, false otherwise.
 */
bool agn_logger_print_all(AgnLogger *logger, FILE *outstream,
                          const char *format, ...);

/**
 * @function Print the error messages associated with this logger to the given
 * file stream. Returns true if errors were printed.
 */
bool agn_logger_print_error(AgnLogger *logger, FILE *outstream,
                            const char *format, ...);

/**
 * @function Print the status messages associated with this logger to the given
 * file stream. Returns true if any messages were printed, false otherwise.
 */
bool agn_logger_print_status(AgnLogger *logger, FILE *outstream,
                             const char *format, ...);

/**
 * @function Print the warning messages associated with this logger to the given
 * file stream. Returns true if any warnings were printed, false otherwise.
 */
bool agn_logger_print_warning(AgnLogger *logger, FILE *outstream,
                              const char *format, ...);

/**
 * @function Reset this logger object.
 */
void agn_logger_unset(AgnLogger *logger);

#endif
