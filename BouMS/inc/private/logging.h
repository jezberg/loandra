/**
 * @file logging.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef LOGGING_H
#define LOGGING_H

#define BOUMS_LOG_LEVEL_OFF     0
#define BOUMS_LOG_LEVEL_WARN    1
#define BOUMS_LOG_LEVEL_INFO    2
#define BOUMS_LOG_LEVEL_VERBOSE 3
#define BOUMS_LOG_LEVEL_TRACE   4

#if BOUMS_LOG_LEVEL
#  include <stdio.h>
#  include <time.h>
#  define INIT_DURATION_MEAS() \
    clock_t startTime;         \
    ((void)startTime)
#  define START_DURATION_MEAS() startTime = clock()
#  define DURATION_STR          " in %.3f sec.\n"
#  define DURATION_VAL          (((double)(clock() - startTime)) / CLOCKS_PER_SEC)
#  if BOUMS_LOG_LEVEL == BOUMS_LOG_LEVEL_WARN
#    define LOG_WARN(...) printf(__VA_ARGS__)
#    define LOG_INFO(...)
#    define LOG_VERBOSE(...)
#    define LOG_TRACE(...)
#    define LOG_WARN_WITH_DURATION(str) printf(str DURATION_STR, DURATION_VAL)
#    define LOG_INFO_WITH_DURATION(str)
#    define LOG_VERBOSE_WITH_DURATION(str)
#    define LOG_TRACE_WITH_DURATION(str)
#  elif BOUMS_LOG_LEVEL == BOUMS_LOG_LEVEL_INFO
#    define LOG_WARN(...) printf(__VA_ARGS__)
#    define LOG_INFO(...) printf(__VA_ARGS__)
#    define LOG_VERBOSE(...)
#    define LOG_TRACE(...)
#    define LOG_WARN_WITH_DURATION(str) printf(str DURATION_STR, DURATION_VAL)
#    define LOG_INFO_WITH_DURATION(str) printf(str DURATION_STR, DURATION_VAL)
#    define LOG_VERBOSE_WITH_DURATION(str)
#    define LOG_TRACE_WITH_DURATION(str)
#  elif BOUMS_LOG_LEVEL == BOUMS_LOG_LEVEL_VERBOSE
#    define LOG_WARN(...)    printf(__VA_ARGS__)
#    define LOG_INFO(...)    printf(__VA_ARGS__)
#    define LOG_VERBOSE(...) printf(__VA_ARGS__)
#    define LOG_TRACE(...)
#    define LOG_WARN_WITH_DURATION(str)    printf(str DURATION_STR, DURATION_VAL)
#    define LOG_INFO_WITH_DURATION(str)    printf(str DURATION_STR, DURATION_VAL)
#    define LOG_VERBOSE_WITH_DURATION(str) printf(str DURATION_STR, DURATION_VAL)
#    define LOG_TRACE_WITH_DURATION(str)
#  elif BOUMS_LOG_LEVEL >= BOUMS_LOG_LEVEL_TRACE
#    define LOG_WARN(...)                  printf(__VA_ARGS__)
#    define LOG_INFO(...)                  printf(__VA_ARGS__)
#    define LOG_VERBOSE(...)               printf(__VA_ARGS__)
#    define LOG_TRACE(...)                 printf(__VA_ARGS__)
#    define LOG_WARN_WITH_DURATION(str)    printf(str DURATION_STR, DURATION_VAL)
#    define LOG_INFO_WITH_DURATION(str)    printf(str DURATION_STR, DURATION_VAL)
#    define LOG_VERBOSE_WITH_DURATION(str) printf(str DURATION_STR, DURATION_VAL)
#    define LOG_TRACE_WITH_DURATION(str)   printf(str DURATION_STR, DURATION_VAL)
#  endif
#else
#  define LOG_WARN(...)
#  define LOG_INFO(...)
#  define LOG_VERBOSE(...)
#  define LOG_TRACE(...)
#  define INIT_DURATION_MEAS()
#  define START_DURATION_MEAS()
#  define LOG_WARN_WITH_DURATION(str)
#  define LOG_INFO_WITH_DURATION(str)
#  define LOG_VERBOSE_WITH_DURATION(str)
#  define LOG_TRACE_WITH_DURATION(str)
#endif

#endif
