/***************************************************************************
 *                                                                         *
 *  utils.h
 *  
 *   
 *
 *                                                                         *
 ***************************************************************************/

#ifndef __UTILS__
#define __UTILS__

#include <stdlib.h>
#include <stdio.h>
//#include <assert.h>
#include <string.h>
#include <time.h>

extern int SHOW_DEBUG;
extern int VERBOSITY;
//extern int VERSION;
//extern char* COMMAND_LINE;

// ===========================================================================
// =                            TRUE FALSE
// ===========================================================================
#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif
// ===========================================================================
// =                            MIN, MAX, ...
// ===========================================================================
#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b)) 
#endif


#define RAND (rand() / (double) RAND_MAX)

#define VERBOSE1(...)         {if (VERBOSITY >= 1) {fprintf(stderr, __VA_ARGS__);}}
#define VERBOSE2(...)         {if (VERBOSITY >= 2) {fprintf(stderr, __VA_ARGS__);}}
#define VERBOSE3(...)         {if (VERBOSITY >= 3) {fprintf(stderr, __VA_ARGS__);}}
#define CHECK_VALUE(a, minval, maxval, msg) {if (a < minval || a > maxval) {ERROR(msg); exit(1);}}

#define ASSERT(assertion, msg)        \
({                                    \
    if (!(assertion))                 \
    {                                 \
        fprintf(stderr, "Internal Error: invalid assertion (%s)", msg); \
        fprintf(stderr, " in %s line %d\n", __FILE__, __LINE__);        \
        exit(1);                      \
    }                                 \
}) 

#define ENSURE(assertion, msg)        \
({                                    \
    if (!(assertion))                 \
    {                                 \
        fprintf(stderr, "Internal Error: %s\n", msg); \
        exit(1);                      \
    }                                 \
})

#define FATAL_ERROR(...)              \
({                                    \
    fprintf(stderr, "Fatal Error: "); \
    fprintf(stderr, __VA_ARGS__);     \
    fprintf(stderr, "\n");            \
    exit(1);                          \
})

#define WARNING(...)                  \
({                                    \
    if (SHOW_WARNING)                 \
    {                                 \
        fprintf(stderr, "Warning: "); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n");        \
    }                                 \
})

#define INFO(...)                     \
({                                    \
    if (1)                            \
    {                                 \
        fprintf(stderr, "[INFO] ");   \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n");        \
        fflush(stderr);               \
    }                                 \
})

#define DEBUG(...)                    \
({                                    \
    if (SHOW_DEBUG)                   \
    {                                 \
        fprintf(stdout, "[DEBUG] ");  \
        fprintf(stdout, __VA_ARGS__); \
        fprintf(stdout, "\n");        \
        fflush(stdout);               \
    }                                 \
})


#endif
