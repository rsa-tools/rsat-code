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

using namespace std;

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

extern int VERBOSITY;
extern int VERSION;
extern char* COMMAND_LINE;

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

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#define RAND (rand() / (double) RAND_MAX)

#define VERBOSE1(...)         {if (VERBOSITY >= 1) {fprintf(stderr, __VA_ARGS__);}}
#define VERBOSE2(...)         {if (VERBOSITY >= 2) {fprintf(stderr, __VA_ARGS__);}}
#define VERBOSE3(...)         {if (VERBOSITY >= 3) {fprintf(stderr, __VA_ARGS__);}}
#define CHECK_VALUE(a, minval, maxval, msg) {if (a < minval || a > maxval){ERROR(msg); exit(1);}}

#define ASSERT(assertion, msg)                          \
{                                                       \
    if (!(assertion)) {                                 \
        fprintf(stderr, "Assertion error (%s)\n", msg); \
        exit(1);                                        \
    }                                                   \
} 

#define WARNING(...)                 \
{                                    \
    fprintf(stderr, "WARNING: ");    \
    fprintf(stderr, __VA_ARGS__);    \
    fprintf(stderr, "\n");           \
}

#define ERROR(...)                   \
{                                    \
    fprintf(stderr, "ERROR: ");      \
    fprintf(stderr, __VA_ARGS__);    \
    fprintf(stderr, "\n");           \
    exit(1);                         \
}

#define DEBUG(...)                    \
{                                     \
    if (VERBOSITY >= 0) {             \
        fprintf(stderr, "[DEBUG] ");  \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n");        \
        fflush(stderr);               \
    }                                 \
}

// TRACE
#define NEW_TRACE(trace, i)   {if (VERBOSITY >= 3) {char name[256]; sprintf(name, "%d.ic", i);trace = fopen(name, "w");}}
#define CLOSE_TRACE(tace)     {if (VERBOSITY >= 3) {fclose(trace);}}
#define TRACE(i, ...)         {if (VERBOSITY >= 3) {fprintf(trace, __VA_ARGS__);}}


#endif
