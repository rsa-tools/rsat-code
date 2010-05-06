#ifndef __UTILS__
#define __UTILS__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// ===========================================================================
// =                            Globals
// ===========================================================================
extern int VERBOSITY;
extern int SHOW_WARNING;
extern int SHOW_DEBUG;

// ===========================================================================
// =                            TRUE FALSE
// ===========================================================================
#define TRUE 1
#define FALSE 0

// ===========================================================================
// =                            MIN, MAX, ...
// ===========================================================================
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#define MAX(a,b) ((a) > (b) ? (a) : (b)) 

// ===========================================================================
// =                            ASSERT, ERROR, ...
// ===========================================================================
#define ASSERT(assertion, msg)        \
({                                    \
    if (!(assertion))                 \
    {                                 \
        fprintf(stderr, "Internal Error: invalid assertion (%s)", msg); \
        fprintf(stderr, " in %s line %d\n", __FILE__, __LINE__);        \
        exit(1);                      \
    }                                 \
})
//#define ASSERT(assertion, msg) {if(!(assertion)) {fprintf(stderr, "InternError invalid assertion (%s)\n", msg); exit(1);}}
#define VERBOSE1(...)          {if (VERBOSITY >= 1) {fprintf(stderr, __VA_ARGS__);}}
#define VERBOSE2(...)          {if (VERBOSITY >= 2) {fprintf(stderr, __VA_ARGS__);}}
#define VERBOSE3(...)          {if (VERBOSITY >= 3) {fprintf(stderr, __VA_ARGS__);}}
#define DEBUG(...)     ({if (SHOW_DEBUG) {fprintf(stdout, "[DEBUG] "); fprintf(stdout, __VA_ARGS__); \
                        fprintf(stdout, "\n"); fflush(stdout);}})
#define ERROR(...)     ({ fprintf(stderr, "Error "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); exit(1);})
#define WARNING(...)   ({if (SHOW_WARNING) {fprintf(stderr, "Warning "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); }})
#define CHECK_ALLOC(ptr) ({if (ptr == NULL) ERROR("Memory error")})

#define CHECK_VALUE(a, minval, maxval, msg) ({if (a < minval || a > maxval){ERROR(msg); exit(1);}})

// ===========================================================================
// =                            string_buffer
// ===========================================================================
typedef struct {
    int allocated_size;
    int current_size;
    char *data;
} string_buffer_t;

#define STRING_BUFFER_REALLOC(buffer, position) \
    if (position >= buffer->allocated_size) { \
        buffer->allocated_size += 1024; \
        buffer->data = (char *) realloc(buffer->data, sizeof(char) * buffer->allocated_size); \
    }

string_buffer_t *new_string_buffer();
void free_string_buffer(string_buffer_t *string_buffer);
inline void string_buffer_set_char(string_buffer_t *buffer, char c, int position);

#endif

