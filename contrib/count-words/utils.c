#include "utils.h"

int VERBOSITY = 0;
int SHOW_WARNING = 0;
int SHOW_DEBUG = 0;

void free_string_buffer(string_buffer_t *string_buffer)
{
    free(string_buffer->data);
}

string_buffer_t *new_string_buffer()
{
    string_buffer_t *new_string_buffer = malloc(sizeof(string_buffer_t));
    new_string_buffer->allocated_size = 1024;
    new_string_buffer->data = malloc(sizeof(char) * new_string_buffer->allocated_size);
    return new_string_buffer;
}

inline void string_buffer_set_char(string_buffer_t *buffer, char c, int position)
{
    if (position >= buffer->allocated_size) 
    {
        buffer->allocated_size += 1024;
        buffer->data = (char *) realloc(buffer->data, sizeof(char) * buffer->allocated_size);
    }
    buffer->data[position] = c;
}
