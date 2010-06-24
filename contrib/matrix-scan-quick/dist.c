#include "dist.h"

values_t *new_values(double min, double max, double e)
{
    values_t *values = (values_t *) malloc(sizeof(values_t));
    values->min = min;
    values->max = max;
    values->e = e;
    values->size = (int) (values->e + (values->max - values->min) / values->e);
    values->data = (int *) malloc(sizeof(int) * values->size);
    int i;
    for (i = 0; i < values->size; i++)
        values->data[i] = 0;

    return values;
}

void free_values(values_t *values)
{
    free(values->data);
    free(values);
}

void values_add(values_t *values, double value)
{
    ASSERT(value > values->min && value < values->max, "score out of range");
    int i = (int) ((value - values->min) / values->e);
    ASSERT(i >= 0 && i < values->size, "score out of table range");
    values->data[i] += 1;
}

void values_print(FILE *fout, values_t *values)
{
    int i;

    // find [a,b]
    i = 0;
    while (i < values->size - 1 && values->data[i] == 0)
        i++;
    int a = i;
    i = values->size;
    while (i > 0 && values->data[i] == 0)
        i--;
    int b = i;
    
    // DEBUG("a=%d b=%d", a, b);
    assert(a >= 0 && b >= 0 && a < values->size && b < values->size);

    // total
    int total = 0;
    for (i = a; i <= b; i++)
        total += values->data[i];
    
    // print
    fprintf(fout, "#score\tocc\tco\tcco\tdCDF\n");
    int cum = 0;
    for (i = a; i <= b; i++)
    {
        cum += values->data[i];
        fprintf(fout, "%G\t%d\t%d\t%d\t%G\n", 
            values->min + i * values->e, values->data[i], cum, total - cum, (total - cum) / (double) total);
    }
}
