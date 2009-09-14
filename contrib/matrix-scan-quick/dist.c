#include "dist.h"

values_t *new_values(int min, int max, double e)
{
    values_t *values = (values_t *) malloc(sizeof(values_t));
    values->size = (int) (max - min) / e;
    values->min = min;
    values->max = max;
    values->e = e;
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
    //value = MAX(MIN(value, values->max), values->min);
    int i = (int) ((value - values->min) / values->e);
    assert(i >= 0 && i < values->size);
    values->data[i] += 1;
}

void values_print(FILE *fout, values_t *values)
{
    int i;

    // find [a,b]
    i = 0;
    while (values->data[i] == 0)
        i++;
    int a = i;
    i = values->size;
    while (values->data[i] == 0)
        i--;
    int b = i;

    // total
    int total = 0;
    for (i = a; i <= b; i++)
    {
        total += values->data[i];
    }
    
    // print
    fprintf(fout, "#score\tocc\tco\tcco\tccdf\n");
    int cum = 0;
    for (i = a; i <= b; i++)
    {
        fprintf(fout, "%G\t%d\t%d\t%d\t%G\n", values->min + i * values->e, values->data[i], cum, total - cum, (total - cum) / (double) total);
        cum += values->data[i];
    }
}
