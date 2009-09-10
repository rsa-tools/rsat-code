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

void values_print(values_t *values)
{
    printf("#score\tocc\n");
    int i;
    for (i = 0; i < values->size; i++)
    {
        if (values->data[i] > 0)
            printf("%.4f\t%d\n", values->min + i * values->e, values->data[i]);
    }
}
