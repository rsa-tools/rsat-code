#include "binomial.h"

// Compute P(x >= n) (binomial)
// n -- number of success
// N -- number of trials
// p -- prob of success
// binomial recursive formula:
//            P[x] * p(N-x)
//   P[x+1] = -------------
//               q(x+1)

//double cache[1024]

// P[x] = P[x-1] * p (N-x+1) / (q x)
double pbinom(int n, int N, double p)
{
    ASSERT(p > 0.0 && p <= 1.0, "invalid probability (p)");
    ASSERT(n <= N && n >= 0, "invalid n value");
    
    if (n == 0 || p == 1.0)
        return 1.0;
    
    double logp = log(p);
    double logq = log(1.0 - p);
    double logbin = N * logq;
    
    double S = 0.0; // sum P(x >= n)
    int x;
    for (x = 1; x <= N; x++) 
    {
        //printf("%d %d\n", N - x + 1, x);
        logbin += log(N - x + 1) - log(x) + logp - logq;
        if (x >= n)
            S += exp(logbin);
        if (x > N*p && exp(logbin) <= 0.0)
            break;
    }
    return MAX(MIN(S, 1.0), 0.0);
}

// int main() 
// {
//     double r = pbinom(2, 10, 0.2);
//     printf("r=%f\n", r);
// }

