/* Estimador Criva — port ligero de vma_methods/criva.py */
#include <math.h>
#include "../include/antipc_native.h"

static const int FIRST_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
};
static const int FIRST_PRIMES_N = sizeof(FIRST_PRIMES) / sizeof(FIRST_PRIMES[0]);

static double euler_product(int layers) {
    double result = 1.0;
    int i;
    if (layers <= 0) return 1.0;
    if (layers > FIRST_PRIMES_N) layers = FIRST_PRIMES_N;
    for (i = 0; i < layers; ++i) {
        result *= (1.0 - 1.0 / (double)FIRST_PRIMES[i]);
    }
    return result;
}

ANTIPC_API double antipc_criva_estimate(double x, int layers, int iterations) {
    double d0, d, t;
    int i, n;

    if (x < 2.0) return 0.0;
    if (layers < 1) layers = 10;
    if (layers > FIRST_PRIMES_N) layers = FIRST_PRIMES_N;
    if (iterations < 1) iterations = 8;

    d0 = euler_product(layers);
    d = d0;

    for (i = 0; i < iterations; ++i) {
        t = d * (1.0 - 1.0 / log(x + 1.0));
        d = (d + t) * 0.5;
    }

    for (n = 1; n < layers; ++n) {
        d += (d0 / (double)(1ULL << n)) * euler_product(n);
    }

    return d;
}