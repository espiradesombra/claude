/* Criba híbrida ascendente/descendente — vma-methods/cribas.py */
#include <stdlib.h>
#include <stdint.h>
#include "../include/antipc_native.h"

ANTIPC_API uint32_t antipc_sieve_hibrida_count(uint32_t limit) {
    uint8_t* is_composite;
    uint32_t mid, p, j, start, r, c;
    uint32_t* small_primes;
    uint32_t n_small, i;

    if (limit < 2) return 0;

    is_composite = (uint8_t*)calloc((size_t)limit + 1, 1);
    if (!is_composite) return 0;
    is_composite[0] = is_composite[1] = 1;

    mid = limit / 2;
    p = 2;
    while (p * p <= mid) {
        if (!is_composite[p]) {
            for (j = p * p; j <= mid; j += p) {
                is_composite[j] = 1;
            }
        }
        ++p;
    }

    n_small = 0;
    for (p = 2; p <= mid; ++p) {
        if (!is_composite[p]) ++n_small;
    }
    small_primes = (uint32_t*)malloc((size_t)n_small * sizeof(uint32_t));
    if (!small_primes) {
        free(is_composite);
        return 0;
    }
    i = 0;
    for (p = 2; p <= mid; ++p) {
        if (!is_composite[p]) small_primes[i++] = p;
    }

    for (i = 0; i < n_small; ++i) {
        p = small_primes[i];
        r = limit % p;
        start = limit - r;
        if (start == p) continue;
        j = start;
        while (j > mid) {
            if (!is_composite[j]) is_composite[j] = 1;
            j -= p;
        }
    }

    for (i = 0; i < n_small; ++i) {
        p = small_primes[i];
        start = p * p;
        if ((mid / p + 1) * p > start) start = (mid / p + 1) * p;
        for (j = start; j <= limit; j += p) {
            is_composite[j] = 1;
        }
    }

    c = 0;
    for (p = 2; p <= limit; ++p) {
        if (!is_composite[p]) ++c;
    }

    free(small_primes);
    free(is_composite);
    return c;
}