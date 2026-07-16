/* Criba Modular 6k±1 — vma-methods/cribas.py CribaModular6k.run() */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/antipc_native.h"

static uint32_t sieve_mod6k_count(uint32_t limit) {
    uint8_t* is_prime;
    uint32_t p, root, count = 0;

    if (limit < 2) {
        return 0;
    }

    is_prime = (uint8_t*)malloc((size_t)limit + 1);
    if (!is_prime) {
        return 0;
    }
    memset(is_prime, 1, (size_t)limit + 1);
    is_prime[0] = is_prime[1] = 0;

    root = 1;
    while ((uint64_t)(root + 1) * (root + 1) <= limit) {
        ++root;
    }

    for (p = 2; p <= root; ++p) {
        uint32_t start;
        if (!is_prime[p]) {
            continue;
        }
        for (start = p * p; start <= limit; start += p) {
            is_prime[start] = 0;
        }
    }

    for (p = 2; p <= limit; ++p) {
        if (is_prime[p]) {
            ++count;
        }
    }
    free(is_prime);
    return count;
}

ANTIPC_API uint32_t antipc_sieve_modular6k_count(uint32_t limit) {
    return sieve_mod6k_count(limit);
}