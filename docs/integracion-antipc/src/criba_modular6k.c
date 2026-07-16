/* Criba Modular 6k±1 — port de vma-methods/cribas.py CribaModular6k.run() */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "../include/antipc_port_v04.h"

static uint8_t* mark_sieve(uint32_t limit, uint32_t* count) {
    uint8_t* is_prime;
    uint32_t p, root;

    if (limit < 2) {
        *count = 0;
        return NULL;
    }
    is_prime = (uint8_t*)malloc((size_t)limit + 1);
    if (!is_prime) {
        *count = 0;
        return NULL;
    }
    memset(is_prime, 1, (size_t)limit + 1);
    is_prime[0] = is_prime[1] = 0;

    root = 1;
    while ((uint64_t)(root + 1) * (root + 1) <= limit) {
        ++root;
    }
    for (p = 2; p <= root; ++p) {
        uint32_t start;
        if (!is_prime[p]) continue;
        for (start = p * p; start <= limit; start += p) {
            is_prime[start] = 0;
        }
    }
    *count = 0;
    for (p = 2; p <= limit; ++p) {
        if (is_prime[p]) ++(*count);
    }
    return is_prime;
}

uint32_t port_sieve_modular6k_count(uint32_t limit) {
    uint8_t* m;
    uint32_t c = 0;
    m = mark_sieve(limit, &c);
    free(m);
    return c;
}