/* Criba Eratóstenes — hot path CPU para vma-methods / AntiPC */
#include <stdlib.h>
#include <string.h>
#include "../include/antipc_native.h"

static uint8_t* sieve_mark(uint32_t limit, uint32_t* count_out) {
    uint8_t* is_prime;
    uint32_t p, root, start;

    if (limit < 2) {
        *count_out = 0;
        return NULL;
    }

    is_prime = (uint8_t*)malloc((size_t)limit + 1);
    if (!is_prime) {
        *count_out = 0;
        return NULL;
    }
    memset(is_prime, 1, (size_t)limit + 1);
    is_prime[0] = is_prime[1] = 0;

    root = 1;
    while ((uint64_t)(root + 1) * (root + 1) <= limit) ++root;

    for (p = 2; p <= root; ++p) {
        if (!is_prime[p]) continue;
        for (start = p * p; start <= limit; start += p) {
            is_prime[start] = 0;
        }
    }

    *count_out = 0;
    for (p = 2; p <= limit; ++p) {
        if (is_prime[p]) ++(*count_out);
    }
    return is_prime;
}

ANTIPC_API uint32_t antipc_sieve_count(uint32_t limit) {
    uint8_t* marks;
    uint32_t count = 0;

    marks = sieve_mark(limit, &count);
    free(marks);
    return count;
}

ANTIPC_API uint32_t antipc_sieve_fill(uint32_t limit, uint32_t* out, uint32_t cap) {
    uint8_t* marks;
    uint32_t total, i, written;

    if (!out || cap == 0) return 0;

    marks = sieve_mark(limit, &total);
    if (!marks) return 0;

    written = 0;
    for (i = 2; i <= limit && written < cap; ++i) {
        if (marks[i]) out[written++] = i;
    }
    free(marks);
    return written;
}