/*
 * k3_suite.c — HASHTOOLCODE + Grafcet: fichero, redundancia, heavy, Hamming
 */
#include "../include/antipc_native.h"
#include "k3hash.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const uint32_t K3_REDUNDANT_SEEDS[3] = {
    0xA5A5A5A5u, 0x5A5A5A5Au, 0x12345678u
};

#define K3_REDUNDANT_MAX 16

static K3HashConfig k3_cfg_with_seed(uint32_t seed) {
    K3HashConfig cfg = k3_config_default();
    if (seed != 0u) {
        cfg.semilla_inicial = seed;
    }
    return cfg;
}

ANTIPC_API int antipc_k3_hash_file(const char* path, uint32_t seed, uint32_t* hash_out) {
    if (!path || !hash_out) {
        return -1;
    }
    K3HashConfig cfg = k3_cfg_with_seed(seed);
    return k3_hash_file(path, &cfg, hash_out);
}

ANTIPC_API int antipc_k3_fingerprint_file(
    const char* path,
    int64_t* size_out,
    uint32_t* hash_out
) {
    if (!path) {
        return -1;
    }

    FILE* f = fopen(path, "rb");
    if (!f) {
        return -1;
    }
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return -1;
    }
    long tam = ftell(f);
    fclose(f);
    if (tam < 0) {
        return -1;
    }

    uint32_t h = 0u;
    if (antipc_k3_hash_file(path, 0u, &h) != 0) {
        return -1;
    }

    if (size_out) {
        *size_out = (int64_t)tam;
    }
    if (hash_out) {
        *hash_out = h;
    }
    return 0;
}

ANTIPC_API uint32_t antipc_k3_heavy_hash(const void* data, size_t len) {
    if (!data && len > 0) {
        return 0u;
    }

    uint32_t digest = k3_hash_buffer(data, len, NULL);
    if (len == 0) {
        return digest;
    }

    /* 4 rondas: hash(digest_le || payload) — grafcet.py _heavy_hash */
    uint8_t* buf = (uint8_t*)malloc(len + 4u);
    if (!buf) {
        return digest;
    }

    for (int round = 0; round < 4; round++) {
        memcpy(buf, &digest, 4);
        memcpy(buf + 4, data, len);
        digest = k3_hash_buffer(buf, len + 4u, NULL);
    }

    free(buf);
    return digest;
}

ANTIPC_API int antipc_k3_redundant_hashes(
    const void* data,
    size_t len,
    int replicas,
    uint32_t* out,
    int out_cap
) {
    if (!out || out_cap <= 0 || replicas <= 0) {
        return -1;
    }
    if (!data && len > 0) {
        return -1;
    }

    int n = replicas;
    if (n > out_cap) {
        n = out_cap;
    }
    if (n > K3_REDUNDANT_MAX) {
        n = K3_REDUNDANT_MAX;
    }

    for (int i = 0; i < n; i++) {
        K3HashConfig cfg = k3_cfg_with_seed(K3_REDUNDANT_SEEDS[i % 3]);
        out[i] = k3_hash_buffer(data, len, &cfg);
    }
    return n;
}

ANTIPC_API int antipc_k3_hamming(uint32_t a, uint32_t b) {
    uint32_t x = a ^ b;
    int count = 0;
    while (x) {
        count += (int)(x & 1u);
        x >>= 1;
    }
    return count;
}

ANTIPC_API double antipc_k3_similarity(uint32_t a, uint32_t b) {
    int dist = antipc_k3_hamming(a, b);
    return 1.0 - ((double)dist / 32.0);
}