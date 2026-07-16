/* Motor acordeón K3 — port Windows-safe de encriptacionGeometrica/k3_core.c */
#include "../include/antipc_native.h"

#ifndef ANTIPC_K3_GOLDEN
#define ANTIPC_K3_GOLDEN 0x9E3779B97F4A7C15ULL
#endif

static void k3_evolve(AntipcK3Motor* m) {
    m->L += (m->v + 1);
    m->v += 2;
    if ((5 * m->L) <= (2 * m->v + 1)) {
        m->L += (m->v * 2);
    }
}

ANTIPC_API void antipc_k3_stream_xor(uint8_t* data, size_t len, uint64_t base, uint64_t rel) {
    AntipcK3Motor motor;
    size_t i;
    uint64_t stream;

    if (!data || len == 0) return;

    motor.v = base;
    motor.L = rel;

    for (i = 0; i < len; ++i) {
        k3_evolve(&motor);
        stream = (motor.L ^ motor.v) * ANTIPC_K3_GOLDEN;
        data[i] ^= (uint8_t)(stream & 0xFF);
    }
}