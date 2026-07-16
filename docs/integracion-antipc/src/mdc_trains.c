/* MDC dos trenes — aritmética entera, sin Fraction (mdc_analyze.py) */
#include <stdint.h>
#include <stdlib.h>
#include "../include/antipc_port_v04.h"

static int64_t diofantic_y_num(int64_t n, int64_t x, int64_t* den_out) {
    int64_t den = 4 * x + 6;
    if (den == 0) return 0;
    *den_out = den;
    return n - 6 * x - 9;
}

static int64_t diofantic_x_num(int64_t n, int64_t y, int64_t* den_out) {
    int64_t den = 4 * y + 6;
    if (den == 0) return 0;
    *den_out = den;
    return n - 6 * y - 9;
}

static int try_collision(int64_t xi, int64_t yi, uint64_t n,
                         PortMdcTrainResult* out, const char* source) {
    uint32_t s, t;
    PortMdcCollision c;
    uint32_t i;

    if (xi < 0 || yi < 0) return 0;
    s = (uint32_t)(2 * xi + 3);
    t = (uint32_t)(2 * yi + 3);
    if (s == 0 || t == 0 || (uint64_t)s * t != n) return 0;

    for (i = 0; i < out->n_collisions; ++i) {
        if (out->hits[i].s == s && out->hits[i].t == t) return 0;
        if (out->hits[i].s == t && out->hits[i].t == s) return 0;
    }
    if (out->n_collisions >= 64) return 0;

    c.x = (int32_t)xi;
    c.y = (int32_t)yi;
    c.s = s < t ? s : t;
    c.t = s < t ? t : s;
    c.k = (c.t - c.s) / 2;
    out->hits[out->n_collisions++] = c;
    (void)source;
    return 1;
}

static void scan_x(uint64_t n, PortMdcTrainResult* out) {
    int64_t n_i = (int64_t)n;
    int64_t x, den, num, y;
    int64_t x_max = (int64_t)(n / 2) + 1;
    if (x_max < 2) x_max = 2;

    out->n_collisions = 0;
    for (x = -3; x <= x_max; ++x) {
        num = diofantic_y_num(n_i, x, &den);
        if (den == 0) continue;
        if (num % den != 0) continue;
        y = num / den;
        try_collision(x, y, n, out, "tren-X");
    }
}

static void scan_y(uint64_t n, PortMdcTrainResult* out) {
    int64_t n_i = (int64_t)n;
    int64_t y, den, num, x;
    int64_t y_max = (int64_t)(n / 2) + 1;
    if (y_max < 2) y_max = 2;

    out->n_collisions = 0;
    for (y = -3; y <= y_max; ++y) {
        num = diofantic_x_num(n_i, y, &den);
        if (den == 0) continue;
        if (num % den != 0) continue;
        x = num / den;
        try_collision(x, y, n, out, "tren-Y");
    }
}

void port_mdc_scan_trains(uint64_t n, PortMdcTrainResult* tx, PortMdcTrainResult* ty) {
    if (!tx || !ty || n < 4) {
        if (tx) tx->n_collisions = 0;
        if (ty) ty->n_collisions = 0;
        return;
    }
    scan_x(n, tx);
    scan_y(n, ty);
}