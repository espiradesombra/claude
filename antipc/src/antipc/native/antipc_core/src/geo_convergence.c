/* Convergencia geométrica — aritmética entera (escala 1000, sin float) */
#include "../include/antipc_native.h"

#define GEO_SCALE 1000ULL

static void lados_figura(int tipo, uint32_t puntos, uint32_t escala, uint64_t* l0, uint64_t* l1) {
    /* base = puntos * 1.5 * GEO_SCALE */
    uint64_t base = (uint64_t)puntos * 3ULL * GEO_SCALE / 2ULL;

    if (tipo == 0) {
        uint64_t lado = base / 3ULL * escala;
        *l0 = *l1 = lado;
        return;
    }
    if (tipo == 1) {
        /* Python: [lado_a, lado_a, lado_b] — 10/11/00 usan lados[0] y lados[1] */
        uint64_t lado_a = base / 4ULL * escala;
        *l0 = lado_a;
        *l1 = lado_a;
        return;
    }
    *l0 = base * 25ULL / 100ULL * escala;
    *l1 = base * 35ULL / 100ULL * escala;
}

ANTIPC_API uint64_t antipc_geo_converge(
    const char* bits,
    size_t bits_len,
    const uint32_t* tales,
    size_t n_tales,
    const uint32_t* puntos,
    size_t n_puntos
) {
    uint64_t perimetro = 0;
    size_t idx = 0;
    size_t iter = 0;

    if (!bits || bits_len == 0 || !tales || n_tales == 0 || !puntos || n_puntos == 0) {
        return 0;
    }

    while (idx < bits_len) {
        char b0 = bits[idx];
        char b1 = (idx + 1 < bits_len) ? bits[idx + 1] : '0';
        uint32_t escala = tales[iter % n_tales];
        int fig = (int)(iter % 3);
        uint32_t pt = puntos[iter % n_puntos];
        uint64_t l0, l1;

        lados_figura(fig, pt, escala, &l0, &l1);

        if (b0 == '1' && b1 == '0') {
            perimetro += l0;
            idx += 1;
        } else if (b0 == '1' && b1 == '1') {
            perimetro += l1;
            idx += 1;
        } else if (b0 == '0' && b1 == '0') {
            perimetro += l0 + l1;
            idx += 1;
        } else if (b0 == '0' && b1 == '1') {
            idx += 2;
        } else {
            break;
        }
        ++iter;
    }

    return perimetro;
}