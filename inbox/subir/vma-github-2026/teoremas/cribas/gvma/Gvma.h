/**
 * @file Gvma.h
 * @brief SunRaman Jet — criba segmentada 6k±1 (VMA / Víctor Manzanares Alberola)
 *
 * v = 2nm + 3n + 3m + 3  <=>  (2n+3)(2m+3) = 2v+3
 * Bloques 1MB con bit-packing; factorizador Flash (pinza).
 */

#ifndef GVMA_H
#define GVMA_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GVMA_BLOC_BYTES 1048576u
#define GVMA_V_PER_BLOC (GVMA_BLOC_BYTES * 8u)

typedef struct {
    bool es_compuesto;
    uint64_t n;
    uint64_t m;
    uint64_t f1;
    uint64_t f2;
} GvmaFlashResult;

static inline uint64_t gvma_v_a_n(uint64_t v) {
    return 2u * v + 3u;
}

static inline void gvma_marca_compuesto(uint8_t *mapa, uint64_t v_local) {
    mapa[v_local >> 3] |= (uint8_t)(1u << (v_local & 7u));
}

static inline bool gvma_es_compuesto(const uint8_t *mapa, uint64_t v_local) {
    return (mapa[v_local >> 3] & (uint8_t)(1u << (v_local & 7u))) != 0;
}

static inline int gvma_rueda_paso(int rueda) {
    return (rueda == 4) ? 1 : 2;
}

static inline int gvma_rueda_toggle(int rueda) {
    return 6 - rueda;
}

/**
 * Criba un segmento [v_inici, v_fi] sobre mapa de 1MB (bit-packed).
 * mapa se inicializa a 0 (todos candidatos a primo).
 */
static inline void gvma_criba_segmentada(uint8_t *mapa, uint64_t v_inici, uint64_t v_fi) {
    if (!mapa || v_fi < v_inici) {
        return;
    }

    memset(mapa, 0, GVMA_BLOC_BYTES);

    uint64_t n_max = (2u * v_fi + 3u) / 5u;
    int rueda_n = 4;

    for (uint64_t n = 1; n <= n_max; ) {
        uint64_t p = 2u * n + 3u;

        int64_t m_calc = (int64_t)((v_inici + 3u * n + 3u + p - 1u) / p);
        uint64_t m = (m_calc < (int64_t)n) ? n : (uint64_t)m_calc;

        int rueda_m = (m == n) ? rueda_n : 4;

        while (1) {
            uint64_t v = 2u * n * m + 3u * n + 3u * m + 3u;
            if (v > v_fi) {
                break;
            }
            if (v >= v_inici) {
                gvma_marca_compuesto(mapa, v - v_inici);
            }

            if (m == 0) {
                break;
            }
            m -= (uint64_t)gvma_rueda_paso(rueda_m);
            rueda_m = gvma_rueda_toggle(rueda_m);
        }

        n += (uint64_t)gvma_rueda_paso(rueda_n);
        rueda_n = gvma_rueda_toggle(rueda_n);
    }
}

/**
 * Variante Jet en memoria densa (1 byte por v) para rangos pequeños.
 */
static inline void gvma_criba_jet(uint8_t *mapa, uint64_t limite_v) {
    if (!mapa) {
        return;
    }

    memset(mapa, 0, (size_t)(limite_v + 1u));

    uint64_t objetivo = 2u * limite_v + 3u;
    uint64_t limite_raiz = (uint64_t)sqrt((double)objetivo);
    bool jet_activado = false;
    int rueda_n = 4;

    for (uint64_t n = 1; (2u * n + 3u) <= objetivo / 5u; ) {
        uint64_t f1 = 2u * n + 3u;
        uint64_t m = n;
        int rueda_m = rueda_n;

        if (f1 > limite_raiz && !jet_activado) {
            uint64_t k = objetivo / f1;
            m = (k >= 3u) ? (k - 3u - ((k - 3u) % 2u)) / 2u : 0u;
            jet_activado = true;
        }

        while (m >= 1) {
            uint64_t v = 2u * n * m + 3u * n + 3u * m + 3u;
            if (v <= limite_v) {
                mapa[v] = 1;
            }
            m -= (uint64_t)gvma_rueda_paso(rueda_m);
            rueda_m = gvma_rueda_toggle(rueda_m);
        }

        n += (uint64_t)gvma_rueda_paso(rueda_n);
        rueda_n = gvma_rueda_toggle(rueda_n);
    }
}

static inline bool gvma_guardar_bloc(const char *dir, uint32_t index_bloc, const uint8_t *mapa) {
    char ruta[512];
    snprintf(ruta, sizeof(ruta), "%s/bloc_%u.bin", dir, index_bloc);
    FILE *f = fopen(ruta, "wb");
    if (!f) {
        return false;
    }
    size_t n = fwrite(mapa, 1, GVMA_BLOC_BYTES, f);
    fclose(f);
    return n == GVMA_BLOC_BYTES;
}

static inline uint32_t gvma_ultimo_bloc_index(const char *dir) {
    uint32_t ultim = 0;
    char ruta[512];
    while (1) {
        snprintf(ruta, sizeof(ruta), "%s/bloc_%u.bin", dir, ultim);
        FILE *test = fopen(ruta, "rb");
        if (!test) {
            break;
        }
        fclose(test);
        ultim++;
    }
    return ultim;
}

static inline bool gvma_afegir_bloc_incremental(const char *dir) {
    uint32_t idx = gvma_ultimo_bloc_index(dir);
    uint8_t *mapa = (uint8_t *)malloc(GVMA_BLOC_BYTES);
    if (!mapa) {
        return false;
    }

    uint64_t v_start = (uint64_t)idx * GVMA_V_PER_BLOC;
    uint64_t v_end = v_start + GVMA_V_PER_BLOC - 1u;

    gvma_criba_segmentada(mapa, v_start, v_end);
    bool ok = gvma_guardar_bloc(dir, idx, mapa);
    free(mapa);
    return ok;
}

static inline GvmaFlashResult gvma_factorizar_flash(uint64_t v) {
    GvmaFlashResult res = {false, 0, 0, 0, 0};
    double raiz_v = sqrt((double)v);

    uint64_t m_max = (uint64_t)(raiz_v * 0.25);
    uint64_t n_ini = (uint64_t)(raiz_v * 0.50);
    uint64_t n_max = (uint64_t)(raiz_v * 0.75);

    for (uint64_t m = 0; m <= m_max; m++) {
        for (uint64_t n = n_ini; n <= n_max; n++) {
            uint64_t v_test = 2u * n * m + 3u * n + 3u * m + 3u;
            if (v_test == v) {
                res.es_compuesto = true;
                res.n = n;
                res.m = m;
                res.f1 = 2u * n + 3u;
                res.f2 = 2u * m + 3u;
                return res;
            }
            if (v_test > v) {
                break;
            }
        }
    }
    return res;
}

#endif /* GVMA_H */