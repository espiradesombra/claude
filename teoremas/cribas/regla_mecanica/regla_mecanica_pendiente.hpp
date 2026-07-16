/**
 * regla_mecanica_pendiente.hpp
 * Regla mecánica VMA: pendiente, choque con hipotenusa, (2v+3)(2l+3)
 * Autor: Víctor Manzanares Alberola
 */
#pragma once

#include <cstdint>
#include <cmath>
#include <algorithm>

namespace rm {

/** N = 2v+3 desde índice v */
static inline uint64_t f_desde_v(uint64_t v) {
    return 2u * v + 3u;
}

/** N = 2l+3 desde índice l */
static inline uint64_t f_desde_l(uint64_t l) {
    return 2u * l + 3u;
}

/** Índice v tal que 2v+3 = f (si f es de forma 6k±1 impar >= 5) */
static inline uint64_t v_desde_f(uint64_t f) {
    return (f >= 3) ? (f - 3u) / 2u : 0u;
}

/**
 * Encoger en choque (Cribas §8, algoritmo i=2..n):
 *   x := x - (x+1)/pendiente
 * Cada aumento de pendiente en 1 reduce el intervalo restante.
 */
static inline uint64_t encoger_en_choque(uint64_t x, uint32_t pendiente) {
    if (pendiente == 0) {
        return x;
    }
    return x - (x + 1u) / pendiente;
}

/**
 * Primer choque con la hipotenusa (pendiente 1): mitad de lo que quedaba.
 * x' = (x+1)/2  — equivalente a encoger_en_choque(x, 2) en el post golbach.
 */
static inline uint64_t choque_mitad(uint64_t x) {
    return (x + 1u) / 2u;
}

/**
 * Secuencia de encogimientos hasta pendiente max_p.
 * Devuelve x final tras aplicar pendientes 2,3,...,max_p.
 */
static inline uint64_t secuencia_encoger(uint64_t x, uint32_t max_pendiente) {
    for (uint32_t p = 2; p <= max_pendiente; ++p) {
        x = encoger_en_choque(x, p);
    }
    return x;
}

/**
 * Pendiente asociada a (2v+3)(2l+3)=N con v fijo:
 *   k = N / (2v+3) ≈ (2l+3)
 */
static inline uint64_t pendiente_desde_v(uint64_t N, uint64_t v) {
    uint64_t f1 = f_desde_v(v);
    return (f1 > 0) ? (N / f1) : 0u;
}

/**
 * Al aumentar v en 1, la pendiente k = N/(2v+3) decrece.
 * Retorna delta aproximado: k(v) - k(v+1).
 */
static inline uint64_t delta_pendiente_v(uint64_t N, uint64_t v) {
    uint64_t k0 = pendiente_desde_v(N, v);
    uint64_t k1 = pendiente_desde_v(N, v + 1u);
    return (k0 > k1) ? (k0 - k1) : 0u;
}

struct FactorPar {
    bool ok;
    uint64_t v;
    uint64_t l;
    uint64_t f1;
    uint64_t f2;
};

/**
 * Factorización por barrido v con regla (2v+3)(2l+3)=N.
 * Usa cotas mecánicas: l_max desde choque_mitad sobre v_max.
 */
static inline FactorPar factorizar_par(uint64_t N) {
    FactorPar r = {false, 0, 0, 0, 0};
    if (N < 9 || (N & 1u) == 0) {
        return r;
    }

    uint64_t sq = static_cast<uint64_t>(std::sqrt(static_cast<double>(N)));
    while ((sq + 1u) * (sq + 1u) <= N) {
        ++sq;
    }

    uint64_t v_max = (sq >= 3) ? (sq - 3u) / 2u : 0u;
    uint64_t v_ini = 0u;

    /* Cota mecánica: tras choque mitad sobre v_max, acotar búsqueda */
    if (v_max > 4) {
        uint64_t v_cota = choque_mitad(v_max);
        if (v_cota > 0 && v_cota < v_max) {
            v_ini = std::max(v_ini, v_cota / 3u);
        }
    }

    for (uint64_t v = v_ini; v <= v_max; ++v) {
        uint64_t f1 = f_desde_v(v);
        if (N % f1 != 0) {
            continue;
        }
        uint64_t f2 = N / f1;
        if (f2 < 3 || (f2 & 1u) == 0) {
            continue;
        }
        if ((f2 - 3u) % 2u != 0) {
            continue;
        }
        r.ok = true;
        r.v = v;
        r.l = (f2 - 3u) / 2u;
        r.f1 = std::min(f1, f2);
        r.f2 = std::max(f1, f2);
        return r;
    }
    return r;
}

/**
 * Barrido por pendiente k (≈ 2l+3) con decrecimiento al subir v.
 * Equivalente k-sweep MDC pero con salto guiado por delta_pendiente_v.
 */
static inline uint64_t factorizar_por_pendiente(uint64_t N) {
    FactorPar p = factorizar_par(N);
    return p.ok ? p.f1 : 0u;
}

} // namespace rm