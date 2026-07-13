/**
 * ============================================================================
 * MÒDUL: factoritzacio_mdc.hpp
 * Nucli del Mètode Diofàntic Cinemàtic: k-sweep i factorizar_mdc complet.
 * Depèn de: roda_p210.hpp
 * ============================================================================
 */
#pragma once
#include <cstdint>
#include <cmath>
#include <algorithm>
#include "roda_p210.hpp"

namespace mdc {

/**
 * k_sweep_mdc(N, m_ini, m_fi): cerca factors de N en el rang D=2m+3
 * iterant sobre les pendents k=N//(2m+3) en lloc dels valors m.
 * Determinista i exhaustiu. Complexitat O(m_fi - m_ini).
 */
static inline uint64_t k_sweep_mdc(uint64_t N, uint64_t m_ini, uint64_t m_fi) {
    if (m_ini < 1) m_ini = 1;
    if (m_fi < m_ini) return 0;

    uint64_t pos_fi  = 2 * m_fi  + 3;
    uint64_t pos_ini = 2 * m_ini + 3;

    uint64_t k_lo = (pos_fi  > 0) ? std::max(uint64_t(1), N / pos_fi ) : 1;
    uint64_t k_hi =                  N / pos_ini;

    for (uint64_t k = k_lo; k <= k_hi; ++k) {
        if (k == 0) continue;
        uint64_t candidat = N / k;
        if (candidat < 3 || (candidat & 1) == 0) continue;
        if (N % candidat == 0 && candidat > 1 && candidat < N)
            return candidat;
    }
    return 0;
}

/**
 * factorizar_mdc(N): factoritza un semiprim senar usant el k-sweep MDC.
 * Estratègia: primers petits → criba roda p210 → k-sweep zona densa
 * (0.859·m_max, derivat de (e-1)/2) → k-sweep zona baixa.
 * Retorna el factor menut p (o 0 si N és primer / fora d'abast).
 */
static inline uint64_t factorizar_mdc(uint64_t N) {
    if (N < 4) return 0;
    if ((N & 1) == 0) return 2;

    static const uint64_t PETITS[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    for (uint64_t p : PETITS) {
        if (N % p == 0 && p < N) return p;
    }

    uint64_t sq = static_cast<uint64_t>(std::sqrt(double(N)));
    while ((sq + 1) * (sq + 1) <= N) ++sq;
    if (sq * sq == N) return sq;

    uint64_t m_max = (sq - 3) / 2;
    if (m_max < 1) return 0;

    const uint64_t LIM_CRIBA = std::min(m_max, uint64_t(500'000));
    {
        // Índex 1 perquè RESIDUS_210[1]==11, el primer candidat de la criba.
        // Salts ja són diferències absolutes entre residus (no *2).
        uint64_t D = 11;
        int idx = 1;
        while (D <= 2 * LIM_CRIBA + 3) {
            if (N % D == 0 && D < N) return D;
            D += SALTS_RODA_210[idx % 48];
            idx = (idx + 1) % 48;
        }
    }

    if (LIM_CRIBA >= m_max) return 0;

    constexpr double E_MINUS_1_OVER_2 = 0.8591409142295227;
    uint64_t m_conv = static_cast<uint64_t>(E_MINUS_1_OVER_2 * double(m_max));
    if (m_conv < LIM_CRIBA) m_conv = LIM_CRIBA;

    uint64_t f = k_sweep_mdc(N, m_conv, m_max);
    if (f) return std::min(f, N / f);

    f = k_sweep_mdc(N, LIM_CRIBA, m_conv);
    if (f) return std::min(f, N / f);

    return 0;
}

} // namespace mdc
