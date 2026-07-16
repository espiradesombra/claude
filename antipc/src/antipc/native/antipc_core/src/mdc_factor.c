/* MDC factorización toy — port de mdc_lib/factoritzacio.py y factoritzacio_mdc.hpp */
#include <math.h>
#include <stdint.h>
#include "../include/antipc_native.h"

static const uint64_t PETITS[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
};
static const unsigned PETITS_N = sizeof(PETITS) / sizeof(PETITS[0]);

static const uint64_t SALTS_RODA_210[48] = {
    2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2,
    6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 2, 4, 6, 2,
    6, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2
};

static uint64_t k_sweep_mdc(uint64_t n, uint64_t m_ini, uint64_t m_fi) {
    uint64_t pos_fi, pos_ini, k_lo, k_hi, k, candidat;

    if (m_ini < 1) m_ini = 1;
    if (m_fi < m_ini) return 0;

    pos_fi = 2 * m_fi + 3;
    pos_ini = 2 * m_ini + 3;
    k_lo = (pos_fi > 0) ? (n / pos_fi) : 1;
    if (k_lo < 1) k_lo = 1;
    k_hi = n / pos_ini;

    for (k = k_lo; k <= k_hi; ++k) {
        if (k == 0) continue;
        candidat = n / k;
        if (candidat < 3 || (candidat & 1) == 0) continue;
        if (n % candidat == 0 && candidat > 1 && candidat < n) return candidat;
    }
    return 0;
}

ANTIPC_API uint64_t antipc_mdc_factor(uint64_t n) {
    uint64_t sq, m_max, lim_criba, m_conv, f, d;
    unsigned idx, i;

    if (n > 9999999999ULL) return 0;
    if (n < 4) return 0;
    if ((n & 1) == 0) return 2;

    for (i = 0; i < PETITS_N; ++i) {
        if (n % PETITS[i] == 0 && PETITS[i] < n) return PETITS[i];
    }

    sq = (uint64_t)sqrt((double)n);
    while ((sq + 1) * (sq + 1) <= n) ++sq;
    if (sq * sq == n) return sq;

    m_max = (sq - 3) / 2;
    if (m_max < 1) return 0;

    lim_criba = m_max < 500000ULL ? m_max : 500000ULL;
    d = 11;
    idx = 1;
    while (d <= 2 * lim_criba + 3) {
        if (n % d == 0 && d < n) return d;
        d += SALTS_RODA_210[idx % 48];
        ++idx;
    }

    if (lim_criba >= m_max) return 0;

    m_conv = (uint64_t)(0.8591409142295227 * (double)m_max);
    if (m_conv < lim_criba) m_conv = lim_criba;

    f = k_sweep_mdc(n, m_conv, m_max);
    if (f) return (f < n / f) ? f : n / f;

    f = k_sweep_mdc(n, lim_criba, m_conv);
    if (f) return (f < n / f) ? f : n / f;

    return 0;
}