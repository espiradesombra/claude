/* K-sweep clásico entero — base para predictivo (ksweep_predictiu.py) */
#include <stdint.h>
#include "../include/antipc_port_v04.h"

uint64_t port_mdc_ksweep_predict(uint64_t N, uint64_t m_ini, uint64_t m_fi) {
    uint64_t pos_fi, pos_ini, k_lo, k_hi, k, candidat;

    if (m_ini < 1) m_ini = 1;
    if (m_fi < m_ini) return 0;

    pos_fi = 2 * m_fi + 3;
    pos_ini = 2 * m_ini + 3;
    k_lo = (pos_fi > 0) ? (N / pos_fi) : 1;
    if (k_lo < 1) k_lo = 1;
    k_hi = N / pos_ini;

    for (k = k_lo; k <= k_hi; ++k) {
        if (k == 0) continue;
        candidat = N / k;
        if (candidat < 3 || (candidat & 1) == 0) continue;
        if (N % candidat == 0 && candidat > 1 && candidat < N) {
            return candidat;
        }
    }
    return 0;
}