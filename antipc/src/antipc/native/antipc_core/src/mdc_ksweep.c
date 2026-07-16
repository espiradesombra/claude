/* K-sweep clásico + predictivo entero — ksweep_predictiu.py */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/antipc_native.h"

#define KS_RADI_PRED 4

static int64_t iabs64(int64_t x) {
    return x < 0 ? -x : x;
}

static uint64_t isqrt_u64(uint64_t n) {
    uint64_t x = n;
    uint64_t y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (x + n / x) / 2;
    }
    return x;
}

static int isqrt_disc(int64_t disc, int64_t* sq_out) {
    uint64_t sq;
    int64_t try_sq;
    if (disc < 0) return 0;
    sq = isqrt_u64((uint64_t)disc);
    for (try_sq = (int64_t)sq; try_sq <= (int64_t)sq + 1; ++try_sq) {
        if (try_sq >= 0 && try_sq * try_sq <= disc) {
            *sq_out = try_sq;
            return 1;
        }
    }
    return 0;
}

ANTIPC_API uint64_t antipc_mdc_ksweep_classic(uint64_t N, uint64_t m_ini, uint64_t m_fi) {
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
        if (N % candidat == 0 && candidat > 1 && candidat < N) return candidat;
    }
    return 0;
}

ANTIPC_API uint64_t antipc_mdc_ksweep_predict(
    uint64_t N, uint64_t m_ini, uint64_t m_fi, uint64_t* evals_out
) {
    uint64_t n_avals = 0;
    uint64_t m;
    int64_t rs[4], Ds[4];
    size_t i, j, n_cand;
    int64_t candidats[8];

    if (evals_out) *evals_out = 0;
    if (m_ini < 1) m_ini = 1;
    if (m_fi < m_ini || N < 4) return 0;

    m = m_fi;

    while (m >= m_ini) {
        if (m < m_ini + 3) {
            int64_t mm;
            for (mm = (int64_t)m; mm >= (int64_t)m_ini; --mm) {
                uint64_t D;
                ++n_avals;
                D = (uint64_t)(2 * mm + 3);
                if (D > 1 && D < N && N % D == 0) {
                    if (evals_out) *evals_out = n_avals;
                    return D;
                }
            }
            if (evals_out) *evals_out = n_avals;
            return 0;
        }

        memset(rs, 0, sizeof(rs));
        memset(Ds, 0, sizeof(Ds));

        for (i = 0; i < 4; ++i) {
            int64_t mi = (int64_t)m - (int64_t)i;
            uint64_t Di, twoDi;
            ++n_avals;
            Di = (uint64_t)(2 * mi + 3);
            if (Di > 1 && Di < N && N % Di == 0) {
                if (evals_out) *evals_out = n_avals;
                return Di;
            }
            twoDi = 2 * Di;
            rs[i] = (int64_t)(N % twoDi);
            Ds[i] = (int64_t)Di;
        }

        {
            int64_t r0 = rs[0], D0 = Ds[0];
            int64_t distancia = r0 - D0;
            int64_t vs_r[3], accs_r[2];
            int64_t v_r, a_r;
            uint64_t thresh;
            int m_saltat = 0;

            vs_r[0] = rs[1] - rs[0];
            vs_r[1] = rs[2] - rs[1];
            vs_r[2] = rs[3] - rs[2];
            accs_r[0] = vs_r[1] - vs_r[0];
            accs_r[1] = vs_r[2] - vs_r[1];

            thresh = (uint64_t)(D0 / 3);
            if ((uint64_t)iabs64(vs_r[0]) > thresh ||
                (uint64_t)iabs64(vs_r[1]) > thresh ||
                (uint64_t)iabs64(vs_r[2]) > thresh) {
                m -= 4;
                continue;
            }

            v_r = (vs_r[0] + vs_r[1] + vs_r[2]) / 3;
            a_r = (accs_r[0] + accs_r[1]) / 2;

            n_cand = 0;
            if (a_r == 0) {
                if (v_r != 0 && distancia * v_r < 0) {
                    double p_float = -(double)distancia / (double)v_r;
                    if (p_float > 0 && n_cand < 8) {
                        candidats[n_cand++] = (int64_t)(p_float + 0.5);
                    }
                }
            } else {
                int64_t disc = 4 * v_r * v_r - 4 * a_r * (2 * distancia);
                int64_t sq;
                if (isqrt_disc(disc, &sq)) {
                    for (j = 0; j < 2; ++j) {
                        int64_t sgn = (j == 0) ? 1 : -1;
                        int64_t num = -2 * v_r + sgn * sq;
                        if (num * a_r > 0) {
                            double p_float = (double)num / (2.0 * (double)a_r);
                            if (p_float > 0 && n_cand < 8) {
                                candidats[n_cand++] = (int64_t)(p_float + 0.5);
                            }
                        }
                    }
                }
            }

            for (i = 0; i < n_cand; ++i) {
                int64_t p_cand = candidats[i];
                int64_t dm;
                if ((int64_t)m - p_cand - KS_RADI_PRED < 4) continue;

                for (dm = -KS_RADI_PRED; dm <= KS_RADI_PRED; ++dm) {
                    int64_t m_t = (int64_t)m - p_cand + dm;
                    uint64_t D_t;
                    if (m_t < (int64_t)m_ini || m_t > (int64_t)m_fi) continue;
                    ++n_avals;
                    D_t = (uint64_t)(2 * m_t + 3);
                    if (D_t > 1 && D_t < N && N % D_t == 0) {
                        if (evals_out) *evals_out = n_avals;
                        return D_t;
                    }
                }
                if (p_cand + KS_RADI_PRED + 1 < (int64_t)m) {
                    m = (uint64_t)((int64_t)m - p_cand - KS_RADI_PRED - 1);
                } else {
                    m = 0;
                }
                m_saltat = 1;
                break;
            }

            if (!m_saltat) m -= 4;
        }
    }

    if (evals_out) *evals_out = n_avals;
    return 0;
}