/* Newton Rápido + oráculo MEcuation — vma-methods/newton.py */
#include <math.h>
#include "../include/antipc_native.h"

#define NW_FE   1.8
#define NW_CTE (-0.85)
#define NW_SBP  1.3
#define NW_TOL  1e-12
#define NW_MAX  200

ANTIPC_API double antipc_newton_oraculo(
    double E, double b, int familia, int n_exp, double k_known
) {
    double lnE, lnb;

    if (E <= 0 || b <= 0 || b == 1.0) return 1.0;
    lnE = log(E);
    lnb = log(b);

    switch (familia) {
        case 1: return lnE / 2.0 / lnb;
        case 2: return lnE / 3.0 / lnb;
        case 3:
            if (n_exp < 1) n_exp = 2;
            return lnE / (double)n_exp / lnb;
        case 4:
            if (k_known <= 0) k_known = 1.0;
            return log(E / k_known) / lnb;
        case 5:
            return log2(E + 1.0) * log(2.0) / lnb;
        default: return 1.0;
    }
}

ANTIPC_API AntipcNewtonResult antipc_newton_rapido(double E, double b, double j0) {
    AntipcNewtonResult out;
    double j, j_prev, j_nuevo, j_ref;
    double d1 = 0, d2 = 0, d3 = 0, d4 = 0;
    double bpow, diff, denom, exp_val;
    int it;

    out.j = 0;
    out.j_exacto = 0;
    out.error = 0;
    out.iterations = 0;
    out.converged = 0;

    if (E <= 0 || b <= 0 || b == 1.0) return out;
    if (j0 == 0.0) j0 = 1.0;

    j = j0;
    j_prev = j;
    j_ref = log(E) / log(b);
    out.j_exacto = j_ref;

    for (it = 1; it <= NW_MAX; ++it) {
        bpow = pow(b, -j);
        if (!isfinite(bpow)) bpow = 0.0;

        j_nuevo = (E + E * bpow - 1.0) / E;
        if (isfinite(pow(b, j)) && pow(b, j) >= sqrt(E)) {
            j_nuevo /= NW_SBP;
        }

        diff = j_nuevo - j;
        d4 = d3;
        d3 = d2;
        d2 = d1;
        d1 = diff;
        j = j_nuevo;

        denom = d3 - d1;
        if (fabs(denom) > 1e-30 && fabs(d3 - d2) > 1e-30) {
            exp_val = fabs(NW_CTE + NW_FE * (d3 - d2) / denom);
            if (exp_val > 0 && j > 0) {
                double jumped = pow(j, exp_val);
                if (isfinite(jumped)) j = jumped;
            }
        }

        out.iterations = it;
        if (fabs(j - j_prev) < NW_TOL) {
            out.converged = 1;
            break;
        }
        j_prev = j;
    }

    out.j = j;
    out.error = fabs(j - j_ref);
    return out;
}

ANTIPC_API AntipcNewtonResult antipc_newton_log(
    double E, double b, int familia, int n_exp, double k_known
) {
    double j0 = antipc_newton_oraculo(E, b, familia, n_exp, k_known);
    return antipc_newton_rapido(E, b, j0);
}