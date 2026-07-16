/* Newton Rápido — port fiel de vma-methods/newton.py */
#include <math.h>
#include "../include/antipc_port_v04.h"

#ifndef PORT_NEWTON_FE
#define PORT_NEWTON_FE  1.8
#define PORT_NEWTON_CTE (-0.85)
#define PORT_NEWTON_SBP 1.3
#define PORT_NEWTON_TOL 1e-12
#define PORT_NEWTON_MAX 200
#endif

PortNewtonResult port_newton_rapido(double E, double b, double j0) {
    PortNewtonResult out;
    double j, j_prev, j_nuevo, d1 = 0, d2 = 0, d3 = 0, d4 = 0;
    double bpow, diff, denom, exp_val;
    int it;

    out.converged = 0;
    out.iterations = 0;
    out.j = j0;

    if (E <= 0 || b <= 0 || b == 1.0) {
        out.j = 0;
        return out;
    }
    if (j0 == 0.0) j0 = 1.0;
    j = j0;
    j_prev = j;

    for (it = 1; it <= PORT_NEWTON_MAX; ++it) {
        bpow = pow(b, -j);
        if (!isfinite(bpow)) bpow = 0.0;

        j_nuevo = (E + E * bpow - 1.0) / E;

        if (isfinite(pow(b, j)) && pow(b, j) >= sqrt(E)) {
            j_nuevo /= PORT_NEWTON_SBP;
        }

        diff = j_nuevo - j;
        d4 = d3;
        d3 = d2;
        d2 = d1;
        d1 = diff;
        j = j_nuevo;

        denom = d3 - d1;
        if (fabs(denom) > 1e-30 && fabs(d3 - d2) > 1e-30) {
            exp_val = fabs(PORT_NEWTON_CTE + PORT_NEWTON_FE * (d3 - d2) / denom);
            if (exp_val > 0 && j > 0) {
                double jumped = pow(j, exp_val);
                if (isfinite(jumped)) j = jumped;
            }
        }

        out.iterations = it;
        if (fabs(j - j_prev) < PORT_NEWTON_TOL) {
            out.converged = 1;
            break;
        }
        j_prev = j;
    }

    out.j = j;
    return out;
}