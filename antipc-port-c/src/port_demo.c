/* Demo CLI — compilar y ejecutar sin tocar el repo */
#include <stdio.h>
#include <stdint.h>
#include "../include/antipc_port_v04.h"

int main(void) {
    uint32_t lim = 50000;
    uint64_t n = 1147;
    PortMdcTrainResult tx, ty;
    PortNewtonResult nw;
    uint32_t i;

    printf("=== antipc-port-c (Escritorio) ===\n");
    printf("Criba 6k count(50000)  = %u\n", port_sieve_modular6k_count(lim));
    printf("Criba hibrida(50000)   = %u\n", port_sieve_hibrida_count(lim));
    printf("MDC trains n=%llu\n", (unsigned long long)n);
    port_mdc_scan_trains(n, &tx, &ty);
    printf("  Tren X colisiones: %u\n", tx.n_collisions);
    for (i = 0; i < tx.n_collisions; ++i) {
        printf("    (%d,%d) S=%u T=%u\n", tx.hits[i].x, tx.hits[i].y, tx.hits[i].s, tx.hits[i].t);
    }
    printf("  Tren Y colisiones: %u\n", ty.n_collisions);
    /* j0 oráculo cuadrados: ln(E)/2/ln(b) — newton.py familia cuadrados */
    nw = port_newton_rapido(121.0, 10.0, 1.04381229068);
    printf("Newton log10(121) j=%.12f iter=%d ok=%d\n", nw.j, nw.iterations, nw.converged);
    printf("K-sweep(1147,1,50)     = %llu\n",
           (unsigned long long)port_mdc_ksweep_predict(1147, 1, 50));
    return 0;
}