/*
 * API propuesta antipc_native v0.4 — prototipos en Escritorio
 * No sustituye repo\antipc\...\antipc_native.h hasta integración.
 */
#ifndef ANTIPC_PORT_V04_H
#define ANTIPC_PORT_V04_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --- Cribas VMA (vma-methods/cribas.py) --- */
uint32_t port_sieve_modular6k_count(uint32_t limit);
uint32_t port_sieve_hibrida_count(uint32_t limit);

/* --- MDC trenes (mdc_lib/mdc_analyze.py) --- */
typedef struct {
    int32_t  x;
    int32_t  y;
    uint32_t s;
    uint32_t t;
    uint32_t k;
} PortMdcCollision;

typedef struct {
    uint32_t n_collisions;
    PortMdcCollision hits[64];
} PortMdcTrainResult;

void port_mdc_scan_trains(uint64_t n, PortMdcTrainResult* tx, PortMdcTrainResult* ty);

/* --- Newton rápido (vma-methods/newton.py) --- */
typedef struct {
    double j;
    int    iterations;
    int    converged;
} PortNewtonResult;

PortNewtonResult port_newton_rapido(double E, double b, double j0);

/* --- K-sweep predictivo (filestot l5/ksweep_predictiu.py) --- */
uint64_t port_mdc_ksweep_predict(uint64_t N, uint64_t m_ini, uint64_t m_fi);

#ifdef __cplusplus
}
#endif

#endif