/*
 * antipc_native.h — Núcleo C unificado AntiPC (VMA)
 * K3 hash + MDC + cribas + Criva + stream geométrico
 */
#ifndef ANTIPC_NATIVE_H
#define ANTIPC_NATIVE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(_WIN64)
  #ifdef ANTIPC_BUILD_SHARED
    #define ANTIPC_API __declspec(dllexport)
  #elif defined(ANTIPC_SHARED)
    #define ANTIPC_API __declspec(dllimport)
  #else
    #define ANTIPC_API
  #endif
#else
  #if defined(ANTIPC_BUILD_SHARED) && __GNUC__ >= 4
    #define ANTIPC_API __attribute__((visibility("default")))
  #else
    #define ANTIPC_API
  #endif
#endif

#define ANTIPC_NATIVE_VERSION "0.10.0-c"
#define ANTIPC_GEO_MAX_PRIMOS 4

#define ANTIPC_NEWTON_GENERAL   0
#define ANTIPC_NEWTON_CUADRADOS 1
#define ANTIPC_NEWTON_CUBOS     2
#define ANTIPC_NEWTON_POTENCIA  3
#define ANTIPC_NEWTON_KP        4
#define ANTIPC_NEWTON_MERSENNE  5
#define ANTIPC_MDC_MAX_HITS 64

typedef struct {
    int32_t  x;
    int32_t  y;
    uint32_t s;
    uint32_t t;
    uint32_t k;
} AntipcMdcCollision;

typedef struct {
    uint32_t steps;
    uint32_t n_collisions;
    AntipcMdcCollision hits[ANTIPC_MDC_MAX_HITS];
} AntipcMdcTrainResult;

ANTIPC_API const char* antipc_native_version(void);

/* MDC toy factorization (≤ 10 dígitos) */
ANTIPC_API uint64_t antipc_mdc_factor(uint64_t n);

/* MDC dos trenes — colisiones enteras (2x+3)(2y+3)=n */
ANTIPC_API void antipc_mdc_scan_trains(
    uint64_t n,
    AntipcMdcTrainResult* train_x,
    AntipcMdcTrainResult* train_y
);

/* Eratóstenes: cuenta primos ≤ limit */
ANTIPC_API uint32_t antipc_sieve_count(uint32_t limit);

/* Criba híbrida VMA (ascendente + descendente) */
ANTIPC_API uint32_t antipc_sieve_hibrida_count(uint32_t limit);

/* Criba modular 6k±1 (vma-methods CribaModular6k) */
ANTIPC_API uint32_t antipc_sieve_modular6k_count(uint32_t limit);

/* Rellena out[] con primos ≤ limit; retorna cantidad (trunca si cap insuficiente) */
ANTIPC_API uint32_t antipc_sieve_fill(uint32_t limit, uint32_t* out, uint32_t cap);

/* Estimador Criva π(x)/x */
ANTIPC_API double antipc_criva_estimate(double x, int layers, int iterations);

/* Newton Rápido + oráculo MEcuation (vma-methods/newton.py) */
typedef struct {
    double j;
    double j_exacto;
    double error;
    int    iterations;
    int    converged;
} AntipcNewtonResult;

ANTIPC_API double antipc_newton_oraculo(
    double E, double b, int familia, int n_exp, double k_known
);
ANTIPC_API AntipcNewtonResult antipc_newton_rapido(double E, double b, double j0);
ANTIPC_API AntipcNewtonResult antipc_newton_log(
    double E, double b, int familia, int n_exp, double k_known
);

/* K-sweep MDC (ksweep_predictiu.py) */
ANTIPC_API uint64_t antipc_mdc_ksweep_classic(uint64_t N, uint64_t m_ini, uint64_t m_fi);
ANTIPC_API uint64_t antipc_mdc_ksweep_predict(
    uint64_t N, uint64_t m_ini, uint64_t m_fi, uint64_t* evals_out
);

/* Motor acordeón K3 — XOR in-place (33×1: base=33, rel=1) */
typedef struct {
    uint64_t v;
    uint64_t L;
} AntipcK3Motor;

ANTIPC_API void antipc_k3_stream_xor(uint8_t* data, size_t len, uint64_t base, uint64_t rel);

/* HASHTOOLCODE suite — fichero, redundancia Grafcet, heavy, similitud Hamming */
ANTIPC_API int antipc_k3_hash_file(const char* path, uint32_t seed, uint32_t* hash_out);
ANTIPC_API int antipc_k3_fingerprint_file(
    const char* path, int64_t* size_out, uint32_t* hash_out
);
ANTIPC_API uint32_t antipc_k3_heavy_hash(const void* data, size_t len);
ANTIPC_API int antipc_k3_redundant_hashes(
    const void* data, size_t len, int replicas, uint32_t* out, int out_cap
);
ANTIPC_API int antipc_k3_hamming(uint32_t a, uint32_t b);
ANTIPC_API double antipc_k3_similarity(uint32_t a, uint32_t b);

/* Clave geométrica Aleatorovix (gemini-code-1784158392232.py) */
typedef struct {
    uint32_t tales_count;
    uint32_t tales[8];
    uint32_t figuras_count;
    int32_t  figuras[8];   /* 0=equilatero 1=isosceles 2=escaleno */
    uint32_t puntos_count;
    uint32_t puntos[8];
    uint32_t saltos_count;
    int32_t  saltos[8];
    uint32_t primos_count;
    int32_t  primos_p1[ANTIPC_GEO_MAX_PRIMOS];
    int32_t  primos_p2[ANTIPC_GEO_MAX_PRIMOS];
    int32_t  porcentajes[ANTIPC_GEO_MAX_PRIMOS];
    int32_t  iteraciones_pi;
} AntipcGeoClave;

ANTIPC_API void antipc_geo_clave_default(AntipcGeoClave* clave);
ANTIPC_API double antipc_geo_fase_encrypt(uint16_t semilla, const AntipcGeoClave* clave);
ANTIPC_API int antipc_geo_fase_decrypt(
    double hash_fase, const AntipcGeoClave* clave, uint16_t* semilla_out
);
ANTIPC_API void antipc_aleatorovix_xor(
    const uint8_t* in, uint8_t* out, size_t len, uint16_t semilla
);
ANTIPC_API int antipc_geo_masivo_crypt(
    const uint8_t* in,
    uint8_t* out,
    size_t len,
    uint16_t semilla,
    const AntipcGeoClave* clave,
    double* hash_fase_out,
    int decrypt
);

/* Convergencia geométrica fija (escala 1000, sin float) */
ANTIPC_API uint64_t antipc_geo_converge(
    const char* bits,
    size_t bits_len,
    const uint32_t* tales,
    size_t n_tales,
    const uint32_t* puntos,
    size_t n_puntos
);

#ifdef __cplusplus
}
#endif

#endif /* ANTIPC_NATIVE_H */