/*
 * k3_audit.h — Sistema de auditoría K3 industrial (Theorem K3 v0)
 * VMA — Víctor Manzanares Alberola
 *
 * NO criptográfico. Integridad industrial / telemetría.
 */

#ifndef K3_AUDIT_H
#define K3_AUDIT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define K3_AUDIT_MAX_REGISTROS 10
#define K3_AUDIT_MAX_MARCAS    4
#define K3_AUDIT_TELEMETRIA_MAX 512

typedef struct {
    int      bits_ancho;
    int      bytes_ancho;
    int      desfase_stride;
    int      num_registros;
    uint32_t semilla_inicial;
    uint32_t marca_inicial;
    uint32_t constante_disp;
    int      total_marcas_encoladas;
    struct {
        uint32_t valor_bits;
        uint8_t  banderita;
        int      activa;
    } cola_marcas[K3_AUDIT_MAX_MARCAS];
} K3AuditConfig;

typedef struct {
    uint32_t hash_final;
    size_t   puntos_telemetria;
    uint32_t estado[K3_AUDIT_MAX_REGISTROS];
} K3AuditResult;

extern const char K3_FIRMA_LOGICA[];

K3AuditConfig k3_audit_config_default(void);

K3AuditConfig k3_audit_crear_config(int bits_lectura, int desfase_trozos,
                                    int n_registros, uint32_t semilla,
                                    uint32_t marca_inicial);

void k3_audit_encolar_marca(K3AuditConfig* config, uint32_t valor_bits, uint8_t banderita);

void k3_audit_encolar_marcas_default(K3AuditConfig* config);

uint32_t k3_audit_hash_buffer(const uint8_t* datos, size_t longitud,
                              const K3AuditConfig* config,
                              uint8_t* telemetria_out, size_t telemetria_cap,
                              size_t* puntos_telemetria,
                              uint32_t* estado_out);

#ifdef __cplusplus
}
#endif

#endif /* K3_AUDIT_H */