#ifndef HASHK3_256_H
#define HASHK3_256_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Estado de 256 bits = 4 x uint64_t
typedef struct {
    uint64_t v[4];       // Registros principales (el "acordeón")
    uint64_t sig;        // Firma acumulada (para anidación de seguridad)
    uint64_t len;        // Longitud procesada
    uint8_t buffer[64];  // Buffer de absorción
    int buf_pos;         // Posición en el buffer
} K3_CTX_256;

// Inicializa el contexto
void k3_256_init(K3_CTX_256 *ctx);

// Actualiza con datos (absorbe)
void k3_256_update(K3_CTX_256 *ctx, const void *data, size_t len);

// Finaliza y escribe el hash de 256 bits (32 bytes)
void k3_256_final(K3_CTX_256 *ctx, uint8_t *out);

// Función de un solo paso (para checksums rápidos)
void k3_256_hash(const void *data, size_t len, uint8_t *out);

#ifdef __cplusplus
}
#endif

#endif // HASHK3_256_H