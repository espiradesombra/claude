#include "hashk3_256.h"
#include <string.h>

// Constantes mágicas (números primos/áureos para dispersión)
#define ROT1 0x9E3779B97F4A7C15ULL  // 1.618...
#define ROT2 0xBF58476D1CE4E5B9ULL  // Otro primo
#define ROT3 0x94D049BB133111EBULL  // Tercera constante

// Secuencia de anidación: rotaciones para cada uno de los 4 registros
// (va cambiando según el "lado" de la iteración)
static inline uint64_t rotar(uint64_t x, int r) {
    return (x >> r) | (x << (64 - r));
}

// Mezcla de firmas (Guardado de vez en cuando)
static void k3_256_mezclar(K3_CTX_256 *ctx) {
    // Anidación al lado según secuencia fija (como un Feistel)
    ctx->v[0] += ctx->v[1] ^ ctx->sig;
    ctx->v[1] += ctx->v[2] ^ ROT1;
    ctx->v[2] += ctx->v[3] ^ ROT2;
    ctx->v[3] += ctx->v[0] ^ ROT3;

    // Rotaciones diferenciales (cada registro se mueve a su lado)
    ctx->v[0] = rotar(ctx->v[0], 7);
    ctx->v[1] = rotar(ctx->v[1], 13);
    ctx->v[2] = rotar(ctx->v[2], 23);
    ctx->v[3] = rotar(ctx->v[3], 37);

    // Actualizar la firma global (introducir más firmas en el proceso)
    ctx->sig ^= (ctx->v[0] & ctx->v[1]) | (ctx->v[2] ^ ctx->v[3]);
    ctx->sig = rotar(ctx->sig, 17) + ROT1;
}

// Inicialización
void k3_256_init(K3_CTX_256 *ctx) {
    memset(ctx->v, 0, sizeof(ctx->v));
    ctx->sig = 0x9E3779B97F4A7C15ULL; // Semilla inicial
    ctx->len = 0;
    ctx->buf_pos = 0;
    memset(ctx->buffer, 0, 64);
}

// Procesar un bloque completo de 64 bytes (con guardado intermedio)
static void k3_256_procesar_bloque(K3_CTX_256 *ctx, const uint8_t *bloque) {
    // Absorber bloque en el estado de 256 bits
    // (Agrupamiento de bits: 4 palabras de 64 bits)
    for (int i = 0; i < 4; i++) {
        uint64_t palabra;
        memcpy(&palabra, bloque + (i * 8), 8);
        ctx->v[i] ^= palabra;
    }

    // Anidar al lado según secuencia (guardado de vez en cuando)
    k3_256_mezclar(ctx);

    // Introducir más firmas: rotar el buffer de firmas
    ctx->sig ^= (ctx->v[0] + ctx->v[1]) ^ (ctx->v[2] - ctx->v[3]);
    ctx->sig = rotar(ctx->sig, 31);
}

// Actualización (modo streaming)
void k3_256_update(K3_CTX_256 *ctx, const void *data, size_t len) {
    const uint8_t *d = (const uint8_t *)data;
    ctx->len += len;

    for (size_t i = 0; i < len; i++) {
        ctx->buffer[ctx->buf_pos++] = d[i];
        if (ctx->buf_pos == 64) {
            k3_256_procesar_bloque(ctx, ctx->buffer);
            ctx->buf_pos = 0;
        }
    }
}

// Finalización (padding y extracción de 256 bits)
void k3_256_final(K3_CTX_256 *ctx, uint8_t *out) {
    // Padding: 0x80 + ceros + longitud (como en SHA-1 pero adaptado)
    ctx->buffer[ctx->buf_pos++] = 0x80;
    if (ctx->buf_pos > 56) {
        while (ctx->buf_pos < 64) ctx->buffer[ctx->buf_pos++] = 0x00;
        k3_256_procesar_bloque(ctx, ctx->buffer);
        ctx->buf_pos = 0;
        memset(ctx->buffer, 0, 64);
    }
    while (ctx->buf_pos < 56) ctx->buffer[ctx->buf_pos++] = 0x00;

    // Escribir la longitud (64 bits) en los últimos 8 bytes
    uint64_t len_bits = ctx->len * 8;
    memcpy(ctx->buffer + 56, &len_bits, 8);

    k3_256_procesar_bloque(ctx, ctx->buffer);

    // Una última mezcla de seguridad (más firmas)
    for (int i = 0; i < 4; i++) k3_256_mezclar(ctx);

    // Extraer los 32 bytes
    for (int i = 0; i < 4; i++) {
        memcpy(out + (i * 8), &ctx->v[i], 8);
    }
}

// Hash de un solo paso
void k3_256_hash(const void *data, size_t len, uint8_t *out) {
    K3_CTX_256 ctx;
    k3_256_init(&ctx);
    k3_256_update(&ctx, data, len);
    k3_256_final(&ctx, out);
}