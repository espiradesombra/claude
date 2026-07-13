#include "k3hash.h"
#include <stdio.h>
#include <string.h>

static const int OFFSETS[8] = {5, 7, 9, 11, 13, 17, 19, 23};
static const uint32_t MAGICO = 0x9E3779B1u;

static inline uint32_t rotar_izquierda_32bits(uint32_t valor, int posiciones) {
    posiciones &= 31;
    if (posiciones == 0) return valor;
    return (valor << posiciones) | (valor >> (32 - posiciones));
}

/* Un paso de mezcla del algoritmo K3 original, sin cambios de comportamiento. */
static void ejecutar_compresion_k3(uint32_t* estado, uint32_t valor_bloque, int N) {
    uint32_t B = valor_bloque;
    uint32_t x = estado[0] ^ (B & estado[1 % N]);

    for (int j = 2; j < N; j++) {
        x ^= rotar_izquierda_32bits(estado[j], OFFSETS[j % 8]);
    }
    x ^= rotar_izquierda_32bits(estado[1 % N], OFFSETS[4]);
    x += rotar_izquierda_32bits(estado[0], OFFSETS[5]);
    x ^= (B * MAGICO);

    x ^= B;
    x += rotar_izquierda_32bits(B, OFFSETS[6]);
    x ^= rotar_izquierda_32bits(B, OFFSETS[7]);
    x ^= rotar_izquierda_32bits(x, OFFSETS[2]);
    x *= MAGICO;
    x ^= (x >> 15);

    uint32_t estado_anterior_0 = estado[0];
    estado[0] = x;

    uint32_t temp_prev = estado_anterior_0;
    for (int i = 1; i < N; i++) {
        uint32_t temp_actual = estado[i];
        estado[i] = rotar_izquierda_32bits(estado[i], OFFSETS[0] + i) ^ temp_prev;
        temp_prev = temp_actual;
        if (i == N - 1) {
            estado[i] ^= B;
        }
    }
}

static uint32_t extraer_hash_final(const uint32_t* registros, int N) {
    uint32_t acumulador = registros[0];
    for (int i = 1; i < N; i++) {
        acumulador ^= rotar_izquierda_32bits(registros[i], 5 + i);
    }
    acumulador ^= acumulador >> 16;
    acumulador *= 0x85EBCA6Bu;
    acumulador ^= acumulador >> 13;
    acumulador *= 0xC2B2AE35u;
    acumulador ^= acumulador >> 16;
    return acumulador;
}

K3HashConfig k3_config_default(void) {
    K3HashConfig cfg;
    cfg.bits_ancho = 32;
    cfg.num_registros = 4;
    cfg.semilla_inicial = 0x1F2E3D4Cu;
    return cfg;
}

static int normalizar_config(K3HashConfig* cfg) {
    if (cfg->bits_ancho != 8 && cfg->bits_ancho != 16 && cfg->bits_ancho != 32) {
        cfg->bits_ancho = 32;
    }
    if (cfg->num_registros < 2) cfg->num_registros = 2;
    if (cfg->num_registros > K3HASH_MAX_REGISTROS) cfg->num_registros = K3HASH_MAX_REGISTROS;
    return cfg->bits_ancho / 8;
}

int k3_hash_init(K3HashCtx* ctx, const K3HashConfig* config) {
    if (!ctx) return -1;

    ctx->config = config ? *config : k3_config_default();
    ctx->ancho_bytes = normalizar_config(&ctx->config);
    ctx->bytes_en_buffer = 0;
    memset(ctx->buffer_parcial, 0, sizeof(ctx->buffer_parcial));

    for (int i = 0; i < ctx->config.num_registros; i++) {
        ctx->registros[i] = ctx->config.semilla_inicial ^ ((uint32_t)i * MAGICO);
    }
    return 0;
}

int k3_hash_update(K3HashCtx* ctx, const void* datos, size_t longitud) {
    if (!ctx || (!datos && longitud > 0)) return -1;

    const uint8_t* p = (const uint8_t*)datos;
    int N = ctx->config.num_registros;
    int ancho = ctx->ancho_bytes;

    size_t i = 0;
    while (i < longitud) {
        /* rellenar el buffer parcial hasta completar un bloque */
        while (ctx->bytes_en_buffer < ancho && i < longitud) {
            ctx->buffer_parcial[ctx->bytes_en_buffer++] = p[i++];
        }
        if (ctx->bytes_en_buffer == ancho) {
            uint32_t bloque = 0;
            for (int b = 0; b < ancho; b++) {
                bloque = (bloque << 8) | ctx->buffer_parcial[b];
            }
            ejecutar_compresion_k3(ctx->registros, bloque, N);
            ctx->bytes_en_buffer = 0;
        }
    }
    return 0;
}

uint32_t k3_hash_final(K3HashCtx* ctx) {
    if (!ctx) return 0;

    /* procesar el último bloque parcial, si lo hay, rellenando con ceros */
    if (ctx->bytes_en_buffer > 0) {
        uint32_t bloque = 0;
        for (int b = 0; b < ctx->ancho_bytes; b++) {
            uint8_t byte = (b < ctx->bytes_en_buffer) ? ctx->buffer_parcial[b] : 0;
            bloque = (bloque << 8) | byte;
        }
        ejecutar_compresion_k3(ctx->registros, bloque, ctx->config.num_registros);
        ctx->bytes_en_buffer = 0;
    }

    return extraer_hash_final(ctx->registros, ctx->config.num_registros);
}

uint32_t k3_hash_buffer(const void* datos, size_t longitud, const K3HashConfig* config) {
    K3HashCtx ctx;
    k3_hash_init(&ctx, config);
    k3_hash_update(&ctx, datos, longitud);
    return k3_hash_final(&ctx);
}

int k3_hash_file(const char* ruta, const K3HashConfig* config, uint32_t* hash_salida) {
    if (!ruta || !hash_salida) return -1;

    FILE* f = fopen(ruta, "rb");
    if (!f) return -1;

    K3HashCtx ctx;
    k3_hash_init(&ctx, config);

    uint8_t buffer[8192];
    size_t leidos;
    while ((leidos = fread(buffer, 1, sizeof(buffer), f)) > 0) {
        k3_hash_update(&ctx, buffer, leidos);
    }

    int error_lectura = ferror(f);
    fclose(f);
    if (error_lectura) return -1;

    *hash_salida = k3_hash_final(&ctx);
    return 0;
}
