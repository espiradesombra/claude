/*
 * ============================================================================
 * k3hash.h — Librería de hash no criptográfico K3 (portable)
 * ============================================================================
 * Basada en el algoritmo de mezcla original (rotaciones + XOR + constante
 * de dispersión 0x9E3779B1) aportado por el autor. Uso previsto: checksums
 * de integridad, deduplicación, tablas hash, fingerprinting, sharding.
 *
 * NO es criptográfico: no usar para contraseñas, autenticación, ni firma.
 * Para eso hay otra librería (HMAC-SHA256) ya entregada.
 *
 * Compatible con Windows (MSVC/MinGW), macOS (clang) y Linux (gcc/clang).
 * API en C puro (C99), sin dependencias externas.
 * ============================================================================
 */

#ifndef K3HASH_H
#define K3HASH_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------------
 * Macro de exportación multiplataforma para compilar como librería dinámica
 * (.dll en Windows, .dylib en macOS, .so en Linux) o estática.
 * Definir K3HASH_BUILD_SHARED al compilar la librería como .dll/.so/.dylib.
 * Definir K3HASH_SHARED (sin _BUILD) al consumirla como librería dinámica.
 * -------------------------------------------------------------------------- */
#if defined(_WIN32) || defined(_WIN64)
  #ifdef K3HASH_BUILD_SHARED
    #define K3HASH_API __declspec(dllexport)
  #elif defined(K3HASH_SHARED)
    #define K3HASH_API __declspec(dllimport)
  #else
    #define K3HASH_API
  #endif
#else
  #if defined(K3HASH_BUILD_SHARED) && __GNUC__ >= 4
    #define K3HASH_API __attribute__((visibility("default")))
  #else
    #define K3HASH_API
  #endif
#endif

#define K3HASH_MAX_REGISTROS 10

/* Configuración del hash. Usa k3_config_default() para valores por defecto
 * razonables, o personaliza si necesitas reproducir un resultado concreto. */
typedef struct {
    int      bits_ancho;      /* 8, 16 o 32 (ancho de bloque de lectura) */
    int      num_registros;   /* 2..K3HASH_MAX_REGISTROS, tamaño del anillo interno */
    uint32_t semilla_inicial; /* semilla base */
} K3HashConfig;

/* Contexto de hashing en streaming (para archivos grandes / datos por partes) */
typedef struct {
    K3HashConfig config;
    uint32_t     registros[K3HASH_MAX_REGISTROS];
    uint8_t      buffer_parcial[4];
    int          bytes_en_buffer;
    int          ancho_bytes;
} K3HashCtx;

/* Devuelve una configuración por defecto: 32 bits, 4 registros, semilla fija. */
K3HASH_API K3HashConfig k3_config_default(void);

/* --- API de un solo golpe (para bloques de memoria ya completos) --- */
K3HASH_API uint32_t k3_hash_buffer(const void* datos, size_t longitud,
                                    const K3HashConfig* config);

/* --- API en streaming (para archivos grandes o datos que llegan por partes) --- */
K3HASH_API int      k3_hash_init(K3HashCtx* ctx, const K3HashConfig* config);
K3HASH_API int      k3_hash_update(K3HashCtx* ctx, const void* datos, size_t longitud);
K3HASH_API uint32_t k3_hash_final(K3HashCtx* ctx);

/* Conveniencia: hashea un fichero completo por su ruta. Devuelve 0 en éxito
 * y escribe el hash en *hash_salida; devuelve -1 si no pudo abrir el fichero. */
K3HASH_API int      k3_hash_file(const char* ruta, const K3HashConfig* config,
                                  uint32_t* hash_salida);

#ifdef __cplusplus
}
#endif

#endif /* K3HASH_H */
