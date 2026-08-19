#include "server_auth.h"

#include <string.h>
#include <openssl/rand.h>
#include <openssl/hmac.h>
#include <openssl/evp.h>

int server_auth_init(ServerAuthCtx* ctx) {
    if (!ctx) return -1;
    if (RAND_bytes(ctx->clave_maestra, AUTH_KEY_LEN) != 1) {
        return -1;
    }
    return 0;
}

int server_auth_init_with_key(ServerAuthCtx* ctx, const uint8_t clave[AUTH_KEY_LEN]) {
    if (!ctx || !clave) return -1;
    memcpy(ctx->clave_maestra, clave, AUTH_KEY_LEN);
    return 0;
}

int server_derive_user_token(const ServerAuthCtx* ctx,
                              const char* user_id,
                              uint8_t token_salida[AUTH_TOKEN_LEN]) {
    if (!ctx || !user_id || !token_salida) return -1;

    unsigned int len_salida = 0;
    unsigned char* resultado = HMAC(EVP_sha256(),
                                     ctx->clave_maestra, AUTH_KEY_LEN,
                                     (const unsigned char*)user_id, strlen(user_id),
                                     token_salida, &len_salida);
    if (resultado == NULL || len_salida != AUTH_TOKEN_LEN) {
        return -1;
    }
    return 0;
}

int server_generate_challenge(uint8_t challenge_salida[AUTH_CHALLENGE_LEN]) {
    if (!challenge_salida) return -1;
    if (RAND_bytes(challenge_salida, AUTH_CHALLENGE_LEN) != 1) {
        return -1;
    }
    return 0;
}

/* Comparación en tiempo constante para no filtrar información por timing. */
static int comparar_tiempo_constante(const uint8_t* a, const uint8_t* b, size_t len) {
    uint8_t diff = 0;
    for (size_t i = 0; i < len; i++) {
        diff |= a[i] ^ b[i];
    }
    return diff == 0;
}

int server_verify_response(const ServerAuthCtx* ctx,
                            const char* user_id,
                            const uint8_t challenge[AUTH_CHALLENGE_LEN],
                            const uint8_t respuesta_cliente[AUTH_RESPONSE_LEN]) {
    if (!ctx || !user_id || !challenge || !respuesta_cliente) return -1;

    uint8_t token_usuario[AUTH_TOKEN_LEN];
    if (server_derive_user_token(ctx, user_id, token_usuario) != 0) {
        return -1;
    }

    unsigned int len_salida = 0;
    uint8_t respuesta_esperada[AUTH_RESPONSE_LEN];
    unsigned char* resultado = HMAC(EVP_sha256(),
                                     token_usuario, AUTH_TOKEN_LEN,
                                     challenge, AUTH_CHALLENGE_LEN,
                                     respuesta_esperada, &len_salida);

    /* Borra el token de la pila antes de salir por cualquier camino. */
    if (resultado == NULL || len_salida != AUTH_RESPONSE_LEN) {
        OPENSSL_cleanse(token_usuario, sizeof(token_usuario));
        return -1;
    }

    int valido = comparar_tiempo_constante(respuesta_esperada, respuesta_cliente, AUTH_RESPONSE_LEN);

    OPENSSL_cleanse(token_usuario, sizeof(token_usuario));
    OPENSSL_cleanse(respuesta_esperada, sizeof(respuesta_esperada));

    return valido ? 1 : 0;
}

void server_auth_wipe(ServerAuthCtx* ctx) {
    if (!ctx) return;
    OPENSSL_cleanse(ctx->clave_maestra, sizeof(ctx->clave_maestra));
}
