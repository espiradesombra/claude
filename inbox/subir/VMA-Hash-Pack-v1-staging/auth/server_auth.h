/*
 * server_auth.h
 * ============================================================================
 * Librería de autenticación challenge-response (servidor).
 *
 * MODELO DE SEGURIDAD:
 *  - El servidor posee una CLAVE MAESTRA (32 bytes) que jamás sale de él.
 *  - Cada usuario recibe, en el registro, un TOKEN DERIVADO:
 *        token_usuario = HMAC-SHA256(clave_maestra, user_id)
 *    El usuario solo conoce SU token, nunca la clave maestra.
 *  - Para autenticar:
 *      1) el servidor genera un challenge aleatorio (nonce) y se lo envía
 *      2) el cliente calcula  respuesta = HMAC-SHA256(token_usuario, challenge)
 *      3) el servidor recalcula el token del usuario a partir de la clave
 *         maestra + user_id, recalcula la respuesta esperada y compara en
 *         tiempo constante.
 *
 * Esto es AUTENTICACIÓN (probar que alguien conoce un secreto) usando HMAC,
 * no es cifrado de confidencialidad. Si además necesitas cifrar el tráfico
 * (que nadie pueda leerlo, no solo verificarlo), eso es un paso aparte
 * (AES-256-GCM / ChaCha20-Poly1305) — decíamelo y lo añadimos.
 *
 * Requiere OpenSSL 3.x (libcrypto). Compilar con: -lcrypto
 * ============================================================================
 */

#ifndef SERVER_AUTH_H
#define SERVER_AUTH_H

#include <stddef.h>
#include <stdint.h>

#define AUTH_KEY_LEN     32   /* clave maestra: 256 bits                */
#define AUTH_TOKEN_LEN   32   /* token derivado por usuario: 256 bits   */
#define AUTH_CHALLENGE_LEN 32 /* nonce de challenge: 256 bits           */
#define AUTH_RESPONSE_LEN  32 /* respuesta HMAC: 256 bits               */

typedef struct {
    uint8_t clave_maestra[AUTH_KEY_LEN];
} ServerAuthCtx;

/* Inicializa el contexto del servidor generando una clave maestra
 * criptográficamente aleatoria (usa RAND_bytes de OpenSSL). */
int server_auth_init(ServerAuthCtx* ctx);

/* Inicializa el contexto del servidor con una clave maestra ya existente
 * (por ejemplo, cargada de un almacén seguro / HSM / variable de entorno). */
int server_auth_init_with_key(ServerAuthCtx* ctx, const uint8_t clave[AUTH_KEY_LEN]);

/* Deriva el token de un usuario a partir de su identificador.
 * El servidor puede llamar a esto en cualquier momento sin guardar
 * los tokens de cada usuario: solo necesita user_id + clave maestra. */
int server_derive_user_token(const ServerAuthCtx* ctx,
                              const char* user_id,
                              uint8_t token_salida[AUTH_TOKEN_LEN]);

/* Genera un challenge aleatorio para iniciar una autenticación. */
int server_generate_challenge(uint8_t challenge_salida[AUTH_CHALLENGE_LEN]);

/* Verifica la respuesta del cliente para un user_id y challenge dados.
 * Devuelve 1 si es válida, 0 si no lo es, -1 en error. */
int server_verify_response(const ServerAuthCtx* ctx,
                            const char* user_id,
                            const uint8_t challenge[AUTH_CHALLENGE_LEN],
                            const uint8_t respuesta_cliente[AUTH_RESPONSE_LEN]);

/* Borra la clave maestra de memoria de forma segura. Llamar siempre
 * antes de que el contexto salga de ámbito. */
void server_auth_wipe(ServerAuthCtx* ctx);

#endif /* SERVER_AUTH_H */
