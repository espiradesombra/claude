/*
 * client_auth.h
 * ============================================================================
 * Librería de autenticación challenge-response (cliente / usuario).
 *
 * El cliente SOLO conoce su token derivado (recibido una vez en el registro,
 * por un canal seguro — TLS, entrega manual, etc.). Nunca conoce la clave
 * maestra del servidor. Con el token y un challenge que le manda el
 * servidor, calcula una respuesta HMAC-SHA256 que demuestra que posee el
 * token sin revelarlo.
 * ============================================================================
 */

#ifndef CLIENT_AUTH_H
#define CLIENT_AUTH_H

#include <stdint.h>

#define AUTH_TOKEN_LEN      32
#define AUTH_CHALLENGE_LEN  32
#define AUTH_RESPONSE_LEN   32

/* Calcula la respuesta al challenge usando el token del usuario.
 * Devuelve 0 en éxito, -1 en error. */
int client_compute_response(const uint8_t token_usuario[AUTH_TOKEN_LEN],
                             const uint8_t challenge[AUTH_CHALLENGE_LEN],
                             uint8_t respuesta_salida[AUTH_RESPONSE_LEN]);

#endif /* CLIENT_AUTH_H */
