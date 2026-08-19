/*
 * demo_main.c
 * Simula:
 *   1) Registro: el servidor deriva el token de "victor" y se lo "entrega"
 *      (en la vida real esto viaja por TLS una sola vez).
 *   2) Login: servidor manda challenge -> cliente responde con su token ->
 *      servidor verifica sin haber almacenado el token en ningún sitio.
 *   3) Caso negativo: un atacante sin el token intenta responder y falla.
 */
#include <stdio.h>
#include <string.h>
#include "server_auth.h"
#include "client_auth.h"

static void imprimir_hex(const char* etiqueta, const uint8_t* datos, size_t len) {
    printf("%s: ", etiqueta);
    for (size_t i = 0; i < len; i++) printf("%02x", datos[i]);
    printf("\n");
}

int main(void) {
    printf("=== Demo autenticacion challenge-response (HMAC-SHA256) ===\n\n");

    /* --- Servidor arranca con clave maestra aleatoria --- */
    ServerAuthCtx servidor;
    if (server_auth_init(&servidor) != 0) {
        fprintf(stderr, "Error inicializando servidor\n");
        return 1;
    }

    /* --- Registro de usuario: servidor deriva su token y se lo entrega --- */
    const char* user_id = "victor";
    uint8_t token_usuario[AUTH_TOKEN_LEN];
    server_derive_user_token(&servidor, user_id, token_usuario);
    printf("[registro] usuario: %s\n", user_id);
    imprimir_hex("[registro] token entregado al usuario", token_usuario, AUTH_TOKEN_LEN);
    printf("\n");

    /* --- Login legitimo --- */
    uint8_t challenge[AUTH_CHALLENGE_LEN];
    server_generate_challenge(challenge);
    imprimir_hex("[servidor] challenge generado", challenge, AUTH_CHALLENGE_LEN);

    uint8_t respuesta[AUTH_RESPONSE_LEN];
    client_compute_response(token_usuario, challenge, respuesta);
    imprimir_hex("[cliente]  respuesta calculada", respuesta, AUTH_RESPONSE_LEN);

    int ok = server_verify_response(&servidor, user_id, challenge, respuesta);
    printf("[servidor] verificacion: %s\n\n", ok == 1 ? "VALIDA" : "RECHAZADA");

    /* --- Intento de un atacante sin el token correcto --- */
    uint8_t token_falso[AUTH_TOKEN_LEN];
    memset(token_falso, 0xAA, sizeof(token_falso)); /* token inventado */
    uint8_t respuesta_falsa[AUTH_RESPONSE_LEN];
    client_compute_response(token_falso, challenge, respuesta_falsa);

    int ok_falso = server_verify_response(&servidor, user_id, challenge, respuesta_falsa);
    printf("[atacante] intenta con token inventado\n");
    printf("[servidor] verificacion: %s\n", ok_falso == 1 ? "VALIDA (FALLO DE SEGURIDAD)" : "RECHAZADA (correcto)");

    server_auth_wipe(&servidor);
    return 0;
}
