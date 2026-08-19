#include "client_auth.h"

#include <openssl/hmac.h>
#include <openssl/evp.h>

int client_compute_response(const uint8_t token_usuario[AUTH_TOKEN_LEN],
                             const uint8_t challenge[AUTH_CHALLENGE_LEN],
                             uint8_t respuesta_salida[AUTH_RESPONSE_LEN]) {
    if (!token_usuario || !challenge || !respuesta_salida) return -1;

    unsigned int len_salida = 0;
    unsigned char* resultado = HMAC(EVP_sha256(),
                                     token_usuario, AUTH_TOKEN_LEN,
                                     challenge, AUTH_CHALLENGE_LEN,
                                     respuesta_salida, &len_salida);
    if (resultado == NULL || len_salida != AUTH_RESPONSE_LEN) {
        return -1;
    }
    return 0;
}
