# Librería de autenticación challenge-response (HMAC-SHA256)

## Qué es
Autenticación usuario-servidor mediante HMAC-SHA256 (OpenSSL). El servidor
guarda una única clave maestra; cada usuario tiene un token derivado que
solo él conoce. El servidor nunca almacena tokens: los recalcula al vuelo
a partir de la clave maestra + user_id.

**No es cifrado de confidencialidad.** Es prueba de posesión de un secreto
(autenticación/integridad). Si necesitas que el tráfico sea ilegible para
terceros, hay que añadir AES-256-GCM / ChaCha20-Poly1305 por encima
(típicamente ya lo hace TLS si viaja por HTTPS).

## Archivos
- `server_auth.h/.c` — lado servidor: clave maestra, derivación de tokens,
  generación de challenge, verificación en tiempo constante.
- `client_auth.h/.c` — lado cliente: cálculo de la respuesta al challenge.
- `demo_main.c` — simula registro + login legítimo + intento de ataque.
- `libauth.a` — librería estática ya compilada.
- `demo_auth` — ejecutable de la demo, ya compilado (Linux x86_64).

## Compilar desde cero
```
gcc -Wall -Wextra -c server_auth.c -o server_auth.o
gcc -Wall -Wextra -c client_auth.c -o client_auth.o
ar rcs libauth.a server_auth.o client_auth.o

gcc tu_programa.c server_auth.o client_auth.o -lcrypto -o tu_programa
```
Requiere `libssl-dev` (OpenSSL 3.x).

## Flujo real de despliegue
1. El servidor genera su clave maestra una vez y la guarda en un sitio
   seguro (variable de entorno, KMS, HSM) — NUNCA en el código fuente.
2. En el registro de cada usuario, el servidor deriva su token y se lo
   entrega **por un canal ya cifrado** (HTTPS/TLS). El usuario lo guarda
   localmente (ej. en su app, cifrado con su PIN/biometría).
3. En cada login: servidor manda challenge (nonce nuevo cada vez, para
   evitar replay) → cliente responde con HMAC(token, challenge) →
   servidor verifica.
4. Este protocolo va SIEMPRE dentro de TLS (HTTPS), no lo sustituye:
   protege contra "no sé tu token" pero no oculta el tráfico de un
   atacante que ya esté en la red si no hay TLS por encima.
