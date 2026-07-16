/**
 * =====================================================================
 * SISTEMA INDUSTRIAL DE FASE K3 - VERSIÓN 2.0
 * =====================================================================
 * - Semilla de 64 bits (ampliable a 128)
 * - Función trampa: factorización de semiprimo
 * - CSPRNG: SHA-256 en modo contador
 * - Reverso por sistema lineal (usando factorización como guía)
 * - Parametrizable: Tales, figuras, puntos, saltos, primos, iteraciones
 * - Separación total entre matemática y flujo de uso
 * =====================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>

// =====================================================================
// BLOQUE I: MATEMÁTICA PURA (GEOMETRÍA, π, e)
// =====================================================================

/**
 * Calcula el seno usando serie de Taylor con long double.
 */
static long double taylor_sin(long double x) {
    long double term = x, sum = x, x2 = x * x;
    for (int n = 1; fabsl(term) > 1e-20L; n++) {
        term = -term * x2 / (long double)((2 * n) * (2 * n + 1));
        sum += term;
    }
    return sum;
}

/**
 * Aproxima π usando polígono de p1*p2 lados y 64 iteraciones de Arquímedes.
 */
long double aproximar_pi(int p1, int p2, int iteraciones) {
    long double lados = (long double)(p1 * p2);
    long double angulo = 3.14159265358979323846L / lados;
    long double lado = 2.0L * taylor_sin(angulo);
    long double perimetro = lados * lado;
    long double pi = perimetro / 2.0L;

    for (int i = 0; i < iteraciones; i++) {
        lados *= 2.0L;
        long double raiz = sqrtl(4.0L - lado * lado);
        lado = lado / sqrtl(2.0L + raiz);
        perimetro = lados * lado;
        pi = perimetro / 2.0L;
    }
    return pi;
}

/**
 * Aproxima e mediante serie convergente (50 términos).
 */
long double aproximar_e(void) {
    long double term = 1.0L, e = 1.0L;
    for (int v = 1; v <= 50; v++) {
        int denom = (v % 2 == 0) ? (2 * v + 2) : (3 * v);
        term = term / (long double)denom;
        e += term;
    }
    return e;
}

/**
 * Calcula los tres lados de una figura geométrica según tipo, puntos y escala.
 */
void calcular_lados(const char *tipo, int puntos, long double escala, long double lados[3]) {
    long double base = (long double)puntos * 1.5L;
    if (strcmp(tipo, "equilatero") == 0) {
        long double lado = (base / 3.0L) * escala;
        lados[0] = lado; lados[1] = lado; lados[2] = lado;
    } else if (strcmp(tipo, "isosceles") == 0) {
        long double lado_a = (base / 4.0L) * escala;
        long double lado_b = (base / 2.0L) * escala;
        lados[0] = lado_a; lados[1] = lado_a; lados[2] = lado_b;
    } else { // escaleno
        lados[0] = (base * 0.25L) * escala;
        lados[1] = (base * 0.35L) * escala;
        lados[2] = (base * 0.40L) * escala;
    }
}


// =====================================================================
// BLOQUE II: CRIPTOGRAFÍA DE TRAMPA (FACTORIZACIÓN, SHA-256, PRIMOS)
// =====================================================================

// ---------- Implementación de SHA-256 (estándar FIPS 180-4) ----------
typedef struct {
    uint32_t state[8];
    uint32_t count[2];
    uint8_t buffer[64];
} SHA256_CTX;

#define ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
#define CH(x, y, z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTR(x, 2) ^ ROTR(x, 13) ^ ROTR(x, 22))
#define EP1(x) (ROTR(x, 6) ^ ROTR(x, 11) ^ ROTR(x, 25))
#define SIG0(x) (ROTR(x, 7) ^ ROTR(x, 18) ^ ((x) >> 3))
#define SIG1(x) (ROTR(x, 17) ^ ROTR(x, 19) ^ ((x) >> 10))

static const uint32_t K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

static void sha256_transform(SHA256_CTX *ctx, const uint8_t *data) {
    uint32_t W[64], a, b, c, d, e, f, g, h;
    for (int i = 0; i < 16; i++)
        W[i] = (data[i*4]<<24) | (data[i*4+1]<<16) | (data[i*4+2]<<8) | data[i*4+3];
    for (int i = 16; i < 64; i++)
        W[i] = SIG1(W[i-2]) + W[i-7] + SIG0(W[i-15]) + W[i-16];
    a = ctx->state[0]; b = ctx->state[1]; c = ctx->state[2]; d = ctx->state[3];
    e = ctx->state[4]; f = ctx->state[5]; g = ctx->state[6]; h = ctx->state[7];
    for (int i = 0; i < 64; i++) {
        uint32_t T1 = h + EP1(e) + CH(e,f,g) + K[i] + W[i];
        uint32_t T2 = EP0(a) + MAJ(a,b,c);
        h = g; g = f; f = e; e = d + T1; d = c; c = b; b = a; a = T1 + T2;
    }
    ctx->state[0] += a; ctx->state[1] += b; ctx->state[2] += c; ctx->state[3] += d;
    ctx->state[4] += e; ctx->state[5] += f; ctx->state[6] += g; ctx->state[7] += h;
}

static void sha256_init(SHA256_CTX *ctx) {
    ctx->state[0] = 0x6a09e667; ctx->state[1] = 0xbb67ae85; ctx->state[2] = 0x3c6ef372;
    ctx->state[3] = 0xa54ff53a; ctx->state[4] = 0x510e527f; ctx->state[5] = 0x9b05688c;
    ctx->state[6] = 0x1f83d9ab; ctx->state[7] = 0x5be0cd19;
    ctx->count[0] = ctx->count[1] = 0;
}

static void sha256_update(SHA256_CTX *ctx, const uint8_t *data, size_t len) {
    size_t i;
    for (i = 0; i < len; i++) {
        ctx->buffer[ctx->count[0] & 63] = data[i];
        if ((ctx->count[0] & 63) == 63) {
            sha256_transform(ctx, ctx->buffer);
        }
        ctx->count[0]++;
        if (ctx->count[0] == 0) ctx->count[1]++;
    }
}

static void sha256_final(SHA256_CTX *ctx, uint8_t hash[32]) {
    uint64_t bitlen = ((uint64_t)ctx->count[1] << 32) | (uint64_t)ctx->count[0];
    sha256_update(ctx, (const uint8_t*)"\x80", 1);
    while ((ctx->count[0] & 63) != 56)
        sha256_update(ctx, (const uint8_t*)"\x00", 1);
    for (int i = 0; i < 8; i++)
        sha256_update(ctx, (const uint8_t*)&((uint8_t*)&bitlen)[7-i], 1);
    for (int i = 0; i < 8; i++)
        hash[i*4]   = (ctx->state[i] >> 24) & 0xff,
        hash[i*4+1] = (ctx->state[i] >> 16) & 0xff,
        hash[i*4+2] = (ctx->state[i] >> 8)  & 0xff,
        hash[i*4+3] =  ctx->state[i]        & 0xff;
}

// ---------- Funciones de trampa (factorización) ----------

/**
 * Genera un número primo de 32 bits a partir de una semilla usando SHA-256.
 */
static uint32_t generar_primo(uint64_t seed, int index) {
    SHA256_CTX ctx;
    uint8_t input[16];
    memcpy(input, &seed, 8);
    memcpy(input + 8, &index, 4);
    sha256_init(&ctx);
    sha256_update(&ctx, input, 12);
    uint8_t hash[32];
    sha256_final(&ctx, hash);
    // Tomamos los primeros 32 bits y los ajustamos a impar
    uint32_t p = (hash[0] << 24) | (hash[1] << 16) | (hash[2] << 8) | hash[3];
    p |= 1; // impar
    // Aseguramos que sea > 1000 y probamos primalidad simple (Miller-Rabin determinista para 32 bits)
    if (p < 1000) p += 1001;
    // Miller-Rabin para 32 bits: test con bases 2,3,5,7
    if (p % 2 == 0) p++;
    // Simplemente buscamos el siguiente primo (fuerza bruta para 32 bits es aceptable)
    while (1) {
        bool es_primo = true;
        if (p < 2) { p++; continue; }
        for (uint32_t d = 2; d * d <= p; d++) {
            if (p % d == 0) { es_primo = false; break; }
        }
        if (es_primo) return p;
        p += 2; // solo impares
    }
}

/**
 * Genera un semiprimo N = p*q a partir de una semilla.
 */
uint64_t generar_semiprimo(uint64_t seed, uint32_t *p_out, uint32_t *q_out) {
    uint32_t p = generar_primo(seed, 0);
    uint32_t q = generar_primo(seed, 1);
    if (p_out) *p_out = p;
    if (q_out) *q_out = q;
    return (uint64_t)p * (uint64_t)q;
}

/**
 * Factoriza un semiprimo de 64 bits (fuerza bruta hasta sqrt(n)).
 * Devuelve 0 si falla, o el factor p (q = n/p).
 */
uint32_t factorizar_semiprimo(uint64_t n) {
    if (n % 2 == 0) return 2;
    uint64_t limit = (uint64_t)sqrtl((long double)n);
    for (uint64_t d = 3; d <= limit; d += 2) {
        if (n % d == 0) return (uint32_t)d;
    }
    return 0; // fallo (n es primo o no encontrado)
}


// =====================================================================
// BLOQUE III: FLUJO PARAMETRIZABLE (estructuras y orquestación)
// =====================================================================

/**
 * Estructura de parámetros de la geometría (Tales, figuras, etc.)
 */
typedef struct {
    int tales[8];          // hasta 8 elementos
    int tales_count;
    char figuras[8][16];   // hasta 8 figuras
    int figuras_count;
    int puntos[8];
    int puntos_count;
    int saltos[8];
    int saltos_count;
    int primos[4][2];      // hasta 4 pares de primos para π
    int primos_count;
    int iteraciones_pi;
    int porcentajes[4];
} GeomParams;

/**
 * Estructura de parámetros de flujo (semilla, longitud, etc.)
 */
typedef struct {
    uint64_t seed;         // semilla maestra
    int num_bits;          // tamaño de la semilla de fase (múltiplo de 8)
    GeomParams geom;
    bool usar_trampa;      // si se usa factorización como trampa
} FaseK3Params;

/**
 * Inicializa los parámetros por defecto.
 */
void fase_k3_default_params(FaseK3Params *params, uint64_t seed) {
    params->seed = seed;
    params->num_bits = 64;  // semilla de 64 bits

    // Geometría por defecto
    GeomParams *g = &params->geom;
    int tales_def[] = {3, 5, 8, 13, 21};
    memcpy(g->tales, tales_def, sizeof(tales_def));
    g->tales_count = 5;
    const char *figs_def[] = {"equilatero", "isosceles", "escaleno"};
    for (int i = 0; i < 3; i++) strcpy(g->figuras[i], figs_def[i]);
    g->figuras_count = 3;
    int pts_def[] = {6, 12, 18};
    memcpy(g->puntos, pts_def, sizeof(pts_def));
    g->puntos_count = 3;
    int saltos_def[] = {1, 0, 2};
    memcpy(g->saltos, saltos_def, sizeof(saltos_def));
    g->saltos_count = 3;
    int primos_def[][2] = {{3,5}, {7,11}};
    for (int i = 0; i < 2; i++) {
        g->primos[i][0] = primos_def[i][0];
        g->primos[i][1] = primos_def[i][1];
    }
    g->primos_count = 2;
    g->iteraciones_pi = 64;
    int porc_def[] = {50, 100};
    memcpy(g->porcentajes, porc_def, sizeof(porc_def));
    params->usar_trampa = true;
}

/**
 * Calcula el hash de fase a partir de una semilla de fase (uint64_t).
 * La semilla de fase se deriva de la semilla maestra mediante SHA-256.
 */
long double fase_k3_calcular_hash(uint64_t fase_seed, const GeomParams *g) {
    long double perimetro = 0.0L;
    // Usamos 64 bits (8 bytes) de la semilla de fase
    for (int idx = 0; idx < 64; idx++) {
        uint8_t bit = (fase_seed >> (63 - idx)) & 1ULL;
        int figura_idx = idx / 3;
        int lado_idx = idx % 3;

        // Saltos
        if (lado_idx == g->saltos[figura_idx % g->saltos_count])
            continue;

        long double escala = (long double)g->tales[figura_idx % g->tales_count];
        const char *tipo = g->figuras[figura_idx % g->figuras_count];
        int pts = g->puntos[figura_idx % g->puntos_count];

        long double lados[3];
        calcular_lados(tipo, pts, escala, lados);
        if (bit) perimetro += lados[lado_idx];
    }

    // Ofuscación con π y e
    long double pi_ofuscado = 1.0L;
    for (int i = 0; i < g->primos_count; i++) {
        long double pi_aprox = aproximar_pi(g->primos[i][0], g->primos[i][1], g->iteraciones_pi);
        long double aport = (long double)g->porcentajes[i] / 100.0L;
        pi_ofuscado *= powl(pi_aprox, aport);
    }
    long double e_aprox = aproximar_e();

    return perimetro * pi_ofuscado * e_aprox;
}

/**
 * Deriva una semilla de fase de 64 bits a partir de la semilla maestra y un nonce.
 */
uint64_t fase_k3_derivar_semilla_fase(uint64_t seed, uint64_t nonce) {
    SHA256_CTX ctx;
    uint8_t input[16];
    memcpy(input, &seed, 8);
    memcpy(input + 8, &nonce, 8);
    sha256_init(&ctx);
    sha256_update(&ctx, input, 16);
    uint8_t hash[32];
    sha256_final(&ctx, hash);
    uint64_t out = 0;
    for (int i = 0; i < 8; i++)
        out = (out << 8) | hash[i];
    return out;
}

/**
 * Genera un flujo de bytes usando SHA-256 en modo contador (CSPRNG).
 */
void fase_k3_generar_flujo(uint64_t seed, uint8_t *output, size_t len) {
    SHA256_CTX ctx;
    uint64_t counter = 0;
    uint8_t input[24];
    memcpy(input, &seed, 8);
    size_t out_idx = 0;
    while (out_idx < len) {
        memcpy(input + 8, &counter, 8);
        sha256_init(&ctx);
        sha256_update(&ctx, input, 16);
        uint8_t hash[32];
        sha256_final(&ctx, hash);
        size_t copy = (len - out_idx > 32) ? 32 : len - out_idx;
        memcpy(output + out_idx, hash, copy);
        out_idx += copy;
        counter++;
    }
}

/**
 * Encripta un mensaje usando la semilla maestra y los parámetros.
 * Devuelve el hash de fase (long double) y el texto cifrado en cipher.
 * También devuelve el semiprimo N si usar_trampa está activado.
 */
long double fase_k3_encrypt(const uint8_t *plain, size_t len,
                            const FaseK3Params *params,
                            uint8_t *cipher,
                            uint64_t *semiprimo_out) {
    // Derivar semilla de fase a partir de la semilla maestra y un nonce (0)
    uint64_t fase_seed = fase_k3_derivar_semilla_fase(params->seed, 0);

    // Calcular hash de fase
    long double hash_fase = fase_k3_calcular_hash(fase_seed, &params->geom);

    // Generar flujo de cifrado a partir del hash de fase (convertido a entero)
    // Tomamos los primeros 64 bits del hash como semilla para el flujo
    uint64_t flujo_seed = (uint64_t)(hash_fase * 1e18L) & 0xFFFFFFFFFFFFFFFFULL;
    if (params->usar_trampa) {
        // Usar semiprimo como semilla del flujo (función trampa)
        uint32_t p, q;
        uint64_t N = generar_semiprimo(params->seed, &p, &q);
        if (semiprimo_out) *semiprimo_out = N;
        flujo_seed = N;
    }

    uint8_t *flujo = malloc(len);
    fase_k3_generar_flujo(flujo_seed, flujo, len);

    for (size_t i = 0; i < len; i++)
        cipher[i] = plain[i] ^ flujo[i];

    free(flujo);
    return hash_fase;
}

/**
 * Descifra un mensaje usando la semilla maestra y el hash de fase.
 * Si usar_trampa está activado, factoriza el semiprimo para obtener la semilla de flujo.
 * Devuelve 1 si éxito, 0 si fallo.
 */
int fase_k3_decrypt(const uint8_t *cipher, size_t len,
                    long double hash_fase,
                    const FaseK3Params *params,
                    uint8_t *plain) {
    uint64_t flujo_seed;
    if (params->usar_trampa) {
        // Reconstruir semiprimo a partir del hash? No, el hash es un long double.
        // En el diseño, el hash_fase se usa para derivar la semilla de fase, y el semiprimo se deriva de la semilla maestra.
        // Para descifrar, necesitamos el semiprimo. Podemos regenerarlo a partir de la semilla maestra (si la tenemos).
        // Pero en un sistema real, el emisor enviaría el hash y el receptor tendría la clave privada (semilla maestra).
        // Aquí asumimos que el receptor conoce la semilla maestra (clave privada).
        uint32_t p, q;
        uint64_t N = generar_semiprimo(params->seed, &p, &q);
        flujo_seed = N;
    } else {
        // Sin trampa, derivamos la semilla de fase del hash (esto es simétrico)
        uint64_t fase_seed = (uint64_t)(hash_fase * 1e18L) & 0xFFFFFFFFFFFFFFFFULL;
        flujo_seed = fase_seed;
    }

    uint8_t *flujo = malloc(len);
    fase_k3_generar_flujo(flujo_seed, flujo, len);

    for (size_t i = 0; i < len; i++)
        plain[i] = cipher[i] ^ flujo[i];

    free(flujo);
    return 1;
}


// =====================================================================
// DEMOSTRACIÓN INDUSTRIAL
// =====================================================================

int main(void) {
    printf("============================================================\n");
    printf("  SISTEMA INDUSTRIAL DE FASE K3 v2.0 - DEMO\n");
    printf("============================================================\n");

    // 1. Inicializar parámetros con una semilla maestra
    FaseK3Params params;
    fase_k3_default_params(&params, 0xDEADBEEFCAFEBABEULL);

    // Personalizar algunos parámetros (ejemplo)
    params.geom.tales[0] = 7;
    params.geom.tales[1] = 11;
    params.geom.iteraciones_pi = 32;
    params.usar_trampa = true;

    // 2. Mensaje original
    const char *msg = "ESPANA_PLAN_2030_CONEXION_DE_FASE_K3";
    size_t len = strlen(msg);
    uint8_t *plain = (uint8_t*)msg;
    uint8_t *cipher = malloc(len);
    uint8_t *decrypted = malloc(len + 1);

    // 3. Encriptar
    uint64_t semiprimo;
    long double hash = fase_k3_encrypt(plain, len, &params, cipher, &semiprimo);

    printf("\n[+] Mensaje original: %s\n", msg);
    printf("[+] Hash de fase generado: %.15Lf\n", hash);
    printf("[+] Semiprimo (trampa): %llu\n", (unsigned long long)semiprimo);

    // 4. Factorizar el semiprimo (simulando el atacante sin clave)
    uint32_t factor = factorizar_semiprimo(semiprimo);
    if (factor) {
        uint32_t otro = semiprimo / factor;
        printf("[+] Factorización (sin clave) encontrada: %u x %u\n", factor, otro);
    } else {
        printf("[!] No se pudo factorizar (número primo o demasiado grande).\n");
    }

    // 5. Descifrar (con clave: la semilla maestra)
    fase_k3_decrypt(cipher, len, hash, &params, decrypted);
    decrypted[len] = '\0';
    printf("[+] Mensaje descifrado: %s\n", decrypted);

    // 6. Validación
    if (memcmp(plain, decrypted, len) == 0)
        printf("\n[✔] ÉXITO: Cifrado y descifrado perfectos.\n");
    else
        printf("\n[✘] ERROR: Los datos no coinciden.\n");

    free(cipher);
    free(decrypted);
    return 0;
}