/* π/e ofuscado + Aleatorovix — encriptacionGeometrica/deepseek_c + gemini-code */
#include <math.h>
#include <string.h>
#include "../include/antipc_native.h"

#define GEO_TOL 1e-9

typedef struct {
    double estado;
    double factor_msl;
    double factor_lsl;
} AleatorovixState;

static double taylor_sin(double x) {
    double term = x, sum = x, x2 = x * x;
    int n = 1;
    while (fabs(term) > 1e-14) {
        term = -term * x2 / (double)((2 * n) * (2 * n + 1));
        sum += term;
        ++n;
    }
    return sum;
}

static double aproximar_pi(int p1, int p2, int iteraciones) {
    double lados = (double)(p1 * p2);
    double pi0 = 3.14159265358979323846;
    double angulo = pi0 / lados;
    double lado = 2.0 * taylor_sin(angulo);
    double perimetro = lados * lado;
    double pi = perimetro / 2.0;
    int i;
    for (i = 0; i < iteraciones; ++i) {
        double raiz;
        lados *= 2.0;
        raiz = sqrt(4.0 - lado * lado);
        lado = lado / sqrt(2.0 + raiz);
        perimetro = lados * lado;
        pi = perimetro / 2.0;
    }
    return pi;
}

static double aproximar_e(void) {
    double term = 1.0, e = 1.0;
    int v;
    for (v = 1; v <= 50; ++v) {
        int denom = (v % 2 == 0) ? (2 * v + 2) : (3 * v);
        if (denom == 0) denom = 2 * v + 2;
        term /= (double)denom;
        e += term;
    }
    return e;
}

static void calcular_lados(int tipo, int puntos, double escala, double lados[3]) {
    double base = (double)puntos * 1.5;
    if (tipo == 0) {
        double lado = (base / 3.0) * escala;
        lados[0] = lados[1] = lados[2] = lado;
    } else if (tipo == 1) {
        double a = (base / 4.0) * escala;
        double b = (base / 2.0) * escala;
        lados[0] = a;
        lados[1] = a;
        lados[2] = b;
    } else {
        lados[0] = (base * 0.25) * escala;
        lados[1] = (base * 0.35) * escala;
        lados[2] = (base * 0.40) * escala;
    }
}

static double pi_ofuscado(const AntipcGeoClave* clave) {
    double pi_o = 1.0;
    uint32_t i;
    for (i = 0; i < clave->primos_count && i < ANTIPC_GEO_MAX_PRIMOS; ++i) {
        double pi_a = aproximar_pi(
            clave->primos_p1[i], clave->primos_p2[i], clave->iteraciones_pi
        );
        double aporte = (double)clave->porcentajes[i] / 100.0;
        pi_o *= pow(pi_a, aporte);
    }
    return pi_o;
}

static void aleatorovix_init(AleatorovixState* ax, uint16_t semilla) {
    ax->estado = (semilla != 0) ? (double)semilla : 12345.0;
    ax->factor_msl = 42.123456789;
    ax->factor_lsl = 0.000123456;
    ax->estado -= floor(ax->estado);
    if (ax->estado == 0.0) ax->estado = 0.5;
}

static uint8_t aleatorovix_bit(AleatorovixState* ax) {
    double fase = (ax->estado * ax->factor_msl) + cos(ax->estado * ax->factor_lsl);
    ax->estado = fase - floor(fase);
    return (ax->estado > 0.5) ? 1 : 0;
}

static int backtrack(double residuo, int idx, const AntipcGeoClave* clave, uint16_t* bits) {
    int fig, lado;
    double lados[3], lado_actual;

    if (idx == 16) return (fabs(residuo) < GEO_TOL);

    fig = idx / 3;
    lado = idx % 3;

    if ((uint32_t)lado < clave->saltos_count &&
        lado == clave->saltos[fig % clave->saltos_count]) {
        int bit;
        for (bit = 0; bit <= 1; ++bit) {
            if (backtrack(residuo, idx + 1, clave, bits)) {
                if (bit) *bits |= (uint16_t)(1 << (15 - idx));
                return 1;
            }
        }
        return 0;
    }

    calcular_lados(
        clave->figuras[fig % clave->figuras_count],
        (int)clave->puntos[fig % clave->puntos_count],
        (double)clave->tales[fig % clave->tales_count],
        lados
    );
    lado_actual = lados[lado];

    if (residuo >= lado_actual - GEO_TOL) {
        if (backtrack(residuo - lado_actual, idx + 1, clave, bits)) {
            *bits |= (uint16_t)(1 << (15 - idx));
            return 1;
        }
    }
    if (backtrack(residuo, idx + 1, clave, bits)) {
        *bits &= (uint16_t)~(1 << (15 - idx));
        return 1;
    }
    return 0;
}

ANTIPC_API void antipc_geo_clave_default(AntipcGeoClave* clave) {
    if (!clave) return;
    memset(clave, 0, sizeof(*clave));
    clave->tales_count = 5;
    clave->tales[0] = 3;
    clave->tales[1] = 5;
    clave->tales[2] = 8;
    clave->tales[3] = 13;
    clave->tales[4] = 21;
    clave->figuras_count = 3;
    clave->figuras[0] = 0;
    clave->figuras[1] = 1;
    clave->figuras[2] = 2;
    clave->puntos_count = 3;
    clave->puntos[0] = 6;
    clave->puntos[1] = 12;
    clave->puntos[2] = 18;
    clave->saltos_count = 3;
    clave->saltos[0] = 1;
    clave->saltos[1] = 0;
    clave->saltos[2] = 2;
    clave->primos_count = 2;
    clave->primos_p1[0] = 3;
    clave->primos_p2[0] = 5;
    clave->primos_p1[1] = 7;
    clave->primos_p2[1] = 11;
    clave->porcentajes[0] = 50;
    clave->porcentajes[1] = 100;
    clave->iteraciones_pi = 15;
}

ANTIPC_API double antipc_geo_fase_encrypt(uint16_t semilla, const AntipcGeoClave* clave) {
    double perimetro = 0.0;
    int idx;

    if (!clave) return 0.0;

    for (idx = 0; idx < 16; ++idx) {
        int bit = (semilla >> (15 - idx)) & 1;
        int fig = idx / 3;
        int lado = idx % 3;
        double lados[3];

        if ((uint32_t)lado < clave->saltos_count &&
            lado == clave->saltos[fig % clave->saltos_count]) {
            continue;
        }

        calcular_lados(
            clave->figuras[fig % clave->figuras_count],
            (int)clave->puntos[fig % clave->puntos_count],
            (double)clave->tales[fig % clave->tales_count],
            lados
        );
        if (bit) perimetro += lados[lado];
    }

    return perimetro * pi_ofuscado(clave) * aproximar_e();
}

ANTIPC_API int antipc_geo_fase_decrypt(double hash_fase, const AntipcGeoClave* clave, uint16_t* semilla_out) {
    double perimetro;
    uint16_t bits = 0;

    if (!clave || !semilla_out) return 0;
    perimetro = hash_fase / (pi_ofuscado(clave) * aproximar_e());
    if (!backtrack(perimetro, 0, clave, &bits)) return 0;
    *semilla_out = bits;
    return 1;
}

ANTIPC_API void antipc_aleatorovix_xor(
    const uint8_t* in, uint8_t* out, size_t len, uint16_t semilla
) {
    AleatorovixState ax;
    size_t i;
    if (!in || !out || len == 0) return;

    aleatorovix_init(&ax, semilla);
    for (i = 0; i < len; ++i) {
        uint8_t byte = in[i];
        uint8_t result = 0;
        int bit;
        for (bit = 0; bit < 8; ++bit) {
            uint8_t b_dato = (byte >> (7 - bit)) & 1;
            uint8_t b_chorro = aleatorovix_bit(&ax);
            result |= (uint8_t)((b_dato ^ b_chorro) << (7 - bit));
        }
        out[i] = result;
    }
}

ANTIPC_API int antipc_geo_masivo_crypt(
    const uint8_t* in,
    uint8_t* out,
    size_t len,
    uint16_t semilla,
    const AntipcGeoClave* clave,
    double* hash_fase_out,
    int decrypt
) {
    uint16_t seed = semilla;

    if (!in || !out || !clave || len == 0) return 0;

    if (decrypt) {
        if (!hash_fase_out || !antipc_geo_fase_decrypt(*hash_fase_out, clave, &seed)) {
            return 0;
        }
    } else if (hash_fase_out) {
        *hash_fase_out = antipc_geo_fase_encrypt(semilla, clave);
    }

    antipc_aleatorovix_xor(in, out, len, seed);
    return 1;
}