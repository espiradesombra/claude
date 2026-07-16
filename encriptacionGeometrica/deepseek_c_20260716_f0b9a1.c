#include "k3_engine.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

// =====================================================================
// BLOQUE I: MATEMÁTICA PURA
// =====================================================================

long double k3_taylor_sin(long double x) {
    long double termino = x;
    long double suma = x;
    long double x2 = x * x;
    int n = 1;
    while (fabsl(termino) > 1e-20L) {
        termino = -termino * x2 / (long double)((2 * n) * (2 * n + 1));
        suma += termino;
        n++;
    }
    return suma;
}

long double k3_aproximar_pi(int p1, int p2, int iteraciones) {
    long double lados = (long double)(p1 * p2);
    long double pi0 = 3.14159265358979323846264338327950288L;
    long double angulo = pi0 / lados;
    long double lado = 2.0L * k3_taylor_sin(angulo);
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

long double k3_aproximar_e(void) {
    long double termino = 1.0L;
    long double e = 1.0L;
    for (int v = 1; v <= 50; v++) {
        int denom = (v % 2 == 0) ? (2 * v + 2) : (3 * v);
        if (denom == 0) denom = 2 * v + 2;
        termino /= (long double)denom;
        e += termino;
    }
    return e;
}

void k3_calcular_lados(const char* tipo, int puntos, long double escala, long double lados_out[3]) {
    long double base = (long double)puntos * 1.5L;
    if (strcmp(tipo, "equilatero") == 0) {
        long double lado = (base / 3.0L) * escala;
        lados_out[0] = lado; lados_out[1] = lado; lados_out[2] = lado;
    } else if (strcmp(tipo, "isosceles") == 0) {
        long double a = (base / 4.0L) * escala;
        long double b = (base / 2.0L) * escala;
        lados_out[0] = a; lados_out[1] = a; lados_out[2] = b;
    } else { // escaleno
        lados_out[0] = (base * 0.25L) * escala;
        lados_out[1] = (base * 0.35L) * escala;
        lados_out[2] = (base * 0.40L) * escala;
    }
}

// =====================================================================
// BLOQUE II: ALEATOROVIX (Generador de flujo caótico)
// =====================================================================

void k3_aleatorovix_init(Aleatorovix* ax, uint16_t semilla) {
    ax->estado = (semilla != 0) ? (long double)semilla : 12345.0L;
    ax->factor_msl = 42.123456789L;
    ax->factor_lsl = 0.000123456L;
    ax->estado -= floorl(ax->estado);
    if (ax->estado == 0.0L) ax->estado = 0.5L;
}

uint8_t k3_aleatorovix_siguiente_bit(Aleatorovix* ax) {
    long double fase = (ax->estado * ax->factor_msl) + cosl(ax->estado * ax->factor_lsl);
    ax->estado = fase - floorl(fase);
    return (ax->estado > 0.5L) ? 1 : 0;
}

// =====================================================================
// BLOQUE III: CRIPTOGRAFÍA DE FASE (Backtracking Geométrico)
// =====================================================================

static bool backtrack(long double residuo, int idx, const ClaveK3* clave, uint16_t* bits_out) {
    if (idx == 16) {
        return (fabsl(residuo) < 1e-11L);
    }

    int fig = idx / 3;
    int lado = idx % 3;

    if (lado == clave->saltos[fig % clave->saltos_count]) {
        for (int bit = 0; bit <= 1; bit++) {
            if (backtrack(residuo, idx + 1, clave, bits_out)) {
                *bits_out |= (bit << (15 - idx));
                return true;
            }
        }
        return false;
    }

    long double escala = (long double)clave->tales[fig % clave->tales_count];
    const char* tipo = clave->figuras[fig % clave->figuras_count];
    int pts = clave->puntos[fig % clave->puntos_count];

    long double lados[3];
    k3_calcular_lados(tipo, pts, escala, lados);
    long double lado_actual = lados[lado];

    if (residuo >= lado_actual - 1e-11L) {
        if (backtrack(residuo - lado_actual, idx + 1, clave, bits_out)) {
            *bits_out |= (1 << (15 - idx));
            return true;
        }
    }

    if (backtrack(residuo, idx + 1, clave, bits_out)) {
        *bits_out &= ~(1 << (15 - idx));
        return true;
    }

    return false;
}

long double k3_encriptar_semilla(uint16_t semilla, const ClaveK3* clave) {
    long double perimetro = 0.0L;

    for (int idx = 0; idx < 16; idx++) {
        uint8_t bit = (semilla >> (15 - idx)) & 1;
        int fig = idx / 3;
        int lado = idx % 3;

        if (lado == clave->saltos[fig % clave->saltos_count]) continue;

        long double escala = (long double)clave->tales[fig % clave->tales_count];
        const char* tipo = clave->figuras[fig % clave->figuras_count];
        int pts = clave->puntos[fig % clave->puntos_count];

        long double lados[3];
        k3_calcular_lados(tipo, pts, escala, lados);

        if (bit) perimetro += lados[lado];
    }

    long double pi_ofuscado = 1.0L;
    for (int i = 0; i < 2; i++) {
        long double pi_aprox = k3_aproximar_pi(clave->primos[i][0], clave->primos[i][1], clave->iteraciones_pi);
        long double aporte = (long double)clave->porcentajes_aportacion[i] / 100.0L;
        pi_ofuscado *= powl(pi_aprox, aporte);
    }

    long double e = k3_aproximar_e();
    return perimetro * pi_ofuscado * e;
}

int k3_desencriptar_semilla(long double hash_fase, const ClaveK3* clave) {
    long double pi_ofuscado = 1.0L;
    for (int i = 0; i < 2; i++) {
        long double pi_aprox = k3_aproximar_pi(clave->primos[i][0], clave->primos[i][1], clave->iteraciones_pi);
        long double aporte = (long double)clave->porcentajes_aportacion[i] / 100.0L;
        pi_ofuscado *= powl(pi_aprox, aporte);
    }
    long double e = k3_aproximar_e();
    long double perimetro = hash_fase / (pi_ofuscado * e);

    uint16_t semilla_recuperada = 0;
    if (backtrack(perimetro, 0, clave, &semilla_recuperada)) {
        return (int)semilla_recuperada;
    }
    return -1;
}

// =====================================================================
// BLOQUE IV: MÁSCARA LILA (Cifrado de flujo masivo)
// =====================================================================

void k3_aplicar_mascara(const uint8_t* origen, uint8_t* destino, uint32_t longitud, uint16_t semilla) {
    Aleatorovix ax;
    k3_aleatorovix_init(&ax, semilla);

    for (uint32_t i = 0; i < longitud; i++) {
        uint8_t byte = origen[i];
        uint8_t resultado = 0;
        for (int bit = 0; bit < 8; bit++) {
            uint8_t b_dato = (byte >> (7 - bit)) & 1;
            uint8_t b_chorro = k3_aleatorovix_siguiente_bit(&ax);
            resultado |= ((b_dato ^ b_chorro) << (7 - bit));
        }
        destino[i] = resultado;
    }
}