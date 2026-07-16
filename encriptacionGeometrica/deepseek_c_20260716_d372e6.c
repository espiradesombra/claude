#ifndef K3_ENGINE_H
#define K3_ENGINE_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Estructura de parámetros geométricos (Clave Maestra)
typedef struct {
    int tales[5];               // Secuencia de Tales
    int tales_count;
    char figuras[3][16];        // "equilatero", "isosceles", "escaleno"
    int figuras_count;
    int puntos[3];              // Puntos base
    int puntos_count;
    int saltos[3];              // Lados omitidos
    int saltos_count;
    int primos[2][2];           // Parejas de primos para polígono de Pi
    int porcentajes_aportacion[2];
    int iteraciones_pi;
} ClaveK3;

// Estructura de estado del generador caótico Aleatorovix
typedef struct {
    long double estado;
    long double factor_msl;
    long double factor_lsl;
} Aleatorovix;

// --- Funciones matemáticas puras ---
long double k3_taylor_sin(long double x);
long double k3_aproximar_pi(int p1, int p2, int iteraciones);
long double k3_aproximar_e(void);
void k3_calcular_lados(const char* tipo, int puntos, long double escala, long double lados_out[3]);

// --- Generador Aleatorovix ---
void k3_aleatorovix_init(Aleatorovix* ax, uint16_t semilla);
uint8_t k3_aleatorovix_siguiente_bit(Aleatorovix* ax);

// --- Cifrado de fase (colapso y backtracking) ---
long double k3_encriptar_semilla(uint16_t semilla, const ClaveK3* clave);
int k3_desencriptar_semilla(long double hash_fase, const ClaveK3* clave);

// --- Cifrado de flujo masivo (Máscara Lila) ---
void k3_aplicar_mascara(const uint8_t* origen, uint8_t* destino, uint32_t longitud, uint16_t semilla);

#ifdef __cplusplus
}
#endif

#endif // K3_ENGINE_H