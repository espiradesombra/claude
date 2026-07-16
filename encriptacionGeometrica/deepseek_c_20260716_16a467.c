#ifndef K3_H
#define K3_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Estructura de clave geométrica (33x1 por defecto)
typedef struct {
    int tales[5];
    int tales_count;
    char figuras[3][16];
    int figuras_count;
    int puntos[3];
    int puntos_count;
    int saltos[3];
    int saltos_count;
    int primos[2][2];
    int porcentajes_aportacion[2];
    int iteraciones_pi;
} ClaveK3;

// Estado del generador Aleatorovix
typedef struct {
    long double estado;
    long double factor_msl;
    long double factor_lsl;
} Aleatorovix;

// Matemáticas
long double k3_sin(long double x);
long double k3_pi(int p1, int p2, int iter);
long double k3_e(void);
void k3_lados(const char* tipo, int puntos, long double escala, long double out[3]);

// Aleatorovix
void k3_aleatorovix_init(Aleatorovix* ax, uint16_t semilla);
uint8_t k3_aleatorovix_bit(Aleatorovix* ax);

// Fase
long double k3_encriptar(uint16_t semilla, const ClaveK3* clave);
int k3_desencriptar(long double hash, const ClaveK3* clave);

// Máscara
void k3_mascara(const uint8_t* in, uint8_t* out, uint32_t len, uint16_t semilla);

#ifdef __cplusplus
}
#endif

#endif