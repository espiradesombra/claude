#ifndef FASE_K3_ALEATOROVIX_H
#define FASE_K3_ALEATOROVIX_H

#include <stdint.h>
#include <stdbool.h>

// Estructura de Clave Geométrica K3 (Parámetros del Desafío)
typedef struct {
    int tales[5];               // Secuencia de Tales (escalas)
    int tales_count;
    char figuras[3][16];        // "equilatero", "isosceles", "escaleno"
    int figuras_count;
    int puntos[3];              // Puntos base por nivel
    int puntos_count;
    int saltos[3];              // Lados omitidos deliberadamente
    int saltos_count;
    int primos[2][2];           // Parejas de primos para el polígono de Pi
    int porcentajes_aportacion[2];
    int iteraciones_pi;
} ClaveK3;

// Estructura de Estado de Aleatorovix (Motor Caótico MSL/LSL)
typedef struct {
    long double estado;
    long double factor_msl;
    long double factor_lsl;
} Aleatorovix;

// --- BLOQUE I: MATEMÁTICA PURA ---
long double taylor_sin(long double x);
long double aproximar_pi_puro(int p1, int p2, int iteraciones);
long double aproximar_e_convergente_50(void);
void calcular_lados_geometricos(const char* tipo, int puntos, long double escala, long double lados_out[3]);

// --- BLOQUE II: CRIPTOGRAFÍA Y CHORRO EN ACORDEÓN ---
void aleatorovix_init(Aleatorovix* ax, uint16_t semilla);
uint8_t aleatorovix_siguiente_bit(Aleatorovix* ax);

long double encriptar_semilla_fase(uint16_t semilla, const ClaveK3* clave);
int desencriptar_semilla_fase(long double hash_fase, const ClaveK3* clave);

void cifrar_masivo_xor(const uint8_t* datos_claros, uint8_t* datos_cifrados, uint32_t longitud_bytes, uint16_t semilla);

#endif // FASE_K3_ALEATOROVIX_H