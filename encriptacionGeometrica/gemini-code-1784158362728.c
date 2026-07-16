#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fase_k3_aleatorovix.h"

// =====================================================================
// BLOC I: MATEMÁTICA PURA DE ALTA PRECISIÓN (long double)
// =====================================================================

long double taylor_sin(long double x) {
    long double termino = x;
    long double suma = x;
    long double x_cuadrado = x * x;
    int n = 1;
    // Precisión extrema para el canal de fase
    while (fabsl(termino) > 1e-20L) {
        termino = -termino * x_cuadrado / (long double)((2 * n) * (2 * n + 1));
        suma += termino;
        n++;
    }
    return suma;
}

long double aproximar_pi_puro(int p1, int p2, int iteraciones) {
    long double lados = (long double)(p1 * p2);
    long double pi_inicial = 3.14159265358979323846264338327950288L;
    
    long double angulo = pi_inicial / lados;
    long double lado = 2.0L * taylor_sin(angulo);
    long double perimetro = lados * lado;
    long double pi_actual = perimetro / 2.0L;

    for (int i = 0; i < iteraciones; i++) {
        lados *= 2.0L;
        long double raiz_soporte = sqrtl(4.0L - lado * lado);
        lado = lado / sqrtl(2.0L + raiz_soporte);
        perimetro = lados * lado;
        pi_actual = perimetro / 2.0L;
    }
    return pi_actual;
}

long double aproximar_e_convergente_50(void) {
    long double termino = 1.0L;
    long double acumulado_e = 1.0L;
    
    for (int v = 1; v <= 50; v++) {
        int eleccion = (v % 2 == 0) ? (2 * v + 2) : (3 * v);
        if (eleccion == 0) eleccion = 2 * v + 2;
        
        termino = termino / (long double)eleccion;
        acumulado_e += termino;
    }
    return acumulado_e;
}

void calcular_lados_geometricos(const char* tipo, int puntos, long double escala, long double lados_out[3]) {
    long double base = (long double)puntos * 1.5L;
    if (strcmp(tipo, "equilatero") == 0) {
        long double lado = (base / 3.0L) * escala;
        lados_out[0] = lado; lados_out[1] = lado; lados_out[2] = lado;
    } else if (strcmp(tipo, "isosceles") == 0) {
        long double lado_a = (base / 4.0L) * escala;
        long double lado_b = (base / 2.0L) * escala;
        lados_out[0] = lado_a; lados_out[1] = lado_a; lados_out[2] = lado_b;
    } else { // escaleno
        lados_out[0] = (base * 0.25L) * escala;
        lados_out[1] = (base * 0.35L) * escala;
        lados_out[2] = (base * 0.40L) * escala;
    }
}


// =====================================================================
// BLOC II: MOTOR ALEATOROVIX EN ACORDEÓN (MSL/LSL)
// =====================================================================

void aleatorovix_init(Aleatorovix* ax, uint16_t semilla) {
    ax->estado = (semilla != 0) ? (long double)semilla : 12345.0L;
    ax->factor_msl = 42.123456789L;
    ax->factor_lsl = 0.000123456L;
    
    // Normalizar estado inicial al intervalo [0, 1)
    ax->estado = ax->estado - floorl(ax->estado);
    if (ax->estado == 0.0L) ax->estado = 0.5L;
}

uint8_t aleatorovix_siguiente_bit(Aleatorovix* ax) {
    // Sincronización caótica: x_{n+1} = (x_n * factor_msl + cos(x_n * factor_lsl)) mod 1
    long double fase = (ax->estado * ax->factor_msl) + cosl(ax->estado * ax->factor_lsl);
    ax->estado = fase - floorl(fase);
    
    return (ax->estado > 0.5L) ? 1 : 0;
}


// =====================================================================
// BLOC III: BACKTRACKING GEOMÉTRICO (long double)
// =====================================================================

long double encriptar_semilla_fase(uint16_t semilla, const ClaveK3* clave) {
    long double perimetro = 0.0L;
    
    for (int idx = 0; idx < 16; idx++) {
        uint8_t bit = (semilla >> (15 - idx)) & 1;
        int figura_idx = idx / 3;
        int lado_idx = idx % 3;

        // Saltar lado si coincide con el patrón
        if (lado_idx == clave->saltos[figura_idx % clave->saltos_count]) {
            continue;
        }

        long double escala = (long double)clave->tales[figura_idx % clave->tales_count];
        const char* tipo = clave->figuras[figura_idx % clave->figuras_count];
        int pts = clave->puntos[figura_idx % clave->puntos_count];

        long double lados[3];
        calcular_lados_geometricos(tipo, pts, escala, lados);

        if (bit == 1) {
            perimetro += lados[lado_idx];
        }
    }

    // Ofuscación multiplicativa por constantes
    long double pi_ofuscado = 1.0L;
    for (int i = 0; i < 2; i++) {
        long double pi_aprox = aproximar_pi_puro(clave->primos[i][0], clave->primos[i][1], clave->iteraciones_pi);
        long double aportacion = (long double)clave->porcentajes_aportacion[i] / 100.0L;
        pi_ofuscado *= powl(pi_aprox, aportacion);
    }

    long double e_convergente = aproximar_e_convergente_50();
    
    return perimetro * pi_ofuscado * e_convergente;
}

// Función interna recursiva de backtracking para C
static bool backtrack_recursivo(long double residuo, int idx, const ClaveK3* clave, uint16_t* bits_out) {
    if (idx == 16) {
        return (fabsl(residuo) < 1e-12L); // Tolerancia física del canal
    }

    int figura_idx = idx / 3;
    int lado_idx = idx % 3;

    if (lado_idx == clave->saltos[figura_idx % clave->saltos_count]) {
        // Bifurcación en salto de fase (probamos bit 0 y bit 1)
        for (int bit = 0; bit <= 1; bit++) {
            if (backtrack_recursivo(residuo, idx + 1, clave, bits_out)) {
                *bits_out |= (bit << (15 - idx));
                return true;
            }
        }
        return false;
    }

    long double escala = (long double)clave->tales[figura_idx % clave->tales_count];
    const char* tipo = clave->figuras[figura_idx % clave->figuras_count];
    int pts = clave->puntos[figura_idx % clave->puntos_count];

    long double lados[3];
    calcular_lados_geometricos(tipo, pts, escala, lados);
    long double lado_actual = lados[lado_idx];

    // Poda: Si el residuo permite restar el lado geométrico, exploramos el camino del bit 1
    if (residuo >= lado_actual - 1e-12L) {
        if (backtrack_recursivo(residuo - lado_actual, idx + 1, clave, bits_out)) {
            *bits_out |= (1 << (15 - idx));
            return true;
        }
    }

    // Exploramos el camino del bit 0 (ignora lado)
    if (backtrack_recursivo(residuo, idx + 1, clave, bits_out)) {
        *bits_out &= ~(1 << (15 - idx));
        return true;
    }

    return false;
}

int desencriptar_semilla_fase(long double hash_fase, const ClaveK3* clave) {
    long double pi_ofuscado = 1.0L;
    for (int i = 0; i < 2; i++) {
        long double pi_aprox = aproximar_pi_puro(clave->primos[i][0], clave->primos[i][1], clave->iteraciones_pi);
        long double aportacion = (long double)clave->porcentajes_aportacion[i] / 100.0L;
        pi_ofuscado *= powl(pi_aprox, aportacion);
    }

    long double e_convergente = aproximar_e_convergente_50();
    long double perimetro_original = hash_fase / (pi_ofuscado * e_convergente);

    uint16_t semilla_recuperada = 0;
    if (backtrack_recursivo(perimetro_original, 0, clave, &semilla_recuperada)) {
        return (int)semilla_recuperada;
    }
    
    return -1; // Fallo de sincronización de fase
}

void cifrar_masivo_xor(const uint8_t* datos_claros, uint8_t* datos_cifrados, uint32_t longitud_bytes, uint16_t semilla) {
    Aleatorovix ax;
    aleatorovix_init(&ax, semilla);

    for (uint32_t i = 0; i < longitud_bytes; i++) {
        uint8_t byte_claros = datos_claros[i];
        uint8_t byte_cifrado = 0;

        for (int bit = 0; bit < 8; bit++) {
            uint8_t bit_dato = (byte_claros >> (7 - bit)) & 1;
            uint8_t bit_chorro = aleatorovix_siguiente_bit(&ax);
            uint8_t bit_xor = bit_dato ^ bit_chorro;
            byte_cifrado |= (bit_xor << (7 - bit));
        }
        datos_cifrados[i] = byte_cifrado;
    }
}


// =====================================================================
// DEMO / TEST INDUSTRIAL DE INTEGRACIÓN DE RED
// =====================================================================
int main(void) {
    // Configuración del Desafío K3 del Servidor
    ClaveK3 clave_servidor = {
        .tales = {3, 5, 8, 13, 21},
        .tales_count = 5,
        .figuras = {"equilatero", "isosceles", "escaleno"},
        .figuras_count = 3,
        .puntos = {6, 12, 18},
        .puntos_count = 3,
        .saltos = {1, 0, 2},
        .saltos_count = 3,
        .primos = {{3, 5}, {7, 11}},
        .porcentajes_aportacion = {50, 100},
        .iteraciones_pi = 15
    };

    uint16_t semilla_original = 43210;
    const char* datos = "ESPANA_PLAN_2030_CONEXION_DE_FASE_K3";
    uint32_t tamano = strlen(datos);

    printf("[+] Mensaje Masivo Original: '%s'\n", datos);
    printf("[+] Semilla de fase secreta elegida: %d\n\n", semilla_original);

    // 1. El Cliente encripta su semilla (creando el hash para 'has.txt')
    printf("[*] Generando perímetro de fase en C...\n");
    long double hash_transmitido = encriptar_semilla_fase(semilla_original, &clave_servidor);
    printf("[+] Phase Hash de la Semilla: %.15Lf\n\n", hash_transmitido);

    // 2. El Cliente cifra el bloque de datos usando Aleatorovix nativo
    uint8_t* buffer_cifrado = malloc(tamano);
    cifrar_masivo_xor((const uint8_t*)datos, buffer_cifrado, tamano, semilla_original);
    
    printf("[+] Datos cifrados (en binario de red): ");
    for (uint32_t i = 0; i < tamano && i < 15; i++) {
        printf("%02X ", buffer_cifrado[i]);
    }
    printf("...\n\n");

    // 3. El Servidor recibe el Hash y el Buffer. Aplica Backtracking para recuperar la Semilla
    printf("[*] Servidor procesando Backtracking Geométrico recursivo en C...\n");
    int semilla_recuperada = desencriptar_semilla_fase(hash_transmitido, &clave_servidor);
    printf("[+] Semilla recuperada por el Servidor: %d\n", semilla_recuperada);

    if (semilla_recuperada == (int)semilla_original) {
        printf("[✅ SUCCESS] Handshake de Fase validado.\n");
        
        // El Servidor descifra los datos con la semilla recuperada
        uint8_t* buffer_descifrado = malloc(tamano + 1);
        cifrar_masivo_xor(buffer_cifrado, buffer_descifrado, tamano, (uint16_t)semilla_recuperada);
        buffer_descifrado[tamano] = '\0';

        printf("[<] Datos recuperados: '%s'\n", buffer_descifrado);
        free(buffer_descifrado);
    } else {
        printf("[❌ FAILURE] Desincronización de fase catastrófica.\n");
    }

    free(buffer_cifrado);
    return 0;
}