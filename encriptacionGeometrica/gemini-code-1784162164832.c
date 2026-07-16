#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h> // Para mlock

// Estructura de fase industrial
typedef struct {
    uint64_t v; // Talla (impar)
    uint64_t L; // Longitud de onda (5L > 2v + 1)
} K3Motor;

// Cifra el bloque y destruye el rastro en RAM
void k3_cifrar_bloque(uint8_t *data, size_t len, K3Motor *m) {
    // 1. Bloquear memoria (evita que el SO escriba esto al disco/swap)
    mlock(data, len);

    // 2. Ejecución del motor (Acordeón de fase)
    for(size_t i = 0; i < len; i++) {
        m->L += (m->v + 1);
        m->v += 2;
        // Regla de estabilidad industrial (Pinza de convergencia)
        if ((5 * m->L) <= (2 * m->v + 1)) m->L += (m->v * 2);
        
        // Transformación inyectiva (Constante Aurea)
        uint64_t stream = (m->L ^ m->v) * 0x9E3779B97F4A7C15ULL;
        data[i] ^= (uint8_t)(stream & 0xFF);
    }

    // 3. SUICIDIO DE MEMORIA (Criba Desmemoriada)
    // Borramos físicamente los datos antes de liberar el puntero
    memset(data, 0, len);
    munlock(data, len);
}