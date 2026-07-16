#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h> // Para mlock

typedef struct {
    uint64_t v; // Talla (impar)
    uint64_t L; // Longitud (5L > 2v + 1)
} K3Motor;

// Cifra un buffer y purga inmediatamente
void k3_cifrar_bloque(uint8_t *data, size_t len, K3Motor *m) {
    // Bloquear memoria para evitar swap (anti-forense)
    mlock(data, len);

    for(size_t i = 0; i < len; i++) {
        m->L += (m->v + 1);
        m->v += 2;
        if ((5 * m->L) <= (2 * m->v + 1)) m->L += (m->v * 2);
        
        uint64_t stream = (m->L ^ m->v) * 0x9E3779B97F4A7C15ULL;
        data[i] ^= (uint8_t)(stream & 0xFF);
    }

    // Suicidio de memoria: borrar rastro en RAM
    memset(data, 0, len);
    munlock(data, len);
}