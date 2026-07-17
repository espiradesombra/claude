/* Criba 6k±1 anexoF — port con conteo corregido (candidatos 1 ∪ confirmados 2).
 * Fuente: archivos-vma/codigo/anexoF_criba6kpm1_openmp.c
 * Compilar (con OpenMP timing):
 *   gcc -O3 -fopenmp src/criba_6k_openmp.c -o criba_6k -lm
 * Sin OpenMP:
 *   gcc -O3 src/criba_6k_openmp.c -o criba_6k -lm -DNO_OPENMP
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef NO_OPENMP
#include <omp.h>
#else
#include <time.h>
static double omp_get_wtime(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}
#endif

static inline int is_cand(uint64_t v) {
    return (v % 6 == 1) || (v % 6 == 5);
}
static inline uint64_t value(uint64_t i) { return 2 * i + 3; }

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Uso: %s N\n", argv[0]);
        return 1;
    }
    uint64_t N = strtoull(argv[1], NULL, 10);
    double t0 = omp_get_wtime();

    unsigned char *P = (unsigned char *)calloc(N + 1, 1);
    if (!P) {
        perror("calloc");
        return 1;
    }

    if (N >= 2) P[2] = 2;
    if (N >= 3) P[3] = 2;
    for (uint64_t v = 5; v <= N; ++v)
        if (is_cand(v)) P[v] = 1;

    uint64_t limit = (uint64_t)sqrt((long double)N);

    for (uint64_t n = 0;; ++n) {
        uint64_t pn = value(n);
        if (pn > limit || pn > N) break;
        if (P[pn] == 0) continue;
        P[pn] = 2;

        uint64_t m = n;
        while (1) {
            uint64_t qm = value(m);
            if (qm == 0) break;
            if (qm > N / pn) break;
            if (qm % 3 == 0 || !is_cand(qm)) {
                m++;
                continue;
            }
            uint64_t j = pn * qm;
            uint64_t salto1 = 2 * pn, salto2 = 4 * pn;
            int toggle = 1;
            while (j <= N) {
                if (P[j] == 1) P[j] = 0;
                j += toggle ? salto1 : salto2;
                toggle = !toggle;
            }
            m++;
        }
    }

    /* Corrección: supervivientes (1) también son primos > √N */
    uint64_t count = 0;
    for (uint64_t v = 2; v <= N; ++v)
        if (P[v] >= 1) count++;

    double t1 = omp_get_wtime();
    printf("Primos <= %llu: %llu\n", (unsigned long long)N, (unsigned long long)count);
    printf("Tiempo: %.3f s\n", t1 - t0);
    free(P);
    return 0;
}
