// anexoF_criba6kpm1_openmp.c
// Compilar:  gcc -O3 -fopenmp anexoF_criba6kpm1_openmp.c -o criba
// Ejecutar:  ./criba 10000000
// Salida: contador de primos y tiempo aproximado.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>

static inline int is_cand(uint64_t v){ return (v%6==1) || (v%6==5); }
static inline uint64_t value(uint64_t i){ return 2*i + 3; } // 5,7,11,13,...

int main(int argc, char** argv){
    if(argc<2){ fprintf(stderr, "Uso: %s N\n", argv[0]); return 1; }
    uint64_t N = strtoull(argv[1], NULL, 10);
    double t0 = omp_get_wtime();

    // Array de candidatos plano [0..N], 0=no cand/compuesto, 1=candidato, 2=primo
    unsigned char *P = (unsigned char*)calloc(N+1, 1);
    if(!P){ perror("calloc"); return 1; }

    // Semilla: marcar candidatos 6k±1 y 2,3
    P[2]=2; P[3]=2;
    for(uint64_t v=5; v<=N; ++v) if(is_cand(v)) P[v]=1; // candidatos

    uint64_t limit = (uint64_t)sqrt((long double)N);

    // Criba principal: marcado único visitando triángulo superior pn<=qm
    for(uint64_t n=0; ; ++n){
        uint64_t pn = value(n);
        if(pn>limit) break;
        if(pn>N) break;
        if(P[pn]==0) continue; // no candidato o ya compuesto

        // Confirmar primo
        P[pn]=2;

        // Marcar productos pn*qm en progresiones +2p/+4p
        uint64_t m = n; // pn <= qm garantiza marcado único
        uint64_t Nmax = N;
        while(1){
            uint64_t qm = value(m);
            if(qm==0) break; // overflow guard
            if(qm > Nmax / pn) break; // fuera de rango
            if(qm%3==0){ m++; continue; } // anteriorM: salta no-candidatos
            if(!is_cand(qm)){ m++; continue; }

            uint64_t comp = pn*qm;
            uint64_t j = comp;
            uint64_t salto1 = 2*pn, salto2 = 4*pn;
            int toggle = 1;
            while(j<=Nmax){
                if(P[j]==1) P[j]=0; // marcar compuesto (de candidato a no-candidato)
                j += toggle ? salto1 : salto2;
                toggle = !toggle;
            }
            m++;
        }
    }

    // Contar primos
    uint64_t count=0;
    for(uint64_t v=2; v<=N; ++v) if(P[v]==2) count++;

    double t1 = omp_get_wtime();
    printf("Primos ≤ %llu: %llu\n", (unsigned long long)N, (unsigned long long)count);
    printf("Tiempo: %.3f s\n", t1-t0);

    free(P);
    return 0;
}
