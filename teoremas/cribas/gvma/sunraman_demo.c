/**
 * Demo SunRaman Jet — compilar: gcc sunraman_demo.c -o sunraman_demo -lm
 */
#include <stdio.h>
#include <stdlib.h>

#include "Gvma.h"

static void demo_jet(uint64_t limite_v) {
    uint8_t *mapa = (uint8_t *)calloc((size_t)limite_v + 1u, 1);
    if (!mapa) {
        return;
    }

    gvma_criba_jet(mapa, limite_v);

    printf("SunRaman Jet — primos 2v+3 hasta v=%llu\n", (unsigned long long)limite_v);
    printf("  v=0 -> N=3\n");

    unsigned count = 0;
    for (uint64_t v = 1; v <= limite_v; v++) {
        if (!mapa[v]) {
            printf("  v=%llu -> N=%llu\n", (unsigned long long)v,
                   (unsigned long long)gvma_v_a_n(v));
            count++;
        }
    }
    printf("Total candidatos primos (sin 2): %u\n", count);
    free(mapa);
}

static void demo_flash(void) {
    uint64_t tests[] = {14, 30, 60};
    for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
        uint64_t v = tests[i];
        GvmaFlashResult r = gvma_factorizar_flash(v);
        printf("Flash v=%llu (N=%llu): ", (unsigned long long)v,
               (unsigned long long)gvma_v_a_n(v));
        if (r.es_compuesto) {
            printf("compuesto %llu * %llu\n", (unsigned long long)r.f1,
                   (unsigned long long)r.f2);
        } else {
            printf("sin factor en pinza\n");
        }
    }
}

int main(void) {
    demo_jet(120);
    demo_flash();
    return 0;
}