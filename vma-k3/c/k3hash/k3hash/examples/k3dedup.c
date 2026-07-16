/*
 * k3dedup.c — Detector de archivos duplicados exactos.
 *
 * Uso:
 *   Linux/macOS:  find /ruta -type f | k3dedup
 *   Windows:      dir /s /b C:\ruta | k3dedup
 *
 * Agrupa por (tamaño, hash K3). Dos archivos en el mismo grupo son
 * duplicados con altísima probabilidad; para CERTEZA absoluta en archivos
 * críticos, compara además byte a byte antes de borrar nada.
 */
#include "k3hash.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char* ruta;
    long  tamano;
    uint32_t hash;
} EntradaArchivo;

static int comparar_entradas(const void* a, const void* b) {
    const EntradaArchivo* ea = (const EntradaArchivo*)a;
    const EntradaArchivo* eb = (const EntradaArchivo*)b;
    if (ea->tamano != eb->tamano) return (ea->tamano > eb->tamano) - (ea->tamano < eb->tamano);
    if (ea->hash != eb->hash) return (ea->hash > eb->hash) - (ea->hash < eb->hash);
    return 0;
}

int main(void) {
    char linea[4096];
    EntradaArchivo* entradas = NULL;
    size_t n = 0, cap = 0;
    K3HashConfig cfg = k3_config_default();

    while (fgets(linea, sizeof(linea), stdin)) {
        size_t len = strlen(linea);
        while (len > 0 && (linea[len - 1] == '\n' || linea[len - 1] == '\r')) linea[--len] = '\0';
        if (len == 0) continue;

        FILE* f = fopen(linea, "rb");
        if (!f) { fprintf(stderr, "[aviso] no se pudo abrir: %s\n", linea); continue; }
        fseek(f, 0, SEEK_END);
        long tam = ftell(f);
        fclose(f);

        uint32_t h;
        if (k3_hash_file(linea, &cfg, &h) != 0) {
            fprintf(stderr, "[aviso] no se pudo leer: %s\n", linea);
            continue;
        }

        if (n == cap) {
            cap = cap ? cap * 2 : 64;
            entradas = (EntradaArchivo*)realloc(entradas, cap * sizeof(EntradaArchivo));
        }
        entradas[n].ruta = strdup(linea);
        entradas[n].tamano = tam;
        entradas[n].hash = h;
        n++;
    }

    if (n == 0) {
        fprintf(stderr, "[k3dedup] no se recibieron rutas de archivo por stdin.\n");
        return 1;
    }

    qsort(entradas, n, sizeof(EntradaArchivo), comparar_entradas);

    size_t i = 0;
    int num_grupos = 0;
    long long bytes_recuperables = 0;

    while (i < n) {
        size_t j = i + 1;
        while (j < n && entradas[j].tamano == entradas[i].tamano && entradas[j].hash == entradas[i].hash) j++;

        if (j - i > 1) {
            num_grupos++;
            printf("=== Grupo %d (tam=%ld bytes, hash=0x%08X, %zu copias) ===\n",
                   num_grupos, entradas[i].tamano, entradas[i].hash, j - i);
            for (size_t k = i; k < j; k++) printf("  %s\n", entradas[k].ruta);
            bytes_recuperables += (long long)entradas[i].tamano * (long long)(j - i - 1);
        }
        i = j;
    }

    fprintf(stderr, "\n[k3dedup] %zu archivos analizados, %d grupos de duplicados, ~%lld bytes recuperables.\n",
            n, num_grupos, bytes_recuperables);
    fprintf(stderr, "[k3dedup] Coincidencia por tamaño+hash de 32 bits. Verifica byte a byte antes de borrar archivos críticos.\n");

    for (size_t k = 0; k < n; k++) free(entradas[k].ruta);
    free(entradas);
    return 0;
}
