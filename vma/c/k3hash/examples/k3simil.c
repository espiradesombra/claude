/*
 * k3simil.c — Detector de similitud entre textos (plagio parcial).
 *
 * Trocea cada documento en "shingles" (frases solapadas de K_SHINGLE
 * palabras), hashea cada shingle con K3, y compara los conjuntos de
 * hashes entre documentos con el índice de Jaccard:
 *
 *     similitud = |interseccion| / |union|
 *
 * Esto SÍ detecta reformulaciones parciales (frases movidas, párrafos
 * copiados con cambios menores), a diferencia de un hash de archivo
 * completo que solo detecta copia exacta byte a byte.
 *
 * Uso:
 *   find /ruta -name "*.txt" | k3simil [umbral=0.30]
 *
 * Nota: usa hash de 32 bits para los shingles; en corpus muy grandes
 * (decenas de miles de shingles únicos) la probabilidad de colisión deja
 * de ser despreciable. Para corpus grandes, pide upgrade a hash de 64 bits.
 */
#include "k3hash.h"
#include "textutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define K_SHINGLE 5

typedef struct {
    char* ruta;
    uint32_t* shingles; /* ordenados, sin duplicados */
    int num_shingles;
} DocSimil;

static int cmp_u32(const void* a, const void* b) {
    uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

static int construir_shingles(const char* ruta, const K3HashConfig* cfg, DocSimil* doc) {
    size_t len;
    char* texto = leer_fichero_completo(ruta, &len);
    if (!texto) return -1;

    int num_palabras;
    char** palabras = tokenizar_palabras(texto, &num_palabras);
    free(texto);

    doc->ruta = strdup(ruta);

    if (num_palabras < K_SHINGLE) {
        liberar_tokens(palabras, num_palabras);
        doc->shingles = NULL;
        doc->num_shingles = 0;
        return 0;
    }

    int num_bruto = num_palabras - K_SHINGLE + 1;
    uint32_t* hashes = (uint32_t*)malloc((size_t)num_bruto * sizeof(uint32_t));
    char buffer_shingle[2048];

    for (int i = 0; i < num_bruto; i++) {
        size_t pos = 0;
        for (int k = 0; k < K_SHINGLE; k++) {
            size_t l = strlen(palabras[i + k]);
            if (pos + l + 1 >= sizeof(buffer_shingle)) break;
            memcpy(buffer_shingle + pos, palabras[i + k], l);
            pos += l;
            buffer_shingle[pos++] = ' ';
        }
        buffer_shingle[pos] = '\0';
        hashes[i] = k3_hash_buffer(buffer_shingle, pos, cfg);
    }

    qsort(hashes, (size_t)num_bruto, sizeof(uint32_t), cmp_u32);
    int num_unicos = 0;
    for (int i = 0; i < num_bruto; i++) {
        if (i == 0 || hashes[i] != hashes[i - 1]) hashes[num_unicos++] = hashes[i];
    }

    doc->shingles = hashes;
    doc->num_shingles = num_unicos;

    liberar_tokens(palabras, num_palabras);
    return 0;
}

static double jaccard(const DocSimil* a, const DocSimil* b) {
    if (a->num_shingles == 0 || b->num_shingles == 0) return 0.0;
    int i = 0, j = 0, interseccion = 0;
    while (i < a->num_shingles && j < b->num_shingles) {
        if (a->shingles[i] == b->shingles[j]) { interseccion++; i++; j++; }
        else if (a->shingles[i] < b->shingles[j]) i++;
        else j++;
    }
    int union_total = a->num_shingles + b->num_shingles - interseccion;
    return union_total > 0 ? (double)interseccion / (double)union_total : 0.0;
}

int main(int argc, char** argv) {
    double umbral = 0.30;
    if (argc > 1) umbral = atof(argv[1]);

    K3HashConfig cfg = k3_config_default();
    DocSimil* docs = NULL;
    int n = 0, cap = 0;
    char linea[4096];

    while (fgets(linea, sizeof(linea), stdin)) {
        size_t len = strlen(linea);
        while (len > 0 && (linea[len - 1] == '\n' || linea[len - 1] == '\r')) linea[--len] = '\0';
        if (len == 0) continue;

        if (n == cap) {
            cap = cap ? cap * 2 : 16;
            docs = (DocSimil*)realloc(docs, (size_t)cap * sizeof(DocSimil));
        }
        if (construir_shingles(linea, &cfg, &docs[n]) == 0) {
            n++;
        } else {
            fprintf(stderr, "[aviso] no se pudo leer: %s\n", linea);
        }
    }

    if (n == 0) {
        fprintf(stderr, "[k3simil] no se recibieron rutas de archivo por stdin.\n");
        return 1;
    }

    fprintf(stderr, "[k3simil] %d documentos indexados (shingles de %d palabras). Umbral: %.0f%%\n\n",
            n, K_SHINGLE, umbral * 100);

    int num_pares_similares = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double sim = jaccard(&docs[i], &docs[j]);
            if (sim >= umbral) {
                printf("%5.1f%%  %s  <->  %s\n", sim * 100.0, docs[i].ruta, docs[j].ruta);
                num_pares_similares++;
            }
        }
    }

    if (num_pares_similares == 0) {
        fprintf(stderr, "[k3simil] ningun par supero el umbral de %.0f%%.\n", umbral * 100);
    }

    for (int i = 0; i < n; i++) { free(docs[i].ruta); free(docs[i].shingles); }
    free(docs);
    return 0;
}
