/*
 * k3search.c — Motor de búsqueda de contenido (índice invertido).
 *
 * Subcomandos:
 *   k3search index <indice.k3idx>          (lee rutas de archivos por stdin)
 *   k3search query <indice.k3idx> palabra1 palabra2 ...
 *
 * Uso típico:
 *   find /mi/biblioteca -name "*.txt" | k3search index biblioteca.k3idx
 *   k3search query biblioteca.k3idx factorizacion primos
 *
 * Internamente cada término se hashea con K3 para ubicarlo en una tabla
 * hash de acceso O(1) (con resolución de colisiones por comparación de
 * string real, no solo por hash). El índice guardado es texto plano,
 * portable entre Windows/macOS/Linux sin conversiones.
 */
#include "k3hash.h"
#include "textutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define TABLA_TAM 65536 /* potencia de 2 */

typedef struct { int doc_id; int contador; } Posting;

typedef struct {
    char* termino;
    Posting* postings;
    int num_postings;
    int cap_postings;
} EntradaTermino;

static uint32_t hash_termino(const char* termino, const K3HashConfig* cfg) {
    return k3_hash_buffer(termino, strlen(termino), cfg);
}

/* --------------------------- construcción del índice --------------------------- */

typedef struct {
    EntradaTermino* terminos;
    int num_terminos;
    int cap_terminos;
    int* tabla_hash;
} IndiceConstruccion;

static int buscar_o_crear_termino(IndiceConstruccion* idx, const char* termino, const K3HashConfig* cfg) {
    uint32_t h = hash_termino(termino, cfg);
    unsigned int pos = h & (TABLA_TAM - 1);
    while (idx->tabla_hash[pos] != -1) {
        int existente = idx->tabla_hash[pos];
        if (strcmp(idx->terminos[existente].termino, termino) == 0) return existente;
        pos = (pos + 1) & (TABLA_TAM - 1);
    }
    if (idx->num_terminos == idx->cap_terminos) {
        idx->cap_terminos = idx->cap_terminos ? idx->cap_terminos * 2 : 1024;
        idx->terminos = (EntradaTermino*)realloc(idx->terminos, (size_t)idx->cap_terminos * sizeof(EntradaTermino));
    }
    int nuevo_id = idx->num_terminos++;
    idx->terminos[nuevo_id].termino = strdup(termino);
    idx->terminos[nuevo_id].postings = NULL;
    idx->terminos[nuevo_id].num_postings = 0;
    idx->terminos[nuevo_id].cap_postings = 0;
    idx->tabla_hash[pos] = nuevo_id;
    return nuevo_id;
}

static void anadir_posting(EntradaTermino* entrada, int doc_id) {
    if (entrada->num_postings > 0 && entrada->postings[entrada->num_postings - 1].doc_id == doc_id) {
        entrada->postings[entrada->num_postings - 1].contador++;
        return;
    }
    if (entrada->num_postings == entrada->cap_postings) {
        entrada->cap_postings = entrada->cap_postings ? entrada->cap_postings * 2 : 4;
        entrada->postings = (Posting*)realloc(entrada->postings, (size_t)entrada->cap_postings * sizeof(Posting));
    }
    entrada->postings[entrada->num_postings].doc_id = doc_id;
    entrada->postings[entrada->num_postings].contador = 1;
    entrada->num_postings++;
}

static int comando_index(const char* ruta_salida) {
    K3HashConfig cfg = k3_config_default();

    IndiceConstruccion idx;
    idx.terminos = NULL;
    idx.num_terminos = 0;
    idx.cap_terminos = 0;
    idx.tabla_hash = (int*)malloc(TABLA_TAM * sizeof(int));
    for (int i = 0; i < TABLA_TAM; i++) idx.tabla_hash[i] = -1;

    char** rutas_docs = NULL;
    int num_docs = 0, cap_docs = 0;
    char linea[4096];

    while (fgets(linea, sizeof(linea), stdin)) {
        size_t len = strlen(linea);
        while (len > 0 && (linea[len - 1] == '\n' || linea[len - 1] == '\r')) linea[--len] = '\0';
        if (len == 0) continue;

        size_t longitud_texto;
        char* texto = leer_fichero_completo(linea, &longitud_texto);
        if (!texto) { fprintf(stderr, "[aviso] no se pudo leer: %s\n", linea); continue; }

        if (num_docs == cap_docs) {
            cap_docs = cap_docs ? cap_docs * 2 : 64;
            rutas_docs = (char**)realloc(rutas_docs, (size_t)cap_docs * sizeof(char*));
        }
        int doc_id = num_docs;
        rutas_docs[doc_id] = strdup(linea);
        num_docs++;

        int num_palabras;
        char** palabras = tokenizar_palabras(texto, &num_palabras);
        free(texto);

        for (int p = 0; p < num_palabras; p++) {
            int term_id = buscar_o_crear_termino(&idx, palabras[p], &cfg);
            anadir_posting(&idx.terminos[term_id], doc_id);
        }
        liberar_tokens(palabras, num_palabras);
    }

    if (num_docs == 0) {
        fprintf(stderr, "[k3search] no se recibieron rutas de archivo por stdin.\n");
        return 1;
    }

    FILE* f = fopen(ruta_salida, "w");
    if (!f) { fprintf(stderr, "[ERROR] no se pudo crear el indice: %s\n", ruta_salida); return 1; }

    fprintf(f, "K3SEARCH_INDEX_V1\n");
    fprintf(f, "DOCS %d\n", num_docs);
    for (int i = 0; i < num_docs; i++) fprintf(f, "%d\t%s\n", i, rutas_docs[i]);

    fprintf(f, "TERMS %d\n", idx.num_terminos);
    for (int t = 0; t < idx.num_terminos; t++) {
        EntradaTermino* e = &idx.terminos[t];
        fprintf(f, "%s\t%d", e->termino, e->num_postings);
        for (int p = 0; p < e->num_postings; p++) fprintf(f, "\t%d:%d", e->postings[p].doc_id, e->postings[p].contador);
        fprintf(f, "\n");
    }
    fclose(f);

    fprintf(stderr, "[k3search] indice creado: %d documentos, %d terminos unicos -> %s\n",
            num_docs, idx.num_terminos, ruta_salida);
    return 0;
}

/* --------------------------- consulta del índice --------------------------- */

typedef struct {
    char** rutas_docs;
    int num_docs;
    EntradaTermino* terminos;
    int num_terminos;
    int* tabla_hash;
} IndiceCargado;

static int cargar_indice(const char* ruta, const K3HashConfig* cfg, IndiceCargado* out) {
    FILE* f = fopen(ruta, "r");
    if (!f) return -1;

    char linea[8192];
    if (!fgets(linea, sizeof(linea), f) || strncmp(linea, "K3SEARCH_INDEX_V1", 17) != 0) {
        fclose(f);
        return -1;
    }

    int num_docs = 0;
    if (!fgets(linea, sizeof(linea), f) || sscanf(linea, "DOCS %d", &num_docs) != 1) { fclose(f); return -1; }
    char** rutas = (char**)malloc((size_t)num_docs * sizeof(char*));
    for (int i = 0; i < num_docs; i++) {
        if (!fgets(linea, sizeof(linea), f)) { fclose(f); return -1; }
        int id;
        char ruta_doc[4096];
        sscanf(linea, "%d\t%4095[^\n]", &id, ruta_doc);
        rutas[id] = strdup(ruta_doc);
    }

    int num_terminos = 0;
    if (!fgets(linea, sizeof(linea), f) || sscanf(linea, "TERMS %d", &num_terminos) != 1) { fclose(f); return -1; }
    EntradaTermino* terminos = (EntradaTermino*)malloc((size_t)num_terminos * sizeof(EntradaTermino));

    for (int t = 0; t < num_terminos; t++) {
        if (!fgets(linea, sizeof(linea), f)) { fclose(f); return -1; }
        char* tok = strtok(linea, "\t\n");
        terminos[t].termino = strdup(tok ? tok : "");
        tok = strtok(NULL, "\t\n");
        int npost = tok ? atoi(tok) : 0;
        terminos[t].num_postings = npost;
        terminos[t].cap_postings = npost;
        terminos[t].postings = (Posting*)malloc((size_t)(npost > 0 ? npost : 1) * sizeof(Posting));
        for (int p = 0; p < npost; p++) {
            tok = strtok(NULL, "\t\n");
            int doc_id = 0, contador = 0;
            if (tok) sscanf(tok, "%d:%d", &doc_id, &contador);
            terminos[t].postings[p].doc_id = doc_id;
            terminos[t].postings[p].contador = contador;
        }
    }
    fclose(f);

    int* tabla = (int*)malloc(TABLA_TAM * sizeof(int));
    for (int i = 0; i < TABLA_TAM; i++) tabla[i] = -1;
    for (int t = 0; t < num_terminos; t++) {
        uint32_t h = hash_termino(terminos[t].termino, cfg);
        unsigned int pos = h & (TABLA_TAM - 1);
        while (tabla[pos] != -1) pos = (pos + 1) & (TABLA_TAM - 1);
        tabla[pos] = t;
    }

    out->rutas_docs = rutas;
    out->num_docs = num_docs;
    out->terminos = terminos;
    out->num_terminos = num_terminos;
    out->tabla_hash = tabla;
    return 0;
}

static int buscar_termino(const IndiceCargado* idx, const char* termino, const K3HashConfig* cfg) {
    uint32_t h = hash_termino(termino, cfg);
    unsigned int pos = h & (TABLA_TAM - 1);
    while (idx->tabla_hash[pos] != -1) {
        int candidato = idx->tabla_hash[pos];
        if (strcmp(idx->terminos[candidato].termino, termino) == 0) return candidato;
        pos = (pos + 1) & (TABLA_TAM - 1);
    }
    return -1;
}

static int comando_query(const char* ruta_indice, char** palabras_consulta, int num_palabras) {
    K3HashConfig cfg = k3_config_default();
    IndiceCargado idx;
    if (cargar_indice(ruta_indice, &cfg, &idx) != 0) {
        fprintf(stderr, "[ERROR] no se pudo cargar el indice: %s\n", ruta_indice);
        return 1;
    }

    double* puntuaciones = (double*)calloc((size_t)idx.num_docs, sizeof(double));

    for (int q = 0; q < num_palabras; q++) {
        char normalizada[512];
        int k = 0;
        for (const char* c = palabras_consulta[q]; *c && k < 511; c++) {
            if (isalnum((unsigned char)*c)) normalizada[k++] = (char)tolower((unsigned char)*c);
        }
        normalizada[k] = '\0';

        int t = buscar_termino(&idx, normalizada, &cfg);
        if (t < 0) {
            fprintf(stderr, "[k3search] sin resultados para: %s\n", normalizada);
            continue;
        }
        for (int p = 0; p < idx.terminos[t].num_postings; p++) {
            puntuaciones[idx.terminos[t].postings[p].doc_id] += idx.terminos[t].postings[p].contador;
        }
    }

    /* orden simple por puntuación descendente */
    int* orden = (int*)malloc((size_t)idx.num_docs * sizeof(int));
    for (int i = 0; i < idx.num_docs; i++) orden[i] = i;
    for (int i = 0; i < idx.num_docs; i++) {
        for (int j = i + 1; j < idx.num_docs; j++) {
            if (puntuaciones[orden[j]] > puntuaciones[orden[i]]) {
                int tmp = orden[i]; orden[i] = orden[j]; orden[j] = tmp;
            }
        }
    }

    int mostrados = 0;
    for (int i = 0; i < idx.num_docs && mostrados < 20; i++) {
        if (puntuaciones[orden[i]] <= 0) break;
        printf("%6.1f  %s\n", puntuaciones[orden[i]], idx.rutas_docs[orden[i]]);
        mostrados++;
    }
    if (mostrados == 0) fprintf(stderr, "[k3search] ningun documento coincide con la consulta.\n");

    free(orden);
    free(puntuaciones);
    return 0;
}

static void imprimir_uso(const char* prog) {
    fprintf(stderr, "Uso:\n");
    fprintf(stderr, "  %s index <indice.k3idx>          (rutas por stdin)\n", prog);
    fprintf(stderr, "  %s query <indice.k3idx> palabra1 [palabra2 ...]\n", prog);
}

int main(int argc, char** argv) {
    if (argc < 3) { imprimir_uso(argv[0]); return 2; }

    if (strcmp(argv[1], "index") == 0) {
        return comando_index(argv[2]);
    } else if (strcmp(argv[1], "query") == 0) {
        if (argc < 4) { imprimir_uso(argv[0]); return 2; }
        return comando_query(argv[2], &argv[3], argc - 3);
    }

    imprimir_uso(argv[0]);
    return 2;
}
