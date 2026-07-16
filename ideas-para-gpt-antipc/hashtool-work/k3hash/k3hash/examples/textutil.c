#include "textutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char* leer_fichero_completo(const char* ruta, size_t* longitud_salida) {
    FILE* f = fopen(ruta, "rb");
    if (!f) return NULL;

    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    long tam = ftell(f);
    if (tam < 0) { fclose(f); return NULL; }
    if (fseek(f, 0, SEEK_SET) != 0) { fclose(f); return NULL; }

    char* buffer = (char*)malloc((size_t)tam + 1);
    if (!buffer) { fclose(f); return NULL; }

    size_t leidos = fread(buffer, 1, (size_t)tam, f);
    fclose(f);
    buffer[leidos] = '\0';

    if (longitud_salida) *longitud_salida = leidos;
    return buffer;
}

char** tokenizar_palabras(const char* texto, int* num_palabras_salida) {
    int cap = 256, n = 0;
    char** tokens = (char**)malloc((size_t)cap * sizeof(char*));
    size_t len = strlen(texto);
    size_t i = 0;

    while (i < len) {
        while (i < len && !isalnum((unsigned char)texto[i])) i++;
        size_t inicio = i;
        while (i < len && isalnum((unsigned char)texto[i])) i++;
        size_t palabra_len = i - inicio;
        if (palabra_len == 0) continue;

        if (n == cap) {
            cap *= 2;
            tokens = (char**)realloc(tokens, (size_t)cap * sizeof(char*));
        }
        char* palabra = (char*)malloc(palabra_len + 1);
        for (size_t k = 0; k < palabra_len; k++) {
            palabra[k] = (char)tolower((unsigned char)texto[inicio + k]);
        }
        palabra[palabra_len] = '\0';
        tokens[n++] = palabra;
    }

    *num_palabras_salida = n;
    return tokens;
}

void liberar_tokens(char** tokens, int num_palabras) {
    for (int i = 0; i < num_palabras; i++) free(tokens[i]);
    free(tokens);
}
