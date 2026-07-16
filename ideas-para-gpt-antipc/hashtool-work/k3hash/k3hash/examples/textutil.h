#ifndef TEXTUTIL_H
#define TEXTUTIL_H

#include <stddef.h>

/* Lee un fichero completo en memoria (reservada con malloc, el llamador
 * debe liberarla). Devuelve NULL si no se pudo abrir/leer. */
char* leer_fichero_completo(const char* ruta, size_t* longitud_salida);

/* Trocea un texto en palabras (solo alfanuméricas ASCII, en minúsculas).
 * Devuelve un array de strings reservado con malloc; liberar con
 * liberar_tokens(). Es una tokenización simple, no consciente de Unicode. */
char** tokenizar_palabras(const char* texto, int* num_palabras_salida);

void liberar_tokens(char** tokens, int num_palabras);

#endif /* TEXTUTIL_H */
