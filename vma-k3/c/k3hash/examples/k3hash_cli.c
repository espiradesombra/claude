/*
 * k3hash_cli.c — Ejemplo de uso de la librería k3hash.
 *
 * Uso:
 *   k3hash_cli --text "algo de texto"
 *   k3hash_cli --file ruta/al/fichero
 *
 * Compilar (si no usas CMake), por ejemplo en Linux/macOS:
 *   gcc -Iinclude src/k3hash.c examples/k3hash_cli.c -o k3hash_cli
 *
 * En Windows con MSVC (desde "Developer Command Prompt"):
 *   cl /I include src\k3hash.c examples\k3hash_cli.c /Fe:k3hash_cli.exe
 */
#include <stdio.h>
#include <string.h>
#include "k3hash.h"

int main(int argc, char** argv) {
    if (argc < 3 || (strcmp(argv[1], "--text") != 0 && strcmp(argv[1], "--file") != 0)) {
        fprintf(stderr, "Uso:\n  %s --text \"cadena\"\n  %s --file ruta\n", argv[0], argv[0]);
        return 2;
    }

    K3HashConfig config = k3_config_default();
    uint32_t hash_resultado = 0;

    if (strcmp(argv[1], "--text") == 0) {
        const char* texto = argv[2];
        hash_resultado = k3_hash_buffer(texto, strlen(texto), &config);
        printf("K3(\"%s\") = 0x%08X\n", texto, hash_resultado);
    } else {
        if (k3_hash_file(argv[2], &config, &hash_resultado) != 0) {
            fprintf(stderr, "[ERROR] No se pudo leer el fichero: %s\n", argv[2]);
            return 1;
        }
        printf("K3(%s) = 0x%08X\n", argv[2], hash_resultado);
    }

    return 0;
}
