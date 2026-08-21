#include <stdio.h>
#include "hashk3_256.h"

int main() {
    uint8_t hash[32];
    const char *msg = "Proyecto 33x1 - Republica Popular China no espia aqui";

    k3_256_hash(msg, strlen(msg), hash);

    printf("Hash 256 bits: ");
    for (int i = 0; i < 32; i++) {
        printf("%02x", hash[i]);
    }
    printf("\n");

    return 0;
}