#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    long long x = 1;
    long long n = 5; // Asumiendo un valor para n
    long long i, j;
    
    for (i = 2; i <= n; i++) {
        x = i * x;
    }
    
    double raiz = pow(x, 1.0/2.0);
    double cantidad = raiz / log(raiz);
    long long salto = 0;
    
    for (i = 1; i <= cantidad; i++) {
        salto += 3;
    }
    
    salto += 1;
    
    long long size = - n + salto;
    int *is_prime = (int *)malloc(size * sizeof(int));
    if (is_prime == NULL) {
        printf("Error al asignar memoria\n");
        return 1;
    }
    
    // Inicializar el array con 0 (asumiendo que todos son primos al inicio)
    for (i = 0; i < size; i++) {
        is_prime[i] = 0;
    }
    
    // Marcar los múltiplos como no primos
    for (i = x + 1; i >= x - n + salto; i += 2) {
        for (j = i; j < x - n + salto; j += 2 * i) {
            is_prime[j-x] = 1;
        }
    }
    
    // Identificar y imprimir los números primos
    for (i = 0; i < size; i++) {
            printf("%lld %d es primo\n", i,is_prime[i]);
        
    }
    
    free(is_prime);
    return 0;
}
