#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <math.h>
#include <time.h>

bool* allocateMemory(int size) {
    bool* array = (bool*)malloc(size * sizeof(bool));
    if (array == NULL) {
        printf("Error al asignar memoria dinámica.\n");
        return NULL;
    }
    return array;
}

void initializePrimes(bool* array, int size) {
    array[2] = true;
    array[3] = true;
    array[5] = true;
    array[7] = true;

    for (int i = 6; i <= size; i = i + 6) {
        if (i + 5 < size) {
            array[5 + i] = true;
        }
        if (i + 7 < size) {
            array[7 + i] = true;
        }
        if (i + 2 < size) {
            array[2 + i] = false;
        }
        if (i + 3 < size) {
            array[3 + i] = false;
        }
        if (i + 4 < size) {
            array[4 + i] = false;
        }
        if (i < size) {
            array[i] = false;
        }
    }
}

void sieveOfEratosthenes(bool* array, int size) {
    for (int i = 5; i * i <= size; i = i + 4) {
        if (array[i] == true) {
            int jump = 2 * i;
            int jump2 = 4 * i;
            for (int j = i * i; j <= size; j += jump2) {
                if (j < size) {
                    if (array[j] == true) {
                        array[j] = false;
                    }
                }
                j += jump;
                if (j < size) {
                    if (array[j] == true) {
                        array[j] = false;
                    }
                }
            }
        }
        i += 2;
        if (array[i] == true) {
            int jump = 2 * i;
            int jump2 = 4 * i;
            for (int j = i * i; j <= size; j += jump) {
                if (j < size) {
                    if (array[j] == true) {
                        array[j] = false;
                    }
                }
                j += jump2;
                if (j < size) {
                    if (array[j] == true) {
                        array[j] = false;
                    }
                }
            }
        }
    }
}

void generatePrimes(bool* array, int size, long long* sequence, int* count) {
    int q = 5;
    for (long long pos = 0; pos <= size; pos++) {
        sequence[pos] = q;
        q += 2;
        while (array[q] == false) {
            q += 2;
        }
    }

    *count = 0;
    for (long long cant = 5; cant <= size; cant += 4) {
        if (array[cant] == true) {
            (*count)++;
        }
        cant += 2;
        if (array[cant] == true) {
            (*count)++;
        }
    }
}

void calculateMultiples(long long* array, int size, long long limit, long long* sequence, int count) {
    for (long long i = 0; i <= count; i++) {
        long long mult = limit % sequence[i];
        if (mult % 2 == 1) {
            mult += sequence[i];
        }
        array[i] = mult;
    }
}

void processPrimes(bool* array, int size, long long* multiples, long long limit, long long* sequence, int count, int* failCount, int* iterationCount) {
    for (long long j = 0; j <= count; j++) {
        long long pos;
        long long prime;
        long long mult;
        mult = multiples[j];
        pos = j;
        prime = limit - mult;
        if ((prime - 3) % 6 == 0) {
            mult += sequence[pos] * 2;
            prime = limit - mult;
        }
        if (array[prime] == true) {
            (*failCount)++;
        }
        array[prime] = false;

        if ((prime - 1) % 6 == 0) {
            if ((sequence[pos] - 1) % 6 == 0) {
                while (prime >= sequence[pos] * sequence[pos] && prime >= size) {
                    mult += sequence[pos] * 2;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    mult += sequence[pos] * 4;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    (*iterationCount)++;
                }
            }
            else {
                while (prime >= sequence[pos] * sequence[pos] && prime >= size) {
                    mult += sequence[pos] * 4;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    mult += sequence[pos] * 2;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    (*iterationCount)++;
                }
            }
        }
        else {
            if ((sequence[pos] - 1) % 6 == 0) {
                while (prime >= sequence[pos] * sequence[pos] && prime >= size) {
                    mult += sequence[pos] * 4;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    mult += sequence[pos] * 2;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    (*iterationCount)++;
                }
            }
            else {
                while (prime >= sequence[pos] * sequence[pos] && prime >= size) {
                    mult += sequence[pos] * 2;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    mult += sequence[pos] * 4;
                    prime = limit - mult;
                    if (array[prime] == true) {
                        (*failCount)++;
                    }
                    array[prime] = false;
                    (*iterationCount)++;
                }
            }
        }
    }
}

void printPrimes(bool* array, int size) {
    for (long long i = 0; i <= size; i++) {
        printf(" %llu es %d\n", i, array[i]);
    }
}

void printTimeTable(clock_t start, clock_t end, const char* label) {
    int time = end - start;
    printf("%s\n", label);
    printf("Inicio: %llu\n", start);
    printf("Fin: %llu\n", end);
    printf("Tiempo: %llu\n", time);
}

int main() {
    long long n = 0;
    printf("Ingrese el valor de n: que es el rango en el que se calculan primos\n");
    scanf("%llu", &n);

    clock_t inicio = clock();

    bool* P6 = allocateMemory(2 * n + 1);
    if (P6 == NULL) {
        return 1;
    }

    initializePrimes(P6, 2 * n);

    sieveOfEratosthenes(P6, 2 * n);

    long long* Secuencia;
    int count;
    generatePrimes(P6, 2 * n, Secuencia, &count);

    long long* array = (long long*)malloc(count * sizeof(long long));
    if (array == NULL) {
        printf("Error al asignar memoria dinámica en array.\n");
        return 1;
    }

    long long limit = Secuencia[count - 1] * Secuencia[count - 1];
    calculateMultiples(array, count, limit, Secuencia, count);

    int failCount = 0;
    int iterationCount = 0;
    processPrimes(P6, 2 * n, array, limit, Secuencia, count, &failCount, &iterationCount);

    clock_t fiEr = clock();

    clock_t inicio2 = clock();

    bool* P7 = allocateMemory(2 * n + 1);
    if (P7 == NULL) {
        return 1;
    }

    initializePrimes(P7, 2 * n);

    sieveOfEratosthenes(P7, 2 * n);

    clock_t fiEr2 = clock();

    clock_t total = fiEr - inicio;
    clock_t erastotes = fiEr2 - inicio2;
    clock_t diferenciP1P2 = total - erastotes;

    printf("Tabla de tiempos:\n");
    printf("Inicio híbrido: %llu\n", inicio);
    printf("Fin eratóstenes híbrido: %llu\n", fiEr);
    printf("Tiempo eratóstenes híbrido: %llu\n", total);
    printf("Inicio eratóstenes: %llu\n", inicio2);
    printf("Fin eratóstenes: %llu\n", fiEr2);
    printf("Tiempo eratóstenes: %llu\n", erastotes);
    printf("Diferencia de tiempo entre eratóstenes y híbrido: %llu\n", diferenciP1P2);

    free(P6);
    free(Secuencia);
    free(array);
    free(P7);

    return 0;
}
