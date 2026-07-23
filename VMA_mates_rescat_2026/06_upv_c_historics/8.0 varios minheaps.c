#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> // Add this line
#include <math.h>

typedef struct {
    MinHeap estr; // Múltiplo primo representa cuanto falta para que el numero sea multiplo de el primo de pos
    long long prim; // Posición en la secuencia representa la posicion del primo
    int vacio; // Posición en la secuencia representa la posicion del primo
} Array;

// Estructura de un nodo del minheap
typedef struct {
    long long mult; // Múltiplo primo representa cuanto falta para que el numero sea multiplo de el primo de pos
    long long pos; // Posición en la secuencia representa la posicion del primo
} HeapNode;

// Estructura de minheap
typedef struct {
    HeapNode *arr; // Array de nodos del heap
    long long size; // Tamaño actual del heap
    long long capacity; // Capacidad máxima del heap
} MinHeap;

// Función para crear un nuevo minheap
MinHeap *createMinHeap(long long capacity) {
    MinHeap *minHeap = (MinHeap *)malloc(sizeof(MinHeap)); //reserva para nodo
    minHeap->arr = (HeapNode *)malloc(capacity * sizeof(HeapNode)); //revserva para array de nodos
    minHeap->size = 0; // vacia
    minHeap->capacity = capacity; // maximo
    return minHeap; // devolvemos el minheap
}

// Función para intercambiar dos elementos
void swap(HeapNode *a, HeapNode *b) {
    HeapNode temp = *a;
    *a = *b;
    *b = temp;
}

// Función para mantener la propiedad del minheap
void minHeapify(MinHeap *minHeap, long long i) {
    long long smallest = i;
    long long left = 2 * i + 1;
    long long right = 2 * i + 2;
    if (left < minHeap->size && minHeap->arr[left].mult < minHeap->arr[smallest].mult)
        smallest = left;
    if (right < minHeap->size && minHeap->arr[right].mult < minHeap->arr[smallest].mult)
        smallest = right;
    if (smallest != i) {
        swap(&minHeap->arr[i], &minHeap->arr[smallest]);
        minHeapify(minHeap, smallest);
    }
}

// Función para eliminar el mínimo del minheap
void eliminarMin(MinHeap *minHeap) {
    if (minHeap->size <= 0) {
        printf("El minheap está vacío.\n");
        return;
    }
    // Reemplazar el nodo raíz con el último nodo
    minHeap->arr[0] = minHeap->arr[minHeap->size - 1];
    minHeap->size--;
    // Reajustar el minheap
    minHeapify(minHeap, 0);
}

// Función para insertar un nodo en el minheap
void insert(MinHeap *minHeap, long long mult, long long pos) {
    if (minHeap->size == minHeap->capacity) {
        printf("El minheap está lleno.\n");
        return;
    }
    minHeap->size++;
    long long i = minHeap->size - 1;
    minHeap->arr[i].mult = mult;
    minHeap->arr[i].pos = pos;
    // Proceso de heapify up para mantener la propiedad del minheap
    while (i != 0 && minHeap->arr[(i - 1) / 2].mult > minHeap->arr[i].mult) {
        swap(&minHeap->arr[i], &minHeap->arr[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

// Función para extraer el mínimo del minheap
HeapNode extractMin(MinHeap *minHeap) {
    if (minHeap->size <= 0) {
        printf("El minheap está vacío.\n");
        HeapNode emptyNode = {0, 0}; // Nodo vacío
        return emptyNode;
    }
    // Extraer el nodo mínimo (que está en la raíz)
    HeapNode root = minHeap->arr[0];
    // Proceso de heapify down para mantener la propiedad del minheap
    return root;
}

int main() {
    long long n = 10000;
    // Inicialización del arreglo P6
    bool *P6 = (bool *)malloc((n + 1) * sizeof(bool));
    if (P6 == NULL) {
        printf("Error al asignar memoria dinámica en primos.\n");
        free(P6);
        return 1;
    }
    // marco a mano los primeros primos
    P6[2] = 1;
    P6[3] = 1;
    P6[5] = 1;
    P6[7] = 1;
    // Como los primos son de forma 6k+-1 marco todos los que no son así a no primo
    for (int i = 6; i <= n; i = i + 6) {
        if (i + 5 < n) {
            P6[5 + i] = 1;
        }
        if (i + 7 < n) {
            P6[7 + i] = 1;
        }
        if (i + 2 < n) {
            P6[2 + i] = 0;
        }
        if (i + 3 < n) {
            P6[3 + i] = 0;
        }
        if (i + 4 < n) {
            P6[4 + i] = 0;
        }
    }
    for (int i = 5; i <= n ^ (1 / 2); i = i + 4) {
        int salto = 2 * i;
        int salto2 = 4 * i;
        if (P6[i] == 1) {
            for (int j = i * i; j <= n; j += salto2) {
                if (j < n) {
                    P6[j] = 0;
                }
                j += salto;
                if (j < n) {
                    P6[j] = 0;
                }
            }
        }
        i += 2;
        if (P6[i] == 1) {
            for (int j = i * i; j <= n; j += salto) {
                if (j < n) {
                    P6[j] = 0;
                }
                j += salto2;
                if (j < n) {
                    P6[j] = 0;
                }
            }
        }
    }
    // se imprimen los primos calculados
    for (int i = 0; i <= n; i++) {
        printf(" %d es %d\n", i, P6[i]);
    }
    // Calcula el último número primo en el rango calculado este numero es muy importante. marca la raiz del rango que vamos a calcular ahora
    int c = 0;
    while (P6[n - c] == 0) {
        c++;
    }
    long long ultimo = (n - c);
    // Contar la cantidad de números primos en el rango y crear un arreglo Secuencia , esto esta doble pero no me apetece pensar, y eso que hay que optimizar
    // creo ya no esta doble
    long long cantidad = 0;
    for (long long cant = 5; cant <= n; cant++) {
        if (P6[cant] == 1) {
            cantidad++;
        }
    }
    // secuencia es una lista de primos hasta ultimo, para con pos saber a que primo nos referimos, la posicion de seciencia esta ligada a la posicion del minheap
    long long *Secuencia = (long long *)malloc(cantidad * sizeof(long long));
    if (Secuencia == NULL) {
        printf("Error al asignar memoria dinámica en Secuencia.\n");
        return 1;
    }
    // llenamos secuencia
    for (long long pos = 5; pos <= n; pos++) {
        if (P6[pos] == 1) {
            Secuencia[pos] = pos;
        }
    }
    long long posult = 0;
    // Crear minheap y llenarlo el resto de limite entre la multiplicacion de el primo y el resto de limite y primo
    // numMult nos dice desde limite cuanto falta para ser multiplo de primo [pos],
    long long tamanyo = cantidad / (cantidad - cantidad ^ (-100));
    long long j;
    long long Z;
    long long limite = ultimo * ultimo;
    Array *array = (Array *)malloc(tamanyo * sizeof(Array));
    if (array == NULL) {
        printf("Error al asignar memoria dinámica en array.\n");
        free(array);
        return 1;
    }
    while (posult < ultimo) {
        // crear un array para almacenar los minheaps
        MinHeap *minHeap4 = createMinHeap(tamanyo);
        Z = posult;
        for (long long i = posult; i - posult < tamanyo; i++) {
            long long mult = ultimo * ultimo % Secuencia[posult];
            if (mult % 2 == 1) {
                mult += Secuencia[i];
            }
            insert(minHeap4, mult, i);
            posult = i;
        }
        array[j].estr = minHeap4;
        array[j].prim = Z;
        array[j].vacio = 0;
        j++;
    }
    // Inicializar arreglo de primos y marcar los números primos
    bool *Primos = (bool *)malloc(limite * sizeof(bool));
    if (Primos == NULL) {
        printf("Error al asignar memoria dinámica en primos.\n");
        free(Primos);
        return 1;
    }
    // marcamos a mano los primeros primos
    Primos[2] = 1;
    Primos[3] = 1;
    Primos[5] = 1;
    Primos[7] = 1;
    // marcamos no primos los que no son 6k+-1
    for (int i = 6; i <= limite; i = i + 6) {
        if (5 + i < limite)
            Primos[5 + i] = 1;
        if (7 + i < limite)
            Primos[7 + i] = 1;
        if (2 + i < limite)
            Primos[2 + i] = 0;
        if (3 + i < limite)
            Primos[3 + i] = 0;
        if (4 + i < limite)
            Primos[4 + i] = 0;
    }
    j = 0;
    int vacioR = 0;
    while (vacioR <= tamanyo) {
        for (j = 0; j < tamanyo; j++) {
            if (array[j].vacio == 0) {
                HeapNode nodo = extractMin(array[j].estr); // obtengo minimo de principal
                if (nodo.mult == 0) {
                    vacioR++;
                    array[j].vacio = 1;
                } else { // si el minimo es 0 incremento el contador de vacios
                    long long pos;
                    long long primo;
                    long long mult;
                    mult = nodo.mult;
                    pos = nodo.pos;
                    do {
                        primo = limite - mult; // calculo  no primo
                        mult += Secuencia[pos] * 2;
                    } while (Primos[primo] == 0); // calculo el siguiente multiplo
                    Primos[primo] = 0; //guardo en memoria
                    eliminarMin(array[j].estr); // borro minimo de principal
                    if (mult <= limite - ultimo) {
                        insert(array[j].estr, mult, pos); //si procede inserto en el auxiliar
                    }
                }
            }
        }
    }

    //imprimimos los primos calculados
    for (long long i = ultimo; i < limite; i++) {
        printf(" %llu es %d\n", i, Primos[i]);
    }
    scanf("%d", &n);
    // Liberar memoria
    free(P6);
    free(Secuencia);
    free(Primos);
    for (long long i = 0; i < tamanyo; i++) {
        free(array[i].estr.arr);
    }
    free(array);

    return 0;
}
