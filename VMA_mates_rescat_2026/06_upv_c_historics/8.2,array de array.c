#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> // Add this line
#include <math.h>

typedef struct {
    Est estr; // Múltiplo primo representa cuanto falta para que el numero sea multiplo de el primo de pos
    long long prim;
    long long fallo;
    int vacio;
   
} Array;

typedef struct {
    
     long long mult; // Múltiplo primo representa cuanto falta para que el numero sea multiplo de el primo de pos
    long long pos; // Posición en la secuencia representa la posicion del primo
    int vacio;
} Est;



// Función para crear un nuevo minheap
Array *create(long long capacity,long long ultimo) {
    Array *array = (Array *)malloc(capacity * sizeof(Array)); //reserva para nodo
    if (array == NULL) {
        printf("Error al asignar memoria dinámica en array.\n");
        free(array);
        return NULL;
    }   
    array->estr = (Est *)malloc((ultimo+1) * sizeof(Est)); //revserva para array de nodos
    if (array->estr == NULL) {
        printf("Error al asignar memoria dinámica en array->estr.\n");
        free(array->estr);
        return NULL;
    }
 
    return array; // devolvemos el minheap
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
    long long j=0;
    long long Z;
    long long limite = ultimo * ultimo;
Array *array = create(tamanyo,cantidad);
    while (posult < cantidad) {
        Z = posult;
        for (long long i = 0; i + posult < tamanyo; i++) {
            long long mult = ultimo * ultimo % Secuencia[posult];
            if (mult % 2 == 1) {
                mult += Secuencia[i];
            }
            array[j].estr[i].mult= mult;
            array[j].estr[i].pos = i;
            posult = i;
        }
       
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
    int i;
    int vacioR = 0;
    while (vacioR < ultimo) {
        
        for (j = 0; j < tamanyo; j++) {
            
            
            if (array[i].estr[j].vacio == 0) {
                if (array[i].estr[j].mult == 0) {
                    vacioR++;
                    array[i].estr[j].vacio = 1;
                } else { // si el minimo es 0 incremento el contador de vacios
                    long long pos;
                    long long primo;
                    long long mult;
                    mult = array[i].estr[j].mult;
                    pos = array[i].estr[j].pos;
                    do {
                        primo = limite - mult; // calculo  no primo
                        mult += Secuencia[pos] * 2;
                    } while (Primos[primo] == 0); // calculo el siguiente multiplo
                    Primos[primo] = 0; //guardo en memoria
                   
                    if (mult <= limite - ultimo) {
                        array[i].estr[j].mult = mult;                        
                    }else{
                        array[i].estr[j]=array[i].estr[j-array[i].fallo-1];
                        array[i].estr[j-array[i].fallo-1].mult=0;
                        array[i].fallo--;
                    }
                        i++;
                    if(i>cantidad){
                        i=0;
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
