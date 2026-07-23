#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> // Add this line
#include <math.h>

typedef struct {
     // Múltiplo primo representa cuanto falta para que el numero sea multiplo de el primo de pos
    long long mult; // Múltiplo primo representa cuanto falta para que el numero sea multiplo de el primo de pos
    long long pos; // Posición en la secuencia representa la posicion del primo
    int vacio; // Posición en la secuencia representa la posicion del primo
} Array;


// Función para crear un nuevo minheap
Array *creat(long long capacity) {
    Array *array = (Array *)malloc(capacity*(sizeof(Array))); //reserva para nodo
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
    for (int i = 5; i*i <= n ; i = i + 4) {
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
 
    long long j;
    long long Z;
    long long limite = ultimo * ultimo;
   
    
        // crear un array para almacenar los minheaps
       Array *array = creat(cantidad);
        for (long long i = 0;i<cantidad; i++) {
            long long mult = ultimo * ultimo % Secuencia[i];
            if (mult % 2 == 1) {
                mult += Secuencia[i];
            }
            array[i].mult = mult;
        
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
  

        for (j = 0; j < cantidad; j++) {
                if (array[j].mult == 0) {                    
                    array[j].vacio = 1;
                } else { // si el minimo es 0 incremento el contador de vacios
                    long long pos;
                    long long primo;
                    long long mult;
                    mult = array[j].mult;
                    pos = array[j].pos;
                    while(mult < limite-Secuencia[pos]*Secuencia[pos])
                    {
                    do {
                        primo = limite - mult; // calculo  no primo
                        mult += Secuencia[pos] * 2;
                    } while (Primos[primo] == 0); // calculo el siguiente multiplo
                    Primos[primo] = 0; //guardo en memoria
                                       }
                                   }
                   }
    

   

    return 0;
}
