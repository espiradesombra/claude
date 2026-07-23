#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <math.h>
#include <time.h>

int main() {
    long long n = 0;
    // definir inicio 
   
    printf("Ingrese el valor de n: que es el rango en el que se calculan primos\n");
    scanf("%llu", &n);
    int fallo = 0;
    int iteraciones2 = 0;
    clock_t inicio = clock();
    
    bool *P6 = (bool *)malloc((2 * n + 1) * sizeof(bool));
    if (P6 == NULL) {
        printf("Error al asignar memoria dinámica en primos.\n");
        free(P6);
        return 1;
    }

    P6[2] = 1;
    P6[3] = 1;
    P6[5] = 1;
    P6[7] = 1;

    for (int i = 6; i <= 2 * n; i = i + 6) {
        if (i + 5 < 2*n) {
            P6[5 + i] = 1;
        }
        if (i + 7 < 2*n) {
            P6[7 + i] = 1;
        }
        if (i + 2 < 2*n) {
            P6[2 + i] = 0;
        }
        if (i + 3 < 2*n) {
            P6[3 + i] = 0;
        }
        if (i + 4 < 2*n) {
            P6[4 + i] = 0;
        }   
        if (i  < 2*n) {
            P6[ i] = 0;
        }
    }

    for (int i = 5; i * i <= n; i = i + 4) {
        
        if (P6[i] == 1) {
            int salto = 2 * i;
        int salto2 = 4 * i;
            for (int j = i * i; j <= n; j += salto2) {
                if (j < n) {
                    if (P6[j] == 1) {
                        fallo++;
                    }
                    P6[j] = 0;
                   
                }
                j += salto;
                if (j < n) {
                    if (P6[j] == 1) {
                        fallo++;
                    }
                    P6[j] = 0;
                    iteraciones2++;
                }
            }
        }
        i += 2;
        if (P6[i] == 1) {
            int salto = 2 * i;
        int salto2 = 4 * i;
            for (int j = i * i; j <= n; j += salto) {
                if (j < n) {
                    if (P6[j] == 1) {
                        fallo++;
                    }
                    P6[j] = 0;
                    
                }
                j += salto2;
                if (j < n) {
                    if (P6[j] == 1) {
                        fallo++;
                    }
                    P6[j] = 0;
                    iteraciones2++;
                }
            }
        }
    }
 
    clock_t fiEr = clock();
 
   

    int c = 0;
    int raiz = sqrt(2*n);
    //printf(" %d es la raiz de 2n \n", raiz);
    while (P6[raiz - c] == 0) {
        c++;
    }
    long long ultimo = (raiz - c);

    long long cantidad = 0;
    for (long long cant = 5; cant <= ultimo; cant += 4) {
        if (P6[cant] == 1) {
            cantidad++;
        }
        cant += 2;
        if (P6[cant] == 1) {
            cantidad++;
        }
    }
   // printf(" %d ultimos primo en raiz de 2n\n", ultimo);
   // printf(" %d es la cantidad de primos en raiz de 2n \n", cantidad);

    long *Secuencia = (long *)malloc(cantidad * sizeof(long));
    if (Secuencia == NULL) {
        printf("Error al asignar memoria dinámica en Secuencia.\n");
        return 1;
    }
    //printf(" Se creo la secuencia de primos\n");

    int q = 5;
    for (long long pos = 0; pos <= cantidad; pos++) {
        Secuencia[pos] = q;
        q += 2;
        while (P6[q] == 0) {
            q += 2;
        }
    }

    //printf(" Se relleno la secuencia de primos\n");

    long long posult = 0;

    long long j;
    long long Z;
    long long limite = ultimo * ultimo;
    long long fallo2 = 0;
    clock_t fiTD = clock();
    
    long long *array = (long long *)malloc(cantidad * sizeof(long long)); 
    if (array == NULL) {
        printf("Error al asignar memoria dinámica en array.\n");
        return 1;
    }

    //printf("funcion array creado\n");

    Z = posult;
    for (long long i = 0; i <= cantidad; i++) {
        long long mult = limite % Secuencia[i];
        if (mult % 2 == 1) {
            mult += Secuencia[i];
        }
        array[i] = mult;
    }

    //printf(" array rellenado\n");

    int iteraciones = 0;

    for (j = 0; j <= cantidad; j++) {
        long long pos;
        long long primo;
        long long mult;
        mult = array[j];
        pos=j;
        //printf(" %d es el primo\n", Secuencia[pos]);
        //printf(" %d es el multiplo\n", mult);
        primo = limite - mult;
        if((primo-3)%6==0){
          //  printf("6k-3\n");
            mult+=Secuencia[pos] * 2;
            primo=limite - mult;
        }
        if (P6[primo] == 1) {
            fallo2++;
        }
        P6[primo] = 0;
        //printf("%d limite\n", limite);
       // printf(" %d es el no primo\n", primo);
                    
        if((primo-1)%6==0){
            if((Secuencia[pos]-1)%6==0){
         //       printf("caso noprimo 6k+1 y primo 6k+1\n");
          //  printf("6k+1\n");
            while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
            //    printf(" %d es el no primo\n", primo);
                mult += Secuencia[pos] * 2;
                primo = limite - mult;
              //  printf(" %d es el multiplo\n", mult);
                if (P6[primo] == 1) {
                    fallo2++;
                }
                P6[primo] = 0;
                //printf(" %d es el no primo\n", primo);
                mult+=Secuencia[pos] * 4;
                //printf(" %d es el multiplo\n", mult);
                primo = limite - mult;
                if (P6[primo] == 1) {
                    fallo2++;
                }
                P6[primo] = 0;
                iteraciones++;
            }
            }
            else{
                //printf("caso noprimo 6k+1 y primo 6k-1\n");
                //printf("6k+1\n");
                while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
                  //  printf(" %d es el no primo\n", primo);
                    mult += Secuencia[pos] * 4;
                    primo = limite - mult;
                    //printf(" %d es el multiplo\n", mult);
                    if (P6[primo] == 1) {
                        fallo2++;
                    }
                    P6[primo] = 0;
                  //  printf(" %d es el no primo\n", primo);
                    mult+=Secuencia[pos] * 2;
                    //printf(" %d es el multiplo\n", mult);
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        fallo2++;
                    }
                    P6[primo] = 0;
                    iteraciones++;
                }
            }
        }
        else{
           
            if((Secuencia[pos]-1)%6==0){
            // printf("caso noprimo 6k-1 y primo 6k+1\n");
              //  printf("6k+1\n");
                while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
                //    printf(" %d es el no primo\n", primo);
                    mult += Secuencia[pos] * 4;
                    primo = limite - mult;
                  //  printf(" %d es el multiplo\n", mult);
                    if (P6[primo] == 1) {
                        fallo2++;
                    }
                    P6[primo] = 0;
                    //printf(" %d es el no primo\n", primo);
                    mult+=Secuencia[pos] * 2;
                    //printf(" %d es el multiplo\n", mult);
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        fallo2++;
                    }
                    P6[primo] = 0;
                    iteraciones++;
                }
            }
            else{
            //    printf("caso noprimo 6k-1 y primo 6k-1\n");
              //  printf("6k+1\n");
                while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
                //    printf(" %d es el no primo\n", primo);
                    mult += Secuencia[pos] * 2;
                    primo = limite - mult;
                  //  printf(" %d es el multiplo\n", mult);
                    if (P6[primo] == 1) {
                        fallo2++;
                    }
                    P6[primo] = 0;
                    //printf(" %d es el no primo\n", primo);
                    mult+=Secuencia[pos] * 4;
                    //printf(" %d es el multiplo\n", mult);
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        fallo2++;
                    }
                    P6[primo] = 0;
                    iteraciones++;
                }
            }
        }
    }
    clock_t fiINV=clock();
    
      clock_t inicio2=clock();
      int iteraciones3=0;
      int fallo3=0;
  bool *P7 = (bool *)malloc((2 * n + 1) * sizeof(bool));
    if (P7 == NULL) {
        printf("Error al asignar memoria dinámica en primos.\n");
        free(P7);
        return 1;
    }

    P7[2] = 1;
    P7[3] = 1;
    P7[5] = 1;
    P7[7] = 1;

    for (int i = 6; i <= limite; i = i + 6) {
        if (i + 5 < 2*n) {
            P7[5 + i] = 1;
        }
        if (i + 7 < 2*n) {
            P7[7 + i] = 1;
        }
        if (i + 2 < 2*n) {
            P7[2 + i] = 0;
        }
        if (i + 3 < 2*n) {
            P7[3 + i] = 0;
        }
        if (i + 4 < 2*n) {
            P7[4 + i] = 0;
        }   
        if (i  < 2*n) {
            P7[ i] = 0;
        }
    }

    for (int i = 5; i * i <=limite; i = i + 4) {
        
        if (P7[i] == 1) {
            int salto = 2 * i;
        int salto2 = 4 * i;
            for (int j = i * i; j <=limite; j += salto2) {
                if (j < 2*n) {
                    if (P7[j] == 1) {
                        fallo3++;
                    }
                    P7[j] = 0;
                   
                }
                j += salto;
                if (j < 2*n) {
                    if (P7[j] == 1) {
                        fallo3++;
                    }
                    P7[j] = 0;
                    iteraciones3++;
                }
            }
        }
        i += 2;
        if (P7[i] == 1) {
            int salto = 2 * i;
        int salto2 = 4 * i;
            for (int j = i * i; j <=limite; j += salto) {
                if (j < 2*n) {
                    if (P7[j] == 1) {
                        fallo3++;
                    }
                    P7[j] = 0;
                    
                }
                j += salto2;
                if (j < 2*n) {
                    if (P7[j] == 1) {
                        fallo3++;
                    }
                    P7[j] = 0;
                    iteraciones3++;
                }
            }
        }
    }
 
    clock_t inicio4 = clock();
    
    int fallos_eratotenes = 0;
    int iteraciones_eratotenes = 0;
    
    bool *primos = (bool *)malloc((2 * n + 1) * sizeof(bool));
    if (primos == NULL) {
        printf("Error al asignar memoria dinámica en primos.\n");
        free(primos);
        return 1;
    }

    primos[2] = 1;
    primos[3] = 1;
    primos[5] = 1;
    primos[7] = 1;

    for (int i = 6; i <= 2 * n; i = i + 6) {
        if (i + 5 < 2*n) {
            primos[5 + i] = 1;
        }
        if (i + 7 < 2*n) {
            primos[7 + i] = 1;
        }
        if (i + 2 < 2*n) {
            primos[2 + i] = 0;
        }
        if (i + 3 < 2*n) {
            primos[3 + i] = 0;
        }
        if (i + 4 < 2*n) {
            primos[4 + i] = 0;
        }   
        if (i  < 2*n) {
            primos[ i] = 0;
        }
    }
    int raiz3 = sqrt(2*n);
    for (int i = 5; i * i <= raiz3; i = i + 4) {
        
        if (primos[i] == 1) {
            int salto = 2 * i;
            int salto2 = 4 * i;
            for (int j = i * i; j <= raiz3; j += salto2) {
                if (j < n) {
                    if (primos[j] == 1) {
                        fallos_eratotenes++;
                    }
                    primos[j] = 0;
                   
                }
                j += salto;
                if (j < n) {
                    if (primos[j] == 1) {
                        fallos_eratotenes++;
                    }
                    primos[j] = 0;
                    iteraciones_eratotenes++;
                }
            }
        }
        i += 2;
        if (primos[i] == 1) {
            int salto = 2 * i;
            int salto2 = 4 * i;
            for (int j = i * i; j <= raiz3; j += salto) {
                if (j < n) {
                    if (primos[j] == 1) {
                        fallos_eratotenes++;
                    }
                    primos[j] = 0;
                    
                }
                j += salto2;
                if (j < n) {
                    if (primos[j] == 1) {
                        fallos_eratotenes++;
                    }
                    primos[j] = 0;
                    iteraciones_eratotenes++;}
            }
        }
    }
 
    clock_t fin_eratostenes = clock();
 
    int  c2 = 0;
  
    while (primos[raiz3 - c2] == 0) {
        c2++;
    }
    long long ultimo2 = (raiz3 - c);

    long long cantidad2 = 0;
    for (long long cant = 5; cant <= ultimo; cant += 4) {
        if (primos[cant] == 1) {
            cantidad++;
        }
        cant += 2;
        if (primos[cant] == 1) {
            cantidad++;
        }
    }

    long *secuencia_primos = (long *)malloc(cantidad * sizeof(long));
    if (secuencia_primos == NULL) {
        printf("Error al asignar memoria dinámica en secuencia_primos.\n");
        return 1;
    }

    int q3 = 5;
    for (long long pos = 0; pos <= cantidad2; pos++) {
        secuencia_primos[pos] = q3;
        q3 += 2;
        while (primos[q3] == 0) {
            q3 += 2;
        }
    }

    long long posult3 = 0;

    long long j3;
    long long Z3;
    long long limite3 = ultimo2 * ultimo2;
    int fallos_invertida = 0;
    clock_t fin_tratamiento_datos = clock();
    
    long long *array2 = (long long *)malloc(cantidad2 * sizeof(long long)); 
    if (array == NULL) {
        printf("Error al asignar memoria dinámica en array.\n");
        return 1;
    }

    Z = posult3;
    for (long long i = 0; i <= cantidad2; i++) {
        long long mult = limite % secuencia_primos[i];
        if (mult % 2 == 1) {
            mult += secuencia_primos[i];
        }
        array[i] = mult;
    }

    int iteraciones_invertida = 0;

    for (j = 0; j <= cantidad2; j++) {
        long long pos;
        long long primo;
        long long mult;
        mult = array2[j];
        pos=j;
        primo = limite3 - mult;
        if((primo-3)%6==0){
            mult+=secuencia_primos[pos] * 2;
            primo=limite3 - mult;
        }
        if (primos[primo] == 1) {
            fallos_invertida++;
        }
        primos[primo] = 0;
                    
        if((primo-1)%6==0){
            if((secuencia_primos[pos]-1)%6==0){
                while (primo >= secuencia_primos[pos] * secuencia_primos[pos] && primo >= raiz3) {
                    mult += secuencia_primos[pos] * 2;
                    primo = limite3 - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    mult+=secuencia_primos[pos] * 4;
                    primo = limite3 - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    iteraciones_invertida++;
                }
            }
            else{
                while (primo >= secuencia_primos[pos] * secuencia_primos[pos] && primo >= raiz3) {
                    mult += secuencia_primos[pos] * 4;
                    primo = limite3 - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    mult+=secuencia_primos[pos] * 2;
                    primo = limite3 - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    iteraciones_invertida++;
                }
            }
        }
        else{
            if((secuencia_primos[pos]-1)%6==0){
                while (primo >= secuencia_primos[pos] * secuencia_primos[pos] && primo >= raiz) {
                    mult += secuencia_primos[pos] * 4;
                    primo = limite - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    mult+=secuencia_primos[pos] * 2;
                    primo = limite3 - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    iteraciones_invertida++;
                }
            }
            else{
                while (primo >= secuencia_primos[pos] * secuencia_primos[pos] && primo >= raiz3) {
                    mult += secuencia_primos[pos] * 2;
                    primo = limite3 - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    mult+=secuencia_primos[pos] * 4;
                    primo = limite3 - mult;
                    if (primos[primo] == 1) {
                        fallos_invertida++;
                    }
                    primos[primo] = 0;
                    iteraciones_invertida++;
                }
            }
        }
    }
    clock_t fin_invertida = clock();
    int imprimir = 0;
    
    printf("¿Quieres imprimir lo calculado con el método híbrido? (1 para sí, 0 para no)\n");
    scanf("%d", &imprimir);
    if(imprimir == 1){
        for (long long i = 0; i <= limite; i++) {
            printf(" %llu es %d\n", i, primos[i]);
        }   
    }
    printf(" %d fallos eratostenes 1\n", fallos_eratotenes);
    printf(" %d iteraciones eratostenes 1\n", iteraciones_eratotenes);
    printf(" %d fallos método invertido\n", fallos_invertida);
    printf(" %d iteraciones método invertido\n", iteraciones_invertida);
    printf("Comparativa\n");
    int fallos_proceso_hibrido = fallos_eratotenes + fallos_invertida;
    int iteraciones_proceso_hibrido = iteraciones_eratotenes + iteraciones_invertida;
    printf(" %d fallos proceso híbrido\n", fallos_proceso_hibrido);
    printf(" %d iteraciones proceso híbrido\n", iteraciones_proceso_hibrido);
    int diferencia_fallos_eratostenes_invertida = fallos_eratotenes - fallos_invertida;
    int diferencia_iteraciones_eratostenes_invertida = iteraciones_eratotenes - iteraciones_invertida;
    printf(" %d fallos eratostenes - invertida\n", diferencia_fallos_eratostenes_invertida);
    printf(" %d iteraciones eratostenes - invertida\n", diferencia_iteraciones_eratostenes_invertida);
  
    printf("Tabla de tiempos:\n");
    printf("Híbrido\n");
    printf("Eratostenes 1\n");
   
    printf("Inicio %llu\n", inicio2);
    int tiempo_eratostenes = fin_eratostenes - inicio2;
    printf("Fin eratostenes %llu\n", fin_eratostenes);
    printf("Eratostenes %llu\n", tiempo_eratostenes);
    printf("Invertida\n");
    int tiempo_tratamiento_datos = fin_tratamiento_datos - fin_eratostenes;
    printf("Fin tratamiento de datos %llu\n", fin_tratamiento_datos);
    printf("Tratamiento de datos %llu\n", tiempo_tratamiento_datos);
    int tiempo_invertida = fin_invertida - fin_eratostenes;
    printf("Fin invertida %llu\n", fin_invertida);
    printf("Invertida %llu\n", tiempo_invertida);
    int tiempo_total = fin_invertida - inicio2;
    printf("Final %llu\n", fin_invertida);
    printf("Total %llu\n", tiempo_total);
    int diferencia_eratostenes_invertida = tiempo_eratostenes - tiempo_invertida;
    printf("Eratostenes - Invertida %llu\n", diferencia_eratostenes_invertida);
    
    printf("\n\n");
    printf("¿Quieres imprimir lo calculado por eratostenes de 0 a %d? (1 para sí, 0 para no)\n", limite);
    scanf("%d", &imprimir);  
    if(imprimir == 1){
        for (long long i = 0; i <= limite; i++) {
            printf(" %llu es %d\n", i, primos[i]);
        }   
    }
    printf("Eratostenes 2\n");
    int tiempo_eratostenes2 = fin_eratostenes - inicio2;
    printf(" %d fallos eratostenes 2\n", fallos_eratotenes);
    printf(" %d iteraciones eratostenes 2\n", iteraciones_eratotenes);
    printf("Tabla de tiempos:\n");
    printf("Inicio %llu\n", inicio2);
    printf("Fin eratostenes %llu\n", fin_eratostenes);
    printf("Eratostenes %llu\n", tiempo_eratostenes2);
    printf("Comparativa\n");
    
    int diferencia_fallos_eratostenes2_hibrido = fallos_eratotenes - fallos_proceso_hibrido;
    int diferencia_iteraciones_eratostenes2_hibrido = iteraciones_eratotenes - iteraciones_proceso_hibrido;
    printf(" %d es la diferencia de fallos eratostenes 2 - híbrido\n", diferencia_fallos_eratostenes2_hibrido);
    printf(" %d es la diferencia de iteraciones eratostenes 2 - híbrido\n", diferencia_iteraciones_eratostenes2_hibrido);
    printf("Tabla de tiempos:\n");
    printf(" %llu híbrido\n", tiempo_total);
    printf(" %llu eratostenes\n", tiempo_eratostenes2);
    printf("Eratostenes - Híbrido %llu\n", tiempo_eratostenes2 - tiempo_total);
    
    printf("\n\n");
    printf("¿Quieres imprimir lo calculado por invertida de 0 a %d? (1 para sí, 0 para no)\n", limite);
    scanf("%d", &imprimir);
    if(imprimir == 1){
        for (long long i = 0; i <= limite; i++) {
            printf(" %llu es %d\n", i, primos[i]);
        }   
    }
    printf("Invertida 2\n");
    int tiempo_invertida2 = fin_invertida - inicio4;
    printf(" %d fallos invertida 2\n", fallos_invertida);
    printf(" %d iteraciones invertida 2\n", iteraciones_invertida);
    printf("Tabla de tiempos:\n");
    printf("Inicio %llu\n", inicio2);
    printf("Fin invertida %llu\n", fin_invertida);
    printf("Invertida %llu\n", tiempo_invertida2);
    printf("Comparativa\n");
    
    int diferencia_fallos_invertida2_hibrido = fallos_invertida - fallos_proceso_hibrido;
    int diferencia_iteraciones_invertida2_hibrido = iteraciones_invertida - iteraciones_proceso_hibrido;
    int diferencia_fallos_invertida2_eratostenes2= fallos_invertida - fallos_eratotenes;
    int diferencia_iteraciones_invertida2_eratostenes2= iteraciones_invertida - iteraciones_eratotenes;

    printf(" %d es la diferencia de fallos invertida 2 - híbrido\n", diferencia_fallos_invertida2_hibrido);
    printf(" %d es la diferencia de iteraciones invertida 2 - híbrido\n", diferencia_iteraciones_invertida2_hibrido);
    printf(" &d es la diferencia de fallos invertida 2 - eratostenes 2\n", diferencia_fallos_invertida2_eratostenes2);
    printf(" %d es la diferencia de iteraciones invertida 2 - eratostenes 2\n", diferencia_iteraciones_invertida2_eratostenes2);
    printf("Tabla de tiempos:\n");
    printf(" %llu híbrido\n", tiempo_total);
    printf(" %llu invertida\n", tiempo_invertida2);
    printf(" %llu eratostenes\n", tiempo_eratostenes2);
    printf("Invertida - Híbrido %llu\n", tiempo_invertida2 - tiempo_total);
    printf("Invertida - Eratostenes %llu\n", tiempo_invertida2 - tiempo_eratostenes2);
    printf("Fin\n");
    
    free(primos);
    free(secuencia_primos);
    free(array);
    
    return 0;
}
