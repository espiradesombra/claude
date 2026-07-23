#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <math.h>
#include <time.h>
#include <omp.h>

int main() {
    long long n = 0;
   
    printf("Ingrese el valor de n: que es el rango en el que se calculan primos\n");
    scanf("%llu", &n);
    int fallo = 0;
    int iteraciones2 = 0;
    clock_t inicio = clock();
    
    bool *P6 = (bool *)malloc((2 * n + 1) * sizeof(bool));
    if (P6 == NULL) {
        printf("Error al asignar memoria dinámica en primos P6.\n");
        free(P6);
        return 1;
    }

    P6[2] = 1;
    P6[3] = 1;
    P6[5] = 1;
    P6[7] = 1;

    #pragma omp parallel for
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

    #pragma omp parallel for
    for (int i = 5; i * i <= n; i = i + 4) {
        
        if (P6[i] == 1) {
            int salto = 2 * i;
            int salto2 = 4 * i;
            for (int j = i * i; j <= n; j += salto2) {
                if (j < n) {
                    if (P6[j] == 1) {
                        #pragma omp atomic
                        fallo++;
                    }
                    P6[j] = 0;
                   
                }
                j += salto;
                if (j < n) {
                    if (P6[j] == 1) {
                        #pragma omp atomic
                        fallo++;
                    }
                    P6[j] = 0;
                    #pragma omp atomic
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
                        #pragma omp atomic
                        fallo++;
                    }
                    P6[j] = 0;
                    
                }
                j += salto2;
                if (j < n) {
                    if (P6[j] == 1) {
                        #pragma omp atomic
                        fallo++;
                    }
                    P6[j] = 0;
                    #pragma omp atomic
                    iteraciones2++;
                }
            }
        }
    }
 
    clock_t fiEr = clock();
 
    int c = 0;
    int raiz = sqrt(2*n);
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

    long *Secuencia = (long *)malloc(cantidad * sizeof(long));
    if (Secuencia == NULL) {
        printf("Error al asignar memoria dinámica en Secuencia.\n");
        return 1;
    }

    int q = 5;
    for (long long pos = 0; pos <= cantidad; pos++) {
        Secuencia[pos] = q;
        q += 2;
        while (P6[q] == 0) {
            q += 2;
        }
    }

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

    Z = posult;
    for (long long i = 0; i <= cantidad; i++) {
        long long mult = limite % Secuencia[i];
        if (mult % 2 == 1) {
            mult += Secuencia[i];
        }
        array[i] = mult;
    }

    int iteraciones = 0;

    #pragma omp parallel for
    for (j = 0; j <= cantidad; j++) {
        long long pos;
        long long primo;
        long long mult;
        mult = array[j];
        pos=j;
        primo = limite - mult;
        if((primo-3)%6==0){
            mult+=Secuencia[pos] * 2;
            primo=limite - mult;
        }
        if (P6[primo] == 1) {
            #pragma omp atomic
            fallo2++;
        }
        P6[primo] = 0;
                    
        if((primo-1)%6==0){
            if((Secuencia[pos]-1)%6==0){
                while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
                    mult += Secuencia[pos] * 2;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    mult+=Secuencia[pos] * 4;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    #pragma omp atomic
                    iteraciones++;
                }
            }
            else{
                while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
                    mult += Secuencia[pos] * 4;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    mult+=Secuencia[pos] * 2;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    #pragma omp atomic
                    iteraciones++;
                }
            }
        }
        else{
            if((Secuencia[pos]-1)%6==0){
                while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
                    mult += Secuencia[pos] * 4;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    mult+=Secuencia[pos] * 2;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    #pragma omp atomic
                    iteraciones++;
                }
            }
            else{
                while (primo >= Secuencia[pos] * Secuencia[pos] && primo >= n) {
                    mult += Secuencia[pos] * 2;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    mult+=Secuencia[pos] * 4;
                    primo = limite - mult;
                    if (P6[primo] == 1) {
                        #pragma omp atomic
                        fallo2++;
                    }
                    P6[primo] = 0;
                    #pragma omp atomic
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

    #pragma omp parallel for
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

    #pragma omp parallel for
    for (int i = 5; i * i <=limite; i = i + 4) {
        
        if (P7[i] == 1) {
            int salto = 2 * i;
            int salto2 = 4 * i;
            for (int j = i * i; j <=limite; j += salto2) {
                if (j < 2*n) {
                    if (P7[j] == 1) {
                        #pragma omp atomic
                        fallo3++;
                    }
                    P7[j] = 0;
                   
                }
                j += salto;
                if (j < 2*n) {
                    if (P7[j] == 1) {
                        #pragma omp atomic
                        fallo3++;
                    }
                    P7[j] = 0;
                    #pragma omp atomic
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
                        #pragma omp atomic
                        fallo3++;
                    }
                    P7[j] = 0;
                    
                }
                j += salto2;
                if (j < 2*n) {
                    if (P7[j] == 1) {
                        #pragma omp atomic
                        fallo3++;
                    }
                    P7[j] = 0;
                    #pragma omp atomic
                    iteraciones3++;
                }
            }
        }
    }
 
    clock_t fiEr2 = clock();
    
    int imprimir = 0;
    
    printf("quieres imprimir lo calculado  con el metodo hibrido 1 si 0 no \n");
    scanf("%d", &imprimir);
    if(imprimir == 1){
        for (long long i = 0; i <= limite; i++) {
            printf(" %llu es %d\n", i, P6[i]);
        }   
    }
    printf(" %d fallo erastotenes 1  \n", fallo);
    printf(" %d iteraciones erast 1 \n", iteraciones2);
    printf(" %d fallo met invertido\n", fallo2);
    printf(" %d iteraciones inverti\n", iteraciones);
    printf("comparativa\n");
    int difF1= fallo-fallo2;
    int difI1=iteraciones2-iteraciones;
    int falloP1=fallo+fallo2;
    int itP1=iteraciones+iteraciones2;
    printf(" %d fallo procesohibrido\n", falloP1);
    printf(" %d iteraciones proceso hibrido\n", itP1);
    printf(" %d fallo erastotes1-inversa\n", difF1);
    printf(" %d iteraciones erastotes1-inversa\n", difI1);
  
    printf("tabla de tiempos:\n");
    printf("hibrido\n");
    printf("erastotenes1\n");
   
    printf(" inicio %llu\n", inicio);
    int erastotes = fiEr - inicio;
    printf(" fin erastotes %llu\n", fiEr);
    printf("erastotes %llu\n", erastotes);
    printf("inversa\n");
    int tratamientodedatos = fiTD - fiEr;
    printf(" fin tratamiento de datos %llu\n", fiTD);
    printf(" tratamiento de datos %llu\n", tratamientodedatos);
    int inversa= fiINV - fiEr;
    printf(" fin inversa %llu\n", fiINV);
    printf(" inversa %llu\n", inversa);
    int total= fiINV - inicio;
    printf(" final %llu\n", fiINV);
    printf(" total %llu\n", total);
    int dif= erastotes-inversa;
    printf(" erastotes-inversa %llu\n", dif);
    printf("\n\n");
    printf("quieres imprimir lo calculado por erastotes de 0 a %d 1 si 0 no \n",limite);
    scanf("%d", &imprimir);  
    if(imprimir == 1){
        for (long long i = 0; i <= limite; i++) {
            printf(" %llu es %d\n", i, P7[i]);
        }   
    }
    printf("erastotenes  \n");
    int Ter2 =  fiEr2 - inicio2;
    printf(" %d fallo erastotenes 2 \n", fallo3);
    printf(" %d iteraciones erast 2 \n", iteraciones3);
    printf("tabla de tiempos:\n");
    printf(" inicio %llu\n", inicio2);
    printf(" fin erastotes %llu\n", fiEr2);
    printf(" erastotes %llu\n", Ter2);
    printf("comparativa\n");
    
    int diferencia=fallo3-falloP1;
    int diferenciP1P2=Ter2-total;  
    int diferencia2=iteraciones3-itP1;
    printf(" %d es la diferencia de fallos erastotes2 - hibrido\n", diferencia);
    printf(" %d es la diferencia de iteraciones erastotes -hibrido\n", diferencia2);
    printf("tabla de tiempos:\n");
    printf(" %llu hibrido\n", total);
    printf(" %llu erastotenes\n",Ter2);
    printf(" erastotes-hibrido %llu\n", diferenciP1P2);
    
    printf("fin\n");
    
    free(P6);
    free(Secuencia);
    free(array);
    free(P7);
    
    return 0;
}
