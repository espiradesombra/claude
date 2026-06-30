#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <sched.h>
#include <stdbool.h>
//#include <enmintrin.h>

int main() {
    long long n =12;
    long long max=n*n;
    long long m =0;
    int num_threads=10;

    double start_time6, end_time6;
    bool *P6 = (bool *) malloc((n + 1-m) * sizeof(bool));
    if (P6 == NULL) {
        printf("Error al asignar memoria dinámica en primos.\n");
        free(P6);
        return 1;
    }

    P6[2] = 1;
    P6[3] = 1;
    P6[5] = 1;
    P6[7] = 1;

    for (int i = 6; i <= n; i = i + 6) {
        P6[5 + i] = 1;
        P6[7 + i] = 1;
        P6[2 + i] = 0;
        P6[3 + i] = 0;
        P6[4 + i] = 0;
    }

    for (int i = 5; i <= n; i = i + 2) {
        int salto=2 * i;
        int salto2=4 * i;
        if(P6[i]==1){
            for (int j = i * i; j <= n; j += salto2) {
                P6[j] = 0;
                j += salto;
                P6[j] = 0;
            }
        }
        i += 4;
        if(P6[i]==1){
            for (int j = i * i; j <= n; j += salto) {
                P6[j] = 0;
                j += salto2;
                P6[j] = 0;
            }
        }
    }
    printf("lista de primos hasta raiz de max, metodo viejo.");

    for (long long i = m; i <= n; i++) {
        printf("%d es %lld\n", P6[i-m], i);
    }
    printf("comienzo del calculo\n");


    int c =0;
    printf("busco el ultimo primo en raiz\n");
    while (P6[n-c]==0){
        c++;
    }
    long long ultimo=(n-c);
    printf("%lld",ultimo);
    printf("cuento los primos");
    long long cantidad=0;
    for (long long cant=5;cant<=n;cant++){
        if(P6[cant]==1){
            cantidad++;
        }
    }
    long long *Secuencia= (long long*) malloc((cantidad)* sizeof(long long));
 if (Secuencia== NULL) {
        printf("Error al asignar memoria dinámica en Secuencia.\n");
        free(Secuencia);
        return 1;
    }
    long long cantid=0;
   for (long long pos=5;pos<=n;pos++){
    if(P6[pos]==1){
        Secuencia[cantid]=pos;
        cantid++;
   }
    }
    printf("Secuencia de primos\n");
    for (long long i = 0; i < cantid; i++) {
        printf("%lld\n", Secuencia[i]);
    }
    long long  *SecuenciaMP= (long long*) malloc((cantid-1)* sizeof(long long));
if (SecuenciaMP== NULL) {
        printf("Error al asignar memoria dinámica en Secuencia.\n");
        free(SecuenciaMP);
        return 1;
    }
    long long numMult;

    for (long long pos=0;pos<cantid-1;pos++){
        numMult=ultimo*ultimo/Secuencia[pos];
        SecuenciaMP[pos]=(ultimo*ultimo)%(  numMult  *Secuencia[pos] );
        if (SecuenciaMP[pos ]%2==1){SecuenciaMP[pos]+=Secuencia[pos];}      
    }
    printf("SecuenciaMP \n");
    for (long long i = 0; i < cantid-1; i++) {
        printf("%lld\n", SecuenciaMP[i]);
    }

    long long obj=0;
    long long s;
    long long o;
    long long primo;
    long long recal;
    long long limite = ultimo*ultimo;
    for (long long sec=0;sec<=limite-ultimo;sec+=2){

  for(long long sx=0;sx< cantid-1;sx++){
                if (SecuenciaMP[sx]<=SecuenciaMP[obj]){
                    obj=sx;
                               printf(" busco el minimo en secuenciaMP y obj es %lld\n",obj);

                     }

            }

        if (SecuenciaMP[obj]==sec){
            primo=limite-sec;
            printf("acabo de entrar a un if y sec es %lld, onj es %lld y SecuenciaMP[obj] es %lld, el numero que en teoria no es primo es %llu\n",sec,obj,SecuenciaMP[obj],primo);
       
              SecuenciaMP[obj]+=Secuencia[obj]*2;
        }else{  
            primo=limite-sec;
            printf("%lld es primo\n",primo);
            if (SecuenciaMP[obj]<sec){
                printf("entro a un if y sec es %lld, obj es %lld y SecuenciaMP[obj] es %lld\n",sec,obj,SecuenciaMP[obj]);
                SecuenciaMP[obj]+=Secuencia[obj]*2;
                printf("SecuenciaMP se ha actualizado a %lld\n",SecuenciaMP[obj]);


        }
    }
       
    }

    free(P6);
    free(Secuencia);
    free(SecuenciaMP);

    return 0;
}
