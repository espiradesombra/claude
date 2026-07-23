#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <math.h>


typedef struct {
    bool primo;
   }Valor;

int main() {

   
    long long int n = 0;
   
    printf("Ingrese el valor de n: que es el rango en el que se calculan primos\n");
    scanf("%llu", &n);
 
    
    Valor *P0 = (  Valor *)malloc((6 * 5) * sizeof(  Valor));
    if (P0 == NULL) {
        printf("Error al asignar memoria dinámica en primos P6.\n");
        free(P0);
        return 1;
    }
    Valor *P6 = P0;

   long long int indice = 5;
   long long int primo = indice;



for(long long int i = 5; i < 6*primo; i+=6) {
    P6[i].primo = 1;
    P6[i+2].primo = 1;
}
//N^(1/2)/ln N^(1/2)
  while(primo*primo<=n){  
    //P6 contiene los primos menores a primo incluyendolo a el y al siguiente.
    //P contiene el patron simetrico de P1
    //P6 sera P1 muchas veces. 
    
     Valor *P = (  Valor *)malloc((6 *primo) * sizeof(  Valor));
    if (P == NULL) {
        printf("Error al asignar memoria dinámica en primos P6.\n");
        free(P);
        return 1;
    }
   //copiamos en P P0 primo veces
       for (long long int j = 0; j < primo; j++) { 
            memcpy(&P[j*(6*primo)], P0, 6*primo * sizeof(Valor));
        }

    //eliminamos los multiplos
    P[6*primo-primo].primo = 0;
    P[primo].primo = 0;


    
//P1 almacena los multiplos de primo que estan en P
     Valor *P1 = (  Valor *)malloc((6 *primo*6*primo) * sizeof(  Valor));
    if (P1 == NULL) {
        printf("Error al asignar memoria dinámica en primos P6.\n");
        free(P1);
        return 1;
    }
    //se copia P en P1, 6*primo veces
    for (long long int j = 0; j < 6*primo; j++) { 
        memcpy(&P1[j*(6*primo)], P, 6*primo * sizeof(Valor));
    }
    //Patron de multiplos generado
    //p1 contiene muchos noprimos a 1 
    //p0 no varia 
    //p6 ahora se actualiza, va creciendo la informacion en P6
    P1[primo].primo = 1;//primo devuelto a 1
    free(P6);
   
   Valor  *P6 = P1; //actualizo P6

    free(P1);
    free(P);//borro P1 y P

// indice apunta de p6 al siguiente primo 
while( P6[indice].primo==0){
    indice++;}
    primo = indice;
  
}
//imprimo los primos
for (long long int i = 0; i < 6*n; i++) {
    if (P6[i].primo == 1) {
        printf(" %llu es primo\n", i);
    }
}
//pregunta si cerrar el programa
    printf("Presione una tecla para cerrar el programa\n");
    getchar();
   
 return 0;
 }