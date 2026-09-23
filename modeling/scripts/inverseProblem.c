#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 
#include <locale.h>
#include <string.h>
#include <math.h>

void main(){
    
    float X=4.5, y_pred=0, y_real=45, erro=9999999,Xi=0; 
    int i=0;
    float tol=0.0001; 
    FILE *erros_file; 
    while(erro>tol){
        i++; 
        y_pred = 6*pow(X, 3) -2*pow(X, 2) + X + 3; // f(X)
        erro = (y_real - y_pred)*(y_real - y_pred); 
        printf("erro em %d:%lf", i, erro); 
        printf("\n");
        Xi = sqrt(X+2); // regra provisoria de erro 
        X=Xi; 
        printf("novo X: %lf", Xi);
        printf("\n"); 
        erros_file = fopen("saida_erros_busca.txt", "w+");
        //fprintf(erros_file, "novo X:%lf, erro:%lf \n", Xi, erro); 
        fprintf(erros_file, "novo X:%lf, erro:%lf \n", Xi, erro);


    }



}