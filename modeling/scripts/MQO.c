#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 

void Transpose(int matriz[200][6], int matriz_destino[200][6]);
void MultiplicadorMatriz(int matriz1[200][6], int matriz2[200][6], int matrizDestino[200][6]);
void MultiplicadorVetor(int vetor[200], int matriz[200][6], int destino[200]);
void main(){
    printf("as funçoes de transpose e multiplicaçao de matriz e vetor ja estao prontas"); 

}
// design das funçoes 
// funçao para transposiçao de matriz 
void Transpose(int matriz[200][6], int matriz_destino[200][6]){
    int i, j; 
    for(i=0; i<2; i++){
        for(j=0; j<2; j++){
            matriz_destino[i][j] = matriz[j][i];
        }
    }
}

// funçao multiplicadora 
void MultiplicadorMatriz(int matriz1[200][6],int matriz2[200][6],int matrizDestino[200][6]){
    int i, j, k;
    for(i = 0; i < 200; i++){
        for(j = 0; j < 6; j++)
        {
        matrizDestino[i][j] = 0;
        for(k = 0; k < 200; k++){
            matrizDestino[i][j] += matriz1[i][k] * matriz2[k][j];
        }
        }
    }
}

void MultiplicadorVetor(int vetor[200], int matriz[200][6], int destino[200]){
    int i, j; 
    for(i=0; i<2; i++){
        for(j=0; j<6; j++){
           destino[j] =  matriz[i][j] * vetor[j]; 
        }
    }
}