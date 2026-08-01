#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 
#include <locale.h>

#define COLUNAS 6
#define LINHAS 171

// esboço das funçoes
void Transpose(double matriz[LINHAS][COLUNAS], double matriz_destino[COLUNAS][LINHAS]);
void MultiplicadorMatriz(double matriz1[COLUNAS][LINHAS],double matriz2[LINHAS][COLUNAS],double matrizDestino[COLUNAS][COLUNAS]);
void MultiplicadorVetor(double vetor[LINHAS], double matriz[COLUNAS][LINHAS], double destino[COLUNAS]);
void Inversa(double matriz[LINHAS][COLUNAS], double matrizDestino[LINHAS][COLUNAS]);
void main(){
    setlocale(LC_NUMERIC, "C");

    FILE *entrada, *saida, *vetorResposta; 
    entrada = fopen("data_modelagem.txt", "rt"); 
    if (entrada == NULL) {
    perror("Erro ao abrir o arquivo");
    return 1;
}
    // setando tamanho da matriz e preenchendo com os dados 
    double Matriz[LINHAS][COLUNAS], vetor[LINHAS], r; 
    int i, j; 
    // preechendo matriz com zeros 
    for(i=0; i<LINHAS; i++){
        for(j=0; j<COLUNAS; j++){
            Matriz[i][j]=0;
        }
    }
    // preenchendo matriz com zeros
    for (i = 0; i < LINHAS; i++) {
        for (j = 0; j < COLUNAS; j++) {
        if (j < COLUNAS - 1)
            r = fscanf(entrada, "%lf,", &Matriz[i][j]);
        else
            r = fscanf(entrada, "%lf", &Matriz[i][j]);
        if (r != 1) {
            printf("Erro em (%d,%d)\n", i, j);
            return 1;
        }
    }
}
    /*
    for(i=0; i<LINHAS; i++){
            for(j=0; j<COLUNAS; j++){
            printf("%lf", Matriz[i][j]);
        }
    }
        */
    
    // preechendo vetor resposta
    vetorResposta = fopen("vetor_resposta.txt", "rt");
    if(vetorResposta==NULL){
        printf("o vetor esta nulo.");
    }
    // preenchendo vetor com zeros 
    for(i=0; i<LINHAS; i++){
        vetor[i]=0;
        
    }
    //preechendo
    for(i=0; i<LINHAS; i++){
        fscanf(vetorResposta, "%lf", &vetor[i]);
    }

    /*
    for(i=0; i<LINHAS; i++){
        printf("%lf \n", vetor[i]);
    }
    */
    fclose(entrada);
    fclose(vetorResposta);
    
    double Xt[COLUNAS][LINHAS], XtX[COLUNAS][COLUNAS], Xty[COLUNAS];

    // calculando Xt 
    Transpose(Matriz, Xt);
    printf("Matriz transposta calculada com Sucesso!");
    printf("\n");
    // calculando XtX
    MultiplicadorMatriz(Xt, Matriz, XtX);
    printf("XtX calculada com sucesso!");
    printf("\n");
    //calculando Xty 
    MultiplicadorVetor(vetor, Xt, Xty); 
    printf("Xty calculado com sucesso!");
    printf("\n"); 
    printf("calculando a matriz aumentada");
    // calculando matriz aumentada 
    double A[COLUNAS][COLUNAS+1]; 

    for(i=0; i<COLUNAS; i++){
        for(j=0; j<COLUNAS; j++){
            A[i][j] = XtX[i][j]; 
        }
    }
    for(i=0; i<COLUNAS; i++){
            A[i][COLUNAS] = Xty[i]; 
        }
        /*
        printf("\n");
        printf("matriz aumentada:");
    
    for(i=0; i<LINHAS; i++){
        for(j=0; j<COLUNAS; j++){
            printf("%lf", A[i][j]); 
        }printf("\n");
    }
        */
        printf("\n");
    // calculando gauss jordan 
    double x[COLUNAS], fator;
    int k; 
    // Eliminação de Gauss
    for (k = 0; k < COLUNAS - 1; k++){
        for (i = k + 1; i < COLUNAS; i++){
            fator = A[i][k] / A[k][k];
                for (j = k; j <= COLUNAS; j++){
                    A[i][j] = A[i][j] - fator * A[k][j];
                }
        }
    }
    // Substituição regressiva
    for (i = COLUNAS - 1; i >= 0; i--){
        x[i] = A[i][COLUNAS];
        for (j = i + 1; j < COLUNAS; j++){
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    printf("\n");
    printf("Solução - Sistema por Gauss Jordan:\n");
    for (i = 0; i < COLUNAS; i++){
        printf("x%d = %.6lf\n", i + 1, x[i]);
    }
    
    // printando elementos matriciais 
    printf("\n");
    printf("XtX:");
    printf("\n");
    for(i=0; i<COLUNAS; i++){
        for(j=0; j<COLUNAS; j++){
            printf(" %lf", XtX[i][j]);
        }printf("\n");
    }
    printf("\n");
    printf("Xty:");
    printf("\n");
    for(i=0; i<COLUNAS; i++){
            printf("%lf", Xty[i]);
            printf("\n");
        }
// printando coeficientes 

}//fim da main()
    


// design das funçoes 
// funçao para transposiçao de matriz 
void Transpose(double matriz[LINHAS][COLUNAS], double matriz_destino[COLUNAS][LINHAS]){
    int i, j; 
    for(i=0; i<LINHAS; i++){
        for(j=0; j<COLUNAS; j++){
            matriz_destino[j][i] = matriz[i][j];
        }
    }
}

// funçao multiplicadora 
void MultiplicadorMatriz(double matriz1[COLUNAS][LINHAS],double matriz2[LINHAS][COLUNAS],double matrizDestino[COLUNAS][COLUNAS]){
    int i,j,k;
    for(i=0;i<COLUNAS;i++){
        for(j=0;j<COLUNAS;j++){
            matrizDestino[i][j]=0;
            for(k=0;k<LINHAS;k++){
                matrizDestino[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }
}


void MultiplicadorVetor(double vetor[LINHAS], double matriz[COLUNAS][LINHAS], double destino[COLUNAS]){
    int i, j; 
    for(j=0;j<COLUNAS;j++){
        destino[j]=0;
        for(i=0;i<LINHAS;i++){
            destino[j] += matriz[j][i] * vetor[i];
        }
    }
}
