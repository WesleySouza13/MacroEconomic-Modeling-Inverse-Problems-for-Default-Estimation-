#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 
#include <locale.h>
#include <string.h>> 
#include <math.h>
#define COLUNAS 7
#define LINHAS 165

// esboço das funçoes
void Transpose(double matriz[LINHAS][COLUNAS], double matriz_destino[COLUNAS][LINHAS]);
void MultiplicadorMatriz(double matriz1[COLUNAS][LINHAS],double matriz2[LINHAS][COLUNAS],double matrizDestino[COLUNAS][COLUNAS]);
void MultiplicadorVetor(double vetor[LINHAS], double matriz[COLUNAS][LINHAS], double destino[COLUNAS]);
void Inversa(double matriz[LINHAS][COLUNAS], double matrizDestino[LINHAS][COLUNAS]);
void main(){
    setlocale(LC_NUMERIC, "C");

    FILE *entrada, *saida_coef, *vetorResposta, *nomes_features; 
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
    // substituiçao regressiva
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
    saida_coef = fopen("resultados.txt", "w+"); 
    // pegando nomes das colunas 
    nomes_features = fopen("nome_variaveis.txt", "rt");
    char features[50];
    for(i=0; i<COLUNAS; i++){
        fscanf(nomes_features, "%s", features);
        fprintf(saida_coef, "%s: %lf \n", features, x[i]);
    }

    // soma dos erros 
    double prev[LINHAS], erro[LINHAS]; 
    for(i=0; i<LINHAS; i++){
        prev[i] = 0;
        for(j=0; j<COLUNAS; j++){
            prev[i] += Matriz[i][j] * x[j];

        }
    }
    printf("\n");
    printf("erros:");
    // calculando termo de erro 
    for(i=0; i<LINHAS; i++){
        erro[i] = 0;
        erro[i] += vetor[i] - prev[i];
        printf("%lf", fabs(erro[i]));
        erro[i] = fabs(erro[i]);
        printf("\n");
    }
    // soma dos erros 
    double somaErros; 
    for(i=0; i<LINHAS; i++){
        somaErros += erro[i];
    }
    fprintf(saida_coef, "\n====== Métricas =======");
    printf("\n");
    fprintf(saida_coef, "\nSOMA DOS ERROS: %lf\n", somaErros);

    double mediaReal=0; 
    for(i=0; i<LINHAS; i++){
        mediaReal += vetor[i];

    }
    mediaReal = mediaReal/LINHAS;
    printf("\n");
    printf("MÉDIA VALORES REAIS: %lf", mediaReal);
    printf("\n");

    double SomaQuadrados=0; 
    for(i=0; i<LINHAS; i++){
        SomaQuadrados += pow(vetor[i] - mediaReal, 2);
    }
    printf("\n");
    printf("SOMA DOS QUADRADOS: %lf", SomaQuadrados);
    printf("\n");

    double SomaQuadradosResid=0;
    for(i=0; i<LINHAS; i++){
        SomaQuadradosResid += pow(erro[i], 2);
    }

    printf("\n");
    printf("SOMA DOS QUADRADOS RESIDUAIS: %lf", SomaQuadradosResid);
    printf("\n");

    // calculo r2
    double r2=0; 
    r2 = 1-(SomaQuadradosResid/SomaQuadrados);
    printf("\n");
    printf("R²: %lf", r2);
    printf("\n");
    fprintf(saida_coef, "\nR2: %lf\n", r2);

    double MAE; 
    
    // calculo MAE 
    for(i=0; i<LINHAS; i++){
        MAE += fabs(vetor[i]-prev[i]);
    }
    MAE /=LINHAS; 
    printf("\n");
    printf("MAE: \n%lf\n", MAE);
    printf("\n");

    // calculo do sigma 
    double sigma; 

    sigma = SomaQuadradosResid/(LINHAS - COLUNAS);
    printf("\n");
    printf("SIGMA: \n%lf\n", sigma);
    fprintf(saida_coef, "\nσ^2: %lf\n", sigma);
    printf("\n");

    fprintf(saida_coef, "\nMAE: %lf\n", MAE);
    
    // calculo de F 
    double SQM=0, QMM=0, SQE=0, QME=0; 
    for(i=0; i<LINHAS; i++){
        SQM += pow(prev[i]-mediaReal, 2);
    }
    QMM = SQM/COLUNAS; 
    
    for(i=0; i<LINHAS; i++){
        SQE += pow(vetor[i]-prev[i],2);
    }
    QME = SQE/(LINHAS - COLUNAS -1);
    double F=0; 
    F = QMM/QME; 
    printf("\n");
    printf("F static: %lf", F);
    fprintf(saida_coef, "F statistic: %lf\n", F);
    
    // calculo de coeficiente de subestimaçao
    double CoefSubestimacao=0, T=0; 
    for(i=0; i<LINHAS; i++){
        if(prev[i]>mediaReal){
            T++;
        }
        CoefSubestimacao = T/LINHAS; 
    }
    printf("\n");
    printf("Coefiente de subestimaçao: %lf", CoefSubestimacao);
    fprintf(saida_coef, "coef. subestimacao: %lf\n", CoefSubestimacao);

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

// funçao multiplicadora de matriz
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
// multiplicaçao de vetor 
void MultiplicadorVetor(double vetor[LINHAS], double matriz[COLUNAS][LINHAS], double destino[COLUNAS]){
    int i, j; 
    for(j=0;j<COLUNAS;j++){
        destino[j]=0;
        for(i=0;i<LINHAS;i++){
            destino[j] += matriz[j][i] * vetor[i];
        }
    }
}
