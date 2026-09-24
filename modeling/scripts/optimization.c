#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 
#include <locale.h>
#include <string.h>
#include <math.h>

#define COLUNAS 6
#define LINHAS 165

void MultiplicadorMatrizVetor(double matriz1[LINHAS][LINHAS],double vetor[COLUNAS],double resultado[COLUNAS]);
void InverteOmega(double rho,double OmegaInv[1][1]);
void Transpose(double matriz[1][COLUNAS], double matriz_destino[COLUNAS][1]);
void InverteHx(double matriz[COLUNAS][COLUNAS], double inversa[COLUNAS][COLUNAS]); 
int main(){
    int i, j; 
    // assumindo X como chutes iniciais 
    double X[1][COLUNAS]={{1,0.0012157578675766612,1.12269006497541,0.0058070313939042,0.332757862491234,0.0024855445478073612}}, y, erro,tol,alpha,  beta[COLUNAS][1] = {-0.017306,-0.255698,0.027675,-0.586581,-0.011552,0.753279};
    double y_; 
/*
    for(i=0; i<1; i++){
        for(j=0; j<COLUNAS; j++){
            y_ += X[0][j] * beta[j][i];
        }
    }
    printf("real: %lf\n", y);
    printf("previsao (chute): %f", y_); 
    erro = y - y_; 
    printf("\n"); 
    printf("erro: %lf", erro); 
    printf("\n"); 
*/

    /* 
    ---------------------------------------------------------------
    farei a transposiçao do vetor beta para 6x1 para uma matriz 1 x 6, para satisfazer a preposiçao b^t 
    -----------------------------------------------------------------
    */
    y = 0.4953393512227677;  // alvo 
    erro = 1; 
    tol = 0.00001; 
    alpha=0.0001; 
    int k=0, maxIter=10000000; 
    double acc[100000]={0}; // limtando pela nao convergencia em um vetor de 10mi de tamanho
while(fabs(erro)>tol){
    y_=0; 
    for(i=0; i<1; i++){
        for(j=0; j<COLUNAS; j++){
            y_ += X[0][j] * beta[j][i];
        }
    }
    erro = y - y_;
    acc[k] = erro; 
    printf("rodada %d: erro: %lf", k, erro); 
    double betaT[1][COLUNAS];

    for(j=0;j<COLUNAS;j++){
    betaT[0][j]=beta[j][0];
                    }
    for(i=0; i<COLUNAS; i++){
        for(j=0; j<1; j++){
            betaT[0][i] = beta[i][0]; 
        }
    }

    // criando matriz omega inversa
    double OmegaInv[1][1]; 
    InverteOmega(0.065508, OmegaInv); 

    
    // multiplicando Omega por 2 
    double _OmegaInv[1][1]= {0}; 
    for(i=0;i<1; i++){
        for(j=0; j<1; j++){
            _OmegaInv[i][j] += 2 * OmegaInv[i][j]; 
        }
    } 
    // 2*OmegaInv * erro
    double resultado1[1][1]={0}, derivada1[1][6]={0} ; 

    resultado1[0][0] = _OmegaInv[0][0] * erro;  
    
    // resultado1 * beta^t 
        for(j=0; j<COLUNAS; j++){
            derivada1[0][j] = resultado1[0][0] *betaT[0][j];

        // gradiente Xn+1 = Xn + alpha*derivada1

            X[0][j]=X[0][j]+alpha*derivada1[0][j];
        }
        printf("Resultado final:\n");

    for(j=0;j<COLUNAS;j++){
        printf("X[%d]=%lf\n",j,X[0][j]);
    }
    

    k++;
    if(k>=maxIter){
        printf("maximo de iteracoes atingido\n");
    break;
}
    if(fabs(erro)<=tol){
        printf("convergencia atendida\n");
    break;
}
    }
    FILE *saida; 
    saida = fopen("errosOtimizacao.csv", "w+"); 
    for(i=0; i<100000; i++){
        fprintf(saida, "%lf\n", acc[i]);
    }
    fclose(saida); 
    
return 0;
} 

void InverteOmega(double rho,double OmegaInv[1][1]){
    double fator;

    fator = 1.0 / (1.0 - rho * rho);

    OmegaInv[0][0] = fator;
}

void Transpose(double matriz[1][COLUNAS], double matriz_destino[COLUNAS][1]){
    int i, j; 
    for(i=0; i<1; i++){
        for(j=0; j<COLUNAS; j++){
            matriz_destino[j][i] = matriz[i][j];
        }
    }
}
void MultiplicadorMatrizVetor(double matriz1[LINHAS][LINHAS],double vetor[COLUNAS],double resultado[COLUNAS]){
    int i, k;
    for(i = 0; i < COLUNAS; i++){
        resultado[i] = 0;
        for(k = 0; k < LINHAS; k++){
            resultado[i] += matriz1[i][k] * vetor[k];
        }
    }
}
