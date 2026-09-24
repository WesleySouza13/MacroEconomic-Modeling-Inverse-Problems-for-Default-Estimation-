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

int main(){
    int i, j; 
    // assumindo X como chutes iniciais 
    double X[1][COLUNAS]={{1,0.0012157578675766612,1.12269006497541,0.0058070313939042,0.332757862491234,0.0024855445478073612}}, y=-0.006666691358189203, erro=0, beta[COLUNAS][1] = {-0.017306,-0.255698,0.027675,-0.586581,-0.011552,0.753279};
    double y_ = 0; 

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


    /* 
    ---------------------------------------------------------------
    farei a transposiçao do vetor beta para 6x1 para uma matriz 1 x 6, para satisfazer a preposiçao b^t 
    -----------------------------------------------------------------
    */
    double betaT[1][COLUNAS];
    for(i=0;i<1; i++){
        for(j=0; j<COLUNAS; j++){
            betaT[0][i] = 0; 
        }
    }
    printf("betaT:");
    printf("\n"); 
    for(i=0; i<COLUNAS; i++){
        for(j=0; j<1; j++){
            betaT[0][i] = beta[i][0]; 
            printf("%lf\t", betaT[0][i]);
        }
    }

    // criando matriz omega inversa
    double OmegaInv[1][1]; 
    InverteOmega(0.065508, OmegaInv); 

    for(i=0; i<1; i++){
        for(j=0; j<1; j++){
            printf("\t %lf", OmegaInv[i][j]); 
        }printf("\n"); 
    }
    // multiplicando Omega por 2 
    double _OmegaInv[1][1]= {0}; 
    for(i=0;i<1; i++){
        for(j=0; j<1; j++){
            _OmegaInv[i][j] += 2 * OmegaInv[i][j]; 
        }
    } printf("multiplicaçao por 2 de OmegaInv feita com sucesso!");

    // 2*OmegaInv * erro
    double resultado1[1][1]={0}, derivada1[1][6]={0} ; 

    resultado1[0][0] = _OmegaInv[0][0] * erro; 
    printf("%lf", resultado1[0][0]); 
    
    // resultado1 * beta^t 
    printf("\n"); 
    printf("derivada primeira:");
        for(j=0; j<COLUNAS; j++){
            derivada1[0][j] = resultado1[0][0] *betaT[0][j];
            printf("derivada1: %lf\t", derivada1[0][j]);
        }
        
        /*
        calculando segunda derivada para Heussiana: 2OmegaInv*Beta*Beta^t = Hx
        */
        double Hx[COLUNAS][COLUNAS] = {0};
        printf("\nHessiana Hx:\n");
    for(i = 0; i < COLUNAS; i++){
        for(j = 0; j < COLUNAS; j++){
            Hx[i][j] =_OmegaInv[0][0] *beta[i][0] * betaT[0][j];
            printf("%lf\t", Hx[i][j]);
    }
    printf("\n");
}





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