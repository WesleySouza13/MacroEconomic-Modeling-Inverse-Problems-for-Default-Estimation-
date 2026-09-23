#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 
#include <locale.h>
#include <string.h>
#include <math.h>
#define COLUNAS 6
#define LINHAS 165

// esboço das funçoes
void Transpose(double matriz[LINHAS][COLUNAS], double matriz_destino[COLUNAS][LINHAS]);
void MultiplicadorMatriz(double matriz1[COLUNAS][LINHAS],double matriz2[LINHAS][COLUNAS],double matrizDestino[COLUNAS][COLUNAS]);
void MultiplicadorVetor(double vetor[LINHAS], double matriz[COLUNAS][LINHAS], double destino[COLUNAS]);
double Var(double resid[LINHAS], double media);
double ACF(double resid[LINHAS], double mediaResid, int lag);
double LjingBox(double acf[12]); 
void CalculoOmega(double matrizOmega[LINHAS][LINHAS], double rho, double destino[LINHAS][LINHAS]);
void InverteOmega(double rho,double OmegaInv[LINHAS][LINHAS]);
void MultiplicadorXOMega(double Xt[COLUNAS][LINHAS], double OmegaInv[LINHAS][LINHAS], double destino[COLUNAS][LINHAS]);
void InvertePt1(double Omega[COLUNAS][COLUNAS], double OmegaInv[COLUNAS][COLUNAS]);
void MultiplicadorMatrizVetor(double matriz1[COLUNAS][LINHAS],double vetor[LINHAS],double resultado[COLUNAS]);
void MultiplicadorVetorMatrizQuadrada(double matriz[COLUNAS][COLUNAS], double vetor[COLUNAS], double resultado[COLUNAS]);
void MultplicadorMatrizQuadrada(double matriz1[COLUNAS][COLUNAS], double matriz2[COLUNAS][COLUNAS], double resuldado[COLUNAS][COLUNAS]);
double F_statistic(double prev[LINHAS], double mediaReal, double vetor[LINHAS]);
// funçao para calculo direto de MQG
void MQGFunc(double rho_MQO, double Omega[LINHAS][LINHAS], double X[LINHAS][COLUNAS], double vetor_y[LINHAS], int maxIter, double tol, double BETAS[COLUNAS]); 

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
    fprintf(saida_coef, "coef. subestimacao: %lf\n \n", CoefSubestimacao);


    // calculo de variancia 
    double var=0; 
    var = Var(erro, mediaReal);

    printf("\n");
    printf("Var: %lf", var);
    
    // calculo de phi para acf 
    double phi, mediaResid=0;
    for(i=0; i<LINHAS; i++){
        mediaResid +=erro[i]; 
    }
    mediaResid /=LINHAS; 
    printf("\n");
    printf("Média Residuos: %lf", mediaResid);
    double acf=0; 
    acf = ACF(erro, mediaResid, 12); 
    printf("\n");
    printf("ACF para lag12: %lf", acf); 
    
    // calculando varios rho para acf
    double rho[12];
    for(i=0; i<12; i++){
        rho[i] = ACF(erro, mediaResid, i+1);
        printf("rho: %lf, lag:%d",rho[i],i+1);
        printf("\n");
        fprintf(saida_coef, "rho:%lf, lag:%d\n", rho[i], i+1); 
        
    }
    // calculando ljingbox
    double pt1=0, pt2=0, Q=0;
    double k_=0; 
    pt1=LINHAS*(LINHAS+2);  
    for(i=0; i<12; i++){
        k_ = i+1;
        pt2 += pow(rho[i],2)/(LINHAS-k_); 

    }
    Q = pt1*pt2;
    printf("\n");
    printf("ljung-box: %lf", Q);
    fprintf(saida_coef, "\nLjing-Box: %lf\n", Q);
    
    /// montando GLS para correçao de autocorrelaçao e heterocedasticidade do modelo
    double rho_new;
    rho_new = ACF(erro, mediaResid, 1);
    // calculo de omega 
    double Omega[LINHAS][LINHAS],destinoOmega[LINHAS][LINHAS]; 

    // preenchendo matrizes com zeros 
    for(i=0; i<LINHAS; i++){
        for(j=0; j<LINHAS; j++){
            Omega[i][j] =0;
            destinoOmega[i][j] = 0; 
        }
    }
    CalculoOmega(Omega, 0.792944736, destinoOmega);
    printf("\n");

    // calculando inversa do Omega 
    double OmegaInv[LINHAS][LINHAS]; 
    for(i=0; i<LINHAS; i++){
        for(j=0; j<LINHAS; j++){
            OmegaInv[i][j] = 0; 
        }
    }
    InverteOmega(0.792944736, OmegaInv); 

    // calculando Xt * OmegaINv 
    double XtOmegaInv[COLUNAS][LINHAS]; 
    MultiplicadorXOMega(Xt, OmegaInv, XtOmegaInv);
    
    //calculando XtOmegaInvX
    double XtOmegaInvX[COLUNAS][COLUNAS];
    MultiplicadorMatriz(XtOmegaInv, Matriz,XtOmegaInvX);
    for(i=0; i<COLUNAS; i++){
        for(j=0; j<COLUNAS; j++){
            printf("%lf", XtOmegaInvX[i][j]); 
        }printf("\n");
    }
    // segunda parte (XtOmegaInvY)
    double XtOmegaInvY[COLUNAS]; 
    MultiplicadorMatrizVetor(XtOmegaInv, vetor, XtOmegaInvY); 
    printf("\n");
    for(i=0; i<COLUNAS; i++){
        printf("%lf", XtOmegaInvY[i]);
        printf("\n");
    }
    // invertendo a primeira parte
    double InvPt1[COLUNAS][COLUNAS]; 
    InvertePt1(XtOmegaInvX, InvPt1);
    printf("\n");
    for(i=0; i<COLUNAS; i++){
        for(j=0; j<COLUNAS; j++){
            printf("%lf", InvPt1[i][j]);
        }
    }
    printf("\n");
    double betas[COLUNAS]; 
    MultiplicadorVetorMatrizQuadrada(InvPt1, XtOmegaInvY, betas);
    for(i=0; i<COLUNAS; i++){
        printf("%lf", betas[i]);
        printf("\n");
    }
    // printando coeficientes no arquivo de saida
    FILE *saida_GLS, *features_gls; 
    saida_GLS = fopen("resultados_MQG.txt", "w"); 
    features_gls = fopen("nome_variaveis.txt", "rt");

    char gls_features[50]; 
    for(i=0; i<COLUNAS; i++){
        fscanf(features_gls, "%s", gls_features);
        fprintf(saida_GLS, "%s: %lf\n", gls_features, betas[i]);
    }

    // fazendo previsoes para MQG
    double prev_MQG[LINHAS]; 
    for(i=0; i<LINHAS; i++){
        prev_MQG[i]=0;
    }
    for(i=0; i<LINHAS; i++){
        for(j=0; j<COLUNAS; j++){
            prev_MQG[i] += Matriz[i][j] * betas[j]; 
        }
    }
    // calculando erros para MQG 
    double resid_MQG[LINHAS],  resid_MQG_media=0; 
    for(i=0; i<LINHAS; i++){
        resid_MQG[i] = vetor[i] - prev_MQG[i];
        resid_MQG_media +=resid_MQG[i];
    }
    double rho_mqg=0;
    resid_MQG_media = (resid_MQG_media/ LINHAS); 
    rho_mqg = ACF(resid_MQG, resid_MQG_media, 1);
    
    // ----------------------------------------------
    // otimizaçao de rho para minimizar erros
    // ----------------------------------------------
    double tol, BETAS_FGLS[COLUNAS]; 
    FILE *saidaFGLS; 
    saidaFGLS = fopen("Saida_FGLS_.txt", "w"); 
    for(i=0; i<COLUNAS; i++){
        BETAS_FGLS[COLUNAS] = 0;
    }
    printf("\n");
    tol = pow(10, -6);
    printf("fazendo otimizaçao do rho por FGLS");
    printf("\n");
    MQGFunc(rho_new, Omega, Matriz, vetor, 10, tol, BETAS_FGLS);

    // fazendo previsoes com FGLS
    double prev_FGLS[LINHAS], resid_FGLS[LINHAS], mediaResidFGLS=0, mediaprevFGLS=0;  
    for(i=0; i<LINHAS; i++){
        for(j=0; j<COLUNAS; j++){
            prev_FGLS[i] += Matriz[i][j] * BETAS_FGLS[j];
            resid_FGLS[i] += vetor[i] - prev_FGLS[i]; 
            mediaprevFGLS += prev_FGLS[i]; 
        }
    }
    mediaprevFGLS /=LINHAS; // media das previsoes de FGLS
    // calculando média de residuos FGLS 
    for(i=0; i<LINHAS; i++){
        mediaResidFGLS += resid_FGLS[i]; 
    }
    mediaResidFGLS /= LINHAS; 

    // CALCULANDO LJING-BOX PARA FGLS 
    double acfFGLS[12] ;
    for(i = 0; i < 12; i++){
        acfFGLS[i] = ACF(resid_FGLS, mediaResidFGLS, i + 1);
        printf("Lag %d: %lf\n", i + 1, acfFGLS[i]);
        
    }
    double ljingBox_FGLS=0;
    ljingBox_FGLS = LjingBox(acfFGLS); 
    FILE *FGLS_features; 
    FGLS_features = fopen("nome_variaveis.txt", "rt"); 
    char features_FGLS[50]; 
    for(i=0; i<6; i++){
        fscanf(FGLS_features, "%s", features_FGLS);
        fprintf(saidaFGLS, "%s: %lf\n", features_FGLS, BETAS_FGLS[i]);
    }
    for(i=0; i<12; i++){
        fprintf(saidaFGLS, "\nrho: %lf, lag: %d", acfFGLS[i], i+1); 
    }
    fprintf(saidaFGLS, "\n\n(Q) Ljing-Box: %lf", ljingBox_FGLS); 

    double F_fgls=0, sub_coef_FGLS=0; 
    F_fgls = F_statistic(prev_FGLS, mediaprevFGLS, vetor); 
    fprintf(saidaFGLS, "\nF_statistic: %lf\n", F_fgls); 
    

    // exportando matriz omega 
    FILE *saida_omega; 
    saida_omega = fopen("matriz_omega.txt", "w+"); 
    for(i=0; i<LINHAS; i++){
        for(j=0; j<LINHAS; j++){
            fprintf(saida_omega, "%lf\n", Omega[i][j]);
        }
    }
    fclose(saida_omega); 
}
//fim da main()
    


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
double Var(double resid[LINHAS], double media){
    double var=0;
    int i;  
    for(i=0; i<LINHAS; i++){
        var +=pow(resid[i]-media,2); 
    }
    return (var)/(LINHAS-1); 
}
double ACF(double resid[LINHAS], double mediaResid, int lag){
    int i;
    double numerador=0, denominador=0; 

    for(i=lag; i<LINHAS; i++){
        numerador += (resid[i]-mediaResid)*(resid[i-lag]-mediaResid);
        denominador+=pow(resid[i]-mediaResid, 2);

    }
    
    return numerador/denominador; 

}
double LjingBox(double acf[12]){
    double soma = 0.0;
    for(int k = 1; k <= 12; k++){
        soma += (acf[k-1] * acf[k-1]) / (LINHAS - k);
    }
    return LINHAS * (LINHAS + 2) * soma;
}
// calculo de matriz omega 
void CalculoOmega(double matrizOmega[LINHAS][LINHAS], double rho, double destino[LINHAS][LINHAS]){
    int i, j;
    double distancia;  
    for(i=0; i<LINHAS; i++){
        for(j=0; j<LINHAS; j++){
            distancia = fabs(i-j);
            matrizOmega[i][j] = pow(rho, (double)distancia); 
            destino[i][j] = matrizOmega[i][j];
        }
    }
}

void InverteOmega(double rho,double OmegaInv[LINHAS][LINHAS]){
    int i, j;

    double fator = 1.0 / (1.0 - rho * rho);
    // zera a matriz
    for(i = 0; i < LINHAS; i++){
        for(j = 0; j < LINHAS; j++){
            OmegaInv[i][j] = 0.0;
        }
    }
    OmegaInv[0][0] = fator;
    OmegaInv[LINHAS-1][LINHAS-1] = fator;

    // Diagonal interna
    for(i=1; i<LINHAS-1; i++){
        OmegaInv[i][i] =
            fator * (1.0 + rho * rho);
    }
    for(i = 0; i < LINHAS-1; i++){
        OmegaInv[i][i+1] = -fator * rho;
        OmegaInv[i+1][i] = -fator * rho;
    }
}
void MultiplicadorXOMega(double Xt[COLUNAS][LINHAS], double OmegaInv[LINHAS][LINHAS], double destino[COLUNAS][LINHAS]){
    int i, j, k; 
    for(i=0; i<COLUNAS; i++){
        for(j=0; j<LINHAS; j++){
            destino[i][j] = 0;
            for(k=0; k<LINHAS; k++){
                destino[i][j] +=Xt[i][k] * OmegaInv[k][j];
            }
        }
    }
}
// funçao para inverter pt1 do GLS
void InvertePt1(
    double matriz[COLUNAS][COLUNAS],
    double inversa[COLUNAS][COLUNAS])
{
    int i, j, k, linhaPivo;
    double fator, maior, temp;
    double A[COLUNAS][2 * COLUNAS];
    // [A | I]
    for(i = 0; i < COLUNAS; i++){
        for(j = 0; j < COLUNAS; j++){
            A[i][j] = matriz[i][j];
        }
        for(j = 0; j < COLUNAS; j++){
            A[i][j + COLUNAS] = (i == j) ? 1.0 : 0.0;
        }
    }
    // Gauss-Jordan com pivoteamento parcial
    for(k = 0; k < COLUNAS; k++){
        // procura maior pivô
        linhaPivo = k;
        maior = fabs(A[k][k]);
        for(i = k + 1; i < COLUNAS; i++){
            if(fabs(A[i][k]) > maior){
                maior = fabs(A[i][k]);
                linhaPivo = i;
            }
        }
        // matriz singular ou quase singular
        if(maior < 1e-12){
            printf("Matriz singular ou quase singular.\n");
            return;
        }
        // troca de linhas
        if(linhaPivo != k){
            for(j = 0; j < 2 * COLUNAS; j++){
                temp = A[k][j];
                A[k][j] = A[linhaPivo][j];
                A[linhaPivo][j] = temp;
            }
        }
        // normaliza pivô
        fator = A[k][k];
        for(j = 0; j < 2 * COLUNAS; j++){
            A[k][j] /= fator;
        }
        // zera coluna
        for(i = 0; i < COLUNAS; i++){
            if(i != k){
                fator = A[i][k];
                for(j = 0; j < 2 * COLUNAS; j++){
                    A[i][j] -= fator * A[k][j];
                }
            }
        }
    }
    // extrai inversa
    for(i=0; i < COLUNAS; i++){
        for(j = 0; j < COLUNAS; j++){
            inversa[i][j] = A[i][j + COLUNAS];
        }
    }
}

void MultiplicadorMatrizVetor(double matriz1[COLUNAS][LINHAS],double vetor[LINHAS],double resultado[COLUNAS]){
    int i, k;
    for(i = 0; i < COLUNAS; i++){
        resultado[i] = 0;
        for(k = 0; k < LINHAS; k++){
            resultado[i] += matriz1[i][k] * vetor[k];
        }
    }
}
void MultiplicadorVetorMatrizQuadrada(double matriz[COLUNAS][COLUNAS], double vetor[COLUNAS], double resultado[COLUNAS]){
    int i, j; 
    for(i=0; i<COLUNAS; i++){
        resultado[i] = 0; 
        for(j=0; j<COLUNAS; j++){
            resultado[i] += matriz[i][j] * vetor[j];
        }
    }
}
void MultplicadorMatrizQuadrada(double matriz1[COLUNAS][COLUNAS], double matriz2[COLUNAS][COLUNAS], double resuldado[COLUNAS][COLUNAS]){
    int i, j, k; 
    for(i=0; i<COLUNAS; i++){
        for(j=0; j<COLUNAS; j++){
            resuldado[i][j] =0; 
            for(k=0; k<COLUNAS; k++){
                resuldado[i][j] +=matriz1[i][k] *matriz2[k][j];
            }
        }
    }
}

void MQGFunc(double rho_MQO, double Omega[LINHAS][LINHAS], double X[LINHAS][COLUNAS], double vetor_y[LINHAS], int maxIter, double tol, double BETAS[COLUNAS]){
    int i, j, k; 

    for(k=0; k< maxIter; k++){
    // inicializaçao de omega
    double Xt[COLUNAS][LINHAS], OmegaDestino[LINHAS][LINHAS], OmegaInv[LINHAS][LINHAS], XtOmegaInv[COLUNAS][LINHAS], XtOmegaInvX[COLUNAS][COLUNAS], pt1[COLUNAS][COLUNAS], pt2[COLUNAS], betas[COLUNAS];
    // fazendo matriz transposta de X 
    Transpose(X, Xt); 
    // calculando Omega
    CalculoOmega(Omega, rho_MQO, OmegaDestino); 
    // inverte omega
    InverteOmega(rho_MQO,OmegaInv);
    // Xt * OmegaInv
    MultiplicadorXOMega(Xt, OmegaInv,XtOmegaInv); 
    // Xt * OmegaInv * X 
    MultiplicadorMatriz(XtOmegaInv, X, XtOmegaInvX); 
    // invertendo primeira parte da forma matricial 
    InvertePt1(XtOmegaInvX, pt1); 

    // ----------
    // pt2 -> Xt * OmegaInv * Y 
    // ---------
    //XtOmegaInv * Y
    MultiplicadorMatrizVetor(XtOmegaInv, vetor_y, pt2); 

    // resolvendo matrizes finais 
    // pt1 * pt2
    MultiplicadorVetorMatrizQuadrada(pt1, pt2, betas); 
    
    double prev[LINHAS], resid[LINHAS], mediaResid=0, rho_new;
    // inicializando prev e resid com 0
    for(i=0; i<LINHAS; i++){
        prev[i]=0; 
        resid[i]=0; 
    }
    for(i=0; i<LINHAS; i++){
        for(j=0;j<COLUNAS; j++){
            prev[i] += X[i][j] * betas[j];
        }
        resid[i] = vetor_y[i] - prev[i]; 
        mediaResid += resid[i];
    }
    mediaResid/=LINHAS; 

    // calculando novo rho 
    rho_new = ACF(resid, mediaResid, 1); 
    printf("rho inicial: %lf", rho_MQO);
    printf("\n");
    printf("novo rho: %lf", rho_new);
    printf("\n");
    if(fabs(rho_MQO - rho_new) < tol){
        printf("FGLS convergido com sucesso");
        rho_MQO = rho_new; 
        printf("\n");
        printf("betas:");
        FILE *saida, *nomes; 
        //saida = fopen("Saida_FGLS.txt", "w");
        nomes = fopen("nome_variaveis.txt", "rt");
        char features[50];
        for(i=0; i<COLUNAS; i++){
            printf("%lf", betas[i]);
            fscanf(nomes, "%s", features);
            //fprintf(saida, "%s: %lf\n\n", features, betas[i]); 
            BETAS[i] = betas[i];
            printf("\n");
        }
        break;
    }
    rho_MQO = rho_new;  
}

}
double F_statistic(double prev[LINHAS], double mediaReal, double vetor[LINHAS]){
    // calculo de F 
    int i; 
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
    return QMM/QME; 
}
