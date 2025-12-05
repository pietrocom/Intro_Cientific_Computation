/*
Autor: Pietro Comin
GRR:   20241955

Observações: 
 1- O trabalho foi pensado e projetado não para ser o mais eficiente possível, pois
    aprimoramentos de eficiência caberão ao trabalho 2, continuação deste. Isso im-
    plica que algumas decisões de implementação foram tomadas não levando como obje-
    tivo principal tornar o código ótimo, mas sim funcional.
 2- O gerenciamento de memória fica resposável pelas funções geradoras da estruturas 
    de dados, e não pela main.
*/

#include <stdio.h>
#include <stdlib.h>

#include "sislin.h"
#include "utils.h"

// Define uma tolerância/range que o valor pode ter para não ser considerado zero
#define ZERO_TOLERANCE 1e-12

// Altera o modo de visualização dos resultados finais
//  - 1 para visualização clean e legível
//  - 0 para visualização como solicitada no enunciado
#define OUTPUT_VISUALIZATION_MODE 0
// Limite de tamanho para print do vetor solução
#define MAX_VECTOR_PRINT 15

// Prorótipo da função que contém o loop e lógica principal.
real_t resolveSL(real_t **A, real_t *b, real_t **X, int n, int maxit, real_t epsilon, real_t omega, rtime_t *tempo_iter, DLU_matrices_t *dlu, int *out_total_iterations);

int main() {
    // Parâmetros de entrada
    int n, k, maxit;
    real_t omega, epsilon;

    // Ponteiros para as matrizes e vetores
    real_t **A = NULL, **ASP = NULL;
    DLU_matrices_t *dlu_mats = NULL;
    real_t *B = NULL, *bsp = NULL, *X = NULL;
    
    // Variáveis para medição de tempo
    rtime_t t_pc = 0.0, t_iter = 0.0, t_residuo = 0.0;
    
    // Variáveis para os resultados finais
    real_t norma_erro = 0.0, norma_residuo = 0.0;
    int total_iterations = 0;
    
    // Lê os 5 valores da entrada padrão (stdin)
    scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon);
    // OBS: A validação dos parâmetros é feita dentro de criaKDiagonal.
    
    srandom(20252);

    // Gera a matriz original A (k-diagonal) e o vetor B
    criaKDiagonal(n, k, &A, &B);

    // Gera a matriz Simétrica Positiva Definida (ASP = At * A) e o vetor (bsp = At * b)
    // Esta etapa corresponde ao "Tempo PC" (tempo de pré-condicionamento).
    genSimetricaPositiva(A, B, n, &ASP, &bsp, &t_pc);

    // Para os pré-condicionadores Gauss-Seidel e SSOR
    if (omega >= 1.0 && omega < 2.0) {
        rtime_t t_dlu = 0.0;
        geraDLU(ASP, n, &dlu_mats, &t_dlu); 
        t_pc += t_dlu;
    }

    // ---- RESOLUÇÃO DO SISTEMA ----
    
    // A função resolveSL alocará a memória para o vetor solução X e retornará a norma do erro final
    norma_erro = resolveSL(ASP, bsp, &X, n, maxit, epsilon, omega, &t_iter, dlu_mats, &total_iterations);

    
    // ---- CÁLCULO DO RESÍDUO ----
    
    // Calcula a norma do resíduo b - Ax usando a matriz A e o vetor B originais
    norma_residuo = calcResiduoSL(A, B, X, n, &t_residuo);


    // ---- IMPRESSÃO ----
    if (OUTPUT_VISUALIZATION_MODE == 1) {
        printf("Dimensão do Sistema (n): %d\n", n);
        
        if (n <= MAX_VECTOR_PRINT) {
            printf("Vetor Solução X:\n");
            for (int i = 0; i < n; ++i) {
                printf("  x[%d] = %.16g\n", i, X[i]);
            }
        } else {
            printf("Vetor Solução X: (suprimido para n > %d)\n", MAX_VECTOR_PRINT);
        }

        printf("\n--- Métricas de Convergência e Desempenho ---\n");
        printf("Iterações para Convergência:               %d de %d\n", total_iterations, maxit);
        printf("Norma Máx. do Erro (||x_k - x_{k-1}||∞):   %.8g\n", norma_erro);
        printf("Norma Euclidiana do Resíduo (||b - Ax||₂): %.8g\n", norma_residuo);
        printf("Tempo de Pré-Cálculo (ms):                 %.8g\n", t_pc);
        printf("Tempo Médio por Iteração (ms):             %.8g\n", t_iter);
        printf("Tempo de Cálculo do Resíduo (ms):          %.8g\n", t_residuo);
        printf("--------------------------------------------\n");
    }
    else {
        printf("%d\n", n);
        
        // Imprime o vetor solução X
        for (int i = 0; i < n; ++i) {
            printf("%.16g ", X[i]);
        }
        printf("\n");

        // Imprime as normas e os tempos
        printf("%.8g\n", norma_erro);
        printf("%.8g\n", norma_residuo);
        printf("%.8g\n", t_pc);
        printf("%.8g\n", t_iter);
        printf("%.8g\n", t_residuo);
    }

    // Libera toda memória alocada
    destroiKDiagonal(n, A, B);
    destroiKDiagonal(n, ASP, bsp);
    destroiDLU(dlu_mats, n);
    if (X) free(X); // Libera o vetor solução

    return 0;
}

real_t resolveSL(real_t **A, real_t *b, real_t **X, int n, int maxit, real_t epsilon, real_t omega, rtime_t *tempo_iter, DLU_matrices_t *dlu, int *out_total_iterations) {
    
    // Aloca o vetor solução X já com zeros.
    *X = calloc(n, sizeof(real_t));
    if(!(*X)) {
        fprintf(stderr, "ERRO: Durante alocação de memória!\n");
        exit(EXIT_FAILURE);
    }

    // Aloca os vetores de trabalho necessários para o algoritmo
    real_t *r = malloc(sizeof(real_t) * n);
    real_t *p = malloc(sizeof(real_t) * n);
    real_t *z = malloc(sizeof(real_t) * n);
    real_t *Ap = malloc(sizeof(real_t) * n);
    real_t *x_velho = malloc(sizeof(real_t) * n); // Usado para o critério de parada

    if (!r || !p || !z || !Ap || !x_velho) {
        fprintf(stderr, "ERRO: Durante alocação de memória!\n");
        exit(EXIT_FAILURE);
    }

    // Medição de tempo total das iterações
    rtime_t t_inicio = timestamp();
    
    // Como x0 é 0, o resíduo inicial r0 = b - A*x0 é simplesmente b.
    copiaVetor(b, r, n);

    // Aplica o pré-condicionador: z0 = M-1 * r0
    aplicaPreCond(A, r, z, n, omega, dlu);

    // A primeira direção de busca é o resíduo pré-condicionado: p0 = z0
    copiaVetor(z, p, n);
    
    // Pré-calcula o produto escalar rk, zk que será usado no loop
    real_t r_z_velho = produtoEscalar(r, z, n);
    real_t norma_erro_final = 0.0;
    int k_final = 0;

    // Loop principal
    for (int k = 0; k < maxit; ++k) {
        k_final = k + 1; // Guarda o número de iterações executadas
        
        // Se o resíduo é suficientemente pequeno, pára para não causar instabilidade na divisão
        if (fabs(r_z_velho) < ZERO_TOLERANCE) { 
            break; 
        }

        // Guarda o X da iteração anterior para calcular o erro no final
        copiaVetor(*X, x_velho, n);

        // Calcula o produto Matriz-Vetor: Ap = A * pk
        multMatVet(A, p, Ap, n);

        // Calcula o tamanho do passo: alpha_k = <rk, zk> / <pk, A*pk>
        real_t p_Ap = produtoEscalar(p, Ap, n);
        if (fabs(p_Ap) < ZERO_TOLERANCE) { // Evita divisão por zero se o método estagnar
            fprintf(stderr, "AVISO: Método estagnou (denominador de alfa próximo de zero).\n");
            break;
        }
        real_t alpha = r_z_velho / p_Ap;

        // Atualiza a solução e o resíduo:
        // x{k+1} = xk + alpha * pk
        // r{k+1} = rk - alpha * (A*pk)
        for (int i = 0; i < n; ++i) {
            (*X)[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        // Verifica o critério de parada: max(|xi - xi_velho|) < epsilon
        norma_erro_final = 0.0;
        for (int i = 0; i < n; ++i) {
            real_t erro_i = fabs((*X)[i] - x_velho[i]);
            if (erro_i > norma_erro_final) {
                norma_erro_final = erro_i;
            }
        }
        // Condição que leva em consideração se a solução parou de andar significativamente
        if (norma_erro_final < epsilon) {
            break; // Convergiu
        }
        
        // Aplica o pré-condicionador novamente: z{k+1} = M-1 * r{k+1}
        aplicaPreCond(A, r, z, n, omega, dlu);

        // Calcula o fator de correção: beta_k = <r{k+1}, z{k+1}> / <rk, zk>
        real_t r_z_novo = produtoEscalar(r, z, n);
        real_t beta = r_z_novo / r_z_velho;
        
        // Atualiza a nova direção de busca: p{k+1} = z{k+1} + beta * pk
        for (int i = 0; i < n; ++i) {
            p[i] = z[i] + beta * p[i];
        }
        
        // Atualiza o produto escalar para a próxima iteração
        r_z_velho = r_z_novo;
    }
    
    rtime_t t_fim = timestamp();
    
    // Caso o loop terminou por atingir o máximo de iterações
    if (k_final == maxit && norma_erro_final >= epsilon) {
        fprintf(stderr, "AVISO: Método não convergiu em %d iterações.\n", maxit);
    }
    
    // Calcula o tempo médio por iteração
    if (k_final > 0) {
        *tempo_iter = (t_fim - t_inicio) / k_final;
    } else {
        *tempo_iter = 0.0;
    }
    
    // Passa total de iterações
    if (out_total_iterations) {
        *out_total_iterations = k_final;
    }

    // Desaloca memória
    free(r);
    free(p);
    free(z);
    free(Ap);
    free(x_velho);
    
    // Retorna a norma do erro da última iteração calculada
    return norma_erro_final; 
}