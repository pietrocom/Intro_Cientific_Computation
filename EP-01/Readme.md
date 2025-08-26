EP-01 - Solução de Equações Não-lineares

Enunciado

1. Implemente uma biblioteca ZeroFuncao contendo as funções abaixo, especializadas em encontrar os zeros de polinômios, onde o parâmetro p é um polinômio do qual se deseja a raiz, os parâmetros a, b, x0 são os chutes iniciais para os métodos, it é o número de iterações do método que foram necessários para chegar à raiz, cujo valor deve ser devolvido pelo parâmetro raiz. O valor de retorno da função é o erro obtido para o valor final calculado para a raiz.

    double newtonRaphson (Polinomio p, double x0, int criterioParada, int *it, double *raiz)
    double bisseccao (Polinomio p, double a, double b, int criterioParada, int *it, double *raiz)

Também devem ser criadas as funções abaixo para o cálculo do valor px de um polinômio p e sua primeira derivada dpx em um ponto x.
Uma função é para o cálculo rápido como visto nos slides e notas de aula, e a outra é a implementação com o
cálculo direto do polinômio e sua derivada, usando a função pow().

    void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx)
    void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx)

Considere a seguinte definição para o tipo Polinomio:

typedef struct { double *p; int grau } Polinomio

onde p é um vetor de tamanho grau+1 contendo os coeficientes do polinômio.

Use os arquivos ZeroFuncao.c e ZeroFuncao.h que estão ao final deste enunciado como ponto de partida para a codificação desta biblioteca.

2. Faça o programa labZero que, utilizando as funções da biblioteca ZeroFuncao, calcule as raízes de um polinômio, cujo grau n, coeficientes p e intervalo [a,b] contendo as raízes são informados pela entrada padrão (stdin). 

Para o cálculo do valor de um polinômio (e de suas derivadas) em um ponto x, use a função calcPolinomio_rapido() e a função calcPolinomio_lento().

O programa deve calcular a raiz para cada um dos métodos na biblioteca ZeroFuncao, mostrando ao final na saída padrão (stdout) (para cada tipo de cálculo de polinômio)  3 (três)  linhas com as colunas

metodo, raiz, erro, iterações e tempo,

onde tempo é o tempo gasto pelo método. Cada linha está associada a um dos critérios de parada indicados abaixo.

 

Devem ser implementados 3 (três) critérios de parada (cada um gerando uma linha de resultados):

    Criterio-01:

;
Criterio-02:
;
Criterio-03: ULP's entre
e

Cada um dos valores destes critérios devem ser exibidos na linha de resultado na coluna erro.

Para efeitos deste programa, considere o número máximo de iterações como sendo 600.

O tempo de execução deve ser medido em milisegundos, considerando tempo de relógio, utilizando a função timestamp(), disponível no módulo utils. O tempo decorrido é medido pela diferença dos “timestamp” medidos antes e depois da região de interesse (veja descrição da função em utils.c).

ATENÇÃO: Os valores reais devem ser representados em ponto flutuante de precisão dupla (tipo double).

Todos os cálculos de ponto flutuante devem ser efetuados com arredondamento para baixo (use a função  fesetround ).

Os valores de saída devem ser mostrados com 15  casas decimais (formato "%.15e"). O valor do tempo de cálculo de deve ser exibido com 8 casas decimais (formato "%.8e").

ATENÇÃO: nada deve ser impresso além dos valores citados, pois os resultados serão verificados via script. Implementações que não respeitarem essa restrição receberão nota zero.

Exemplo de Execução do programa

Entrada: (2x⁴ + 4x³ + 3x² − 10x − 15 = 0, ξ ∈ [0, 3])

4
2 4 3 -10 -15
0 3

Saída:

RAPIDO

bissec  +1.492878690361976e+00 5.988897183511489e-08   26  1.51919983e-02
bissec  +1.492878708663603e+00 1.065814103640150e-14  501  2.18512002e-01
bissec  +1.492878708663603e+00 2.000000000000000e+00   53  2.51480005e-02
newton  +1.492878708663603e+00 1.280976475499041e-09    3  1.31199881e-03
newton  +1.492878708663603e+00 1.065814103640150e-14  500  1.09547000e-01
newton  +1.492878708663603e+00 0.000000000000000e+00    4  1.38799846e-03

LENTO

bissec  +1.492878690361976e+00 5.988897183511489e-08   26  6.01199977e-02
bissec  +1.492878708663603e+00 3.552713678800500e-15  501  1.51478400e+00
bissec  +1.492878708663603e+00 2.000000000000000e+00   53  2.03621000e-01
newton  +1.492878708663603e+00 1.280976475499041e-09    3  7.83599913e-03
newton  +1.492878708663603e+00 3.552713678800500e-15  500  9.57649998e-01
newton  +1.492878708663603e+00 0.000000000000000e+00    4  1.11569985e-02

 

O que deve ser entregue

Deve ser entregue um pacote de software completo contendo os fontes em linguagem C. O pacote deve ser arquivado e compactado com tar(1) e gzip(1) em um arquivo chamado login1.tgz , onde login1 é o login do aluno no DINF.

O pacote deve ter a seguinte estrutura de diretório e arquivos:

    ./login1/: diretório principal contendo o código-fonte;
    ./login1/Makefile;

Note que a extração dos arquivos de login1.tgz deve criar o diretório login1 contendo todos os arquivos acima. Os arquivos fonte também devem estar contidos no diretório, ou em algum sub-diretório, desde que o Makefile funcione.

Arquivos auxiliares 

O arquivo polinomios.dat apresenta vários exemplos de matrizes para inversão.  Observe  que o programa labZero deve ler APENAS 1 (UM) polinômio e  apresentar o resultado. O arquivo polinomios.dat  é apenas um conjunto de exemplos, não a entrada do programa.

Os arquivos utils.* contém a função timestamp(), a ser usada neste trabalho e nos próximos.

O arquivo Makefile é um exemplo que pode ser usado e já contém as regras all, testeFormato, clean, purge, e adicionalmente a regra dist, que pode ser usada para gerar o arquivo tgz a ser entregue.