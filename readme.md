TRABALHO 2: Otimização de Desempenho para Sistemas Lineares Esparsos
Disciplina: Introdução à Computação Científica
Autor: Pietro Comin (GRR: 20241955)

ESTRUTURA DE DIRETÓRIOS E ARQUIVOS

Este pacote contém dois diretórios principais:

Trabalho1_v2/

Contém o código original do Trabalho 1 (Matrizes Densas).

NOTA IMPORTANTE: A função de geração da matriz foi levemente otimizada
(mantendo a estrutura densa O(N^2)) apenas para permitir que a inicialização
do sistema concluísse em tempo hábil para N=20.000. O algoritmo do Solver
permaneceu inalterado (ingênuo).

Contém instrumentação LIKWID.

Trabalho2/

Contém a versão otimizada (Formato DIA + Vetorização AVX).

Contém os scripts de automação para coleta de dados e geração de gráficos.

Contém os resultados brutos (CSV) e os gráficos gerados (PNG).

Contém o Relatório em PDF.


COMO COMPILAR

Para compilar ambos os projetos, entre no diretório 'Trabalho2' e execute o
script de comparação ou compile manualmente cada pasta com 'make'.

Requisitos: GCC com suporte a AVX, Python 3, LIKWID.

COMO REPRODUZIR OS RESULTADOS

Todos os scripts devem ser executados a partir do diretório 'Trabalho2'.

Coleta de Tempos (Speedup):
Execute: ./comparar_versoes.sh

Compila T1 e T2.

Executa para N = {32..20000}.

Gera o arquivo: resultados_comparacao.csv

Coleta de Métricas de Hardware (LIKWID):
Execute: python3 coleta_final.py

Mede FLOPS_DP (AVX) e Cache Miss Ratio.

Gera o arquivo: resultados_likwid_final.csv

Nota: Nomes de grupos ajustados para arquitetura AMD EPYC (Zen).

Geração de Gráficos:
Execute: python3 gera_graficos.py

Lê os dois CSVs acima.

Gera as imagens na pasta 'graficos/'.


RESULTADOS INCLUÍDOS

Os testes já foram executados e os dados estão salvos em 'Trabalho2':

resultados_comparacao.csv (Dados de Tempo e Speedup)

resultados_likwid_final.csv (Dados de Performance e Cache)

graficos/*.png (6 gráficos utilizados no relatório)

relatorio.pdf (Documento final com a análise)


OBSERVAÇÕES

Para N >= 8000, a versão V1 (Trabalho1_v2) falha ao rodar sob instrumentação
LIKWID devido ao consumo excessivo de memória da representação densa (> 9GB),
resultando em métricas nulas nos gráficos para esses pontos. Isso é um
comportamento esperado que demonstra a limitação da abordagem não otimizada.