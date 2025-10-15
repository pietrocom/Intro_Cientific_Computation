# EP-02: Resolução de EDOs com o Método de Gauss-Seidel

- **Autor:** Pietro Comin
- **GRR:** GRR20241955

## Descrição do Projeto

Este trabalho implementa uma solução em linguagem C para resolver uma família de Equações Diferenciais Ordinárias (EDOs) de segunda ordem com condições de 
contorno. A EDO genérica é discretizada utilizando o método de diferenças finitas, o que a transforma em um sistema linear tridiagonal. Este sistema é então 
resolvido iterativamente pelo método de Gauss-Seidel.

O projeto também inclui a instrumentação do código com a biblioteca LIKWID para permitir a medição de desempenho da rotina de solução numérica, conforme solicitado pelo enunciado.

## Estrutura do Projeto

O código foi organizado de forma modular para garantir clareza e manutenibilidade:

- `resolveEDO.c`: Módulo principal que contém a função `main`. É responsável por orquestrar a leitura dos dados de entrada, a chamada das rotinas de cálculo e a 
impressão dos resultados.
- `edo.c` / `edo.h`: Módulo que encapsula toda a lógica relacionada à EDO, incluindo a geração do sistema linear tridiagonal e a implementação do método de 
Gauss-Seidel.
- `utils.c` / `utils.h`: Módulo com funções utilitárias de propósito geral, como a medição de tempo (`timestamp`).
- `Makefile`: Arquivo de compilação com regras para gerar o executável (`all`), limpar arquivos temporários (`clean`) e remover todos os arquivos gerados 
(`purge`).
- `executa.sh`: Script final para compilar e executar o programa, gerando uma saída formatada idêntica à do enunciado.
- `teste.dat` / `edos.dat`: Arquivos com dados de entrada para teste.

## Compilação e Execução

O projeto pode ser compilado utilizando o `Makefile` fornecido.

**1. Para compilar o programa:**
```bash
make
```

**2. Para executar a análise completa e obter a saída formatada:**
O script `executa.sh` foi criado para automatizar todo o processo, gerando a saída numérica seguida pelas métricas de desempenho.
```bash
# Dar permissão de execução (apenas uma vez)
chmod +x executa.sh

# Executar o script
./executa.sh
```

## Desafios Encontrados e Metodologia

A implementação da solução numérica e a estruturação do código ocorreram conforme o planejado. O principal desafio do projeto surgiu durante a etapa de medição 
de desempenho no ambiente de desenvolvimento (CPU Intel Core i7 de 11ª Geração).

Após a correta instrumentação do código com a Marker API do LIKWID, a ferramenta `likwid-perfctr` se mostrou incapaz de gerar relatórios para as regiões de 
código marcadas (`-m`). Uma longa e aprofundada investigação foi realizada para diagnosticar a causa raiz. Este processo de depuração incluiu:

1.  A reinstalação completa da biblioteca LIKWID para corrigir uma aparente corrupção nos arquivos de configuração de grupos de desempenho.
2.  Testes com múltiplos grupos de medição, incluindo `FLOPS_DP` e o universal `CLOCK`.
3.  A exploração de diversas técnicas de script para contornar possíveis problemas de permissão e ambiente de execução com `sudo`, como o uso da flag `-E` e o 
encapsulamento de comandos.

A conclusão final foi que existe uma **incompatibilidade fundamental da Marker API do LIKWID** com a combinação específica de hardware, kernel e/ou BIOS do 
sistema, que impede a medição de regiões específicas de código.

Este processo de depuração para esgotar todas as possíveis causas do problema com a ferramenta LIKWID foi complexo e consumiu um tempo considerável, sendo a 
razão de eu demorar tanto para entregar o trabalho final.

O script `executa.sh` submetido representa a minha versão de implementação teoricamente correta e ideal, que funcionaria em um ambiente onde a Marker API do 
LIKWID é 
totalmente compatível.