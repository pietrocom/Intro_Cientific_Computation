#!/bin/bash

# Autor: Pietro Comin (pc24 - GRR20241955)

# NOTA: Este script representa a solução teoricamente correta para o problema.
# Em sistemas com incompatibilidade na Marker API (-m) do LIKWID, a segunda
# parte da saída (medição de desempenho) pode não ser gerada.

# Lembre-se de dar permissão de execução: chmod +x executa.sh

make 

# Mude o arquivo de entrada para ter mais casos de teste
./resolveEDO < teste.dat


# Roda e coleta dados de desempenho
sudo -E likwid-perfctr -C 0 -g FLOPS_DP -m -O ./resolveEDO < teste.dat 1>/dev/null 2>&1 | \
grep "FP_ARITH_INST_RETIRED_SCALAR_DOUBLE" | \
awk -F',' '{print $1","$3}'