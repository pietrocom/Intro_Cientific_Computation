#!/bin/bash

# ==============================================================================
# Script de Comparação de TEMPO (T1 vs T2)
# Gera: resultados_comparacao.csv
# ==============================================================================

# Força ponto decimal
export LC_NUMERIC=C

# --- Configurações ---
# Apontando para a pasta correta (V2)
DIR_T1="../Trabalho1_v2"
DIR_T2="."
EXEC_NAME="cgSolver"
RESULT_FILE="resultados_comparacao.csv"

# Tamanhos de N (Lista completa)
SIZES=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# Parâmetros
K=7
MAXIT=25
EPSILON="1e-8"
OMEGA="0.0" # Jacobi

# Cores
GREEN='\033[0;32m'
NC='\033[0m'

echo -e "${GREEN}===> Iniciando Comparação de TEMPO <===${NC}"

# 1. Compilação
echo "Compilando T1 ($DIR_T1)..."
make -C $DIR_T1 clean > /dev/null && make -C $DIR_T1 > /dev/null

echo "Compilando T2 ($DIR_T2)..."
make -C $DIR_T2 clean > /dev/null && make -C $DIR_T2 > /dev/null

# 2. Prepara CSV
echo "N,T1_Iter_ms,T2_Iter_ms,Speedup_Iter,T1_Res_ms,T2_Res_ms,Speedup_Res" > $RESULT_FILE
echo -e "N\t| T1 Iter\t| T2 Iter\t| Speedup"
echo "--------------------------------------------------------"

# 3. Loop
for N in "${SIZES[@]}"; do
    INPUT="$N $K $OMEGA $MAXIT $EPSILON"

    # Executa T1 (Pega as ultimas linhas para achar os tempos)
    OUT_T1=$(echo "$INPUT" | $DIR_T1/$EXEC_NAME | tail -n 5)
    T1_ITER=$(echo "$OUT_T1" | sed -n '4p')
    T1_RES=$(echo "$OUT_T1" | sed -n '5p')

    # Executa T2
    OUT_T2=$(echo "$INPUT" | $DIR_T2/$EXEC_NAME | tail -n 5)
    T2_ITER=$(echo "$OUT_T2" | sed -n '4p')
    T2_RES=$(echo "$OUT_T2" | sed -n '5p')

    # Calcula Speedup
    SPEEDUP_ITER=$(awk -v t1="$T1_ITER" -v t2="$T2_ITER" 'BEGIN { if(t2==0) print 0; else print t1/t2 }')
    SPEEDUP_RES=$(awk -v t1="$T1_RES" -v t2="$T2_RES" 'BEGIN { if(t2==0) print 0; else print t1/t2 }')

    # Imprime e Salva
    printf "%d\t| %.4f\t| %.4f\t| ${GREEN}%.2fx${NC}\n" $N $T1_ITER $T2_ITER $SPEEDUP_ITER
    echo "$N,$T1_ITER,$T2_ITER,$SPEEDUP_ITER,$T1_RES,$T2_RES,$SPEEDUP_RES" >> $RESULT_FILE
done

echo "--------------------------------------------------------"
echo "CSV Gerado: $RESULT_FILE"