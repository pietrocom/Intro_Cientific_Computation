#!/bin/bash

# ==============================================================================
# Script de Comparação de Desempenho (T1 vs T2)
# Autor: Auxiliar de IA (para Pietro Comin)
# ==============================================================================

# --- CORREÇÃO DE LOCALE ---
# Força o shell a usar '.' como separador decimal para evitar erros no printf
export LC_NUMERIC=C

# --- Configurações ---
DIR_T1="../Trabalho1"
DIR_T2="."
EXEC_NAME="cgSolver"
RESULT_FILE="resultados_comparacao.csv"

# Tamanhos de N conforme enunciado
SIZES=(32 64 128 256 512 1000 2000 4000 8000 9000 10000 20000)

# Parâmetros fixos do enunciado
K=7
MAXIT=25
EPSILON="1e-8"
OMEGA="0.0" # Jacobi

# Cores para saída no terminal
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}===> Iniciando Comparação T1 (Naive) vs T2 (Otimizado) <===${NC}"

# 1. Compilação
echo "------------------------------------------------------------"
echo "Compilando Trabalho 1..."
make -C $DIR_T1 clean > /dev/null
make -C $DIR_T1 > /dev/null
if [ $? -ne 0 ]; then
    echo -e "${RED}ERRO: Falha ao compilar Trabalho 1${NC}"
    exit 1
fi

echo "Compilando Trabalho 2..."
make -C $DIR_T2 clean > /dev/null
make -C $DIR_T2 > /dev/null
if [ $? -ne 0 ]; then
    echo -e "${RED}ERRO: Falha ao compilar Trabalho 2${NC}"
    exit 1
fi
echo "------------------------------------------------------------"

# 2. Prepara CSV
# Colunas: N, TempoIter_T1, TempoIter_T2, Speedup_Iter, TempoRes_T1, TempoRes_T2, Speedup_Res
echo "N,T1_Iter_ms,T2_Iter_ms,Speedup_Iter,T1_Res_ms,T2_Res_ms,Speedup_Res" > $RESULT_FILE
echo -e "N\t| T1 Iter\t| T2 Iter\t| Speedup\t| T1 Res\t| T2 Res\t| Speedup"
echo "------------------------------------------------------------------------------------------------"

# 3. Loop de Testes
for N in "${SIZES[@]}"; do
    # Monta a string de entrada
    INPUT="$N $K $OMEGA $MAXIT $EPSILON"

    # --- Executa T1 ---
    # Captura as últimas 5 linhas da saída (onde estão os tempos) para evitar processar o vetor gigante X
    OUT_T1=$(echo "$INPUT" | $DIR_T1/$EXEC_NAME | tail -n 5)
    
    # Extrai T_Iter (linha 4 do tail) e T_Res (linha 5 do tail)
    T1_ITER=$(echo "$OUT_T1" | sed -n '4p')
    T1_RES=$(echo "$OUT_T1" | sed -n '5p')

    # --- Executa T2 ---
    OUT_T2=$(echo "$INPUT" | $DIR_T2/$EXEC_NAME | tail -n 5)
    
    T2_ITER=$(echo "$OUT_T2" | sed -n '4p')
    T2_RES=$(echo "$OUT_T2" | sed -n '5p')

    # --- Cálculos (usando awk para float) ---
    # Speedup = Tempo Antigo / Tempo Novo
    # Proteção contra divisão por zero se o tempo for muito pequeno (0.0)
    # A precisão do awk é mantida internamente
    
    SPEEDUP_ITER=$(awk -v t1="$T1_ITER" -v t2="$T2_ITER" 'BEGIN { if(t2==0 || t2=="") print 0; else print t1/t2 }')
    SPEEDUP_RES=$(awk -v t1="$T1_RES" -v t2="$T2_RES" 'BEGIN { if(t2==0 || t2=="") print 0; else print t1/t2 }')

    # Formatação para o terminal (tabs)
    # Agora que LC_NUMERIC=C, o printf entenderá os pontos decimais.
    printf "%d\t| %.4f\t| %.4f\t| ${GREEN}%.2fx${NC}\t\t| %.4f\t| %.4f\t| ${GREEN}%.2fx${NC}\n" \
           "$N" "$T1_ITER" "$T2_ITER" "$SPEEDUP_ITER" "$T1_RES" "$T2_RES" "$SPEEDUP_RES"

    # Salva no CSV
    echo "$N,$T1_ITER,$T2_ITER,$SPEEDUP_ITER,$T1_RES,$T2_RES,$SPEEDUP_RES" >> $RESULT_FILE

done

echo "------------------------------------------------------------"
echo -e "${GREEN}Concluído! Resultados salvos em: $RESULT_FILE${NC}"