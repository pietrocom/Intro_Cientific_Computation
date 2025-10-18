#!/bin/bash

# ==============================================================================
# Script de Teste para o Trabalho 1 - cgSolver
# Autor: Pietro Comin
# GRR:   20241955
# ==============================================================================

# --- Configurações ---
EXECUTABLE="./cgSolver"
OUTPUT_DIR="test_results"

# --- Início do Script ---

# 1. Compila o projeto
echo "--> Compilando o projeto com 'make'..."
make

# 2. Prepara o ambiente de teste
echo "--> Preparando o ambiente de testes..."
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

echo "--> Iniciando a bateria de testes..."
echo "Resultados serão salvos no diretório: $OUTPUT_DIR/"
echo "----------------------------------------------------"

# --- Bateria de Testes ---

# Para cada teste executado, verifica-se se houve alguma saída de erro.

echo "Executando Teste 1: Sanity Check (n=11, k=3, sem PC)..."
echo "11 3 -1.0 100 0.000001" | $EXECUTABLE > "$OUTPUT_DIR/1_SanityCheck.out" 2> "$OUTPUT_DIR/1_SanityCheck.err"
if [ -s "$OUTPUT_DIR/1_SanityCheck.err" ]; then
    echo "  Resultado: FALHA (AVISO detectado)"
else
    echo "  Resultado: SUCESSO"
fi

echo "Executando Teste 2a: Difícil (n=500, k=5, sem PC)..."
echo "500 5 -1.0 1000 0.00000001" | $EXECUTABLE > "$OUTPUT_DIR/2a_HardCase.out" 2> "$OUTPUT_DIR/2a_HardCase.err"
if [ -s "$OUTPUT_DIR/2a_HardCase.err" ]; then
    echo "  Resultado: FALHA (AVISO detectado)"
else
    echo "  Resultado: SUCESSO"
fi

echo "Executando Teste 2b: Difícil (n=500, k=5, com Jacobi)..."
echo "500 5 0.0 1000 0.00000001" | $EXECUTABLE > "$OUTPUT_DIR/2b_HardCase_Jacobi.out" 2> "$OUTPUT_DIR/2b_HardCase_Jacobi.err"
if [ -s "$OUTPUT_DIR/2b_HardCase_Jacobi.err" ]; then
    echo "  Resultado: FALHA (AVISO detectado)"
else
    echo "  Resultado: SUCESSO"
fi

echo "Executando Teste 2c: Difícil (n=500, k=5, com Gauss-Seidel)..."
echo "500 5 1.0 1000 0.00000001" | $EXECUTABLE > "$OUTPUT_DIR/2c_HardCase_GaussSeidel.out" 2> "$OUTPUT_DIR/2c_HardCase_GaussSeidel.err"
if [ -s "$OUTPUT_DIR/2c_HardCase_GaussSeidel.err" ]; then
    echo "  Resultado: FALHA (AVISO detectado)"
else
    echo "  Resultado: SUCESSO"
fi

echo "Executando Teste 2d: Difícil (n=500, k=5, com SSOR)..."
echo "500 5 1.5 1000 0.00000001" | $EXECUTABLE > "$OUTPUT_DIR/2d_HardCase_SSOR.out" 2> "$OUTPUT_DIR/2d_HardCase_SSOR.err"
if [ -s "$OUTPUT_DIR/2d_HardCase_SSOR.err" ]; then
    echo "  Resultado: FALHA (AVISO detectado)"
else
    echo "  Resultado: SUCESSO"
fi

echo "Executando Teste 3: Stress Test (n=1500, k=7, com Jacobi)..."
echo "1500 7 0.0 3000 0.00000001" | $EXECUTABLE > "$OUTPUT_DIR/3_StressTest_Jacobi.out" 2> "$OUTPUT_DIR/3_StressTest_Jacobi.err"
if [ -s "$OUTPUT_DIR/3_StressTest_Jacobi.err" ]; then
    echo "  Resultado: FALHA (AVISO detectado)"
else
    echo "  Resultado: SUCESSO"
fi

echo "Executando Teste 4: Falha por Limite de Iterações (n=200)..."
echo "200 5 0.0 5 0.0000000001" | $EXECUTABLE > "$OUTPUT_DIR/4_Failure_Maxit.out" 2> "$OUTPUT_DIR/4_Failure_Maxit.err"
if [ -s "$OUTPUT_DIR/4_Failure_Maxit.err" ]; then
    echo "  Resultado: FALHA (AVISO detectado)"
else
    echo "  Resultado: SUCESSO"
fi

echo "----------------------------------------------------"
echo "Bateria de testes concluída."