import subprocess
import csv
import os
import sys

# --- CONFIGURAÇÕES ---
EXEC_T1 = "../Trabalho1_v2/cgSolver"
EXEC_T2 = "./cgSolver"

# Todos os tamanhos solicitados no enunciado
SIZES = [32, 64, 128, 256, 512, 1000, 2000, 4000, 8000, 9000, 10000, 20000]

# Métricas disponíveis na máquina (AMD EPYC)
METRICS = {
    "FLOPS_DP": "DP MFLOP/s (AVX assumed)", 
    "CACHE": "miss ratio"
}

# Parâmetros do Trabalho
K = 7
OMEGA = 0.0
MAXIT = 25
EPSILON = 1e-8
OUTPUT_FILE = "resultados_likwid_final.csv"

def run_likwid(executable, n, group):
    input_str = f"{n} {K} {OMEGA} {MAXIT} {EPSILON}\n"
    # Executa com flag -O para saída em CSV
    cmd = ["likwid-perfctr", "-C", "0", "-g", group, "-O", "-m", executable]
    try:
        result = subprocess.run(cmd, input=input_str, capture_output=True, text=True, errors='ignore')
        return result.stdout
    except Exception:
        return ""

def parse_likwid_csv(output, region, keyword):
    if not output: return "0.0"

    # Localiza a tabela da região
    marker = f"TABLE,Region {region}"
    start_idx = output.find(marker)
    if start_idx == -1: return "0.0"

    block = output[start_idx:]
    next_idx = block.find("Region ", 10)
    if next_idx != -1:
        block = block[:next_idx]

    # Busca o valor na segunda coluna do CSV
    for line in block.splitlines():
        if keyword in line:
            parts = line.split(",")
            if len(parts) >= 2:
                try:
                    return str(float(parts[1]))
                except ValueError:
                    continue
    return "0.0"

def main():
    print("==> Iniciando Coleta LIKWID (Lista Completa)...")
    
    # Recompilação silenciosa
    os.system("make -C ../Trabalho1_v2 clean > /dev/null && make -C ../Trabalho1_v2 > /dev/null")
    os.system("make -C . clean > /dev/null && make -C . > /dev/null")

    # Prepara CSV
    header = ["N", "Grupo", "Metrica", "T1_OP1", "T1_OP2", "T2_OP1", "T2_OP2"]
    results = []

    for n in SIZES:
        print(f"Processando N = {n}...")
        for group, keyword in METRICS.items():
            
            # Coleta T1
            out_t1 = run_likwid(EXEC_T1, n, group)
            v_t1_op1 = parse_likwid_csv(out_t1, "OP1_ITERACAO", keyword)
            v_t1_op2 = parse_likwid_csv(out_t1, "OP2_RESIDUO", keyword)
            
            # Coleta T2
            out_t2 = run_likwid(EXEC_T2, n, group)
            v_t2_op1 = parse_likwid_csv(out_t2, "OP1_ITERACAO", keyword)
            v_t2_op2 = parse_likwid_csv(out_t2, "OP2_RESIDUO", keyword)
            
            results.append([n, group, group, v_t1_op1, v_t1_op2, v_t2_op1, v_t2_op2])

    # Salva
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(results)
    
    print(f"==> Coleta finalizada. Arquivo: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
