import subprocess
import re
import csv
import os
import sys

# --- CONFIGURAÇÕES ---
EXEC_T1 = "../Trabalho1_v2/cgSolver"
EXEC_T2 = "./cgSolver"

# Se o T1 demorar demais, reduza essa lista
SIZES = [32, 64, 128, 256, 512, 1000, 2000, 4000, 8000, 10000] 

METRICS = {
    "FLOPS_DP": "MFLOP/s",          # Ajuste conforme likwid-perfctr -a
    "L2CACHE": "L2 miss ratio",
    "L3": "Memory bandwidth"
}

K = 7
OMEGA = 0.0
MAXIT = 25
EPSILON = 1e-8
OUTPUT_FILE = "resultados_likwid_final.csv"

def run_likwid(executable, n, group):
    input_str = f"{n} {K} {OMEGA} {MAXIT} {EPSILON}\n"
    cmd = ["likwid-perfctr", "-C", "0", "-g", group, "-m", executable]
    
    try:
        result = subprocess.run(cmd, input=input_str, capture_output=True, text=True, errors='ignore')
        return result.stdout
    except FileNotFoundError:
        print("[ERRO] likwid-perfctr não encontrado.")
        sys.exit(1)

def parse_likwid_output(output, region, metric_keyword):
    if not output: return "0.0"
    regions = output.split("Region ")
    target_region_text = ""
    for r in regions:
        if r.startswith(region):
            target_region_text = r
            break
    if not target_region_text: return "0.0"

    for line in target_region_text.splitlines():
        if metric_keyword in line:
            matches = re.findall(r"[\d]+\.[\d]+", line)
            if matches: return matches[0]
    return "0.0"

def main():
    print(f"==> Iniciando Coleta LIKWID...")
    # Limpa e Recompila
    os.system("make -C ../Trabalho1_v2 clean > /dev/null && make -C ../Trabalho1_v2 > /dev/null")
    os.system("make -C . clean > /dev/null && make -C . > /dev/null")

    results = [["N", "Grupo", "Metrica", "T1_OP1", "T1_OP2", "T2_OP1", "T2_OP2"]]

    for n in SIZES:
        print(f"\n--- N = {n} ---")
        for group, metric_name in METRICS.items():
            print(f"  > Grupo {group}...", end=" ", flush=True)
            
            out_t1 = run_likwid(EXEC_T1, n, group)
            v_t1_op1 = parse_likwid_output(out_t1, "OP1_ITERACAO", metric_name)
            v_t1_op2 = parse_likwid_output(out_t1, "OP2_RESIDUO", metric_name)
            
            out_t2 = run_likwid(EXEC_T2, n, group)
            v_t2_op1 = parse_likwid_output(out_t2, "OP1_ITERACAO", metric_name)
            v_t2_op2 = parse_likwid_output(out_t2, "OP2_RESIDUO", metric_name)
            
            results.append([n, group, metric_name, v_t1_op1, v_t1_op2, v_t2_op1, v_t2_op2])
            print("OK")

    with open(OUTPUT_FILE, "w", newline="") as f:
        csv.writer(f).writerows(results)
    print(f"\n==> Arquivo salvo: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()