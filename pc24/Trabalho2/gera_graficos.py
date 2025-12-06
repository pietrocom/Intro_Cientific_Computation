import pandas as pd
import matplotlib.pyplot as plt
import os

# Nomes dos arquivos de entrada
FILE_TEMPO = "resultados_comparacao.csv"
FILE_LIKWID = "resultados_likwid_final.csv"

# Cria pasta para salvar imagens
if not os.path.exists("graficos"):
    os.mkdir("graficos")

def plot_graph(df, metric_col_v1, metric_col_v2, ylabel, title, filename, log_y=False):
    plt.figure(figsize=(10, 6))
    
    # Eixo X deve ser logarítmico conforme enunciado
    plt.xscale('log')
    if log_y:
        plt.yscale('log')

    plt.plot(df['N'], df[metric_col_v1], marker='o', linestyle='--', label='V1 (Naive)')
    plt.plot(df['N'], df[metric_col_v2], marker='s', linestyle='-', label='V2 (Otimizado)')
    
    plt.xlabel('Tamanho do Sistema (N) - Escala Log')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)
    
    path = f"graficos/{filename}.png"
    plt.savefig(path)
    print(f"Gerado: {path}")
    plt.close()

def main():
    # --- 1. GRÁFICOS DE TEMPO ---
    if os.path.exists(FILE_TEMPO):
        print("Lendo dados de tempo...")
        try:
            df_tempo = pd.read_csv(FILE_TEMPO, decimal=',')
            # Se falhar conversão, tenta decimal ponto
            if df_tempo['T1_Iter_ms'].dtype == 'object':
                 df_tempo = pd.read_csv(FILE_TEMPO, decimal='.')
        except:
            df_tempo = pd.read_csv(FILE_TEMPO)

        # Gráfico Tempo Iteração (Log-Log)
        plot_graph(df_tempo, 'T1_Iter_ms', 'T2_Iter_ms', 
                   'Tempo Médio (ms)', 
                   'Tempo Médio por Iteração (OP1)', 
                   'tempo_iteracao', log_y=True)

        # Gráfico Tempo Resíduo (Log-Log)
        plot_graph(df_tempo, 'T1_Res_ms', 'T2_Res_ms', 
                   'Tempo Total (ms)', 
                   'Tempo Cálculo Resíduo (OP2)', 
                   'tempo_residuo', log_y=True)
    else:
        print(f"AVISO: {FILE_TEMPO} não encontrado.")

    # --- 2. GRÁFICOS DE LIKWID ---
    if os.path.exists(FILE_LIKWID):
        print("Lendo dados do LIKWID...")
        df_lik = pd.read_csv(FILE_LIKWID)

        # Filtra por Grupo
        flops = df_lik[df_lik['Grupo'] == 'FLOPS_DP']
        cache = df_lik[df_lik['Grupo'] == 'CACHE']

        # Gráfico FLOPS - Iteração
        plot_graph(flops, 'T1_OP1', 'T2_OP1', 
                   'MFLOP/s', 
                   'Performance Aritmética - Iteração (OP1)', 
                   'flops_iteracao')

        # Gráfico FLOPS - Resíduo
        plot_graph(flops, 'T1_OP2', 'T2_OP2', 
                   'MFLOP/s', 
                   'Performance Aritmética - Resíduo (OP2)', 
                   'flops_residuo')

        # Gráfico Cache Miss - Iteração
        plot_graph(cache, 'T1_OP1', 'T2_OP1', 
                   'Miss Ratio', 
                   'Cache Miss Ratio - Iteração (OP1)', 
                   'cache_iteracao')

        # Gráfico Cache Miss - Resíduo
        plot_graph(cache, 'T1_OP2', 'T2_OP2', 
                   'Miss Ratio', 
                   'Cache Miss Ratio - Resíduo (OP2)', 
                   'cache_residuo')
    else:
        print(f"AVISO: {FILE_LIKWID} não encontrado.")

if __name__ == "__main__":
    main()
