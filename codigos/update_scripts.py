import os

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    if 'sys.argv' not in content:
        content = content.replace('M0_EARTH = 0.01', 'M0_EARTH = float(sys.argv[1]) if len(sys.argv) > 1 else 0.01')
    
    if 'analisis_pa3_10myr.py' in filepath or 'analisis_pa3_nuevos_alphas.py' in filepath:
        # Check if already replaced
        if 'FIG_DIR = f"data/figures/10Myr' in content:
            content = content.replace('FIG_DIR = f"data/figures/10Myr_M0_{M0_EARTH}"\n', '')
            content = content.replace('FIG_DIR = f"data/figures/10Myr_nuevos_alphas_M0_{M0_EARTH}"\n', '')
            content = content.replace('os.makedirs(FIG_DIR, exist_ok=True)\n', '')
            
            # Replace function signatures
            content = content.replace('def plot_lines_grouped_by_rgap(data_alpha, alpha_val):', 'def plot_lines_grouped_by_rgap(data_alpha, alpha_val, fig_dir):')
            content = content.replace('def plot_lines_grouped_by_mgap(data_alpha, alpha_val):', 'def plot_lines_grouped_by_mgap(data_alpha, alpha_val, fig_dir):')
            content = content.replace('def plot_heatmaps(data_alpha, alpha_val):', 'def plot_heatmaps(data_alpha, alpha_val, fig_dir):')
            
            # Replace savefig
            content = content.replace('os.path.join(FIG_DIR, ', 'os.path.join(fig_dir, ')
            
            # Replace main
            main_old = '''    print("\\nGenerando gráficos para cada alpha...")
    for alpha_val, data_alpha in data.items():
        print(f" -> Generando plots para alpha = {alpha_val} ({len(data_alpha)} corridas)")
        plot_lines_grouped_by_rgap(data_alpha, alpha_val)
        plot_lines_grouped_by_mgap(data_alpha, alpha_val)
        plot_heatmaps(data_alpha, alpha_val)
        
    print(f"\\n¡Todos los gráficos se han guardado en {FIG_DIR}!")'''
            
            main_new = '''    print("\\nGenerando gráficos para cada alpha...")
    for alpha_val, data_alpha in data.items():
        print(f" -> Generando plots para alpha = {alpha_val} ({len(data_alpha)} corridas)")
        fig_dir = f"data/figures/M_{M0_EARTH}/alpha_{alpha_val}"
        os.makedirs(fig_dir, exist_ok=True)
        plot_lines_grouped_by_rgap(data_alpha, alpha_val, fig_dir)
        plot_lines_grouped_by_mgap(data_alpha, alpha_val, fig_dir)
        plot_heatmaps(data_alpha, alpha_val, fig_dir)
        
    print(f"\\n¡Todos los gráficos se han guardado en data/figures/M_{M0_EARTH}!")'''
            
            content = content.replace(main_old, main_new)
        
    elif 'analisis_retrasado.py' in filepath:
        if '10Myr_retrasado_M0' in content:
            content = content.replace('FIG_DIR = f"data/figures/10Myr_retrasado_M0_{M0_EARTH}"', 'FIG_DIR = f"data/figures/M_{M0_EARTH}/delayed"')
        
    elif 'analisis_sinusoidal.py' in filepath:
        if '10Myr_Sinusoidal_M0' in content:
            content = content.replace('FIG_DIR = f"data/figures/10Myr_Sinusoidal_M0_{M0_EARTH}"', 'FIG_DIR = f"data/figures/M_{M0_EARTH}/sinusoidal"')

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

for f in ['analisis_pa3_10myr.py', 'analisis_pa3_nuevos_alphas.py', 'analisis_retrasado.py', 'analisis_sinusoidal.py']:
    process_file(f)
