with open('C:/astro/Codigos practica + docs + papers/codigos/analisis_interactivo.py', 'r', encoding='utf-8') as f:
    content = f.read()

import re

# We will replace the run_population_grid function
old_func_start = "def run_population_grid():"
old_func_end = "def main():"

parts = content.split(old_func_start)
if len(parts) == 2:
    subparts = parts[1].split(old_func_end)
    if len(subparts) == 2:
        new_function = """def run_population_grid():
    clear_screen()
    print("=====================================================")
    print("        GRILLA CUANTITATIVA (Matriz Heatmap)         ")
    print("=====================================================")
    
    base_data = "data/runs"
    vfrags = ['1', '3', '10']
    
    import os, pickle, glob, re
    import numpy as np
    scenarios_found = set()
    exclude_scenarios = ['rounded']
    for vf in vfrags:
        vf_path = os.path.join(base_data, f"vf_{vf}ms")
        if os.path.isdir(vf_path):
            for s in os.listdir(vf_path):
                if os.path.isdir(os.path.join(vf_path, s)) and s not in exclude_scenarios:
                    scenarios_found.add(s)
                    
    if not scenarios_found: return
    s_list = sorted(list(scenarios_found))
    s_options = {str(i+1): s for i, s in enumerate(s_list)}
    s_choice = get_choice("Seleccione el escenario:", s_options, allow_back=True)
    if s_choice == '0': return
    scenario_target = s_options[s_choice]
    
    global_data = []
    print(f"\\nRecolectando cachés para el escenario: {scenario_target}...")
    for vf in vfrags:
        v_val = float(vf)
        scen_path = os.path.join(base_data, f"vf_{vf}ms", scenario_target)
        if not os.path.isdir(scen_path): continue
        for cf in glob.glob(os.path.join(scen_path, "**/cache_*.pkl"), recursive=True):
            try:
                with open(cf, 'rb') as f:
                    data_list = pickle.load(f)
                    for d in data_list:
                        d['vfrag'] = v_val
                        if 'm0_earth' not in d:
                            m_match = re.search(r'_M0_([0-9\\.e\\-]+)\\.pkl', cf)
                            if m_match: d['m0_earth'] = float(m_match.group(1))
                        if 'alpha' not in d:
                            a_match = re.search(r'_a_([0-9\\.e\\-]+)_', cf)
                            if a_match: d['alpha'] = float(a_match.group(1))
                        global_data.append(d)
            except: pass
                
    if len(global_data) == 0: return
    
    from pipeline_analysis.plotters import PlotterPopulation
    out_dir = "data/figures/benchmarks"
    os.makedirs(out_dir, exist_ok=True)
    
    # Extraemos todos los valores únicos de masa inicial
    m0_values = sorted(list(set([d.get('m0_earth', np.nan) for d in global_data if not np.isnan(d.get('m0_earth', np.nan))])))
    
    # 1. Generar caso General
    print(f"\\nGenerando grilla cuantitativa general (todos los M0)...")
    fig_path_general = os.path.join(out_dir, f"population_grid_{scenario_target}_general.png")
    PlotterPopulation.plot_grilla_cuantitativa(global_data, fig_path_general, scenario_name=f"{scenario_target.capitalize()} (General)")
    
    # 2. Generar caso por cada M0
    for m0 in m0_values:
        print(f"Generando grilla para M0 = {m0:g} M_E...")
        subset_data = [d for d in global_data if 'm0_earth' in d and np.isclose(d['m0_earth'], m0)]
        if not subset_data: continue
        
        m0_str = f"{m0:g}".replace('.', 'p')
        fig_path_m0 = os.path.join(out_dir, f"population_grid_{scenario_target}_M0_{m0_str}.png")
        PlotterPopulation.plot_grilla_cuantitativa(subset_data, fig_path_m0, scenario_name=f"{scenario_target.capitalize()} (M0 = {m0:g} M_E)")
        
    print("\\n¡Todas las grillas generadas con éxito!")
    input("Presione Enter para volver al menú principal...")

"""
        new_content = parts[0] + new_function + "def main():" + subparts[1]
        
        with open('C:/astro/Codigos practica + docs + papers/codigos/analisis_interactivo.py', 'w', encoding='utf-8') as f:
            f.write(new_content)
        print("Updated analisis_interactivo.py successfully.")
    else:
        print("Error subparts")
else:
    print("Error parts")
