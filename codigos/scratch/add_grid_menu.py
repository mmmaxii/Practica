with open('C:/astro/Codigos practica + docs + papers/codigos/analisis_interactivo.py', 'r', encoding='utf-8') as f:
    content = f.read()

new_function = """
def run_population_grid():
    clear_screen()
    print("=====================================================")
    print("        GRILLA CUANTITATIVA (Matriz Heatmap)         ")
    print("=====================================================")
    
    base_data = "data/runs"
    vfrags = ['1', '3', '10']
    
    import os, pickle, glob, re
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
    fig_path = os.path.join(out_dir, f"population_grid_{scenario_target}.png")
    
    print(f"Generando grilla cuantitativa...")
    PlotterPopulation.plot_grilla_cuantitativa(global_data, fig_path, scenario_name=scenario_target.capitalize())
    input("\\nPresione Enter para volver al menú principal...")

def main():
"""

if "def main():" in content:
    content = content.replace("def main():", new_function)

old_menu = """        v_frags_with_all = {
            '0': 'PROCESAR TODO (Lote automático)', 
            **v_frags, 
            '4': 'Money Plot (Población Sintética Global)',
            '5': 'Figura 4 (Física Oculta: f_post para casos estructurados)',
            '6': 'Panel Conceptual (Gráfico Único Representativo)',
            '7': 'Estadísticas Poblacionales (Barras Apiladas)'
        }
        v_choice = get_choice("Seleccione una opción de v_frag:", v_frags_with_all, allow_back=False)
        if v_choice == '4':
            run_population_synthesis(alphas_full, m0_map)
            continue
        elif v_choice == '5':
            run_fpost_analysis()
            continue
        elif v_choice == '6':
            run_population_single()
            continue
        elif v_choice == '7':
            run_population_bars()
            continue"""

new_menu = """        v_frags_with_all = {
            '0': 'PROCESAR TODO (Lote automático)', 
            **v_frags, 
            '4': 'Money Plot (Población Sintética Global)',
            '5': 'Figura 4 (Física Oculta: f_post para casos estructurados)',
            '6': 'Panel Conceptual (Gráfico Único Representativo)',
            '7': 'Estadísticas Poblacionales (Barras Apiladas)',
            '8': 'Grilla Cuantitativa (Matriz Numérica de Porcentajes)'
        }
        v_choice = get_choice("Seleccione una opción de v_frag:", v_frags_with_all, allow_back=False)
        if v_choice == '4':
            run_population_synthesis(alphas_full, m0_map)
            continue
        elif v_choice == '5':
            run_fpost_analysis()
            continue
        elif v_choice == '6':
            run_population_single()
            continue
        elif v_choice == '7':
            run_population_bars()
            continue
        elif v_choice == '8':
            run_population_grid()
            continue"""

if old_menu in content:
    content = content.replace(old_menu, new_menu)

with open('C:/astro/Codigos practica + docs + papers/codigos/analisis_interactivo.py', 'w', encoding='utf-8') as f:
    f.write(content)

print("Appended grid option to analisis_interactivo.py")
