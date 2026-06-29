import os
import sys

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

def print_header(breadcrumbs):
    clear_screen()
    print("=====================================================")
    print("        PIPELINE DE ANÁLISIS DE ACRECIÓN (PA3)       ")
    print("=====================================================")
    if breadcrumbs:
        print(f"Ruta actual: {' > '.join(breadcrumbs)}")
    print("-----------------------------------------------------")

def get_choice(prompt, options, allow_back=True):
    while True:
        print("\n" + prompt)
        for key, val in options.items():
            print(f"[{key}] {val}")
        if allow_back:
            print("[0] Volver atrás")
            
        choice = input("\nSeleccione una opción: ").strip()
        if allow_back and choice == '0':
            return '0'
        if choice in options:
            return choice
        print("Opción inválida. Intente de nuevo.")


def run_population_synthesis(alphas_full, m0_map):
    clear_screen()
    print("=====================================================")
    print("        MONEY PLOT (Población Sintética Global)      ")
    print("=====================================================")
    print("Generando mapa de fase condensado M_final vs f_H2O")
    print("cruzando todas las velocidades de fragmentación.")
    print("-----------------------------------------------------")
    
    base_data = "data/runs"
    vfrags = ['1', '3', '10']
    
    import os
    scenarios_found = set()
    exclude_scenarios = ['rounded']
    for vf in vfrags:
        vf_path = os.path.join(base_data, f"vf_{vf}ms")
        if os.path.isdir(vf_path):
            for s in os.listdir(vf_path):
                if os.path.isdir(os.path.join(vf_path, s)) and s not in exclude_scenarios:
                    scenarios_found.add(s)
                    
    if not scenarios_found:
        print("No se encontraron escenarios.")
        input("Presione Enter...")
        return
        
    s_list = sorted(list(scenarios_found))
    s_options = {str(i+1): s for i, s in enumerate(s_list)}
    
    s_choice = get_choice("Seleccione el escenario base a analizar:", s_options, allow_back=True)
    if s_choice == '0': return
    
    scenario_target = s_options[s_choice]
    
    import pickle
    import glob
    import re
    import pandas as pd
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró la tabla compilada en {csv_path}.")
        input("Presione Enter...")
        return
        
    df_master = pd.read_csv(csv_path, encoding='utf-8')
    global_data = df_master[df_master['scenario'].str.lower() == scenario_target.lower()].copy()
    
    print(f"Total de simulaciones cargadas: {len(global_data)}")
    if global_data.empty:
        input("Presione Enter...")
        return
        
    from pipeline_analysis.plotters import PlotterPopulation
    
    out_dir = "data/figures/benchmarks/money_plot"
    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, f"population_{scenario_target}.png")
    
    print(f"Generando gráfica usando PlotterPopulation...")
    
    print("\n¿Desea usar escala logarítmica para el eje X (Fracción de agua)?")
    log_x_input = input("[S/N] (Por defecto S): ").strip().lower()
    log_x = False if log_x_input == 'n' else True

    PlotterPopulation.plot_poblacion_sintetica(global_data, fig_path, scenario_name=scenario_target.capitalize(), mass_threshold=0.1, water_threshold=0.10, log_x=log_x)
    
    input("\nPresione Enter para volver al menú principal...")

def run_fpost_analysis():
    clear_screen()
    print("=====================================================")
    print("        FIGURA 4 (FÍSICA OCULTA: f_post)             ")
    print("=====================================================")
    
    base_data = "data/runs"
    vfrags = ['1', '3', '10']
    
    import os
    scenarios_found = set()
    exclude_scenarios = ['rounded']
    for vf in vfrags:
        vf_path = os.path.join(base_data, f"vf_{vf}ms")
        if os.path.isdir(vf_path):
            for s in os.listdir(vf_path):
                if os.path.isdir(os.path.join(vf_path, s)) and s not in exclude_scenarios:
                    scenarios_found.add(s)
                    
    if not scenarios_found:
        print("No se encontraron escenarios estructurados.")
        input("Presione Enter...")
        return
        
    s_list = sorted(list(scenarios_found))
    s_options = {str(i+1): s for i, s in enumerate(s_list)}
    
    s_choice = get_choice("Seleccione el escenario estructurado a analizar:", s_options, allow_back=True)
    if s_choice == '0': return
    
    scenario_target = s_options[s_choice]
    
    import pickle
    import glob
    import re
    import pandas as pd
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró la tabla compilada en {csv_path}.")
        input("Presione Enter...")
        return
        
    df_master = pd.read_csv(csv_path, encoding='utf-8')
    global_data = df_master[df_master['scenario'].str.lower() == scenario_target.lower()].copy()
    
    print(f"Total de simulaciones cargadas: {len(global_data)}")
    if global_data.empty:
        input("Presione Enter...")
        return
        
    from pipeline_analysis.plotters import PlotterPopulation
    
    out_dir = "data/figures/benchmarks/f_post"
    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, f"population_fpost_{scenario_target}.png")
    
    print(f"Generando gráfica de f_post usando PlotterPopulation...")
    PlotterPopulation.plot_poblacion_sintetica(global_data, fig_path, scenario_name=scenario_target.capitalize(), mass_threshold=0.1, water_threshold=0.10, color_metric='f_post')
    
    input("\nPresione Enter para volver al menú principal...")



def run_population_single():
    clear_screen()
    print("=====================================================")
    print("        PANEL CONCEPTUAL (Gráfico Único)             ")
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
    
    import pandas as pd
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró la tabla compilada en {csv_path}.")
        input("Presione Enter...")
        return
        
    df_master = pd.read_csv(csv_path, encoding='utf-8')
    global_data = df_master[df_master['scenario'].str.lower() == scenario_target.lower()].copy()
    if global_data.empty:
        print(f"No hay datos compilados para el escenario {scenario_target}.")
        input("Presione Enter...")
        return
    
    from pipeline_analysis.plotters import PlotterPopulation
    out_dir = "data/figures/benchmarks/panel_conceptual"
    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, f"population_single_{scenario_target}.png")
    
    print(f"Generando panel conceptual...")
    PlotterPopulation.plot_poblacion_sintetica_single(global_data, fig_path, scenario_name=scenario_target.capitalize(), mass_threshold=0.1)
    input("\nPresione Enter para volver al menú principal...")

def run_population_facet_scatter():
    clear_screen()
    print("=====================================================")
    print("     DISPERSIÓN EN FACETAS (Masa vs. H2O por Alpha)  ")
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
    
    import pandas as pd
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró la tabla compilada en {csv_path}.")
        input("Presione Enter...")
        return
        
    df_master = pd.read_csv(csv_path, encoding='utf-8')
    global_data = df_master[df_master['scenario'].str.lower() == scenario_target.lower()].copy()
    if global_data.empty:
        print(f"No hay datos compilados para el escenario {scenario_target}.")
        input("Presione Enter...")
        return
        
    m0_options = {
        '1': '1e-5 (0.00001 M_earth)',
        '2': '1e-4 (0.0001 M_earth)',
        '3': '1e-3 (0.001 M_earth)',
        '4': '1e-2 (0.01 M_earth)',
        '5': '1e-1 (0.1 M_earth)'
    }
    m0_map = {'1': 0.00001, '2': 0.0001, '3': 0.001, '4': 0.01, '5': 0.1}

    print("\nOpciones de masa inicial del embrión (M0):")
    for k, v in m0_options.items():
        print(f"[{k}] {v}")
    print("\nIngrese los números de las masas que desea EXCLUIR separados por coma (ej: 1, 2, 3).")
    print("Deje en blanco para incluir todas las masas.")
    excluir_input = input("Excluir: ").strip()
    
    if excluir_input:
        excluir_keys = [k.strip() for k in excluir_input.split(',')]
        excluir_vals = [m0_map[k] for k in excluir_keys if k in m0_map]
        if excluir_vals:
            for val in excluir_vals:
                global_data = global_data[~np.isclose(global_data['M_emb0'], val)]
                
    if global_data.empty:
        print("Todos los datos fueron excluidos. No hay nada que graficar.")
        input("Presione Enter...")
        return
    
    from pipeline_analysis.plotters import PlotterPopulation
    out_dir = "data/figures/benchmarks/facet_scatter"
    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, f"population_facet_scatter_{scenario_target}.png")
    
    print(f"Generando gráfico de dispersión en facetas...")
    PlotterPopulation.plot_poblacion_facet_scatter(global_data, fig_path, scenario_name=scenario_target.capitalize(), mass_threshold=0.1)
    input("\nPresione Enter para volver al menú principal...")


def run_population_bars():
    clear_screen()
    print("=====================================================")
    print("        ESTADÍSTICAS POBLACIONALES (Barras)          ")
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
    
    import pandas as pd
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró la tabla compilada en {csv_path}.")
        input("Presione Enter...")
        return
        
    df_master = pd.read_csv(csv_path, encoding='utf-8')
    global_data = df_master[df_master['scenario'].str.lower() == scenario_target.lower()].copy()
    if global_data.empty:
        print(f"No hay datos compilados para el escenario {scenario_target}.")
        input("Presione Enter...")
        return
    
    from pipeline_analysis.plotters import PlotterPopulation
    out_dir = "data/figures/benchmarks/estadisticas_poblacionales"
    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, f"population_bars_{scenario_target}.png")
    
    print(f"Generando barras poblacionales...")
    PlotterPopulation.plot_poblacion_barras_apiladas(global_data, fig_path, scenario_name=scenario_target.capitalize())
    input("\nPresione Enter para volver al menú principal...")


def run_population_grid():
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
    
    import pandas as pd
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró la tabla compilada en {csv_path}.")
        input("Presione Enter...")
        return
        
    df_master = pd.read_csv(csv_path, encoding='utf-8')
    global_data = df_master[df_master['scenario'].str.lower() == scenario_target.lower()].copy()
    if global_data.empty:
        print(f"No hay datos compilados para el escenario {scenario_target}.")
        input("Presione Enter...")
        return
    
    from pipeline_analysis.plotters import PlotterPopulation
    out_dir = "data/figures/benchmarks/grilla_cuantitativa"
    os.makedirs(out_dir, exist_ok=True)
    
    m0_values = sorted(global_data['M_emb0'].dropna().unique())
    
    # 1. Generar caso General
    print(f"\nGenerando grilla cuantitativa general (todos los M0)...")
    fig_path_general = os.path.join(out_dir, f"population_grid_{scenario_target}_general.png")
    PlotterPopulation.plot_grilla_cuantitativa(global_data, fig_path_general, scenario_name=f"{scenario_target.capitalize()} (General)")
    
    # 2. Generar caso por cada M0
    for m0 in m0_values:
        print(f"Generando grilla para M0 = {m0:g} M_E...")
        subset_data = global_data[np.isclose(global_data['M_emb0'], m0)].copy()
        if subset_data.empty: continue
        
        m0_str = f"{m0:g}".replace('.', 'p')
        fig_path_m0 = os.path.join(out_dir, f"population_grid_{scenario_target}_M0_{m0_str}.png")
        PlotterPopulation.plot_grilla_cuantitativa(subset_data, fig_path_m0, scenario_name=f"{scenario_target.capitalize()} (M0 = {m0:g} M_E)", save_main=False)
        
    print("\n¡Todas las grillas generadas con éxito!")
    input("Presione Enter para volver al menú principal...")

def main():



    v_frags = {'1': '1 m/s', '2': '3 m/s', '3': '10 m/s'}
    v_frag_map = {'1': 1, '2': 3, '3': 10}
    
    alphas_full = {
        '1': '0.0001', '2': '0.0003', '3': '0.0005', '4': '0.0007',
        '5': '0.001', '6': '0.003', '7': '0.005', '8': '0.01'
    }
    
    m0_options = {
        '1': '1e-5 (0.00001 M_earth)',
        '2': '1e-4 (0.0001 M_earth)',
        '3': '1e-3 (0.001 M_earth)',
        '4': '1e-2 (0.01 M_earth)',
        '5': '1e-1 (0.1 M_earth)'
    }
    m0_map = {'1': 0.00001, '2': 0.0001, '3': 0.001, '4': 0.01, '5': 0.1}
    m0_names = {'1': 'Memb_e-5', '2': 'Memb_e-4', '3': 'Memb_e-3', '4': 'Memb_e-2', '5': 'Memb_e-1'}

    while True:
        breadcrumbs = []
        print_header(breadcrumbs)
        v_frags_with_all = {
            '0': 'PROCESAR TODO (Lote automático)', 
            **v_frags, 
            '4': 'Money Plot (Población Sintética Global)',
            '5': 'Figura 4 (Física Oculta: f_post para casos estructurados)',
            '6': 'Panel Conceptual (Gráfico Único Representativo)',
            '7': 'Estadísticas Poblacionales (Barras Apiladas)',
            '8': 'Grilla Cuantitativa (Matriz Numérica de Porcentajes)',
            '9': 'Gráficos de Dispersión en Facetas (Facet Grid Scatter)',
            '10': 'Compilar / Actualizar Base de Datos (CSV)'
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
            continue
        elif v_choice == '9':
            run_population_facet_scatter()
            continue
        elif v_choice == '10':
            from compilar_tabla_global import compilar_tabla_maestra
            print("\nCompilando...")
            compilar_tabla_maestra()
            input("\nPresione Enter para volver al menú principal...")
            continue
            
        if v_choice == '0':
            plot_options = {
                '1': 'Evolución M_embr vs t (por R_gap)',
                '2': 'Evolución M_embr vs t (por M_gap)',
                '3': 'Mosaico Evolutivo (por R_gap)',
                '4': 'Mosaico Evolutivo (por M_gap)',
                '5': 'Heatmaps (Masa y Agua)',
                '6': 'Mapa de Evolución Continente del Éxito',
                '7': 'TODAS LAS GRÁFICAS'
            }
            print_header(["PROCESAMIENTO POR LOTE (Run All)"])
            p_choice = get_choice("¿Qué gráficas desea generar para TODOS los casos?", plot_options)
            if p_choice == '0':
                continue
                
            run_all_cases_batch(v_frag_map, alphas_full, m0_map, m0_names, m0_options, p_choice)
            input("\nProceso por lote completado. Presione Enter para continuar...")
            continue
            
        v_val = v_frag_map.get(v_choice)
        if not v_val: continue
        breadcrumbs.append(f"v_frag: {v_val} m/s")
        
        # Filtrar opciones de alpha según el v_frag seleccionado
        if v_val in [1, 3]:
            alphas = {
                '1': '0.0001', '2': '0.0003', '3': '0.0005',
                '4': '0.0007', '5': '0.001'
            }
        else:
            alphas = alphas_full
        
        while True:
            print_header(breadcrumbs)
            a_choice = get_choice("Elija un alpha:", alphas)
            if a_choice == '0': break
            a_val = float(alphas[a_choice])
            
            # temporal breadcrumbs copy
            bc_a = breadcrumbs + [f"alpha: {a_val}"]
            
            while True:
                print_header(bc_a)
                
                # Check available scenarios for the chosen v_frag
                base_run_dir = f"data/runs/vf_{v_val}ms"
                if not os.path.exists(base_run_dir):
                    print(f"\nNo se encontró la carpeta {base_run_dir}. Asegúrese de que las simulaciones estén ahí.")
                    input("\nPresione Enter para volver...")
                    break
                    
                scenarios_available = [d for d in os.listdir(base_run_dir) if os.path.isdir(os.path.join(base_run_dir, d))]
                
                # Filtrar escenarios que realmente tienen datos para el alpha seleccionado
                valid_scenarios = []
                for s in scenarios_available:
                    s_path = os.path.join(base_run_dir, s)
                    subdirs = [d for d in os.listdir(s_path) if os.path.isdir(os.path.join(s_path, d))]
                    
                    if s == 'delayed':
                        # Los casos delayed mantienen hardcode de a=0.001 por ahora
                        if a_val == 0.001:
                            valid_scenarios.append(s)
                    else:
                        # general o rounded: buscar si existe el flag de alpha en los nombres
                        # Formato: a_val=0.001 -> _a0.001
                        # Evitar confusión entre 0.001 y 0.0015 etc, agregando "_" o chequeando exactitud
                        has_alpha = any(f"_a{a_val}" in d for d in subdirs)
                        if has_alpha:
                            valid_scenarios.append(s)
                            
                if not valid_scenarios:
                    print(f"\nNo se encontraron escenarios con simulaciones para alpha = {a_val}.")
                    input("\nPresione Enter para volver...")
                    break
                    
                s_options = {str(i+1): s for i, s in enumerate(valid_scenarios)}
                s_choice = get_choice("Elija el escenario:", s_options)
                if s_choice == '0': break
                s_name = s_options[s_choice]
                
                if s_name == 'smooth':
                    keys_to_run = list(m0_options.keys())
                    run_smooth_global(v_val, keys_to_run, m0_map, m0_names, alphas_full, batch_mode=False, plot_choice='7')
                    input("\nPresione Enter para volver al menú principal...")
                    break
                
                bc_s = bc_a + [f"Escenario: {s_name}"]
                
                while True:
                    print_header(bc_s)
                    mode_options = {
                        '1': 'Todos los casos (1e-5, 1e-4, 1e-3, 1e-2, 1e-1 M_earth)', 
                        '2': 'Personalizado (Elegir uno)'
                    }
                    mode_choice = get_choice("¿Qué masas iniciales de embrión (M0) desea analizar?", mode_options)
                    if mode_choice == '0': break
                    
                    if mode_choice == '1':
                        keys_to_run = list(m0_options.keys())
                        run_analysis_multiple(v_val, a_val, s_name, keys_to_run, m0_map, m0_names)
                        input("\nPresione Enter para volver al menú principal...")
                        return
                    else:
                        while True:
                            print_header(bc_s + ["M0: Personalizado"])
                            m_choice = get_choice("Seleccione la masa inicial del embrión (M0):", m0_options)
                            if m_choice == '0': break
                            
                            run_analysis_multiple(v_val, a_val, s_name, [m_choice], m0_map, m0_names)
                            input("\nPresione Enter para volver al menú principal...")
                            return

def run_smooth_global(v_val, m_keys, m0_map, m0_names, alphas_full, batch_mode=False, plot_choice='7'):
    if not batch_mode:
        clear_screen()
    runs_path = f"data/runs/vf_{v_val}ms/smooth"
    fig_dir = f"data/figures/vf_{v_val}ms/smooth"
    
    print("=====================================================")
    print("           RESUMEN DEL ANÁLISIS SMOOTH GLOBAL        ")
    print("=====================================================")
    print(f" Velocidad de frag.: {v_val} m/s")
    print(f" Escenario:          Smooth (Agrupado por alpha y m0)")
    print("-----------------------------------------------------")
    print(f" Directorio de runs: {runs_path}")
    print("=====================================================")
    
    if not batch_mode:
        confirm = input("\n¿Desea proceder con el análisis? [S/N]: ").strip().lower()
        if confirm != 's':
            return
            
    print("\nCargando librerías y preparando el pipeline...")
    from pipeline_analysis.analyzer_base import SmoothAnalyzer
    from pipeline_analysis.plotters import PlotterSmooth
    import os
    
    if v_val in [1, 3]:
        alphas_list = ['0.0001', '0.0003', '0.0005', '0.0007', '0.001']
    else:
        alphas_list = list(alphas_full.values())
        
    global_data = []
    
    for a_str in alphas_list:
        a_val = float(a_str)
        if not os.path.exists(runs_path):
            continue
            
        subdirs = [d for d in os.listdir(runs_path) if os.path.isdir(os.path.join(runs_path, d))]
        if not any(f"_a{a_val}" in d for d in subdirs):
            continue
            
        for m_choice in m_keys:
            m0_val = m0_map[m_choice]
            cache_path = os.path.join(runs_path, f"cache_smooth_a_{a_val}_M0_{m0_val}.pkl")
            analyzer = SmoothAnalyzer(runs_path, cache_path, m0_val, alpha_val=a_val)
            data_filtered = analyzer.load_or_extract()
            for d in data_filtered:
                d['m0_earth'] = m0_val
            global_data.extend(data_filtered)
            
    print(f"\n-> Se encontraron {len(global_data)} simulaciones agrupadas para Smooth.")
    if len(global_data) == 0:
        print("No se generarán gráficas ya que no hay datos para smooth.")
        return
        
    print(f"Generando gráficas en: {fig_dir}")
    if plot_choice in ['1', '2', '3', '4', '7']:
        PlotterSmooth.plot_evolution(global_data, fig_dir, alpha_val="Todos")
    if plot_choice in ['5', '6', '7']:
        PlotterSmooth.plot_heatmaps(global_data, fig_dir)
        PlotterSmooth.plot_categorical_map(global_data, fig_dir)
        
    print("\n¡Análisis Smooth Global completado exitosamente!")

def run_analysis_multiple(v_val, a_val, s_name, m_keys, m0_map, m0_names, batch_mode=False, plot_choice='7'):
    if not batch_mode:
        clear_screen()
    runs_path = f"data/runs/vf_{v_val}ms/{s_name}"
    
    print("=====================================================")
    print("                RESUMEN DEL ANÁLISIS                 ")
    print("=====================================================")
    print(f" Velocidad de frag.: {v_val} m/s")
    print(f" Alpha:              {a_val}")
    print(f" Escenario:          {s_name.capitalize()}")
    if len(m_keys) > 1:
        print(f" Masa inicial (M0):  TODOS LOS {len(m_keys)} CASOS")
    else:
        print(f" Masa inicial (M0):  {m0_map[m_keys[0]]} M_earth")
    print("-----------------------------------------------------")
    print(f" Directorio de runs: {runs_path}")
    print("=====================================================")
    
    if not batch_mode:
        confirm = input("\n¿Desea proceder con el análisis? [S/N]: ").strip().lower()
        if confirm != 's':
            return
        
    print("\nCargando librerías y preparando el pipeline...")
    from pipeline_analysis.analyzer_base import GeneralAnalyzer, SinusoidalAnalyzer, DelayedAnalyzer, SmoothAnalyzer
    from pipeline_analysis.plotters import PlotterGeneral, PlotterSinusoidal, PlotterSmooth
    
    for m_choice in m_keys:
        m0_val = m0_map[m_choice]
        m0_str = m0_names[m_choice]
        print(f"\n--- Procesando caso M0 = {m0_val} M_earth ---")
        
        cache_path = os.path.join(runs_path, f"cache_{s_name}_a_{a_val}_M0_{m0_val}.pkl")
        fig_dir = f"data/figures/vf_{v_val}ms/a_{a_val}/{s_name}/{m0_str}"
        
        # Select Analyzer based on scenario
        if s_name == 'sinusoidal':
            analyzer = SinusoidalAnalyzer(runs_path, cache_path, m0_val, alpha_val=a_val)
        elif s_name == 'delayed':
            analyzer = DelayedAnalyzer(runs_path, cache_path, m0_val, alpha_val=a_val)
        elif s_name == 'smooth':
            analyzer = SmoothAnalyzer(runs_path, cache_path, m0_val, alpha_val=a_val)
        else:
            # general or rounded
            analyzer = GeneralAnalyzer(runs_path, cache_path, m0_val, alpha_val=a_val)
            
        data_filtered = analyzer.load_or_extract()
            
        print(f"-> Se encontraron {len(data_filtered)} simulaciones para estos parámetros.")
        if len(data_filtered) == 0:
            print("No se generarán gráficas ya que no hay datos para este alpha.")
            continue
            
        print(f"Generando gráficas en: {fig_dir}")
        if s_name == 'sinusoidal':
            if plot_choice in ['1', '3', '7']:
                PlotterSinusoidal.plot_lines_grouped_by_gap_type(data_filtered, fig_dir)
            if plot_choice in ['5', '6', '7']:
                PlotterSinusoidal.plot_heatmaps(data_filtered, fig_dir)
        elif s_name == 'smooth':
            if plot_choice in ['1', '2', '3', '4', '7']:
                PlotterSmooth.plot_evolution(data_filtered, fig_dir, alpha_val=a_val)
        else:
            if plot_choice in ['1', '7']:
                PlotterGeneral.plot_lines_grouped_by_rgap(data_filtered, fig_dir, alpha_val=a_val)
            if plot_choice in ['2', '7']:
                PlotterGeneral.plot_lines_grouped_by_mgap(data_filtered, fig_dir, alpha_val=a_val)
            if plot_choice in ['3', '7']:
                PlotterGeneral.plot_lines_grouped_by_rgap_mosaic(data_filtered, fig_dir, alpha_val=a_val)
            if plot_choice in ['4', '7']:
                PlotterGeneral.plot_lines_grouped_by_mgap_mosaic(data_filtered, fig_dir, alpha_val=a_val)
            if s_name in ['general', 'rounded']:
                if plot_choice in ['5', '7']:
                    PlotterGeneral.plot_heatmaps(data_filtered, fig_dir, alpha_val=a_val)
                if plot_choice in ['6', '7']:
                    PlotterGeneral.plot_categorical_map(data_filtered, fig_dir, alpha_val=a_val)
                
    print("\n¡Análisis completado exitosamente!")

def run_all_cases_batch(v_frag_map, alphas_full, m0_map, m0_names, m0_options, plot_choice='7'):
    print("\nIniciando PROCESAMIENTO POR LOTE (Run All)...")
    for v_choice, v_val in v_frag_map.items():
        if v_val in [1, 3]:
            alphas = {
                '1': '0.0001', '2': '0.0003', '3': '0.0005',
                '4': '0.0007', '5': '0.001'
            }
        else:
            alphas = alphas_full
            
        base_run_dir = f"data/runs/vf_{v_val}ms"
        if not os.path.exists(base_run_dir):
            continue
            
        scenarios_available = [d for d in os.listdir(base_run_dir) if os.path.isdir(os.path.join(base_run_dir, d))]
        
        smooth_run = False
        for a_choice, a_str in alphas.items():
            a_val = float(a_str)
            valid_scenarios = []
            for s in scenarios_available:
                s_path = os.path.join(base_run_dir, s)
                subdirs = [d for d in os.listdir(s_path) if os.path.isdir(os.path.join(s_path, d))]
                
                if s == 'smooth':
                    if not smooth_run:
                        keys_to_run = list(m0_options.keys())
                        run_smooth_global(v_val, keys_to_run, m0_map, m0_names, alphas_full, batch_mode=True, plot_choice=plot_choice)
                        smooth_run = True
                    continue
                elif s == 'delayed':
                    if a_val == 0.001:
                        valid_scenarios.append(s)
                else:
                    has_alpha = any(f"_a{a_val}" in d for d in subdirs)
                    if has_alpha:
                        valid_scenarios.append(s)
            
            for s_name in valid_scenarios:
                keys_to_run = list(m0_options.keys())
                run_analysis_multiple(v_val, a_val, s_name, keys_to_run, m0_map, m0_names, batch_mode=True, plot_choice=plot_choice)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nSaliendo...")
        sys.exit(0)
