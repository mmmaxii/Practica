import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# path PA3Py
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "PA3Py"))
from PebbleAccretion3 import PebbleAccretionModule3
import dustpy.constants as c

DATAROOT = "data/1e5yr"
SAVEDIR = "figures/1e5yr/gradient_gaps"
os.makedirs(SAVEDIR, exist_ok=True)

# Configuraciones
EMBRYOS_AU = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
CATEGORIES = ["sup_jup", "jup", "sat", "nep", "sup_earth"]
GAPS_AU = [1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

# Colormap donde Rojo = 1 (bajo), Azul = 20 (alto)
# "jet_r", "coolwarm_r", "Spectral"
# Usaremos un cmap interpolado personalizado o 'rainbow_r' (Rojo -> Azul)
cmap = cm.get_cmap("rainbow_r")
norm = mcolors.Normalize(vmin=1.0, vmax=20.0)

def final_metrics(hist):
    if hist is None or len(hist) == 0:
        return np.nan, np.nan
    row = hist[-1]
    M = row[1]
    M_h2o = row[2]
    f_h2o = 100.0 * M_h2o / (M + 1e-30)
    return M / c.M_earth, f_h2o

def plot_for_category(cat):
    print(f"\n=======================================================")
    print(f" Analizando categoría: single_{cat}")
    print(f"=======================================================")
    
    fig, (ax_mass, ax_water) = plt.subplots(1, 2, figsize=(14, 5.5))
    
    runs_found = 0
    for gap in GAPS_AU:
        run_name = f"single_{cat}_{gap}au"
        run_path = os.path.join(DATAROOT, run_name)
        
        if not os.path.isdir(run_path) or not glob.glob(os.path.join(run_path, "data*.hdf5")):
            print(f"  [MISSING] {run_name}")
            continue
        
        runs_found += 1
        print(f"  → Procesando: {run_name}")
        
        try:
            # Inicializamos PA3 desde el t_min apropiado
            pa3 = PebbleAccretionModule3.from_datadir(run_path, t_min_yr=100.0)
        except Exception as e:
            print(f"    [!] Error al cargar datos: {e}")
            continue
            
        # Corremos todos los embriones (sin depleción cruzada)
        results = pa3.run_growth(EMBRYOS_AU)
        
        masses = []
        waters = []
        
        for r_emb in sorted(EMBRYOS_AU):
            hist = results.get(r_emb)
            m_core, f_h2o = final_metrics(hist)
            masses.append(m_core)
            waters.append(f_h2o)
            
        color = cmap(norm(gap))
        
        # Plot Mass vs Embryo pos
        ax_mass.plot(EMBRYOS_AU, masses, marker='o', markersize=4, 
                     color=color, linewidth=2, alpha=0.8,
                     label=f"Gap {gap} AU")
        
        # Plot Water Fraction vs Embryo pos
        ax_water.plot(EMBRYOS_AU, waters, marker='o', markersize=4, 
                      color=color, linewidth=2, alpha=0.8)

    if runs_found == 0:
        print("  [!] No se encontró ninguna simulación para esta categoría.")
        plt.close(fig)
        return

    # ----- ESTILIZADO DE LOS PLOTS -----
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=[ax_mass, ax_water], 
                        fraction=0.03, aspect=30, pad=0.02)
    cbar.set_label("Posición del Gap Planetario [AU]", fontsize=11, fontweight="bold")
    
    # Ax Mass
    ax_mass.set_title(f"Masa Final del Embrión ({cat})", fontsize=12, fontweight="bold")
    ax_mass.set_xlabel("Posición Inicial del Embrión [AU]", fontsize=10)
    ax_mass.set_ylabel(r"Masa Final [$M_\oplus$]", fontsize=10)
    ax_mass.set_yscale("log")
    ax_mass.grid(True, linestyle="--", alpha=0.5)
    
    # Ax Water
    ax_water.set_title(f"Fracción de Agua Final ({cat})", fontsize=12, fontweight="bold")
    ax_water.set_xlabel("Posición Inicial del Embrión [AU]", fontsize=10)
    ax_water.set_ylabel(r"Fracción de H$_2$O Final [%]", fontsize=10)
    ax_water.grid(True, linestyle="--", alpha=0.5)
    
    # Dibuja la línea de Waterworld
    ax_water.axhline(10.0, color="gray", linestyle="--", alpha=0.8, label="Límite Waterworld (10%)")
    ax_water.legend(loc="upper left")

    plt.suptitle(f"Perfil de Acreción de Pebbles bajo el Efecto del Gap: {cat}", 
                 fontsize=14, fontweight="bold", y=0.98)
    
    fig_name = f"gradient_{cat}.pdf"
    fig_path = os.path.join(SAVEDIR, fig_name)
    plt.savefig(fig_path, bbox_inches="tight", dpi=150)
    print(f"  ✓ Figura guardada en: {fig_path}")
    plt.close(fig)

if __name__ == "__main__":
    for cat in CATEGORIES:
        plot_for_category(cat)
    print(f"\n¡Análisis completo! Revisa la carpeta {SAVEDIR}")
