import os
import numpy as np
import dustpy.constants as c
import matplotlib.pyplot as plt
from pipeline_snowlines import WaterworldPipeline

def main():
    print("Test pipeline...")
    pipeline = WaterworldPipeline(
        datadir="data/test_scratch",
        active_species=[],
        grid_rmin=0.5 * c.au,
        grid_rmax=100 * c.au,
        Nr=300,
        M_star_Msun=1.0,
        R_star_Rsun=2.1,
        T_star_K=4000.0,
        gap_positions_au=[5.0, 10.0, 20.0, 40.0, 80.0],
        snowline_au=2.0
    )

    planets_to_add = [
        {"M_planet": 317.8 * c.M_earth, "a_planet_au": 5.0},
        {"M_planet": 317.8 * c.M_earth, "a_planet_au": 10.0},
        {"M_planet": 317.8 * c.M_earth, "a_planet_au": 20.0},
        {"M_planet": 317.8 * c.M_earth, "a_planet_au": 40.0},
        {"M_planet": 317.8 * c.M_earth, "a_planet_au": 80.0}
    ]

    print("Testing Duffell Multi para todos los gaps...")
    pipeline.setup_gap_duffell_multi(planets=planets_to_add, imprint=True)
    
    # ---------------- PLOT ----------------
    print("\nGenerando plot de Sigma...")
    r_au = pipeline.sim.grid.r / c.au
    sigma_gas = pipeline.sim.gas.Sigma
    sigma_dust = pipeline.sim.dust.Sigma.sum(axis=-1)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot de Sigma
    ax.plot(r_au, sigma_gas, lw=2, color='blue', label=r'$\Sigma_{gas}$')
    ax.plot(r_au, sigma_dust, lw=2, color='brown', label=r'$\Sigma_{dust}$')
    
    for p in planets_to_add:
        ax.axvline(p['a_planet_au'], color='gray', linestyle='--', alpha=0.6)
        
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-4)
    ax.set_xlabel('Distancia (AU)', fontsize=12)
    ax.set_ylabel(r'Densidad superficial $\Sigma$ [g/cm²]', fontsize=12)
    ax.set_title('Perfil de Densidad (Gaps Múltiples de Duffell)', fontsize=14)
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    
    # Guardar la imagen
    plt.tight_layout()
    plt.savefig('sigma_duffell.png', dpi=150)
    print("Plot guardado como 'sigma_duffell.png'")
    # --------------------------------------

if __name__ == "__main__":
    main()
