import numpy as np

def refinegrid_tripodpy(ri, r0, num=3):
    """Function to refine the radial grid

    Parameters
    ----------
    ri : array
        Radial interfaces
    r0 : float
        Radial location around which grid should be refined
    num : int, option, default : 3
        Number of refinement iterations

    Returns
    -------
    ri : array
        New refined radial grid"""
    if num == 0:
        return ri
    for _ in range(num):
        ind = np.argmin(r0 > ri) - 1
        ind_left  = ind-num
        ind_right = ind+num+1
        ri_left  = ri[:ind_left]
        ri_right = ri[ind_right:]
        number_refined_cells = (2*num+1)*2
        ri_mid = np.empty(number_refined_cells)
        for i in range(0, number_refined_cells, 2):
            j = ind-num+int(i/2)
            ri_mid[i] = ri[j]
            ri_mid[i+1] = 0.5*(ri[j]+ri[j+1])
        ri = np.concatenate((ri_left, ri_mid, ri_right))
    return ri



# Esta funcion fue extraido de dustpylib, funcion auxiliar para modificar la grilla de una 
# simulacion
def refinegrid_dustpy(ri, r0, num=3):
    """
    Refina la grilla radial logarítmicamente alrededor de r0.
    En cada nivel recursivo inserta el punto medio geométrico
    sqrt(r_i * r_{i+1}) dentro de una ventana centrada en r0.
    La ventana se reduce en cada llamada (num-1), concentrando
    el refinamiento cerca de r0.
    Parameters
    ----------
    ri : np.ndarray, shape (Nr+1,)
        Interfaces de la grilla radial [cm]
    r0 : float
        Radio alrededor del cual refinar [cm]
    num : int, optional, default=3
        Número de niveles de refinamiento
    Returns
    -------
    ri_fine : np.ndarray
        Grilla refinada [cm]
    """
    if num == 0:
        return ri
    # Índice de la interfaz inmediatamente a la izquierda de r0
    i0 = np.argmax(ri > r0) - 1
    # Límites de la ventana de refinamiento
    il = max(0, i0 - num + 1)
    ir = min(i0 + num, ri.shape[0] - 1)
    # Regiones sin modificar
    ril = ri[:il]
    rir = ri[ir:]
    # Región refinada: media geométrica (= punto medio en log-space)
    N = ir - il
    rim = np.empty(2 * N)
    for i in range(N):
        j = il + i
        rim[2 * i]     = ri[j]
        rim[2 * i + 1] = np.sqrt(ri[j] * ri[j + 1])
    ri_fine = np.concatenate((ril, rim, rir))
    # Siguiente nivel de refinamiento (recursivo) con ventana más angosta
    return refinegrid_dustpy(ri_fine, r0, num=num - 1)



def refine_time_local(ti, t0, num=3):
    if num == 0:
        return ti

    # índice a la izquierda de t0
    i0 = np.argmax(ti > t0) - 1

    il = max(0, i0 - num + 1)
    ir = min(i0 + num, ti.shape[0] - 1)

    til = ti[:il]
    tir = ti[ir:]

    N = ir - il
    tim = np.empty(2 * N)

    for i in range(N):
        j = il + i
        tim[2 * i]     = ti[j]
        tim[2 * i + 1] = np.sqrt(ti[j] * ti[j + 1])  # log midpoint

    ti_fine = np.concatenate((til, tim, tir))

    return refine_time_local(ti_fine, t0, num=num - 1)



import matplotlib.pyplot as plt
import numpy as np
import dustpy.constants as c

def plot_disk_diagnostics_multicomp(sim, snap=-1, trim_boundary=True, trim_cells=2):
    """
    Genera un dashboard para visualizar un disco en TriPoDPy con múltiples 
    componentes químicos (snowlines, etc).
    """
    r_au = sim.grid.r / c.au
    
    # --- 1. Cálculo Manual de Densidades Totales sumando componentes ---
    # Inicializamos arrays de ceros con la forma de la snapshot elegida
    shape_gas = sim.components.Default.gas.Sigma[snap].shape if sim.components.Default.gas.Sigma.ndim == 2 else sim.components.Default.gas.Sigma.shape
    shape_dust = sim.components.Default.dust.Sigma[snap].shape if sim.components.Default.dust.Sigma.ndim == 3 else sim.components.Default.dust.Sigma.shape
    
    sigma_gas_total = np.zeros(shape_gas[0]) if len(shape_gas) == 1 else np.zeros(shape_gas)
    sigma_dust_total = np.zeros(shape_gas[0]) if len(shape_gas) == 1 else np.zeros(shape_gas) # Queremos el polvo total sumado en masas

    # Sumar todas las componentes de gas activas
    for name, comp in sim.components.__dict__.items():
        if name.startswith("_"): continue
        if comp.gas._active:
            sigma_comp = comp.gas.Sigma[snap] if comp.gas.Sigma.ndim == 2 else comp.gas.Sigma
            sigma_gas_total += sigma_comp
            
    # Sumar todas las componentes de polvo activas
    for name, comp in sim.components.__dict__.items():
        if name.startswith("_"): continue
        if comp.dust._active:
            sigma_comp = comp.dust.Sigma[snap] if comp.dust.Sigma.ndim == 3 else comp.dust.Sigma
            sigma_dust_total += sigma_comp.sum(axis=-1) # Sumamos sobre las masas
            
    # Recalculamos epsilon manualmente
    eps_recalc = sigma_dust_total / sigma_gas_total

    # --- 2. Extracción de las demás variables ---
    eta = sim.gas.eta[snap] if sim.gas.eta.ndim == 2 else sim.gas.eta
    a_max = sim.dust.s.max[snap] if sim.dust.s.max.ndim == 2 else sim.dust.s.max
    
    if sim.dust.St.ndim == 3:
        st_max = sim.dust.St[snap, :, -1]
    else:
        st_max = sim.dust.St[:, -1]

    # --- 3. Recorte de bordes ---
    if trim_boundary and trim_cells > 0:
        r_au = r_au[:-trim_cells]
        sigma_gas_total = sigma_gas_total[:-trim_cells]
        sigma_dust_total = sigma_dust_total[:-trim_cells]
        eps_recalc = eps_recalc[:-trim_cells]
        eta = eta[:-trim_cells]
        a_max = a_max[:-trim_cells]
        st_max = st_max[:-trim_cells]

    # --- 4. Creación de la figura ---
    plt.style.use('dark_background')
    fig, axs = plt.subplots(3, 2, figsize=(14, 15), sharex=True)
    fig.suptitle(f'Diagnóstico Multicomponente (Snapshot: {snap})', fontsize=16, fontweight='bold', color='white')

    # Plot 1: Densidad Superficial Total
    axs[0, 0].loglog(r_au, sigma_gas_total, label='Gas (Suma total)', color='#1f77b4', linewidth=2)
    axs[0, 0].loglog(r_au, sigma_dust_total, label='Polvo (Suma total)', color='#ff7f0e', linewidth=2)
    axs[0, 0].set_ylabel('Densidad Superficial [g/cm²]')
    axs[0, 0].set_title('Densidad Total Recalculada')
    axs[0, 0].grid(True, which="both", ls="--", alpha=0.3)
    axs[0, 0].legend()

    # Plot 2: Parámetro Eta
    axs[0, 1].semilogx(r_au, eta, label='Gradiente de presión ($\eta$)', color='#2ca02c', linewidth=2)
    axs[0, 1].axhline(0, color='white', linestyle='--', linewidth=1, alpha=0.5)
    axs[0, 1].set_ylabel('Parámetro $\eta$')
    axs[0, 1].set_title('Estructura de Presión del Gas')
    axs[0, 1].grid(True, which="both", ls="--", alpha=0.3)
    axs[0, 1].legend()

    # Plot 3: Tamaño Máximo
    axs[1, 0].loglog(r_au, a_max, label='Tamaño máximo ($a_{max}$)', color='#d62728', linewidth=2)
    axs[1, 0].set_ylabel('Tamaño de Partícula [cm]')
    axs[1, 0].set_title('Dinámica del Polvo')
    axs[1, 0].grid(True, which="both", ls="--", alpha=0.3)
    axs[1, 0].legend()

    # Plot 4: Acumulación Relativa (EPSILON RECALCULADO)
    axs[1, 1].loglog(r_au, eps_recalc, label='Razón Polvo/Gas Recalculada ($\epsilon$)', color='#9467bd', linewidth=2)
    # Aquí puedes ver los "saltos" donde las especies de hielo se suman al polvo
    axs[1, 1].axhline(1.0, color='white', linestyle=':', label='Límite $\epsilon = 1$', alpha=0.5)
    axs[1, 1].set_ylabel('$\epsilon$')
    axs[1, 1].set_title('Acumulación Relativa ($\Sigma_{dust} / \Sigma_{gas}$)')
    axs[1, 1].grid(True, which="both", ls="--", alpha=0.3)
    axs[1, 1].legend()

    # Plot 5: Stokes
    axs[2, 0].loglog(r_au, st_max, label='Stokes de $a_{max}$', color='#e377c2', linewidth=2)
    axs[2, 0].axhline(1.0, color='white', linestyle=':', alpha=0.5)
    axs[2, 0].set_ylabel('Número de Stokes ($St$)')
    axs[2, 0].set_xlabel('Radio [AU]')
    axs[2, 0].set_title('Acoplamiento Aerodinámico')
    axs[2, 0].grid(True, which="both", ls="--", alpha=0.3)
    axs[2, 0].legend()

    axs[2, 1].axis('off')

    plt.tight_layout()
    plt.subplots_adjust(top=0.94) 
    plt.show()
    plt.style.use('default')