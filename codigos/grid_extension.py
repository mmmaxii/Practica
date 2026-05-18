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
