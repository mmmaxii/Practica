"""
plot_diagnostics.py
===================
Herramienta de diagnóstico visual para simulaciones de tripodpy con snowlines.
Lee los archivos HDF5 generados por el pipeline y produce los siguientes gráficos:

  1. Hovmöller (eps / Sigma_dust / Sigma_gas)  →  diagrama espacio-tiempo
  2. Distribución de tamaño de grano           →  Sigma(a, r) en una snapshot
  3. Pebble flux Ṁ_pebble (r, t)              →  flujo másico espacio-temporal
  4. Perfiles de disco (η, St, a_max, Σ)       →  4-panel en un snapshot
  5. Densidades superficiales de gas por componente   → plot_gas_components
  6. Densidades superficiales de polvo por componente → plot_dust_components

Uso desde línea de comandos:
    python plot_diagnostics.py <datadir> [savedir]

Uso como módulo:
    from plot_diagnostics import SnowlineDiagnostics
    d = SnowlineDiagnostics("output_test_pipeline", savedir="figs", r_trim=0.93)
    d.plot_hovmoller()
    d.plot_size_distribution(it=-1)
    d.plot_pebble_flux()
    d.plot_profiles(it=-1)
    d.plot_gas_components(it=-1)
    d.plot_dust_components(it=-1)
    plt.show()

Notas de diseño:
    - Todas las unidades internas están en CGS (tal como las almacena tripodpy).
    - Los ejes de los gráficos usan AU (radio) y años (tiempo) para legibilidad.
    - r_trim  recorta el último (1 - r_trim) del grid para evitar artefactos
      numéricos en el borde exterior (SigmaFloor explosion).
    - t_min_yr recorta la evolución temporal, omitiendo los primeros años
      donde la simulación aún no tiene polvo físicamente relevante.
    - Los archivos se guardan SOLO en formato PDF.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from dustpy import hdf5writer
import dustpy.constants as c

# ──────────────────────────────────────────────────────────────────────────────
# Estética global
# ──────────────────────────────────────────────────────────────────────────────
SNOWLINE_STYLES = {
    "H2O": {"color": "#4FC3F7", "label": r"$\mathrm{H_2O}$ snowline", "T": 150.0},
    "CO2": {"color": "#80CBC4", "label": r"$\mathrm{CO_2}$ snowline", "T": 70.0},
    "CO":  {"color": "#FFCC80", "label": r"$\mathrm{CO}$ snowline",   "T": 25.0},
}

# Colores fijos para componentes conocidos
COMP_COLORS = {
    "H2O":      "#4FC3F7",
    "CO2":      "#80CBC4",
    "CO":       "#FFCC80",
    "silicates":"#EF9A9A",
    "Default":  "#CE93D8",
    "CH4":      "#A5D6A7",
    "C2H6":     "#FFCC80",
    "N2":       "#B0BEC5",
    "NH3":      "#F48FB1",
}

COMP_LABELS = {
    "H2O":      r"$\mathrm{H_2O}$",
    "CO2":      r"$\mathrm{CO_2}$",
    "CO":       r"$\mathrm{CO}$",
    "silicates":r"Silicatos",
    "Default":  r"H$_2$ (Default)",
}

plt.rcParams.update({
    "font.family":    "DejaVu Sans",
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize":9,
    "ytick.labelsize":9,
    "axes.grid":      True,
    "grid.alpha":     0.2,
    "figure.dpi":     130,
})


# ──────────────────────────────────────────────────────────────────────────────
# Clase principal
# ──────────────────────────────────────────────────────────────────────────────
class SnowlineDiagnostics:
    """
    Carga los datos HDF5 de un directorio de salida de tripodpy y expone
    6 métodos de diagnóstico científico.

    Parameters
    ----------
    datadir : str
        Directorio con los archivos data0000.hdf5, data0001.hdf5, ...
    savedir : str, optional
        Directorio donde guardar las figuras en PDF (None = no guarda).
    r_trim : float, optional
        Fracción del grid radial a conservar [0, 1]. Default 0.93 recorta el
        7 % exterior para evitar artefactos numéricos del SigmaFloor.
    t_min_yr : float, optional
        Tiempo mínimo (en años) a mostrar en gráficos temporales.
        Omite los snapshots antes de este valor. Default 100.
    """

    def __init__(self, datadir: str, savedir: str = None,
                 r_trim: float = 0.93, t_min_yr: float = 100.0):
        self.datadir   = datadir
        self.savedir   = savedir
        self.r_trim    = r_trim
        self.t_min_yr  = t_min_yr

        print(f"[SnowlineDiagnostics] Cargando datos desde: {datadir}")
        wrtr = hdf5writer()
        wrtr.datadir = datadir
        self.data = wrtr.read.all()

        # Grid completo (CGS)
        self.r_full   = self.data.grid.r[0, :]   # (Nr,)  [cm]
        self.t_full   = self.data.t               # (Nt,)  [s]
        self.Nr_full  = self.r_full.size
        self.Nt_full  = self.t_full.size

        # Máscara radial (recorte del borde exterior)
        Nr_keep   = max(2, int(self.Nr_full * r_trim))
        self.r_mask = slice(None, Nr_keep)
        self.r    = self.r_full[self.r_mask]
        self.r_au = self.r / c.au
        self.Nr   = self.r.size

        # Máscara temporal (omitir snapshots tempranos sin polvo relevante)
        t_full_yr = self.t_full / c.year
        self.t_mask = t_full_yr >= t_min_yr
        if not self.t_mask.any():
            self.t_mask = np.ones(self.Nt_full, dtype=bool)
        self.t    = self.t_full[self.t_mask]
        self.t_yr = self.t / c.year
        self.Nt   = self.t.size

        print(f"  → {self.Nt_full} snapshots totales | Nr_full = {self.Nr_full}")
        print(f"  → Usando {self.Nt} snapshots  (t ≥ {t_min_yr:.0f} yr)")
        print(f"  → Usando {self.Nr} celdas radiales de {self.Nr_full}  (r_trim={r_trim})")
        print(f"  → r = {self.r_au[0]:.2f} – {self.r_au[-1]:.1f} AU")
        print(f"  → t = {self.t_yr[0]:.1e} – {self.t_yr[-1]:.1e} yr")

        self._detect_components()

    # ── Detección de componentes ─────────────────────────────────────────────
    def _detect_components(self):
        self.components = []
        if hasattr(self.data, "components"):
            for name in vars(self.data.components):
                if not name.startswith("_"):
                    self.components.append(name)
        print(f"  → Componentes detectados: {self.components or 'ninguno'}")

    def debug_components(self):
        """Imprime la estructura detallada de los componentes en el HDF5 para diagnóstico."""
        comp_ns = getattr(self.data, "components", None)
        if comp_ns is None:
            print("[debug] No hay 'components' en los datos.")
            return
        print("[debug] Estructura de componentes leidos del HDF5:")
        for name, comp in vars(comp_ns).items():
            if name.startswith("_"):
                continue
            print(f"  {name}: {type(comp).__name__}")
            if comp is None:
                continue
            try:
                attrs = [k for k in vars(comp) if not k.startswith("_")]
                print(f"    atributos: {attrs}")
                for subname in attrs:
                    sub = getattr(comp, subname)
                    if hasattr(sub, '__dict__'):
                        subattrs = [k for k in vars(sub) if not k.startswith("_")]
                        print(f"    .{subname}: {type(sub).__name__}  attrs={subattrs}")
                        for sa in subattrs:
                            v = getattr(sub, sa, None)
                            shape = getattr(v, 'shape', type(v).__name__)
                            print(f"      .{subname}.{sa}: shape={shape}")
                    else:
                        shape = getattr(sub, 'shape', type(sub).__name__)
                        print(f"    .{subname}: shape={shape}")
            except Exception as e:
                print(f"    [error] {e}")

    # ── Guardar figura (solo PDF) ────────────────────────────────────────────
    def _save(self, fig, name):
        if self.savedir:
            os.makedirs(self.savedir, exist_ok=True)
            path = os.path.join(self.savedir, name + ".pdf")
            fig.savefig(path, bbox_inches="tight")
            print(f"  → Guardado: {path}")

    # ── Snowlines: lectura, series temporales, anotaciones ──────────────────
    def _get_rsnow_at(self, species_name, it_abs=None):
        """
        Devuelve la posición del snowline [AU] para una especie en un snapshot.

        Lee primero desde data.grid.rsnow_<species> (Field guardado por el
        pipeline con save=True). Si no existe, computa desde gas.T como fallback.

        Parameters
        ----------
        species_name : str
            Nombre de la especie: 'H2O', 'CO2', 'CO'.
        it_abs : int or None
            Índice en el subconjunto de snapshots válidos. None = último.
        """
        it_use = it_abs if it_abs is not None else (self.Nt - 1)
        full_indices = np.where(self.t_mask)[0]
        it_full = full_indices[it_use % len(full_indices)]

        # Intentar leer desde el Field guardado en HDF5
        field_name = f"rsnow_{species_name}"
        rsnow_data = getattr(self.data.grid, field_name, None)
        if rsnow_data is not None:
            arr = np.asarray(rsnow_data)
            if arr.ndim >= 1 and arr.shape[0] > it_full:
                r_cm = float(arr[it_full])
                if r_cm > 0:
                    return r_cm / c.au

        # Fallback: calcular desde el perfil de T
        T_sub = SNOWLINE_STYLES.get(species_name, {}).get("T", None)
        if T_sub is None or not hasattr(self.data.gas, "T"):
            return None
        T_radial = self.data.gas.T[it_full, self.r_mask]
        idx = np.where(T_radial < T_sub)[0]
        return float(self.r_au[idx[0]]) if idx.size > 0 else None

    def _get_rsnow_series(self, species_name):
        """
        Devuelve la evolución temporal del snowline [AU] para todos los
        snapshots válidos (shape: (Nt,)).

        Lee desde data.grid.rsnow_<species> si está disponible (HDF5 guardado
        con el campo del pipeline). Si no, computa snapshot a snapshot desde T.
        """
        field_name = f"rsnow_{species_name}"
        rsnow_data = getattr(self.data.grid, field_name, None)
        if rsnow_data is not None:
            arr = np.asarray(rsnow_data)
            if arr.ndim >= 1 and arr.shape[0] >= self.Nt_full:
                r_au = arr[self.t_mask] / c.au

                # Sanity check: si >=80% de los valores estan pegados al borde
                # interior (= ri[0] por bug de T_bind en vez de T_sub),
                # descartar y usar fallback desde gas.T.
                r_inner_cm = float(np.asarray(self.data.grid.r)[0, 0])
                n_stuck = np.sum(np.abs(arr[self.t_mask] - r_inner_cm) < 0.05 * r_inner_cm)
                if n_stuck / len(arr[self.t_mask]) > 0.8:
                    pass  # cae al fallback abajo
                else:
                    r_au = np.where(r_au > self.r_au[0] * 0.5, r_au, np.nan)
                    return r_au

        # Fallback: calcular desde T para cada snapshot
        T_sub = SNOWLINE_STYLES.get(species_name, {}).get("T", None)
        if T_sub is None or not hasattr(self.data.gas, "T"):
            return None
        full_indices = np.where(self.t_mask)[0]
        result = []
        for it_full in full_indices:
            T_radial = self.data.gas.T[it_full, self.r_mask]
            idx = np.where(T_radial < T_sub)[0]
            result.append(float(self.r_au[idx[0]]) if idx.size > 0 else np.nan)
        return np.array(result)

    def _add_snowline_vlines(self, ax, alpha=0.6, it_abs=None):
        """
        Dibuja líneas verticales en la posición del snowline para cada especie.
        Usa el snapshot correcto (it_abs), NO siempre el primero.
        Para plots de snapshot (perfiles, distribución de tamaño, componentes).
        """
        for name, style in SNOWLINE_STYLES.items():
            r_snow = self._get_rsnow_at(name, it_abs)
            # Permitir snowlines muy cerca del borde interior (H2O puede estar a ~1 AU)
            if r_snow is None or r_snow < self.r_au[0] * 0.5 or r_snow >= self.r_au[-1]:
                continue
            ax.axvline(r_snow, color=style["color"], ls="--",
                       lw=1.2, alpha=alpha,
                       label=f'{style["label"]} ({r_snow:.1f} AU)')

    def _add_snowline_hlines(self, ax, alpha=0.6, it_abs=None):
        """
        Fallback estático: dibuja axhline en r_snow para plots con r en Y.
        Solo se usa cuando _add_snowline_curves no tiene datos disponibles.
        """
        for name, style in SNOWLINE_STYLES.items():
            r_snow = self._get_rsnow_at(name, it_abs)
            if r_snow is None or r_snow < self.r_au[0] * 0.5 or r_snow >= self.r_au[-1]:
                continue
            ax.axhline(r_snow, color=style["color"], ls="--",
                       lw=1.2, alpha=alpha,
                       label=f'{style["label"]} ({r_snow:.1f} AU)')

    def _add_snowline_curves(self, ax, alpha=0.85):
        """
        Dibuja curvas r_snow(t) en un diagrama espacio-tiempo (Hovmöller).
        Muestra cómo migran las snowlines con la evolución del disco.
        Para plots 2D: ejes radio (X) vs tiempo (Y).
        """
        any_plotted = False
        for name, style in SNOWLINE_STYLES.items():
            r_series = self._get_rsnow_series(name)
            if r_series is None:
                continue
            # Filtrar puntos fuera del grid graficado
            r_plot = np.where(
                (r_series > self.r_au[0]) & (r_series < self.r_au[-1]),
                r_series, np.nan
            )
            if np.all(np.isnan(r_plot)):
                continue
            # r en Y, t en X: ax.plot(t, r)
            ax.plot(self.t_yr,
                    r_plot,
                    color=style["color"], ls="--", lw=2.0,
                    alpha=alpha, label=style["label"])
            any_plotted = True

        if not any_plotted:
            # Fallback: líneas horizontales estáticas desde T[0]
            self._add_snowline_hlines(ax, alpha=alpha * 0.7)

    # =========================================================================
    # 1. HOVMÖLLER: Sigma_dust (o ε) vs (radio, tiempo)
    # =========================================================================
    def plot_hovmoller(self, quantity="eps", cmap="magma", vmin=None, vmax=None):
        """
        Diagrama espacio-tiempo (Hovmöller).

        Parameters
        ----------
        quantity : "eps" | "Sigma_dust" | "Sigma_gas"
        """
        rm, tm = self.r_mask, self.t_mask

        Sigma_dust = self.data.dust.Sigma[np.ix_(tm, np.arange(self.Nr_full))][
            :, self.r_mask].sum(-1)                              # (Nt, Nr)
        Sigma_gas  = self.data.gas.Sigma[np.ix_(tm, np.arange(self.Nr_full))][
            :, self.r_mask]                                      # (Nt, Nr)

        if quantity == "eps":
            Z       = Sigma_dust / (Sigma_gas + 1e-300)
            label_cb = r"log$_{10}$ Relación polvo/gas $\epsilon$"
            title    = r"Diagrama Espacio-Tiempo: $\epsilon = \Sigma_{dust}/\Sigma_{gas}$"
        elif quantity == "Sigma_dust":
            Z       = Sigma_dust
            label_cb = r"log$_{10}$ $\Sigma_{dust}$ [g cm$^{-2}$]"
            title    = r"Diagrama Espacio-Tiempo: $\Sigma_{dust}$"
        else:
            Z       = Sigma_gas
            label_cb = r"log$_{10}$ $\Sigma_{gas}$ [g cm$^{-2}$]"
            title    = r"Diagrama Espacio-Tiempo: $\Sigma_{gas}$"

        Z = np.clip(Z, 1e-20, None)

        fig, ax = plt.subplots(figsize=(10, 5))
        r_edges = _log_edges(self.r_au)                         # Y edges
        t_edges = _log_edges(np.maximum(self.t_yr, 1.0))       # X edges

        logZ    = np.log10(Z)                                   # (Nt, Nr)
        vlim_lo = vmin if vmin is not None else float(np.percentile(logZ[Z > 1e-20], 2))
        vlim_hi = vmax if vmax is not None else float(np.percentile(logZ[Z > 1e-20], 98))

        # Transponer: (Nt, Nr) -> (Nr, Nt) para r en Y, t en X
        pcm = ax.pcolormesh(t_edges, r_edges, logZ.T,
                            cmap=cmap, vmin=vlim_lo, vmax=vlim_hi, shading="flat")
        fig.colorbar(pcm, ax=ax, pad=0.02, label=label_cb)

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Tiempo [años]")
        ax.set_ylabel("Radio [AU]")
        ax.set_title(title, pad=8)
        self._add_snowline_curves(ax, alpha=0.85)

        fig.tight_layout()
        self._save(fig, "hovmoller_" + quantity)
        return fig, ax

    # =========================================================================
    # 1b. HOVMÖLLER R(t): radio en Y (log), tiempo en X — estilo paper
    # =========================================================================
    def plot_hovmoller_rt(self, quantities=None, cmaps=None,
                          t_unit="kyr", t_log=False,
                          figsize=None, percentile_range=(2, 98)):
        """
        Diagrama espacio-tiempo orientado R(t): radio en Y (log), tiempo en X.
        Produce N paneles apilados verticalmente con eje X compartido.

        Es el formato estándar de artículos científicos (Drążkowska et al.,
        Booth & Ilee, etc.) que permite ver la evolución temporal de estructuras
        radiales como snowlines, bumps y traffic jams de pebbles.

        Parameters
        ----------
        quantities : list of str, optional
            Cantidades a graficar. Opciones:
            "Sigma_gas", "Sigma_dust", "T", "a_max", "eps", "St"
            Default: ["Sigma_gas", "T", "a_max"]
        cmaps : dict, optional
            Colormaps por cantidad. Ej: {"T": "inferno"}.
        t_unit : "yr" | "kyr" | "Myr"
            Unidades del eje temporal.
        t_log : bool
            Si True, eje X en escala logarítmica.
        figsize : tuple, optional
            Default: (10, 3.5 * n_paneles).
        percentile_range : tuple
            (vmin_pct, vmax_pct) para el colorbar. Default: (2, 98).
        """
        if quantities is None:
            quantities = ["Sigma_gas", "T", "a_max"]
        n_panels = len(quantities)

        _default_cmaps = {
            "Sigma_gas":  "viridis",
            "Sigma_dust": "magma",
            "T":          "inferno",
            "a_max":      "plasma",
            "eps":        "RdYlBu_r",
            "St":         "cividis",
        }
        _labels = {
            "Sigma_gas":  r"log $\Sigma_g$ [g cm$^{-2}$]",
            "Sigma_dust": r"log $\Sigma_d$ [g cm$^{-2}$]",
            "T":          r"log $T$ [K]",
            "a_max":      r"log $a_\mathrm{max}$ [cm]",
            "eps":        r"log $\epsilon = \Sigma_d/\Sigma_g$",
            "St":         r"log $St$",
        }
        if cmaps is None:
            cmaps = {}

        full_t_idx = np.where(self.t_mask)[0]
        Nr_idx     = np.arange(self.Nr_full)

        t_factor = {"yr": 1.0, "kyr": 1e-3, "Myr": 1e-6}.get(t_unit, 1e-3)
        t_plot   = self.t_yr * t_factor

        def _load(qty):
            ix = np.ix_(full_t_idx, Nr_idx)
            if qty == "Sigma_gas":
                return self.data.gas.Sigma[ix][:, self.r_mask]
            elif qty == "Sigma_dust":
                return self.data.dust.Sigma[ix][:, self.r_mask].sum(-1)
            elif qty == "T":
                return self.data.gas.T[ix][:, self.r_mask]
            elif qty == "a_max":
                return self.data.dust.s.max[ix][:, self.r_mask]
            elif qty == "eps":
                Sd = self.data.dust.Sigma[ix][:, self.r_mask].sum(-1)
                Sg = self.data.gas.Sigma[ix][:, self.r_mask]
                return Sd / (Sg + 1e-300)
            elif qty == "St":
                return self.data.dust.St[ix][:, self.r_mask, -1]
            return None

        # Bordes del pcolormesh
        r_edges = _log_edges(self.r_au)         # Y edges (Nr+1,)
        if t_log:
            t_edges = _log_edges(np.maximum(t_plot, t_plot[t_plot > 0].min() * 0.5))
        else:
            dt_lo   = t_plot[1] - t_plot[0]
            dt_hi   = t_plot[-1] - t_plot[-2]
            t_edges = np.concatenate([
                [t_plot[0]  - dt_lo / 2],
                0.5 * (t_plot[:-1] + t_plot[1:]),
                [t_plot[-1] + dt_hi / 2],
            ])                                   # X edges (Nt+1,)

        if figsize is None:
            figsize = (11, 3.5 * n_panels)
        fig, axes = plt.subplots(n_panels, 1, figsize=figsize,
                                 sharex=True, squeeze=False)
        axes = axes[:, 0]

        for ax, qty in zip(axes, quantities):
            Z = _load(qty)
            if Z is None:
                ax.set_visible(False)
                continue

            Z    = np.clip(Z, 1e-30, None)
            logZ = np.log10(Z)                   # (Nt, Nr)
            valid = Z > 1e-20
            vlo  = float(np.percentile(logZ[valid], percentile_range[0]))
            vhi  = float(np.percentile(logZ[valid], percentile_range[1]))

            cmap = cmaps.get(qty, _default_cmaps.get(qty, "viridis"))

            # Transponer (Nt, Nr) → (Nr, Nt): r en Y, t en X
            pcm = ax.pcolormesh(t_edges, r_edges, logZ.T,
                                cmap=cmap, vmin=vlo, vmax=vhi,
                                shading="flat", rasterized=True)

            cb = fig.colorbar(pcm, ax=ax, pad=0.015, aspect=20)
            cb.set_label(_labels.get(qty, qty), fontsize=9)

            ax.set_yscale("log")
            if t_log:
                ax.set_xscale("log")
            ax.set_ylabel(r"$R$ [AU]", fontsize=10)

            # Snowlines como líneas HORIZONTALES r_snow(t)
            for sp_name, style in SNOWLINE_STYLES.items():
                r_series = self._get_rsnow_series(sp_name)
                if r_series is None:
                    continue
                r_plot = np.where(
                    (r_series > self.r_au[0] * 0.5) & (r_series < self.r_au[-1]),
                    r_series, np.nan
                )
                if not np.all(np.isnan(r_plot)):
                    ax.plot(t_plot, r_plot,
                            color=style["color"], ls="--", lw=1.5, alpha=0.85,
                            label=style["label"])

        axes[-1].set_xlabel(f"$t$ [{t_unit}]", fontsize=10)

        handles, lbls = axes[-1].get_legend_handles_labels()
        if handles:
            axes[-1].legend(handles, lbls, fontsize=8,
                            loc="upper right", framealpha=0.5,
                            facecolor="#111111", labelcolor="white",
                            edgecolor="gray")

        fig.suptitle(r"Diagrama Espacio-Tiempo $R(t)$",
                     fontsize=13, fontweight="bold", y=1.01)
        fig.tight_layout()
        self._save(fig, "hovmoller_rt")
        return fig, axes

    # =========================================================================
    # 2. DISTRIBUCIÓN DE TAMAÑO vs RADIO — Sigma(a, r) en una snapshot
    # =========================================================================
    def plot_size_distribution(self, it=-1, cmap="viridis",
                               vmin=None, vmax=None, sigma_spread=0.5):
        """
        Mapa 2D de densidad superficial en ejes (radio, tamaño de grano).

        TriPoD usa 2 bins de superficie (pequeño / grande).  Construimos el
        mapa distribuyendo cada bin como gaussiana en log(a) centrada en:
            bin 0 (pequeño)  →  s.min
            bin 1 (grande)   →  s.max

        Parameters
        ----------
        it : int
            Índice temporal dentro de los snapshots válidos (-1 = último).
        sigma_spread : float
            Ancho en décadas de log(a) de la gaussiana por bin.
        """
        it_abs = _resolve_it(it, self.Nt)

        a_small = self.data.dust.s.min[it_abs, self.r_mask]   # (Nr,) [cm]
        a_large = self.data.dust.s.max[it_abs, self.r_mask]   # (Nr,) [cm]

        Sigma_bins = self.data.dust.Sigma[it_abs, self.r_mask, :]  # (Nr, 2)
        Sig0 = Sigma_bins[:, 0]
        Sig1 = Sigma_bins[:, 1]

        a_lo = np.nanmin(a_small[a_small > 0]) * 0.5
        a_hi = np.nanmax(a_large) * 2.0
        Na_plot    = 300
        log_a_plot = np.linspace(np.log10(a_lo), np.log10(a_hi), Na_plot)
        a_plot     = 10.0 ** log_a_plot

        Z = np.zeros((self.Nr, Na_plot))
        for i in range(self.Nr):
            if Sig0[i] > 0 and a_small[i] > 0:
                lc = np.log10(a_small[i])
                g  = np.exp(-0.5 * ((log_a_plot - lc) / sigma_spread) ** 2)
                Z[i] += Sig0[i] * g / (g.sum() + 1e-300)
            if Sig1[i] > 0 and a_large[i] > 0:
                lc = np.log10(a_large[i])
                g  = np.exp(-0.5 * ((log_a_plot - lc) / sigma_spread) ** 2)
                Z[i] += Sig1[i] * g / (g.sum() + 1e-300)

        Z = np.clip(Z, 1e-30, None)
        t_label = f"{self.t_yr[it_abs]:.2e} yr"

        fig, ax = plt.subplots(figsize=(10, 5))
        valid   = Z > 1e-20
        vlim_lo = vmin if vmin is not None else float(np.percentile(np.log10(Z[valid]), 5))
        vlim_hi = vmax if vmax is not None else float(np.percentile(np.log10(Z[valid]), 99))

        pcm = ax.pcolormesh(_log_edges(self.r_au), _log_edges(a_plot),
                            np.log10(Z).T, cmap=cmap,
                            vmin=vlim_lo, vmax=vlim_hi, shading="flat")
        fig.colorbar(pcm, ax=ax, pad=0.02,
                     label=r"log$_{10}$ $\Sigma(a,r)$ [g cm$^{-2}$]")

        ax.plot(self.r_au, a_large, "w--", lw=1.6, label=r"$a_\mathrm{max}$")
        ax.plot(self.r_au, a_small, "w:",  lw=1.2, label=r"$a_\mathrm{min}$")
        ax.legend(fontsize=8, loc="upper right", framealpha=0.4)

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Radio [AU]")
        ax.set_ylabel(r"Tamaño de grano $a$ [cm]")
        ax.set_title(f"Distribución de tamaño de grano vs radio  |  t = {t_label}", pad=8)
        self._add_snowline_vlines(ax, alpha=0.6, it_abs=it_abs)

        fig.tight_layout()
        self._save(fig, f"size_distribution_it{it_abs:04d}")
        return fig, ax

    # =========================================================================
    # 3. PEBBLE FLUX Ṁ vs (radio, tiempo)
    # =========================================================================
    def plot_pebble_flux(self, cmap="RdYlBu_r", vmin=None, vmax=None):
        """
        Mapa espacio-temporal del flujo másico radial del polvo.
        Ṁ = -2π r Σ_dust v_r  [g/s]  → convertido a M⊕/yr

        dust.v.rad tiene shape (Nt, Nr, 5) con velocidades en:
            [a0, fudge*a1, a1, 0.5*amax, amax]
        Mapeamos:  bin 0 (pequeño) → v[..., 2]  (a1)
                   bin 1 (grande)  → v[..., 4]  (amax)
        """
        rm, tm = self.r_mask, self.t_mask
        Nr_idx = np.arange(self.Nr_full)

        vr_full    = self.data.dust.v.rad       # (Nt_full, Nr_full, 5) o (Nt_full, Nr_full)
        Sig_full   = self.data.dust.Sigma        # (Nt_full, Nr_full, 2)

        vr       = vr_full[np.ix_(np.where(tm)[0], Nr_idx)][:, rm]   # (Nt, Nr, 5)
        Sig_bins = Sig_full[np.ix_(np.where(tm)[0], Nr_idx)][:, rm]  # (Nt, Nr, 2)
        Sigma_dust = Sig_bins.sum(-1)                                  # (Nt, Nr)

        if vr.ndim == 3:
            Na_v = vr.shape[-1]
            if Na_v == 5:
                vr_2bin = np.stack([vr[..., 2], vr[..., 4]], axis=-1)
            elif Na_v == 2:
                vr_2bin = vr
            else:
                vr_mean = vr.mean(axis=-1)
                vr_2bin = np.stack([vr_mean, vr_mean], axis=-1)
            with np.errstate(invalid="ignore", divide="ignore"):
                vr_eff = np.sum(vr_2bin * Sig_bins, axis=-1) / (Sigma_dust + 1e-300)
        else:
            vr_eff = vr

        r2d  = self.r[np.newaxis, :]
        Mdot = -2.0 * np.pi * r2d * Sigma_dust * vr_eff          # [g/s]
        Mdot_abs = np.abs(Mdot * c.year / 5.972e27)               # [M⊕/yr]

        fig, ax = plt.subplots(figsize=(10, 5))
        log_flux = np.log10(np.clip(Mdot_abs, 1e-30, None))    # (Nt, Nr)
        vlim_lo  = vmin if vmin is not None else float(np.percentile(log_flux, 5))
        vlim_hi  = vmax if vmax is not None else float(np.percentile(log_flux, 97))

        # Transponer (Nt, Nr) -> (Nr, Nt) para r en Y, t en X
        pcm = ax.pcolormesh(_log_edges(np.maximum(self.t_yr, 1.0)),
                            _log_edges(self.r_au),
                            log_flux.T, cmap=cmap,
                            vmin=vlim_lo, vmax=vlim_hi, shading="flat")
        fig.colorbar(pcm, ax=ax, pad=0.02,
                     label=r"log$_{10}$ $|\dot{M}_{pebble}|$ [$M_\oplus$ yr$^{-1}$]")

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Tiempo [años]")
        ax.set_ylabel("Radio [AU]")
        ax.set_title(r"Flujo másico de Pebbles $\dot{M}_{pebble}(r, t)$", pad=8)
        self._add_snowline_curves(ax, alpha=0.85)

        fig.tight_layout()
        self._save(fig, "pebble_flux_hovmoller")
        return fig, ax

    # =========================================================================
    # 4. PERFILES DE DISCO — η, St, a_max, Σ en un snapshot
    # =========================================================================
    def plot_profiles(self, it=-1):
        """
        4 paneles en un solo snapshot: η, St (bin grande), a_max, Σ gas+dust.
        """
        it_abs  = _resolve_it(it, self.Nt)
        t_label = f"{self.t_yr[it_abs]:.2e} yr"
        r       = self.r_au

        fig, axes = plt.subplots(2, 2, figsize=(12, 7))
        (ax_eta, ax_st), (ax_amax, ax_sig) = axes

        # η
        eta = self.data.gas.eta[it_abs, self.r_mask]
        ax_eta.plot(r, eta, color="#80DEEA", lw=1.8)
        ax_eta.set_xscale("log")
        ax_eta.set_xlabel("Radio [AU]"); ax_eta.set_ylabel(r"$\eta$")
        ax_eta.set_title(r"Parámetro de gradiente de presión $\eta$", pad=6)
        self._add_snowline_vlines(ax_eta,  alpha=0.6, it_abs=it_abs)

        # St
        St = self.data.dust.St[it_abs, self.r_mask, :]   # (Nr, 2 o 5)
        St_large = St[:, -1]   # bin o tamaño más grande
        ax_st.plot(r, St_large, color="#FFAB91", lw=1.8)
        ax_st.set_xscale("log"); ax_st.set_yscale("log")
        ax_st.set_xlabel("Radio [AU]"); ax_st.set_ylabel(r"$St$")
        ax_st.set_title(r"Número de Stokes $St$ (bin grande)", pad=6)
        self._add_snowline_vlines(ax_st,   alpha=0.6, it_abs=it_abs)

        # a_max
        a_max = self.data.dust.s.max[it_abs, self.r_mask]
        a_min = self.data.dust.s.min[it_abs, self.r_mask]
        ax_amax.plot(r, a_max, color="#CE93D8", lw=1.8, label=r"$a_\mathrm{max}$")
        ax_amax.plot(r, a_min, color="#CE93D8", lw=1.0, ls="--",
                     label=r"$a_\mathrm{min}$")
        ax_amax.set_xscale("log"); ax_amax.set_yscale("log")
        ax_amax.set_xlabel("Radio [AU]"); ax_amax.set_ylabel("Tamaño [cm]")
        ax_amax.set_title(r"Tamaño de grano $a$", pad=6)
        ax_amax.legend(fontsize=8)
        self._add_snowline_vlines(ax_amax, alpha=0.6, it_abs=it_abs)

        # Σ gas + dust
        Sig_d = self.data.dust.Sigma[it_abs, self.r_mask, :].sum(-1)
        Sig_g = self.data.gas.Sigma [it_abs, self.r_mask]
        ax_sig.plot(r, Sig_d, color="#EF9A9A", lw=1.8, label=r"$\Sigma_{dust}$")
        ax_sig.plot(r, Sig_g, color="#90CAF9", lw=1.8, label=r"$\Sigma_{gas}$")
        ax_sig.set_xscale("log"); ax_sig.set_yscale("log")
        ax_sig.set_xlabel("Radio [AU]")
        ax_sig.set_ylabel(r"$\Sigma$ [g cm$^{-2}$]")
        ax_sig.set_title("Densidades superficiales", pad=6)
        ax_sig.legend(fontsize=8)
        self._add_snowline_vlines(ax_sig,  alpha=0.5, it_abs=it_abs)

        fig.suptitle(f"Perfiles de disco  |  t = {t_label}",
                     fontsize=14, fontweight="bold", y=1.01)
        fig.tight_layout()
        self._save(fig, f"profiles_it{it_abs:04d}")
        return fig, axes

    # =========================================================================
    # 5. DENSIDADES SUPERFICIALES DE GAS — por componente
    # =========================================================================
    def plot_gas_components(self, it=-1,
                            ylim=(1e-5, 1e3), legend_outside=True,
                            dumpfile=None):
        """
        Perfil de Sigma_gas de cada componente volátil activo + gas total.

        IMPORTANTE: Las densidades superficiales por componente NO se guardan
        en los archivos HDF5 por defecto. Este método requiere:
          a) El dump file frame.dmp (que si tiene Sigma al momento del pause), O
          b) Haber activado save=True en el pipeline antes de correr.

        Parameters
        ----------
        it : int
            Índice del snapshot HDF5 para el tiempo en el título (-1 = último).
        dumpfile : str, optional
            Ruta al frame.dmp. Si None, busca en datadir/frame.dmp.
        """
        it_abs  = _resolve_it(it, self.Nt)
        t_label = f"{self.t_yr[it_abs]:.2e} yr"

        # Intentar leer desde el dump file
        sim = _load_dumpfile(self.datadir, dumpfile)

        fig, ax = plt.subplots(figsize=(9, 5))
        fallback_colors = list(plt.cm.tab20.colors)
        idx_color = 0
        n_plotted = 0

        if sim is None:
            print("  [!] No se encontró frame.dmp. Para graficar componentes individuales")
            print("      añade esto al pipeline ANTES de correr la simulación:")
            print("          sim.components.H2O.gas.Sigma.save = True")
            print("          sim.components.CO2.gas.Sigma.save = True  (etc.)")
        else:
            # Seguir exactamente el patrón de la documentación de tripodpy:
            #   for name, comp in sim.components.__dict__.items():
            #       if comp.gas._active:
            #           plt.loglog(r/au, comp.gas.Sigma, label=name)
            for name, comp in sim.components.__dict__.items():
                if name.startswith("_") or comp is None:
                    continue
                gas = getattr(comp, "gas", None)
                if gas is None:
                    continue
                if not getattr(gas, "_active", True):
                    continue
                sig = getattr(gas, "Sigma", None)
                if sig is None:
                    continue
                sig_r = np.asarray(sig)[self.r_mask]
                if not np.any(sig_r > 0):
                    continue
                color = COMP_COLORS.get(name, fallback_colors[idx_color % len(fallback_colors)])
                idx_color += 1
                label = COMP_LABELS.get(name, name)
                ax.loglog(self.r_au, np.clip(sig_r, 1e-100, None),
                          color=color, lw=1.5, label=label)
                n_plotted += 1

        if n_plotted == 0:
            print("  [!] Ningún componente de gas graficado.")

        # Gas total desde HDF5 (siempre disponible)
        Sig_gas_total = self.data.gas.Sigma[it_abs, self.r_mask]
        ax.loglog(self.r_au, Sig_gas_total, "k--", lw=1.5, label="Gas total")

        ax.set_xlabel("Radio [AU]")
        ax.set_ylabel(r"$\Sigma_{gas}$ [g cm$^{-2}$]")
        ax.set_title(f"Densidades superficiales de gas por componente  |  t = {t_label}")
        ax.set_ylim(ylim)
        if legend_outside:
            ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=8)
            fig.subplots_adjust(right=0.78)
        else:
            ax.legend(fontsize=8)
        self._add_snowline_vlines(ax, alpha=0.5)
        fig.tight_layout()
        self._save(fig, f"gas_components_it{it_abs:04d}")
        return fig, ax

    # =========================================================================
    # 6. DENSIDADES SUPERFICIALES DE POLVO — por componente
    # =========================================================================
    def plot_dust_components(self, it=-1,
                             ylim=(1e-8, 1e3), legend_outside=True,
                             dumpfile=None):
        """
        Perfil de Sigma_dust de cada componente sólido activo.
        Silicatos en negro como referencia de fondo.

        Requiere el dump file frame.dmp (ver plot_gas_components).
        """
        it_abs  = _resolve_it(it, self.Nt)
        t_label = f"{self.t_yr[it_abs]:.2e} yr"

        sim = _load_dumpfile(self.datadir, dumpfile)

        fig, ax = plt.subplots(figsize=(9, 5))
        fallback_colors = list(plt.cm.tab20.colors)
        idx_color = 0
        n_plotted = 0
        sil_data  = None

        if sim is None:
            print("  [!] No se encontró frame.dmp. Para graficar componentes individuales")
            print("      añade esto al pipeline ANTES de correr la simulación:")
            print("          sim.components.H2O.dust.Sigma.save = True")
            print("          sim.components.CO2.dust.Sigma.save = True  (etc.)")
        else:
            # Patrón exacto de la doc de tripodpy:
            #   for name, comp in sim.components.__dict__.items():
            #       if comp.dust._active:
            #           plt.loglog(r/au, comp.dust.Sigma.sum(-1), label=name)
            for name, comp in sim.components.__dict__.items():
                if name.startswith("_") or comp is None:
                    continue
                dust = getattr(comp, "dust", None)
                if dust is None:
                    continue
                if not getattr(dust, "_active", False):
                    continue
                sig = getattr(dust, "Sigma", None)
                if sig is None:
                    continue
                sig = np.asarray(sig)
                # Sigma puede ser (Nr, Na) o (Nr,)
                if sig.ndim == 2:
                    sig_r = sig[self.r_mask, :].sum(-1)
                else:
                    sig_r = sig[self.r_mask]
                if not np.any(sig_r > 0):
                    continue

                if name == "silicates":
                    sil_data = sig_r
                    continue

                color = COMP_COLORS.get(name, fallback_colors[idx_color % len(fallback_colors)])
                idx_color += 1
                label = COMP_LABELS.get(name, name)
                ax.loglog(self.r_au, np.clip(sig_r, 1e-30, None),
                          color=color, lw=1.5, label=label)
                n_plotted += 1

            if sil_data is not None:
                ax.loglog(self.r_au, np.clip(sil_data, 1e-30, None),
                          "k", lw=2.0, label="Silicatos")
                n_plotted += 1

        if n_plotted == 0:
            print("  [!] Ningún componente de polvo graficado.")

        ax.set_xlabel("Radio [AU]")
        ax.set_ylabel(r"$\Sigma_{dust}$ [g cm$^{-2}$]")
        ax.set_title(f"Densidades superficiales de polvo por componente  |  t = {t_label}")
        ax.set_ylim(ylim)
        if legend_outside:
            ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=8)
            fig.subplots_adjust(right=0.78)
        else:
            ax.legend(fontsize=8)
        self._add_snowline_vlines(ax, alpha=0.5)
        fig.tight_layout()
        self._save(fig, f"dust_components_it{it_abs:04d}")
        return fig, ax


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────
def _log_edges(centers):
    """Genera aristas log-espaciadas a partir de centros para pcolormesh."""
    centers = np.asarray(centers, dtype=float)
    lc      = np.log10(np.maximum(centers, 1e-300))
    edges   = np.empty(len(centers) + 1)
    edges[1:-1] = 10 ** (0.5 * (lc[:-1] + lc[1:]))
    edges[0]    = 10 ** (lc[0]  - (lc[1]  - lc[0])  / 2)
    edges[-1]   = 10 ** (lc[-1] + (lc[-1] - lc[-2]) / 2)
    return edges


def _resolve_it(it, Nt):
    """Convierte índice relativo al subset de snapshots válidos."""
    return int(it % Nt)


def _load_dumpfile(datadir, dumpfile=None):
    """
    Carga el dump file frame.dmp y devuelve el objeto sim completo.

    Los archivos HDF5 de tripodpy NO guardan component.gas.Sigma ni
    component.dust.Sigma porque son arrays Python (no simframe Fields).
    El dump file SÍ contiene el estado completo — simframe lo serializa
    con dill (pickle extendido), incluyendo todos los componentes y sus Sigma.

    Returns None si no se puede cargar.
    """
    if dumpfile is None:
        dumpfile = os.path.join(datadir, "frame.dmp")
    if not os.path.isfile(dumpfile):
        return None
    try:
        import dill
        with open(dumpfile, "rb") as f:
            sim = dill.load(f)
        print(f"  → Dump file cargado: {dumpfile}")
        return sim
    except Exception as e:
        print(f"  [!] No se pudo cargar dump file: {e}")
        return None


def _extract_dust_sigma(dust_ns, it_abs, r_mask):
    """
    Extrae la densidad superficial de polvo de un componente de forma defensiva.
    tripodpy puede guardar la data en 'Sigma' o en 'value' dependiendo de si el
    componente fue declarado como dust_active=True o como tracer.

    Returns None si no encuentra datos.
    """
    for attr in ("Sigma", "value"):
        sig = getattr(dust_ns, attr, None)
        if sig is None:
            continue
        try:
            sig = np.asarray(sig)
            if sig.ndim == 3:
                return sig[it_abs, r_mask, :].sum(-1)
            elif sig.ndim == 2:
                return sig[it_abs, r_mask]
            elif sig.ndim == 1:
                return sig[r_mask]
        except (IndexError, TypeError):
            continue
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Modo standalone
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    datadir = sys.argv[1] if len(sys.argv) > 1 else "output_test_pipeline"
    savedir = sys.argv[2] if len(sys.argv) > 2 else "diagnostics_output"

    d = SnowlineDiagnostics(datadir, savedir=savedir, r_trim=0.93, t_min_yr=100.0)

    print("\n[1/6] Hovmöller (eps)…")
    d.plot_hovmoller(quantity="eps")

    print("[2/6] Hovmöller (Sigma_dust)…")
    d.plot_hovmoller(quantity="Sigma_dust")

    print("[3/6] Distribución de tamaño (último snapshot)…")
    d.plot_size_distribution(it=-1)

    print("[4/6] Pebble flux…")
    d.plot_pebble_flux()

    print("[5/6] Perfiles (η, St, a_max, Σ)…")
    d.plot_profiles(it=-1)

    print("[6a/6] Densidades de gas por componente…")
    d.plot_gas_components(it=-1)

    print("[6b/6] Densidades de polvo por componente…")
    d.plot_dust_components(it=-1)

    print("[7/7] Hovmöller R(t) — 3 paneles (Sigma_gas, T, a_max)…")
    d.plot_hovmoller_rt(
        quantities=["Sigma_gas", "T", "a_max"],
        t_unit="kyr",
        t_log=False,
    )

    plt.show()
    print("\n¡Listo!")
