"""
PA3Py/plot_pa3.py
=================
Modulo de diagnosticos visuales para los resultados de PebbleAccretionModule3.

Plots disponibles:
  1. plot_waterworld_map()     -- M_final vs r, color=f_H2O, con M_iso(r)
  2. plot_composition()        -- Composicion quimica apilada vs tiempo
  3. plot_hovmoller_mass()     -- Mapa 2D masa(r_emb, t) y f_H2O(r_emb, t)
  4. plot_with_disk_temperature() -- M_final vs r, color=T(r) del disco

Uso:
    from PA3Py.plot_pa3 import PA3Diagnostics
    from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

    pam     = PebbleAccretionModule3.from_datadir(DATADIR)
    results = pam.run_growth([1, 2, 3, 5, 10, 20])

    diag = PA3Diagnostics(pam, results, savedir="figs_pa3")
    diag.plot_waterworld_map()
    diag.plot_composition()
    diag.plot_hovmoller_mass()
    diag.plot_with_disk_temperature()
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ── Estetica global ────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "axes.titlesize":  13,
    "axes.labelsize":  11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "axes.grid":       True,
    "grid.alpha":      0.2,
    "figure.dpi":      130,
})

COMP_COLORS = {
    "H2O":       "#4FC3F7",
    "CO2":       "#80CBC4",
    "silicates": "#EF9A9A",
    "CO":        "#FFCC80",
}
SNOWLINE_STYLES = {
    "H2O": {"color": "#4FC3F7", "label": "H2O snowline", "T": 150.0},
    "CO2": {"color": "#80CBC4", "label": "CO2 snowline", "T": 70.0},
    "CO":  {"color": "#FFCC80", "label": "CO snowline",  "T": 25.0},
}

_YR = 3.156e7   # s per yr


def _t_yr(t_s):
    return np.asarray(t_s) / _YR


def _has_data(hist):
    return hist is not None and len(hist) > 1


# ──────────────────────────────────────────────────────────────────────────────
class PA3Diagnostics:
    """
    Herramienta de diagnostico para resultados de PebbleAccretionModule3.

    Parameters
    ----------
    pam : PebbleAccretionModule3
        Objeto ya cargado con from_datadir(). Contiene todos los datos
        del disco (gas, polvo, snowlines) en memoria.
    results : dict
        Salida de pam.run_growth(). Formato: {r_au: array(Nt, 7)}.
    savedir : str or None
        Directorio donde guardar las figuras en PDF. None = no guarda.
    r_trim : float
        Fraccion del grid radial a conservar. Default 0.93.
    """

    def __init__(self, pam, results, savedir=None, r_trim=0.93):
        self.pam     = pam
        self.results = results
        self.savedir = savedir
        self.r_trim  = r_trim

        self.r_locs = sorted(results.keys())

        Nr_keep     = max(2, int(len(pam.r) * r_trim))
        self.r      = pam.r[:Nr_keep]
        self.r_au   = self.r / pam.AU
        self.t_yr   = _t_yr(pam.times)

        # Rango de embriones con datos (para limites de ejes)
        locs_with_data = [r for r in self.r_locs if _has_data(results[r])]
        self.r_min = min(locs_with_data) if locs_with_data else self.r_au[0]
        self.r_max = max(locs_with_data) if locs_with_data else self.r_au[-1]
        self.xlim  = (self.r_min * 0.9, self.r_max + 2.0)

        print(f"[PA3Diagnostics] {len(self.r_locs)} embriones | "
              f"{pam.Nt} snapshots | "
              f"r_emb = {self.r_min:.1f} -- {self.r_max:.1f} AU")
        print(f"  -> xlim radial para plots: {self.xlim[0]:.1f} -- {self.xlim[1]:.1f} AU")

        # M_iso(r) del ultimo snapshot para plots radiales
        self._compute_iso_profile()

    # ── Helpers ───────────────────────────────────────────────────────────────
    def _compute_iso_profile(self):
        last = self.pam.Nt - 1
        self.M_iso_profile = np.array([
            self.pam._isolation_mass(r * self.pam.AU, last) / self.pam.M_EARTH
            for r in self.r_au
        ])

    def _active_snowlines(self):
        """
        Devuelve el subconjunto de SNOWLINE_STYLES cuyas claves existen
        en pam.rsnow (es decir, las especies realmente simuladas).
        Evita aparicion de snowlines fantasma de especies no activas.
        """
        return {
            sp: style
            for sp, style in SNOWLINE_STYLES.items()
            if sp in self.pam.rsnow and self.pam.rsnow[sp] is not None
        }

    def _get_rsnow_series(self, sp):
        arr = self.pam.rsnow.get(sp, None)
        if arr is None:
            return None
        return np.where(np.isfinite(arr), arr / self.pam.AU, np.nan)

    def _add_snowline_vlines(self, ax, it=-1, alpha=0.75):
        it_use = it % self.pam.Nt
        for sp, style in self._active_snowlines().items():
            rs = self._get_rsnow_series(sp)
            if rs is None:
                continue
            rv = rs[it_use]
            if np.isnan(rv) or rv < self.r_au[0] or rv > self.r_au[-1]:
                continue
            ax.axvline(rv, color=style["color"], ls="--", lw=1.3,
                       alpha=alpha, label=f"{sp} snowline ({rv:.1f} AU)")

    def _add_snowline_curves(self, ax, swap_axes=False, alpha=0.85):
        for sp, style in self._active_snowlines().items():
            rs = self._get_rsnow_series(sp)
            if rs is None:
                continue
            rs_plot = np.where(
                (rs > self.r_au[0]) & (rs < self.r_au[-1]), rs, np.nan
            )
            if np.all(np.isnan(rs_plot)):
                continue
            if swap_axes:
                ax.plot(rs_plot, self.t_yr,
                        color=style["color"], ls="--", lw=1.5,
                        alpha=alpha, label=style["label"])
            else:
                ax.plot(self.t_yr, rs_plot,
                        color=style["color"], ls="--", lw=1.5,
                        alpha=alpha, label=style["label"])

    def _select_embryos(self, r_selected=None, max_n=6):
        if r_selected is not None:
            return [r for r in r_selected if r in self.results]
        locs = [r for r in self.r_locs if _has_data(self.results[r])]
        if len(locs) <= max_n:
            return locs
        idx = np.round(np.linspace(0, len(locs) - 1, max_n)).astype(int)
        return [locs[i] for i in idx]

    def _save(self, fig, name):
        if self.savedir:
            os.makedirs(self.savedir, exist_ok=True)
            path = os.path.join(self.savedir, name + ".pdf")
            fig.savefig(path, bbox_inches="tight")
            print(f"  -> Guardado: {path}")

    # ==========================================================================
    # 1. Waterworld Map — M_final vs r, color = f_H2O
    # ==========================================================================
    def plot_waterworld_map(self, cmap_name="RdYlBu"):
        """
        Masa final vs radio del embrion, coloreado por fraccion de H2O.

        Plot principal del proyecto. Conecta estructura del disco con
        el resultado planetario final.

        Sobreimprime:
          - M_iso(r): curva de masa de aislamiento del ultimo snapshot
          - Snowlines verticales (H2O, CO2, CO)

        Eje X recortado a: (r_min_emb * 0.9, r_max_emb + 2 AU).
        """
        r_arr    = []
        M_arr    = []
        fH2O_arr = []

        for r_au in self.r_locs:
            hist = self.results[r_au]
            if not _has_data(hist):
                continue
            _, M, H2O, CO2, sil, _, M_iso = hist[-1]
            r_arr.append(r_au)
            M_arr.append(M / self.pam.M_EARTH)
            fH2O_arr.append(100 * H2O / (M + 1e-30))

        if not r_arr:
            print("[plot_waterworld_map] Sin datos.")
            return None, None

        r_arr    = np.array(r_arr)
        M_arr    = np.array(M_arr)
        fH2O_arr = np.array(fH2O_arr)

        fig, ax = plt.subplots(figsize=(11, 6))

        norm = mcolors.Normalize(vmin=0, vmax=min(50, fH2O_arr.max() + 1))
        sc   = ax.scatter(r_arr, M_arr, c=fH2O_arr, cmap=cmap_name,
                          norm=norm, s=70, zorder=5, edgecolors="none")
        cb   = fig.colorbar(sc, ax=ax, pad=0.02)
        cb.set_label("Fraccion H2O final [%]", fontsize=10)

        # Perfil M_iso(r) — mascara al rango de embriones
        r_mask = (self.r_au >= self.xlim[0]) & (self.r_au <= self.xlim[1])
        ax.plot(self.r_au[r_mask], self.M_iso_profile[r_mask],
                color="white", lw=2.2, ls="--", alpha=0.9,
                label="M_iso(r)", zorder=6)

        self._add_snowline_vlines(ax, it=-1)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(self.xlim)
        ax.set_ylim(bottom=1e-4)
        ax.set_xlabel("Radio del embrion [AU]")
        ax.set_ylabel("Masa final [ME]")
        ax.set_title("Waterworld Map: Masa final y fraccion de H2O vs radio", pad=8)
        ax.legend(fontsize=9, framealpha=0.6)

        fig.tight_layout()
        self._save(fig, "pa3_waterworld_map")
        return fig, ax

    # ==========================================================================
    # 2. Composicion quimica apilada vs tiempo
    # ==========================================================================
    def plot_composition(self, r_selected=None, max_n=4):
        """
        Stacked area plot de la composicion quimica acumulada por embrion.

        Para cada embrion seleccionado, un panel con areas apiladas de:
          H2O / CO2 / silicatos (como fraccion del total acretado).

        Marca con una linea vertical cuando la snowline de H2O cruza
        la posicion del embrion.

        Parameters
        ----------
        r_selected : list of float or None
            Posiciones en AU. None = distribucion automatica.
        max_n : int
            Numero maximo de paneles. Default 4.
        """
        locs = self._select_embryos(r_selected, max_n)
        if not locs:
            print("[plot_composition] Sin datos.")
            return None, None

        n    = len(locs)
        fig, axes = plt.subplots(n, 1, figsize=(10, 3.5 * n),
                                 sharex=True, squeeze=False)
        axes = axes[:, 0]

        rsnow_h2o = self._get_rsnow_series("H2O")

        for ax, r_au in zip(axes, locs):
            hist = self.results[r_au]
            if not _has_data(hist):
                ax.set_visible(False)
                continue

            t    = _t_yr(hist[:, 0])
            M    = hist[:, 1]
            H2O  = hist[:, 2] / (M + 1e-30)
            CO2  = hist[:, 3] / (M + 1e-30)
            sil  = hist[:, 4] / (M + 1e-30)
            rest = np.clip(1.0 - H2O - CO2 - sil, 0, 1)

            ax.stackplot(t,
                         sil * 100, CO2 * 100, H2O * 100, rest * 100,
                         labels=["Silicatos", "CO2", "H2O", "Semilla"],
                         colors=[COMP_COLORS["silicates"], COMP_COLORS["CO2"],
                                 COMP_COLORS["H2O"], "#BDBDBD"],
                         alpha=0.85)

            # Tiempo en que la snowline H2O cruza la posicion del embrion
            if rsnow_h2o is not None:
                rsnow_interp = np.interp(t, self.t_yr, rsnow_h2o,
                                         left=rsnow_h2o[0], right=rsnow_h2o[-1])
                # Cruce: de interior (rsnow > r) a exterior (rsnow < r)
                inside  = rsnow_interp > r_au
                entries = np.where(np.diff(inside.astype(int)))[0]
                for ci in entries:
                    if ci < len(t):
                        ax.axvline(t[ci], color=COMP_COLORS["H2O"], lw=1.5,
                                   ls="--", alpha=0.9,
                                   label="Snowline H2O cruza embrion")

            ax.set_ylim(0, 100)
            ax.set_ylabel("Fraccion [%]")
            ax.set_title(f"Composicion acumulada  r = {r_au:.2f} AU", pad=6)
            ax.legend(loc="upper left", fontsize=7, framealpha=0.6, ncol=2)
            ax.set_xscale("log")
            if len(t) > 1:
                ax.set_xlim(float(t[t > 0].min()) if (t > 0).any() else float(t[0]),
                            float(t[-1]))

        axes[-1].set_xlabel("Tiempo [yr]")
        fig.suptitle("Evolucion de composicion quimica por embrion",
                     fontsize=14, fontweight="bold", y=1.01)
        fig.tight_layout()
        self._save(fig, "pa3_composition")
        return fig, axes

    # ==========================================================================
    # 3. Hovmoller de crecimiento — mapa 2D (r_emb, t)
    # ==========================================================================
    def plot_hovmoller_mass(self, cmap_mass="magma", cmap_comp="RdYlBu"):
        """
        Mapa espacio-tiempo 2D del crecimiento planetario.

        Panel superior: log10(M_core) en color vs (t, r_emb).
        Panel inferior: fraccion H2O [%] en color vs (t, r_emb).

        Eje Y limitado al rango de embriones con datos.
        Sobreimprime curvas de snowline r_snow(t).

        Requiere grilla densa de embriones (>= 10) para buen resultado visual.

        Parameters
        ----------
        cmap_mass : str
            Colormap para log M_core. Default "magma".
        cmap_comp : str
            Colormap para f_H2O. Default "RdYlBu".
        """
        r_list = [r for r in self.r_locs if _has_data(self.results[r])]
        if len(r_list) < 3:
            print("[plot_hovmoller_mass] Se necesitan >= 3 embriones con datos.")
            return None, None

        t_common = self.t_yr
        Nr_emb   = len(r_list)
        M_grid    = np.full((Nr_emb, len(t_common)), np.nan)
        fH2O_grid = np.full((Nr_emb, len(t_common)), np.nan)

        for ir, r_au in enumerate(r_list):
            hist = self.results[r_au]
            t_h  = _t_yr(hist[:, 0])
            M_h  = hist[:, 1] / self.pam.M_EARTH
            fH2O = 100 * hist[:, 2] / (hist[:, 1] + 1e-30)

            M_grid[ir]    = np.interp(t_common, t_h, M_h,
                                       left=M_h[0], right=M_h[-1])
            fH2O_grid[ir] = np.interp(t_common, t_h, fH2O,
                                       left=fH2O[0], right=fH2O[-1])

        r_emb_arr = np.array(r_list)
        y_lim     = (r_emb_arr[0] * 0.95, r_emb_arr[-1] * 1.05)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

        # Panel 1: log M_core
        logM = np.log10(np.clip(M_grid, 1e-5, None))
        pcm1 = ax1.pcolormesh(t_common, r_emb_arr, logM,
                               cmap=cmap_mass, shading="auto")
        cb1  = fig.colorbar(pcm1, ax=ax1, pad=0.01)
        cb1.set_label("log$_{10}$ M_core [ME]")
        ax1.set_ylabel("r_emb [AU]")
        ax1.set_title("Evolucion de masa del planeta M(r, t)", pad=6)
        ax1.set_yscale("log")
        ax1.set_ylim(y_lim)
        if len(t_common) > 1:
            t_pos = t_common[t_common > 0]
            ax1.set_xlim(float(t_pos.min()) if t_pos.size else float(t_common[0]),
                         float(t_common[-1]))
        self._add_snowline_curves(ax1, swap_axes=False)
        ax1.legend(fontsize=7, framealpha=0.4, loc="upper right")

        # Panel 2: f_H2O
        fmax = min(40, np.nanpercentile(fH2O_grid, 98))
        pcm2 = ax2.pcolormesh(t_common, r_emb_arr, fH2O_grid,
                               cmap=cmap_comp, vmin=0, vmax=fmax,
                               shading="auto")
        cb2  = fig.colorbar(pcm2, ax=ax2, pad=0.01)
        cb2.set_label("Fraccion H2O [%]")
        ax2.set_ylabel("r_emb [AU]")
        ax2.set_xlabel("Tiempo [yr]")
        ax2.set_title("Evolucion de fraccion de H2O f(r, t)", pad=6)
        ax2.set_yscale("log")
        ax2.set_ylim(y_lim)
        if len(t_common) > 1:
            t_pos = t_common[t_common > 0]
            ax2.set_xlim(float(t_pos.min()) if t_pos.size else float(t_common[0]),
                         float(t_common[-1]))
        self._add_snowline_curves(ax2, swap_axes=False)

        ax1.set_xscale("log")
        ax2.set_xscale("log")

        fig.tight_layout()
        self._save(fig, "pa3_hovmoller")
        return fig, (ax1, ax2)

    # ==========================================================================
    # 4. Temperatura del disco vs masa final — color = T(r)
    # ==========================================================================
    def plot_with_disk_temperature(self, it=-1):
        """
        M_final vs radio, con color = T(r) del disco en el snapshot 'it'.

        Conecta la termodinamica del disco con el resultado planetario.
        Permite ver directamente que zonas termicas favorecen la acrecion.

        Sobreimprime M_iso(r) y snowlines.
        Eje X recortado a: (r_min_emb * 0.9, r_max_emb + 2 AU).
w
        Parameters
        ----------
        it : int
            Snapshot del disco para T(r). Default -1 (ultimo).
        """
        it_use  = it % self.pam.Nt
        Nr_keep = len(self.r)
        T_disk  = self.pam.gas["T"][it_use, :Nr_keep]

        r_arr = []
        M_arr = []
        T_arr = []

        for r_au in self.r_locs:
            hist = self.results[r_au]
            if not _has_data(hist):
                continue
            _, M, *_ = hist[-1]
            r_arr.append(r_au)
            M_arr.append(M / self.pam.M_EARTH)
            T_arr.append(float(np.interp(
                np.log(r_au * self.pam.AU),
                np.log(self.pam.r[:Nr_keep]),
                T_disk
            )))

        if not r_arr:
            print("[plot_with_disk_temperature] Sin datos.")
            return None, None

        r_arr = np.array(r_arr)
        M_arr = np.array(M_arr)
        T_arr = np.array(T_arr)

        fig, ax = plt.subplots(figsize=(11, 6))

        norm = mcolors.LogNorm(vmin=max(T_arr.min(), 10), vmax=T_arr.max())
        sc   = ax.scatter(r_arr, M_arr, c=T_arr, cmap="inferno",
                          norm=norm, s=70, zorder=5, edgecolors="none")
        cb   = fig.colorbar(sc, ax=ax, pad=0.02)
        cb.set_label("T disco [K]", fontsize=10)

        # M_iso(r) en el rango de embriones
        r_mask = (self.r_au >= self.xlim[0]) & (self.r_au <= self.xlim[1])
        ax.plot(self.r_au[r_mask], self.M_iso_profile[r_mask],
                color="white", lw=2.2, ls="--", alpha=0.9,
                label="M_iso(r)", zorder=6)

        self._add_snowline_vlines(ax, it=it_use)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(self.xlim)
        ax.set_ylim(bottom=1e-4)
        ax.set_xlabel("Radio del embrion [AU]")
        ax.set_ylabel("Masa final [ME]")
        ax.set_title(
            f"Masa final vs temperatura del disco  (t = {self.t_yr[it_use]:.2e} yr)",
            pad=8
        )
        ax.legend(fontsize=9, framealpha=0.6)

        fig.tight_layout()
        self._save(fig, "pa3_mass_vs_temperature")
        return fig, ax

    # ==========================================================================
    # Conveniencia: todos los plots activos
    # ==========================================================================
    def plot_all(self, r_selected=None):
        """
        Genera los 4 plots de diagnostico en secuencia.

        Parameters
        ----------
        r_selected : list of float or None
            Embriones para los plots temporales. None = automatico.
        """
        print("[PA3Diagnostics] Generando plots...")
        self.plot_waterworld_map()
        self.plot_composition(r_selected)
        self.plot_hovmoller_mass()
        self.plot_with_disk_temperature()
        print("[PA3Diagnostics] Listo.")
