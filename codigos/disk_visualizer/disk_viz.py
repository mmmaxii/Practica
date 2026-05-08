"""
disk_viz.py — Visualizador 2D del disco protoplanetario por revolución
======================================================================

Toma los snapshots HDF5 de un run de tripodpy y genera una animación 2D
del disco completo (cara frontal) revolviendo los perfiles 1D en azimut.

Campos visualizados en 4 paneles:
  ┌──────────────────┬──────────────────┐
  │  Σ_gas (g/cm²)   │  Σ_dust (g/cm²)  │
  ├──────────────────┼──────────────────┤
  │  Temperatura (K) │  Fracción H₂O    │
  └──────────────────┴──────────────────┘

La snowline H₂O se marca como contorno blanco punteado.
El gap planetario en α(r) se infiere del perfil de alpha.

Uso
---
    python disk_viz.py --datadir ../data_1myr/single_jup_3.0au --output frames/
    python disk_viz.py --datadir ../data_1myr/single_jup_3.0au --gif disco.gif
    python disk_viz.py --datadir ../data_1myr/single_jup_3.0au --interactive

Requisitos
----------
    pip install h5py numpy matplotlib Pillow tqdm
"""

import os
import sys
import glob

# Forzar UTF-8 en la consola de Windows (evita UnicodeEncodeError cp1252)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import argparse
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ═══════════════════════════════════════════════════════════════════════════════
# Constantes
# ═══════════════════════════════════════════════════════════════════════════════
AU  = 1.49598e13   # cm
YR  = 3.15576e7    # s
PC  = 3.085678e18  # cm

# Resolución angular del disco 2D
N_THETA = 720   # número de puntos azimutales (mayor → más suave)


# ═══════════════════════════════════════════════════════════════════════════════
# Lectura de HDF5
# ═══════════════════════════════════════════════════════════════════════════════

def load_snapshot(path):
    """Carga un snapshot HDF5 y devuelve un dict con los arrays relevantes."""
    with h5py.File(path, "r") as f:
        r       = f["grid/r"][:]            # cm, shape (Nr,)
        t       = float(f["t"][()])         # s
        sig_gas = f["gas/Sigma"][:]         # g/cm², (Nr,)
        T       = f["gas/T"][:]             # K,     (Nr,)
        alpha   = f["gas/alpha"][:]         # adim,  (Nr,)

        # Sigma polvo total (suma sobre species × bins)
        sig_dust_2d = f["dust/Sigma"][:]    # (Nr, Nspec)
        sig_dust    = sig_dust_2d.sum(axis=1)   # (Nr,)

        # Sigma H2O dust (del campo personalizado del grid)
        if "grid/SigmaDust_H2O" in f:
            sig_ice = f["grid/SigmaDust_H2O"][:]
        else:
            # fallback: componente H2O si existe
            try:
                sig_ice = f["components/H2O/dust/Sigma"][:].sum(axis=1)
            except KeyError:
                sig_ice = np.zeros_like(sig_dust)

        # Snowline H2O
        if "grid/rsnow_H2O" in f:
            rsnow = float(f["grid/rsnow_H2O"][()])
        else:
            rsnow = np.nan

        # Tamaño máximo de pebbles [cm]
        if "dust/s/max" in f:
            s_max = f["dust/s/max"][:]
        else:
            s_max = None

    # Fracción de hielo (H2O dust / dust total) — clamped a [0,1]
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_ice = np.where(sig_dust > 0, sig_ice / sig_dust, 0.0)
    frac_ice = np.clip(frac_ice, 0.0, 1.0)

    return dict(
        r=r, t=t,
        sig_gas=sig_gas, sig_dust=sig_dust, T=T, alpha=alpha,
        frac_ice=frac_ice, rsnow=rsnow, s_max=s_max,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Revolución 1D → 2D
# ═══════════════════════════════════════════════════════════════════════════════

def revolve(r_1d, values_1d, n_theta=N_THETA):
    """
    Convierte un perfil radial 1D en una imagen 2D por revolución azimutal.

    Devuelve (R, Theta, Z) en coordenadas polares, listos para pcolormesh.
    Z tiene shape (n_theta, Nr).
    """
    theta = np.linspace(0, 2 * np.pi, n_theta + 1)  # bordes de celdas
    Z = np.tile(values_1d, (n_theta, 1))             # (n_theta, Nr)
    return r_1d, theta, Z


# ═══════════════════════════════════════════════════════════════════════════════
# Figura principal
# ═══════════════════════════════════════════════════════════════════════════════

# ─── Colores del tema oscuro ──────────────────────────────────────────────────
BG     = "#08090d"
FG     = "#d8d9e8"
GRID_C = "#2a2d3a"
CMAP_GAS  = "magma"
CMAP_DUST = "YlOrBr_r"


def _polar_fill(ax, r, theta_edges, Z, cmap, norm, add_star=True):
    """
    Dibuja el disco 2D usando pcolormesh en coordenadas cartesianas.
    theta_edges : (n_theta+1,)  radians
    Z           : (n_theta, Nr)
    Retorna el objeto pcolormesh (pcm).
    """
    r_au = r / AU
    r_edges_au = np.concatenate([
        [r_au[0] * np.sqrt(r_au[0] / r_au[1])],
        np.sqrt(r_au[:-1] * r_au[1:]),
        [r_au[-1] * np.sqrt(r_au[-1] / r_au[-2])],
    ])
    T_grid, R_grid = np.meshgrid(theta_edges, r_edges_au, indexing="ij")
    X = R_grid * np.cos(T_grid)
    Y = R_grid * np.sin(T_grid)
    pcm = ax.pcolormesh(X, Y, Z, cmap=cmap, norm=norm,
                        shading="flat", rasterized=True, zorder=2)
    if add_star:
        ax.add_patch(plt.Circle((0, 0), r_edges_au[0] * 0.5,
                                color="#fffbe6", zorder=10))
        ax.add_patch(plt.Circle((0, 0), r_edges_au[0] * 1.2,
                                color="#fffbe6", alpha=0.12, zorder=9))
    return pcm


def _setup_disk_ax(ax, r_max_au):
    """
    Configura un eje cartesiano para la vista del disco:
    - Fondo oscuro
    - Ticks y grid en AU
    - Ejes x/y etiquetados
    """
    ax.set_facecolor(BG)
    ax.set_xlim(-r_max_au, r_max_au)
    ax.set_ylim(-r_max_au, r_max_au)
    ax.set_aspect("equal")

    # Calcular step para tener siempre ~10 ticks en el eje, redondeado a valor "bonito"
    raw_step = 2 * r_max_au / 10.0
    for nice in [1, 2, 5, 10, 20, 25, 50, 100, 200]:
        if nice >= raw_step:
            step = nice
            break
    else:
        step = int(round(raw_step / 100) * 100)
    ticks_all = np.arange(-r_max_au, r_max_au + step, step, dtype=int)
    ticks = ticks_all[(ticks_all >= -r_max_au) & (ticks_all <= r_max_au)]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([str(v) for v in ticks], fontsize=11, color=FG)
    ax.set_yticklabels([str(v) for v in ticks], fontsize=11, color=FG)
    ax.set_xlabel("x  [AU]", color=FG, fontsize=13, labelpad=6)
    ax.set_ylabel("y  [AU]", color=FG, fontsize=13, labelpad=6)

    # Grid cartesiano
    ax.grid(True, color=GRID_C, lw=0.6, linestyle="--", alpha=0.7, zorder=1)

    # Spines
    for sp in ax.spines.values():
        sp.set_edgecolor(GRID_C)
        sp.set_linewidth(0.8)

    # Tick marks
    ax.tick_params(axis="both", colors=FG, length=4, width=0.9)


def _safe_lognorm(arr, vmin_frac=1e-4):
    pos = arr[arr > 0]
    if pos.size == 0:
        return mcolors.LogNorm(vmin=1e-10, vmax=1.0)
    vmax = float(np.nanmax(pos))
    vmin = float(max(np.nanmin(pos), vmax * vmin_frac))
    return mcolors.LogNorm(vmin=vmin, vmax=vmax)


def make_frame(snap, r_max_au=40, fig=None, axes=None, cbars=None, norms=None):
    """
    Genera (o actualiza) un frame de 2 paneles: Sigma_gas | Sigma_dust.

    Parametros
    ----------
    snap           : dict devuelto por load_snapshot()
    r_max_au       : radio maximo del plot [AU]
    fig, axes, cbars, norms : objetos previos para reutilizar en GIF/interactivo.
                              Pasar de vuelta el retorno de cada llamada.

    Retorna
    -------
    (fig, axes, cbars, norms)
    """
    r    = snap["r"]
    mask = r / AU <= r_max_au
    r_m  = r[mask]

    sig_gas  = snap["sig_gas"][mask]
    sig_dust = snap["sig_dust"][mask]
    rsnow    = snap["rsnow"]
    t_yr     = snap["t"] / YR

    n_theta = N_THETA
    theta_e = np.linspace(0, 2 * np.pi, n_theta + 1)

    # Norms: se calculan UNA sola vez (primer frame) y se reutilizan
    if norms is None:
        norms = [
            _safe_lognorm(sig_gas),
            _safe_lognorm(sig_dust, 1e-6),
        ]
    norm_gas, norm_dust = norms

    panels = [
        (sig_gas,  norm_gas,  CMAP_GAS,  r"$\Sigma_{\rm gas}$   [g cm$^{-2}$]"),
        (sig_dust, norm_dust, CMAP_DUST, r"$\Sigma_{\rm dust}$  [g cm$^{-2}$]"),
    ]

    # ── Crear figura la primera vez ───────────────────────────────────────────
    if fig is None:
        fig = plt.figure(figsize=(13, 6.8), facecolor=BG)
        fig.subplots_adjust(left=0.06, right=0.94, top=0.87,
                            bottom=0.10, wspace=0.18)
        axes  = [fig.add_subplot(1, 2, i + 1) for i in range(2)]
        cbars = [None, None]
    else:
        for ax in axes:
            ax.cla()

    # ── Dibujar cada panel ────────────────────────────────────────────────────
    for i, (ax, (vals, norm, cmap, label)) in enumerate(zip(axes, panels)):
        _setup_disk_ax(ax, r_max_au)

        Z   = np.tile(vals, (n_theta, 1))
        pcm = _polar_fill(ax, r_m, theta_e, Z, cmap, norm, add_star=True)

        # Colorbar: crear la primera vez con la norm fija, solo actualizar despues
        if cbars[i] is None:
            cb = fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.015,
                              aspect=22, shrink=0.88)
            cb.ax.tick_params(colors=FG, labelsize=11, length=4)
            cb.outline.set_edgecolor(GRID_C)
            cb.ax.yaxis.set_tick_params(labelcolor=FG, labelsize=11)
            cbars[i] = cb
        else:
            # La norm ya esta fijada en norms; update_normal solo redirige el
            # mappable al nuevo pcm manteniendo el rango del colorbar estatico.
            cbars[i].update_normal(pcm)

        ax.set_title(label, color=FG, fontsize=14, pad=10)

        # Snowline H2O
        if np.isfinite(rsnow) and 0 < rsnow / AU <= r_max_au:
            snow_au = rsnow / AU
            th = np.linspace(0, 2 * np.pi, 500)
            ax.plot(snow_au * np.cos(th), snow_au * np.sin(th),
                    ls="--", lw=1.4, color="white", alpha=0.8, zorder=15)
            ax.text(0, snow_au + r_max_au * 0.04,
                    f"H$_2$O  {snow_au:.1f} AU",
                    color="white", fontsize=9.5, ha="center", va="bottom",
                    alpha=0.9, zorder=16)

    # ── Titulo global ─────────────────────────────────────────────────────────
    if t_yr < 1e4:
        t_str = f"{t_yr:.2e} yr"
    elif t_yr < 1e6:
        t_str = f"{t_yr / 1e3:.1f} kyr"
    else:
        t_str = f"{t_yr / 1e6:.3f} Myr"

    fig.suptitle(f"Disco protoplanetario  —  t = {t_str}",
                 color=FG, fontsize=14, fontweight="bold", y=0.96)

    return fig, axes, cbars, norms


# ═══════════════════════════════════════════════════════════════════════════════
# Modos de ejecución
# ═══════════════════════════════════════════════════════════════════════════════

def get_hdf5_list(datadir):
    files = sorted(glob.glob(os.path.join(datadir, "data*.hdf5")))
    if not files:
        raise FileNotFoundError(f"No se encontraron HDF5 en {datadir}")
    return files


def run_interactive(datadir, r_max_au=40, step=1):
    """Modo interactivo: flechas ← → para navegar snapshots."""
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

    files = get_hdf5_list(datadir)
    N = len(files)
    state = {"idx": 0}

    snap = load_snapshot(files[0])
    fig, axes, cbars, norms = make_frame(snap, r_max_au=r_max_au)
    plt.ion()
    plt.show()

    def on_key(event):
        if event.key in ("right", "d", " "):
            state["idx"] = min(state["idx"] + step, N - 1)
        elif event.key in ("left", "a"):
            state["idx"] = max(state["idx"] - step, 0)
        elif event.key == "q":
            plt.close(fig)
            return
        snap = load_snapshot(files[state["idx"]])
        make_frame(snap, r_max_au=r_max_au, fig=fig, axes=axes,
                   cbars=cbars, norms=norms)
        title = os.path.basename(files[state["idx"]])
        print(f"\r  [{state['idx']+1}/{N}] {title}  "
              f"t = {snap['t']/YR:.3e} yr", end="", flush=True)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)
    print(f"\nCargando {N} snapshots de: {datadir}")
    print("Controles: -> avanzar  |  <- retroceder  |  q salir")
    plt.ioff()
    plt.show()




def run_save_frames(datadir, output_dir, r_max_au=40, every=1, dpi=140):
    """Guarda cada snapshot como PNG en output_dir."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False

    files = get_hdf5_list(datadir)
    files = files[::every]
    os.makedirs(output_dir, exist_ok=True)

    print(f"Guardando {len(files)} frames en {output_dir}/ ...")
    iterable = tqdm(enumerate(files), total=len(files)) if use_tqdm else enumerate(files)

    fig, axes, cbars, norms = None, None, None, None
    for i, fpath in iterable:
        snap = load_snapshot(fpath)
        fig, axes, cbars, norms = make_frame(snap, r_max_au=r_max_au,
                                             fig=fig, axes=axes,
                                             cbars=cbars, norms=norms)
        out = os.path.join(output_dir, f"frame_{i:04d}.png")
        fig.savefig(out, dpi=dpi, facecolor=BG, bbox_inches="tight")

    plt.close(fig)
    print(f"[OK] {len(files)} frames guardados en {output_dir}/")


def run_gif(datadir, gif_path, r_max_au=40, every=1, dpi=100, fps=6):
    """Genera un GIF animado del disco."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from PIL import Image
    except ImportError:
        raise ImportError("Pillow requerido para GIF: pip install Pillow")

    files = get_hdf5_list(datadir)
    files = files[::every]
    print(f"Generando GIF con {len(files)} frames → {gif_path}")

    frames_pil = []
    fig, axes, cbars, norms = None, None, None, None

    try:
        from tqdm import tqdm
        it = tqdm(files)
    except ImportError:
        it = files

    for fpath in it:
        snap = load_snapshot(fpath)
        fig, axes, cbars, norms = make_frame(snap, r_max_au=r_max_au,
                                             fig=fig, axes=axes,
                                             cbars=cbars, norms=norms)
        fig.canvas.draw()
        buf = fig.canvas.buffer_rgba()
        w, h = fig.canvas.get_width_height()
        img = np.frombuffer(buf, dtype=np.uint8).reshape(h, w, 4)
        frames_pil.append(Image.fromarray(img[:, :, :3]))

    plt.close(fig)

    dur = int(1000 / fps)
    frames_pil[0].save(
        gif_path,
        save_all=True,
        append_images=frames_pil[1:],
        optimize=False,
        duration=dur,
        loop=0,
    )
    print(f"[OK] GIF guardado: {gif_path}  ({len(frames_pil)} frames @ {fps} fps)")


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Visualizador 2D de disco protoplanetario (revolución)"
    )
    parser.add_argument(
        "--datadir", "-d",
        required=True,
        help="Directorio con los HDF5 del run (ej: ../data_1myr/single_jup_3.0au)"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Carpeta donde guardar los frames PNG"
    )
    parser.add_argument(
        "--gif", "-g",
        default=None,
        help="Ruta del GIF de salida (ej: disco.gif)"
    )
    parser.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Modo interactivo (flechas para navegar snapshots)"
    )
    parser.add_argument(
        "--rmax", type=float, default=40,
        help="Radio máximo del plot en AU (default: 40)"
    )
    parser.add_argument(
        "--every", type=int, default=1,
        help="Usar 1 de cada N snapshots (default: 1 = todos)"
    )
    parser.add_argument(
        "--dpi", type=int, default=120,
        help="DPI de los frames guardados (default: 120)"
    )
    parser.add_argument(
        "--fps", type=int, default=6,
        help="Frames por segundo del GIF (default: 6)"
    )
    parser.add_argument(
        "--snapshot", "-s", type=int, default=None,
        help="Mostrar solo este índice de snapshot (0-based) y guardarlo"
    )

    args = parser.parse_args()

    if args.snapshot is not None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        files = get_hdf5_list(args.datadir)
        snap  = load_snapshot(files[args.snapshot])
        fig, _, __ = make_frame(snap, r_max_au=args.rmax)
        name = f"snap_{args.snapshot:04d}.png"
        out  = args.output if args.output and args.output.endswith(".png") else name
        fig.savefig(out, dpi=args.dpi, facecolor=BG, bbox_inches="tight")
        print(f"[OK] Snapshot {args.snapshot} guardado -> {out}")
        plt.close()
    elif args.interactive:
        run_interactive(args.datadir, r_max_au=args.rmax)
    elif args.gif:
        run_gif(args.datadir, args.gif,
                r_max_au=args.rmax, every=args.every,
                dpi=args.dpi, fps=args.fps)
    elif args.output:
        run_save_frames(args.datadir, args.output,
                        r_max_au=args.rmax, every=args.every, dpi=args.dpi)
    else:
        # Por defecto: guarda el primer y último snapshot para preview
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        files = get_hdf5_list(args.datadir)
        for idx, label in [(0, "inicio"), (-1, "final")]:
            snap = load_snapshot(files[idx])
            fig, _, __, ___ = make_frame(snap, r_max_au=args.rmax)
            name = f"disco_{label}.png"
            fig.savefig(name, dpi=args.dpi, facecolor=BG, bbox_inches="tight")
            print(f"[OK] {name}")
            plt.close()


if __name__ == "__main__":
    main()
