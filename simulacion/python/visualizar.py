#!/usr/bin/env python3
"""Render static figures from MATLAB-exported THERMOPLATE results.

This program is deliberately a reader/plotter: it never assembles or solves the
physical model. See README.md for the accepted export schema and CLI examples.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCENARIOS = ("escenario_base", "convergencia", "refinamiento", "sensibilidad", "estres_numerico")
SOLVER_LABELS = {
    "gauss-seidel": "Gauss–Seidel",
    "conjugate gradient": "Gradiente Conjugado",
    "gradiente conjugado": "Gradiente Conjugado",
}
ALIASES = {
    "x": ("x", "x_m", "x_coord", "x_coordinate"),
    "y": ("y", "y_m", "y_coord", "y_coordinate"),
    "temperature": ("temperature", "temp", "t", "temperature_c", "t_c", "temperatura", "temperatura_c"),
    "solver": ("solver", "method", "metodo", "algorithm", "algoritmo"),
    "iteration": ("iteration", "iter", "iterations", "iteracion", "k"),
    "residual": ("relative_residual", "rel_residual", "residual_relative", "residuo_relativo", "residual", "residuo"),
    "time": ("time_s", "elapsed_s", "seconds", "tiempo_s", "runtime_s", "time"),
    "time_median": ("time_median_s", "median_time_s", "tiempo_mediano_s"),
    "time_q25": ("time_q25_s", "q25_time_s", "time_p25_s", "tiempo_q25_s"),
    "time_q75": ("time_q75_s", "q75_time_s", "time_p75_s", "tiempo_q75_s"),
    "time_iqr": ("time_iqr_s", "iqr_time_s", "tiempo_iqr_s"),
    "time_std": ("time_std_s", "std_time_s", "tiempo_desv_estandar_s"),
    "nx": ("nx", "n_x", "nodes_x", "nodos_x"),
    "ny": ("ny", "n_y", "nodes_y", "nodos_y"),
    "unknowns": ("unknowns", "n_unknowns", "n_interior", "incognitas"),
    "tmax": ("tmax", "tmax_c", "max_temperature", "maximum_temperature", "t_max", "temperatura_maxima"),
    "cpu_power": ("cpu_power_w", "cpu_power", "power_cpu_w", "potencia_cpu_w", "p_cpu_w"),
    "gpu_power": ("gpu_power_w", "gpu_power", "power_gpu_w", "potencia_gpu_w", "p_gpu_w"),
    "tolerance": ("tolerance", "tol", "tolerancia"),
    "status": ("status", "converged", "state", "estado"),
    "case": ("case", "scenario", "name", "caso", "escenario"),
}


class InputError(Exception):
    """Error de entrada que puede explicarse al usuario desde la CLI."""

    pass


def display_solver(value: Any) -> str:
    """Traduce etiquetas conocidas para figuras, sin cambiar datos exportados."""
    raw = str(value).strip()
    return SOLVER_LABELS.get(raw.casefold(), raw)


def display_stress_case(value: Any) -> str:
    """Acorta y traduce los dos casos de estrés definidos por MATLAB."""
    raw = str(value).strip()
    labels = {
        "extreme_positive_1000x": "Potencias × 1000",
        "nonphysical_negative_cpu": "CPU negativa",
    }
    return labels.get(raw.casefold(), raw.replace("_", " "))


def norm(value: Any) -> str:
    """Normaliza un nombre de columna para comparar alias sin distinguir formato."""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def find_col(df: pd.DataFrame, logical: str) -> str | None:
    """Devuelve la primera columna que coincide con el nombre lógico solicitado."""
    lookup = {norm(c): c for c in df.columns}
    for alias in ALIASES.get(logical, (logical,)):
        if norm(alias) in lookup:
            return lookup[norm(alias)]
    return None


def read_json(path: Path) -> Any:
    """Lee JSON UTF-8 y convierte fallos de archivo o sintaxis en InputError."""
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError(f"No se pudo leer JSON {path}: {exc}") from exc


def read_table(path: Path) -> pd.DataFrame:
    """Lee una tabla exportada por MATLAB; no interpreta ni recalcula el modelo."""
    suffix = path.suffix.lower()
    try:
        if suffix == ".csv":
            return pd.read_csv(path)
        if suffix in (".json", ".jsonl"):
            obj = read_json(path)
            if isinstance(obj, list):
                return pd.DataFrame(obj)
            if isinstance(obj, dict):
                # Accept a direct row array or the common {data: [...]} wrapper.
                rows = obj.get("data", obj.get("rows"))
                if isinstance(rows, list):
                    return pd.DataFrame(rows)
                # Column-oriented JSON arrays are also understood.
                arrays = {k: v for k, v in obj.items() if isinstance(v, list)}
                if arrays:
                    return pd.DataFrame(arrays)
            raise InputError(f"JSON tabular no reconocido en {path}; use lista de filas o columnas.")
        if suffix == ".mat":
            try:
                from scipy.io import loadmat  # optional, only for legacy MATLAB exports
            except ImportError as exc:
                raise InputError(f"{path.name} requiere SciPy para leerse; instale requirements.txt o exporte CSV/JSON.") from exc
            raw = loadmat(path, squeeze_me=True, struct_as_record=False)
            arrays = {k: np.asarray(v).reshape(-1) for k, v in raw.items()
                      if not k.startswith("__") and np.asarray(v).ndim <= 2 and np.asarray(v).size > 1}
            if not arrays:
                raise InputError(f"No se encontraron variables tabulares simples en {path}.")
            lengths = {len(v) for v in arrays.values()}
            if len(lengths) != 1:
                raise InputError(f"Las variables de {path} no tienen longitudes compatibles.")
            return pd.DataFrame(arrays)
    except (OSError, ValueError, NotImplementedError, pd.errors.ParserError) as exc:
        raise InputError(f"No se pudo leer tabla {path}: {exc}") from exc
    raise InputError(f"Formato no soportado: {path.suffix}. Use CSV, JSON o MAT.")


class Results:
    """Carga manifiestos y tablas de los cinco experimentos exportados."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.manifests: dict[str, dict[str, Any]] = {}
        self.datasets: dict[str, dict[str, pd.DataFrame]] = {}
        self.files: dict[str, list[Path]] = {}
        self.warnings: list[str] = []
        for scenario in SCENARIOS:
            folder = root / scenario
            if not folder.is_dir():
                self.warnings.append(f"Falta carpeta de escenario: {folder}")
                continue
            manifest_path = folder / "manifest.json"
            manifest = read_json(manifest_path) if manifest_path.is_file() else {}
            if not isinstance(manifest, dict):
                raise InputError(f"El manifiesto debe ser un objeto JSON: {manifest_path}")
            metadata_path = folder / "metadata.json"
            if metadata_path.is_file():
                metadata = read_json(metadata_path)
                if isinstance(metadata, dict):
                    existing = manifest.get("metadata", {})
                    manifest["metadata"] = {**metadata, **(existing if isinstance(existing, dict) else {})}
            self.manifests[scenario] = manifest
            tables: dict[str, pd.DataFrame] = {}
            named = manifest.get("datasets", {})
            if isinstance(named, dict):
                for key, rel in named.items():
                    if isinstance(rel, dict):
                        rel = rel.get("file")
                    if not isinstance(rel, str):
                        self.warnings.append(f"Dataset {scenario}.{key}: falta campo file en manifest.json")
                        continue
                    p = (folder / rel).resolve()
                    if not p.is_relative_to(folder.resolve()):
                        raise InputError(f"Ruta de dataset fuera del escenario en {manifest_path}: {rel}")
                    if p.is_file():
                        # MAT is an archival MATLAB artifact; the visualization
                        # contract uses its CSV/JSON siblings, not full workspaces.
                        if str(key).lower() == "metadata" or p.suffix.lower() == ".mat":
                            continue
                        tables[str(key)] = read_table(p)
                    else:
                        self.warnings.append(f"No existe dataset declarado {scenario}.{key}: {p}")
            self.files[scenario] = []
            # Discovery fallback: recognized filename stem becomes dataset name.
            for p in sorted(folder.rglob("*")):
                if p.is_file() and p.suffix.lower() in (".csv", ".json", ".mat") and p.name not in ("manifest.json", "metadata.json"):
                    self.files[scenario].append(p)
                    key = p.stem.lower()
                    if key not in tables:
                        try:
                            tables[key] = read_table(p)
                        except InputError as exc:
                            if p.suffix.lower() != ".mat":
                                self.warnings.append(str(exc))
            self.datasets[scenario] = tables

    def table(self, scenario: str, *names: str) -> pd.DataFrame | None:
        """Busca el primer dataset disponible entre los nombres indicados."""
        data = self.datasets.get(scenario, {})
        for name in names:
            if name in data:
                return data[name]
        return None


def finite_numeric(df: pd.DataFrame, col: str) -> pd.Series:
    """Convierte una columna a números y representa entradas inválidas como NaN."""
    return pd.to_numeric(df[col], errors="coerce")


def field_data(results: Results) -> tuple[pd.DataFrame, dict[str, Any]] | None:
    """Prepara el campo exportado en coordenadas x, y y temperatura en °C."""
    df = results.table("escenario_base", "field", "temperature_field", "campo_termico", "temperatures")
    if df is None:
        return None
    xcol, ycol, tcol = (find_col(df, k) for k in ("x", "y", "temperature"))
    if xcol and ycol and tcol:
        out = df[[xcol, ycol, tcol]].copy()
        out.columns = ["x", "y", "temperature"]
    else:
        # Matrix export: JSON object with x/y axes and a temperature matrix.
        manifest = results.manifests.get("escenario_base", {})
        grid = manifest.get("field_grid")
        if not isinstance(grid, dict):
            results.warnings.append("Campo base: se requieren columnas x,y,T o manifest.field_grid con x,y,T.")
            return None
        try:
            xx, yy, tt = np.asarray(grid["x"], float), np.asarray(grid["y"], float), np.asarray(grid["T"], float)
            if tt.shape != (len(yy), len(xx)):
                raise ValueError("T debe tener forma (len(y), len(x))")
            X, Y = np.meshgrid(xx, yy)
            out = pd.DataFrame({"x": X.ravel(), "y": Y.ravel(), "temperature": tt.ravel()})
        except (KeyError, ValueError, TypeError) as exc:
            results.warnings.append(f"field_grid inválido: {exc}")
            return None
    for c in ("x", "y", "temperature"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.dropna()
    manifest = results.manifests.get("escenario_base", {})
    meta = manifest.get("metadata", {}) if isinstance(manifest.get("metadata", {}), dict) else {}
    params = manifest.get("parameters", {}) if isinstance(manifest.get("parameters", {}), dict) else {}
    meta = {**params, **meta}
    return (out, meta) if not out.empty else None


def save(fig: plt.Figure, out: Path, stem: str) -> None:
    """Guarda una figura en PNG de alta resolución y PDF vectorial."""
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_geometry(meta: dict[str, Any], out: Path) -> bool:
    """Dibuja la placa y las fuentes usando únicamente metadatos exportados."""
    geom = meta.get("geometry", {})
    if not isinstance(geom, dict):
        geom = {}
    Lx = geom.get("Lx_m", geom.get("Lx", meta.get("Lx_m")))
    Ly = geom.get("Ly_m", geom.get("Ly", meta.get("Ly_m")))
    sources = geom.get("sources", meta.get("sources", []))
    if not Lx or not Ly or not isinstance(sources, list):
        return False
    fig, ax = plt.subplots(figsize=(6.2, 3.5))
    ax.add_patch(plt.Rectangle((0, 0), float(Lx), float(Ly), fill=False, lw=1.6, color="black"))
    palette = {"cpu": "#d95f02", "gpu": "#1b9e77"}
    for i, src in enumerate(sources):
        if not isinstance(src, dict):
            continue
        name = str(src.get("name", f"Source {i+1}"))
        x = src.get("x_m", src.get("x")); y = src.get("y_m", src.get("y"))
        w = src.get("width_m", src.get("width")); h = src.get("height_m", src.get("height"))
        if None in (x, y, w, h):
            continue
        color = palette.get(name.lower(), "#7570b3")
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, alpha=.30, edgecolor=color, lw=1.5))
        ax.text(float(x)+float(w)/2, float(y)+float(h)/2, name, ha="center", va="center", fontsize=9)
    ax.set(xlim=(0, float(Lx)), ylim=(0, float(Ly)), xlabel="x [m]", ylabel="y [m]", aspect="equal", title="Placa y fuentes térmicas\nTemperatura prescrita en el contorno")
    save(fig, out, "01_geometria")
    return True


def figure_field(field: pd.DataFrame, out: Path) -> None:
    """Genera dos vistas del campo de temperatura exportado por MATLAB."""
    pivot = field.pivot(index="y", columns="x", values="temperature").sort_index().sort_index(axis=1)
    x, y, z = pivot.columns.to_numpy(float), pivot.index.to_numpy(float), pivot.to_numpy(float)
    if not np.isfinite(z).all():
        raise InputError("El campo térmico exportado contiene nodos faltantes o no finitos.")
    X, Y = np.meshgrid(x, y)
    vmin, vmax = float(z.min()), float(z.max())
    fig, ax = plt.subplots(figsize=(6.3, 4.1), layout="constrained")
    im = ax.pcolormesh(X, Y, z, shading="auto", cmap="inferno", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, label="Temperatura [°C]")
    ax.set(xlabel="x [m]", ylabel="y [m]", title="Campo estacionario de temperatura")
    save(fig, out, "02_mapa_termico")

    # Se conservan los nodos exportados; la superficie no interpola ni resuelve el modelo.
    fig = plt.figure(figsize=(3.55, 3.25), layout="constrained")
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(
        X, Y, z, cmap="inferno", vmin=vmin, vmax=vmax,
        rcount=z.shape[0], ccount=z.shape[1], linewidth=0, antialiased=True,
    )
    ax.set_proj_type("ortho")
    ax.view_init(elev=29, azim=-132)
    ax.set_box_aspect((float(x[-1] - x[0]), float(y[-1] - y[0]), 0.58 * float(x[-1] - x[0])))
    ax.set(xlabel="x [m]", ylabel="y [m]", zlabel="Temperatura [°C]")
    ax.tick_params(labelsize=7)
    ax.xaxis.labelpad = 1
    ax.yaxis.labelpad = 1
    ax.zaxis.labelpad = 2
    fig.colorbar(surf, ax=ax, shrink=.72, pad=.08, aspect=18, label="Temperatura [°C]")
    save(fig, out, "03_superficie_3d")


def figure_history(results: Results, out: Path) -> tuple[bool, bool]:
    """Grafica historial de residuos y, si existe, estados espaciales iterativos."""
    hist = results.table("convergencia", "history", "convergence", "convergencia", "residuals")
    if hist is None:
        hist = results.table("escenario_base", "history", "convergence", "convergencia", "residuals")
    if hist is None:
        return False, False
    iter_col, solver_col, res_col = (find_col(hist, k) for k in ("iteration", "solver", "residual"))
    if not iter_col or not solver_col or not res_col:
        results.warnings.append("Convergencia: history necesita columnas iter, solver y relative_residual.")
        return False, False
    h = hist.copy()
    h["_iter"] = finite_numeric(h, iter_col); h["_res"] = finite_numeric(h, res_col)
    h = h.dropna(subset=["_iter", "_res"])
    if h.empty:
        results.warnings.append("Convergencia: history no contiene filas numéricas utilizables.")
        return False, False
    fig, ax = plt.subplots(figsize=(6.1, 3.8))
    for name, group in h.groupby(solver_col):
        group = group.sort_values("_iter")
        ax.plot(group["_iter"], group["_res"], label=display_solver(name))
    ax.set(xlabel="Iteración", ylabel="Residuo relativo", title="Comparación visual de Gauss–Seidel y Gradiente Conjugado")
    ax.grid(True, alpha=.25); ax.legend()
    save(fig, out, "05_comparacion_gs_cg")

    fig, ax = plt.subplots(figsize=(6.1, 3.8))
    for name, group in h.groupby(solver_col):
        group = group.sort_values("_iter")
        positive = group[group["_res"] > 0]
        if not positive.empty:
            ax.semilogy(positive["_iter"], positive["_res"], label=display_solver(name))
    ax.set(xlabel="Iteración", ylabel="Residuo relativo", title="Convergencia de los solucionadores")
    ax.grid(True, which="both", alpha=.25); ax.legend()
    save(fig, out, "06_residuo_iteraciones")

    # State snapshots are only shown if MATLAB exported temperature states.
    state = results.table("convergencia", "states", "state", "snapshots", "estados")
    if state is None:
        state = results.table("escenario_base", "states", "state", "snapshots", "estados")
    if state is None:
        results.warnings.append("No hay estados intermedios exportados: no se genera evolución espacial de aproximaciones.")
        return True, False
    sc = [find_col(state, k) for k in ("iteration", "solver", "x", "y", "temperature")]
    if not all(sc):
        results.warnings.append("Estados intermedios: se requieren iter,solver,x_m,y_m,T_C para evolución espacial.")
        return True, False
    state = state.rename(columns={sc[0]: "iteration", sc[1]: "solver", sc[2]: "x", sc[3]: "y", sc[4]: "temperature"})
    for col in ("iteration", "x", "y", "temperature"):
        state[col] = pd.to_numeric(state[col], errors="coerce")
    state = state.dropna(subset=["iteration", "x", "y", "temperature"])
    solver_values = state["solver"].dropna().astype(str)
    solver_lookup = {value.casefold(): value for value in solver_values.unique()}
    solvers = [solver_lookup.get(name.casefold()) for name in ("Gauss-Seidel", "Conjugate Gradient")]
    if any(name is None for name in solvers):
        results.warnings.append("Estados intermedios: se requieren Gauss-Seidel y Gradiente Conjugado.")
        return True, False
    initial_values = state.loc[state["iteration"] == 0, "temperature"]
    if initial_values.empty or not np.allclose(initial_values.to_numpy(), 25.0, rtol=0, atol=1e-10):
        results.warnings.append("Estados intermedios: el estado inicial MATLAB no es uniforme a 25 °C.")
        return True, False

    # Exigir tres estados por método evita que una exportación incompleta se
    # presente como una comparación entre etapas equivalentes.
    snapshots: dict[str, list[tuple[int, pd.DataFrame]]] = {}
    for solver_name in solvers:
        method = state[state["solver"].astype(str) == solver_name]
        iterations = sorted(int(value) for value in method["iteration"].unique())
        if len(iterations) != 3 or iterations[0] != 0 or iterations[1] <= 0 or iterations[2] <= iterations[1]:
            results.warnings.append(
                f"Estados intermedios: {display_solver(solver_name)} debe incluir k=0, un estado medio y el convergido."
            )
            return True, False
        panels: list[tuple[int, pd.DataFrame]] = []
        for iteration in iterations:
            snap = method[method["iteration"] == iteration]
            panel = snap.pivot(index="y", columns="x", values="temperature").sort_index().sort_index(axis=1)
            if panel.empty or not np.isfinite(panel.to_numpy()).all():
                results.warnings.append(f"Estado k={iteration} de {display_solver(solver_name)} incompleto o no finito.")
                return True, False
            panels.append((iteration, panel))
        snapshots[solver_name] = panels

    vmin, vmax = float(state["temperature"].min()), float(state["temperature"].max())
    fig, axs = plt.subplots(2, 3, figsize=(9.2, 5.6), squeeze=False,
                            sharex=True, sharey=True, layout="constrained")
    im = None
    for row, solver_name in enumerate(solvers):
        for col, (iteration, panel) in enumerate(snapshots[solver_name]):
            ax = axs[row, col]
            im = ax.pcolormesh(
                panel.columns, panel.index, panel.to_numpy(), shading="auto",
                cmap="inferno", vmin=vmin, vmax=vmax,
            )
            stage = ("Estado inicial" if col == 0 else
                     "Aproximación intermedia" if col == 1 else "Estado convergido")
            ax.set_title(f"{stage}\n$k={iteration}$")
            if row == 1:
                ax.set_xlabel("x [m]")
        axs[row, 0].set_ylabel(f"{display_solver(solver_name)}\ny [m]")
    if im is None:
        plt.close(fig); return True, False
    colorbar_ticks = [vmin, 40.0, 55.0, 70.0, vmax]
    fig.colorbar(
        im, ax=axs.ravel().tolist(), label="Temperatura [°C]", shrink=.88,
        pad=.025, ticks=colorbar_ticks, format="%.1f",
    )
    fig.suptitle("Estados de las aproximaciones algebraicas", fontsize=13)
    save(fig, out, "04_evolucion_aproximaciones")
    return True, True


def figure_mesh(results: Results, out: Path) -> bool:
    """Compara tiempos exportados por método y tamaño de malla."""
    df = results.table("refinamiento", "mesh", "refinement", "refinamiento")
    if df is None:
        return False
    time_col, solver_col = find_col(df, "time"), find_col(df, "solver")
    median_col = find_col(df, "time_median")
    nx_col, ny_col = find_col(df, "nx"), find_col(df, "ny")
    if not (time_col or median_col) or not (nx_col or find_col(df, "unknowns")):
        results.warnings.append("Refinamiento: tabla requiere nx (o unknowns) y time_s.")
        return False
    d = df.copy()
    d["_time"] = finite_numeric(df, median_col or time_col)
    q25_col, q75_col = find_col(df, "time_q25"), find_col(df, "time_q75")
    interval_kind: str | None = None
    if q25_col and q75_col:
        d["_q25"] = finite_numeric(df, q25_col)
        d["_q75"] = finite_numeric(df, q75_col)
        interval_kind = "iqr"
    else:
        # When summary quartiles are absent, derive the interval from MATLAB's
        # per-repetition export. The scientific solver is never rerun here.
        repeats = results.table("refinamiento", "timing", "timing_repeats")
        if repeats is not None:
            repeat_time = find_col(repeats, "time")
            repeat_solver = find_col(repeats, "solver")
            repeat_nx, repeat_ny = find_col(repeats, "nx"), find_col(repeats, "ny")
            repeat_unknowns = find_col(repeats, "unknowns")
            left_keys: list[str] = []
            right_keys: list[str] = []
            if solver_col and repeat_solver:
                left_keys.append(solver_col)
                right_keys.append(repeat_solver)
            if nx_col and ny_col and repeat_nx and repeat_ny:
                left_keys.extend([nx_col, ny_col])
                right_keys.extend([repeat_nx, repeat_ny])
            elif not nx_col and find_col(df, "unknowns") and repeat_unknowns:
                left_keys.append(find_col(df, "unknowns"))
                right_keys.append(repeat_unknowns)
            if repeat_time and left_keys and len(left_keys) == len(right_keys):
                samples = repeats[[*right_keys, repeat_time]].copy()
                samples["_sample_time"] = pd.to_numeric(samples[repeat_time], errors="coerce")
                quartiles = samples.dropna(subset=["_sample_time"]).groupby(right_keys)["_sample_time"].agg(
                    _q25=lambda values: values.quantile(.25),
                    _q75=lambda values: values.quantile(.75),
                ).reset_index()
                d = d.merge(quartiles, left_on=left_keys, right_on=right_keys, how="left")
                if d["_q25"].notna().any() and d["_q75"].notna().any():
                    interval_kind = "iqr"
        if interval_kind is None:
            iqr_col = find_col(df, "time_iqr")
            std_col = find_col(df, "time_std")
            if iqr_col:
                d["_spread"] = finite_numeric(df, iqr_col)
                interval_kind = "iqr_width"
            elif std_col:
                d["_spread"] = finite_numeric(df, std_col)
                interval_kind = "std"
    if nx_col and ny_col:
        d["_mesh"] = d[nx_col].astype(str) + "×" + d[ny_col].astype(str)
    else:
        d["_mesh"] = df[find_col(df, "unknowns")].astype(str) + " incógnitas"
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    group_col = solver_col if solver_col else None
    groups = d.groupby(group_col) if group_col else [("", d)]
    for name, g in groups:
        g = g.dropna(subset=["_time"]).copy()
        label = display_solver(name) if name else None
        if interval_kind == "iqr" and {"_q25", "_q75"}.issubset(g.columns):
            valid = g["_q25"].notna() & g["_q75"].notna() & (g["_time"] > 0)
            if valid.any():
                yerr = np.vstack([
                    (g.loc[valid, "_time"] - g.loc[valid, "_q25"]).clip(lower=0),
                    (g.loc[valid, "_q75"] - g.loc[valid, "_time"]).clip(lower=0),
                ])
                ax.errorbar(g.loc[valid, "_mesh"], g.loc[valid, "_time"], yerr=yerr,
                            marker="o", capsize=3, label=label)
            if (~valid).any():
                ax.plot(g.loc[~valid, "_mesh"], g.loc[~valid, "_time"], marker="o", label=label)
        elif interval_kind in ("iqr_width", "std") and "_spread" in g.columns:
            valid = g["_spread"].notna() & (g["_time"] > 0)
            if valid.any():
                half_width = g.loc[valid, "_spread"] / (2 if interval_kind == "iqr_width" else 1)
                lower = np.minimum(half_width, g.loc[valid, "_time"] * .99)
                yerr = np.vstack([lower, half_width])
                ax.errorbar(g.loc[valid, "_mesh"], g.loc[valid, "_time"], yerr=yerr,
                            marker="o", capsize=3, label=label)
            if (~valid).any():
                ax.plot(g.loc[~valid, "_mesh"], g.loc[~valid, "_time"], marker="o", label=label)
        else:
            ax.plot(g["_mesh"], g["_time"], marker="o", label=label)
    ax.set(xlabel="Malla [nodos x × nodos y]", ylabel="Tiempo computacional [s]", title="Costo por tamaño de malla")
    statistic_label = "Mediana del tiempo [s]" if median_col or interval_kind else "Tiempo computacional [s]"
    interval_labels = {
        "iqr": "; barras: intervalo intercuartílico Q1–Q3",
        "iqr_width": "; barras: semiancho del IQR, centrado en la mediana",
        "std": "; barras: ± desviación estándar",
    }
    interval_label = interval_labels.get(interval_kind, "")
    ax.set_ylabel(statistic_label)
    ax.set_title(f"Tiempo por tamaño de malla{interval_label}")
    ax.set_yscale("log")
    ax.tick_params(axis="x", rotation=25); ax.grid(True, alpha=.25)
    if solver_col: ax.legend(title="Solucionador")
    save(fig, out, "07_costo_refinamiento")
    return True


def figure_sensitivity(results: Results, out: Path) -> bool:
    """Grafica Tmax exportada en función de la potencia variada."""
    df = results.table("sensibilidad", "sensitivity", "sensibilidad")
    if df is None: return False
    cpu, gpu, tmax = (find_col(df, k) for k in ("cpu_power", "gpu_power", "tmax"))
    if not tmax or not (cpu or gpu):
        results.warnings.append("Sensibilidad: se necesitan potencia CPU o GPU y temperatura máxima.")
        return False
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    varied = find_col(df, "varied_source")
    if varied:
        source_values = df[varied].astype(str).str.strip().str.lower()
        groups = ((cpu, "cpu", "Variación CPU"), (gpu, "gpu", "Variación GPU"))
        for col, source_name, label in groups:
            if col:
                subset = df[source_values.str.contains(source_name, regex=False)]
                xx = pd.to_numeric(subset[col], errors="coerce")
                yy = pd.to_numeric(subset[tmax], errors="coerce")
                good = xx.notna() & yy.notna()
                if good.any(): ax.plot(xx[good], yy[good], marker="o", label=label)
    else:
        for col, label in ((cpu, "Potencia CPU [W]"), (gpu, "Potencia GPU [W]")):
            if col:
                xx, yy = pd.to_numeric(df[col], errors="coerce"), pd.to_numeric(df[tmax], errors="coerce")
                good = xx.notna() & yy.notna()
                if good.any(): ax.plot(xx[good], yy[good], marker="o", label=label)
    ax.set(xlabel="Potencia de fuente [W]", ylabel="Temperatura máxima [°C]", title="Sensibilidad de Tmax a las fuentes")
    ax.grid(True, alpha=.25); ax.legend()
    save(fig, out, "08_sensibilidad_termica")
    return True


def figure_stress(results: Results, out: Path) -> bool:
    """Resume residuos frente a tolerancia y tiempos de casos de estrés."""
    df = results.table("estres_numerico", "stress", "stress_test", "estres", "estres_numerico")
    if df is None: return False
    tol, res, solver, case = (find_col(df, k) for k in ("tolerance", "residual", "solver", "case"))
    if not tol:
        results.warnings.append("Estrés numérico: tabla requiere columna tol/tolerance.")
        return False
    d = df.copy(); d["_tol"] = pd.to_numeric(df[tol], errors="coerce")
    d["_label"] = df[case].astype(str) if case else [f"Caso {i+1}" for i in range(len(df))]
    fig, (ax, ax_time) = plt.subplots(1, 2, figsize=(10.0, 4.2), layout="constrained")
    plotted = False
    tolerance_cases = d
    if case:
        tolerance_cases = d[d[case].astype(str).str.lower().str.contains("tolerance|tolerancia")]
    grouped = tolerance_cases.groupby(solver) if solver else [("", tolerance_cases)]
    for name, g in grouped:
        g = g.dropna(subset=["_tol"]).copy()
        if res:
            g["_res"] = pd.to_numeric(g[res], errors="coerce")
            good = g["_res"].notna() & (g["_res"] > 0) & (g["_tol"] > 0)
            if good.any():
                g = g.loc[good].sort_values("_tol")
                ax.loglog(g["_tol"], g["_res"], marker="o", label=display_solver(name) if name else "Residuo")
                plotted = True
    if plotted:
        ax.set(xlabel="Tolerancia solicitada", ylabel="Residuo relativo alcanzado", title="Residuo relativo frente a tolerancia")
        ax.legend(title="Solucionador" if solver else None)
    else:
        status = find_col(df, "status")
        if not status:
            results.warnings.append("Estrés: sin residuo válido ni estado/status; no se genera figura.")
            plt.close(fig); return False
        counts = df.groupby([case, status]).size().unstack(fill_value=0) if case else df.groupby(status).size().to_frame("count")
        counts.plot(kind="bar", ax=ax_time)
        ax_time.set(ylabel="Número de ejecuciones", title="Estados de ejecución")
    time_col = find_col(df, "time")
    if time_col and case:
        other = df[~df[case].astype(str).str.lower().str.contains("tolerance|tolerancia")].copy()
        if not other.empty:
            other["_time"] = pd.to_numeric(other[time_col], errors="coerce")
            if solver:
                pivot = other.pivot_table(index=case, columns=solver, values="_time", aggfunc="mean")
            else:
                pivot = other.groupby(case)["_time"].mean().to_frame("Tiempo")
            cases = list(pivot.index)
            positions = np.arange(len(cases), dtype=float)
            solver_names = list(pivot.columns)
            width = 0.22
            offsets = (np.arange(len(solver_names)) - (len(solver_names) - 1) / 2) * width
            for offset, name in zip(offsets, solver_names):
                values = pd.to_numeric(pivot[name], errors="coerce").to_numpy(dtype=float)
                valid = np.isfinite(values) & (values > 0)
                x = positions + offset
                ax_time.plot(x[valid], values[valid], marker="o", linestyle="none",
                             label=display_solver(name))
                for xx, yy in zip(x[valid], values[valid]):
                    ax_time.annotate(f"{yy:.2g}", (xx, yy), xytext=(0, 5),
                                     textcoords="offset points", ha="center", fontsize=8)
            ax_time.set_xticks(positions, [display_stress_case(value) for value in cases])
            ax_time.set_yscale("log")
            ax_time.set(ylabel="Mediana del tiempo computacional [s]",
                        title="Tiempo en casos extremos/no físicos")
            if solver_names:
                ax_time.legend(title="Solucionador")
    ax.grid(True, which="both", alpha=.25); ax_time.grid(True, axis="y", alpha=.25)
    save(fig, out, "09_estres_numerico")
    return True


def main(argv: list[str] | None = None) -> int:
    """Punto de entrada de CLI; en modo estricto exige todos los escenarios."""
    script_dir = Path(__file__).resolve().parent
    sim_dir = script_dir.parent
    parser = argparse.ArgumentParser(description="Genera figuras estáticas a partir de exportaciones MATLAB; no resuelve el modelo físico.")
    parser.add_argument("--results", type=Path, default=sim_dir / "resultados", help="Carpeta de exportaciones MATLAB")
    parser.add_argument("--output", type=Path, default=script_dir / "figuras", help="Carpeta de salida PNG/PDF")
    parser.add_argument("--strict", action="store_true", help="Devuelve error si falta alguno de los cinco escenarios esperados")
    args = parser.parse_args(argv)
    if not args.results.is_dir():
        print(f"ERROR: no existe la carpeta de resultados: {args.results}\nEjecute primero los experimentos MATLAB y exporte los CSV/JSON.", file=sys.stderr)
        return 2
    try:
        results = Results(args.results.resolve())
    except InputError as exc:
        print(f"ERROR: {exc}", file=sys.stderr); return 2
    if args.strict and results.warnings:
        print("ERROR: faltan entradas para modo estricto:", file=sys.stderr)
        for warning in results.warnings: print(f"  - {warning}", file=sys.stderr)
        return 2
    generated: list[str] = []
    base = field_data(results)
    if base:
        field, meta = base
        if figure_geometry(meta, args.output): generated.append("01_geometria")
        figure_field(field, args.output); generated.extend(["02_mapa_termico", "03_superficie_3d"])
    else:
        results.warnings.append("Sin campo térmico base: se omiten mapa y superficie 3D.")
    has_history, has_states = figure_history(results, args.output)
    if has_history: generated.extend(["05_comparacion_gs_cg", "06_residuo_iteraciones"])
    if has_states: generated.append("04_evolucion_aproximaciones")
    if figure_mesh(results, args.output): generated.append("07_costo_refinamiento")
    if figure_sensitivity(results, args.output): generated.append("08_sensibilidad_termica")
    if figure_stress(results, args.output): generated.append("09_estres_numerico")
    if not generated:
        print(f"ERROR: no se generaron figuras. Revise el esquema de entrada en {script_dir / 'README.md'}.", file=sys.stderr)
        for warning in results.warnings: print(f"  - {warning}", file=sys.stderr)
        return 2
    print(f"Figuras generadas en {args.output.resolve()}:")
    for name in generated: print(f"  - {name}.png y {name}.pdf")
    if results.warnings:
        print("\nAvisos de entradas ausentes/incompletas:")
        for warning in results.warnings: print(f"  - {warning}")
        if args.strict:
            print("ERROR: el modo estricto requiere exportaciones completas.", file=sys.stderr)
            return 2
    print("Nota: todas las figuras representan datos exportados; Python no recalcula la física.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
