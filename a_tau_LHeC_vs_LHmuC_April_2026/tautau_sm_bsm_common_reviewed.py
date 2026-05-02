from __future__ import annotations

import gzip
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

PDG_TAU_MINUS = 15
PDG_TAU_PLUS = -15
HEADER_XSEC_RE = re.compile(r"Integrated weight \(pb\)\s*:\s*([+-]?\d*\.?\d+(?:[Ee][+-]?\d+)?)")
RUN_NUMBER_RE = re.compile(r"run_(\d+)")


@dataclass
class RunInfo:
    path: Path
    events: int
    cross_section_pb: float
    invariant_masses: np.ndarray
    pair_rapidities: np.ndarray


@dataclass
class SampleData:
    label: str
    file_paths: list[Path]
    runs: list[RunInfo]
    events: int
    cross_section_pb: float
    invariant_masses: np.ndarray
    pair_rapidities: np.ndarray


@dataclass
class RatioSeries:
    centers: np.ndarray
    edges: np.ndarray
    ratio: np.ndarray
    ratio_err: np.ndarray
    valid: np.ndarray
    num_counts: np.ndarray | None = None
    den_counts: np.ndarray | None = None


@dataclass
class FitResult:
    slope: float
    slope_err: float
    intercept: float
    intercept_err: float
    chi2: float
    dof: int
    chi2_dof: float


def _natural_run_key(path: Path) -> tuple[int, str]:
    m = RUN_NUMBER_RE.search(path.name)
    return (int(m.group(1)) if m else 10**9, path.name)


def open_lhe_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def find_run_lhe_files(events_dir: Path, max_runs: int | None = None) -> list[Path]:
    if not events_dir.exists():
        raise FileNotFoundError(f"Events directory does not exist: {events_dir}")

    run_dirs = sorted(
        [p for p in events_dir.iterdir() if p.is_dir() and p.name.startswith("run_")],
        key=_natural_run_key,
    )
    if max_runs is not None:
        run_dirs = run_dirs[:max_runs]
    if not run_dirs:
        raise FileNotFoundError(f"No run_* directories found in: {events_dir}")

    files: list[Path] = []
    for run_dir in run_dirs:
        candidates = [
            run_dir / "unweighted_events.lhe.gz",
            run_dir / "unweighted_events.lhe",
            run_dir / "events.lhe.gz",
            run_dir / "events.lhe",
        ]
        found = next((c for c in candidates if c.exists()), None)
        if found is None:
            raise FileNotFoundError(f"No LHE/LHE.GZ file found in: {run_dir}")
        files.append(found)
    return files


def read_lhe_cross_section_pb(path: Path) -> float:
    with open_lhe_text(path) as handle:
        for line in handle:
            match = HEADER_XSEC_RE.search(line)
            if match:
                return float(match.group(1))
            if "<event>" in line:
                break
    raise ValueError(f"Could not find 'Integrated weight (pb)' in header of {path}")


def parse_lhe_kinematics(path: Path) -> tuple[np.ndarray, np.ndarray, int]:
    masses: list[float] = []
    rapidities: list[float] = []
    event_count = 0

    with open_lhe_text(path) as handle:
        in_event = False
        first_event_line = False
        tau_plus = None
        tau_minus = None

        for raw_line in handle:
            line = raw_line.strip()
            if line == "<event>":
                in_event = True
                first_event_line = True
                tau_plus = None
                tau_minus = None
                continue

            if line == "</event>":
                in_event = False
                first_event_line = False
                event_count += 1

                if tau_plus is not None and tau_minus is not None:
                    px_pair = tau_plus[0] + tau_minus[0]
                    py_pair = tau_plus[1] + tau_minus[1]
                    pz_pair = tau_plus[2] + tau_minus[2]
                    e_pair = tau_plus[3] + tau_minus[3]

                    m2 = e_pair * e_pair - (px_pair * px_pair + py_pair * py_pair + pz_pair * pz_pair)
                    masses.append(math.sqrt(max(0.0, m2)))

                    numer = e_pair + pz_pair
                    denom = e_pair - pz_pair
                    if numer > 0.0 and denom > 0.0:
                        rapidities.append(0.5 * math.log(numer / denom))
                continue

            if not in_event:
                continue

            if first_event_line:
                first_event_line = False
                continue

            if not line or line[0] not in "-0123456789":
                continue

            parts = line.split()
            if len(parts) < 10:
                continue

            try:
                pdg_id = int(parts[0])
                status = int(parts[1])
                if status != 1:
                    continue
                px = float(parts[6])
                py = float(parts[7])
                pz = float(parts[8])
                energy = float(parts[9])
            except ValueError:
                continue

            if pdg_id == PDG_TAU_PLUS:
                tau_plus = (px, py, pz, energy)
            elif pdg_id == PDG_TAU_MINUS:
                tau_minus = (px, py, pz, energy)

    return np.asarray(masses, dtype=float), np.asarray(rapidities, dtype=float), event_count


def load_sample(events_dir: Path, label: str, max_runs: int | None = None) -> SampleData:
    file_paths = find_run_lhe_files(events_dir, max_runs=max_runs)

    runs: list[RunInfo] = []
    for path in file_paths:
        xsec_pb = read_lhe_cross_section_pb(path)
        masses, rapidities, n_events = parse_lhe_kinematics(path)
        if n_events == 0:
            raise ValueError(f"No events were parsed from {path}")
        runs.append(RunInfo(path, n_events, xsec_pb, masses, rapidities))

    event_counts = np.asarray([r.events for r in runs], dtype=float)
    xsecs = np.asarray([r.cross_section_pb for r in runs], dtype=float)
    total_events = int(np.sum(event_counts))
    weighted_cross_section_pb = float(np.average(xsecs, weights=event_counts))

    return SampleData(
        label=label,
        file_paths=file_paths,
        runs=runs,
        events=total_events,
        cross_section_pb=weighted_cross_section_pb,
        invariant_masses=np.concatenate([r.invariant_masses for r in runs]),
        pair_rapidities=np.concatenate([r.pair_rapidities for r in runs]),
    )


def effective_luminosity_fb(sample: SampleData) -> float:
    return sample.events / (sample.cross_section_pb * 1000.0)


def average_effective_luminosity_fb(sm: SampleData, bsm: SampleData) -> float:
    return 0.5 * (effective_luminosity_fb(sm) + effective_luminosity_fb(bsm))


def histogram_ds_dx(
    values: np.ndarray,
    cross_section_pb: float,
    n_events: int,
    bins: int | Sequence[float],
    value_range: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    counts, edges = np.histogram(values, bins=bins, range=value_range)
    widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    if n_events <= 0:
        raise ValueError("Number of events must be positive.")
    scale = cross_section_pb / n_events
    ds_dx = counts * scale / widths
    ds_dx_err = np.sqrt(counts) * scale / widths
    return centers, edges, ds_dx, ds_dx_err, counts


def build_ratio_series_from_weighted_spectra(
    num: np.ndarray,
    num_err: np.ndarray,
    den: np.ndarray,
    den_err: np.ndarray,
    centers: np.ndarray,
    edges: np.ndarray,
    num_counts: np.ndarray | None = None,
    den_counts: np.ndarray | None = None,
    min_den_counts: int = 1,
    min_num_counts: int = 1,
) -> RatioSeries:
    valid = (num > 0.0) & (den > 0.0) & np.isfinite(num_err) & np.isfinite(den_err)
    if num_counts is not None:
        valid &= num_counts >= min_num_counts
    if den_counts is not None:
        valid &= den_counts >= min_den_counts

    ratio = np.full_like(num, np.nan, dtype=float)
    ratio_err = np.full_like(num, np.nan, dtype=float)
    ratio[valid] = num[valid] / den[valid]
    ratio_err[valid] = ratio[valid] * np.sqrt(
        (num_err[valid] / num[valid]) ** 2 + (den_err[valid] / den[valid]) ** 2
    )
    return RatioSeries(centers=centers, edges=edges, ratio=ratio, ratio_err=ratio_err, valid=valid, num_counts=num_counts, den_counts=den_counts)


def linear_model(x: np.ndarray, a: float, b: float) -> np.ndarray:
    return a * x + b


def fit_ratio_direct(series: RatioSeries, fit_range: tuple[float, float] | None = None) -> FitResult:
    fit_mask = series.valid & np.isfinite(series.ratio) & np.isfinite(series.ratio_err) & (series.ratio_err > 0.0)
    if fit_range is not None:
        lo, hi = fit_range
        fit_mask &= (series.centers >= lo) & (series.centers <= hi)
    if np.count_nonzero(fit_mask) < 3:
        raise ValueError("Not enough valid bins to perform a linear fit.")

    xfit = series.centers[fit_mask]
    yfit = series.ratio[fit_mask]
    yerr = series.ratio_err[fit_mask]

    popt, pcov = curve_fit(
        linear_model,
        xfit,
        yfit,
        sigma=yerr,
        absolute_sigma=True,
        p0=(0.0, 1.0),
        maxfev=10000,
    )
    slope, intercept = popt
    slope_err = float(np.sqrt(pcov[0, 0]))
    intercept_err = float(np.sqrt(pcov[1, 1]))

    residuals = (yfit - linear_model(xfit, *popt)) / yerr
    chi2 = float(np.sum(residuals**2))
    dof = int(len(xfit) - 2)
    chi2_dof = chi2 / dof if dof > 0 else float("nan")

    return FitResult(float(slope), slope_err, float(intercept), intercept_err, chi2, dof, float(chi2_dof))


def configure_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 300,
        "axes.linewidth": 1.4,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 8,
        "ytick.major.size": 8,
        "xtick.minor.size": 4,
        "ytick.minor.size": 4,
        "legend.frameon": False,
        "font.size": 15,
        "axes.labelsize": 18,
        "axes.titlesize": 22,
        "legend.fontsize": 15,
    })


def default_custom_ratio_edges(mass_range: tuple[float, float]) -> np.ndarray:
    lo, hi = mass_range
    if (lo, hi) != (10.0, 500.0):
        return np.linspace(lo, hi, 11)
    return np.array([10, 35, 60, 85, 110, 140, 175, 215, 260, 310, 365, 430, 500], dtype=float)


def _format_lumi_title(collider_label: str, delta_a_tau: float, lumi_fb: float, scientific: bool) -> str:
    lumi_text = f"{lumi_fb:.2e}" if scientific else f"{lumi_fb:.1f}"
    return rf"{collider_label} ($\delta a_\tau = {delta_a_tau:.3f}, \mathcal{{L}} = {lumi_text}$ fb$^{{-1}}$)"


def save_run_table(outdir: Path, sm: SampleData, bsm: SampleData) -> None:
    lines = ["sample,run,path,events,cross_section_pb"]
    for sample in (sm, bsm):
        for i, run in enumerate(sample.runs, 1):
            lines.append(f"{sample.label},{i},{run.path},{run.events},{run.cross_section_pb:.8f}")
    (outdir / "run_table.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_summary(outdir: Path, collider_label: str, sm: SampleData, bsm: SampleData, delta_a_tau: float, lumi_fb: float, physical_fit: FitResult, uniform_fit: FitResult, custom_fit: FitResult) -> None:
    text = f"""Collider: {collider_label}
Delta a_tau (BSM label): {delta_a_tau}

SM sample:
  runs read             : {len(sm.file_paths)}
  total events          : {sm.events}
  cross section [pb]    : {sm.cross_section_pb:.8f}
  effective lumi [fb^-1]: {effective_luminosity_fb(sm):.8f}

BSM sample:
  runs read             : {len(bsm.file_paths)}
  total events          : {bsm.events}
  cross section [pb]    : {bsm.cross_section_pb:.8f}
  effective lumi [fb^-1]: {effective_luminosity_fb(bsm):.8f}

Average effective luminosity [fb^-1]: {lumi_fb:.8f}

Fit to weighted physical ratio, fine mass bins:
  slope                 : {physical_fit.slope:.8f} +/- {physical_fit.slope_err:.8f}
  intercept             : {physical_fit.intercept:.8f} +/- {physical_fit.intercept_err:.8f}
  chi2/dof              : {physical_fit.chi2_dof:.6f}

Fit to weighted ratio, uniform presentation bins:
  slope                 : {uniform_fit.slope:.8f} +/- {uniform_fit.slope_err:.8f}
  intercept             : {uniform_fit.intercept:.8f} +/- {uniform_fit.intercept_err:.8f}
  chi2/dof              : {uniform_fit.chi2_dof:.6f}

Fit to weighted ratio, custom presentation bins:
  slope                 : {custom_fit.slope:.8f} +/- {custom_fit.slope_err:.8f}
  intercept             : {custom_fit.intercept:.8f} +/- {custom_fit.intercept_err:.8f}
  chi2/dof              : {custom_fit.chi2_dof:.6f}
"""
    (outdir / "analysis_summary.txt").write_text(text, encoding="utf-8")


def plot_mass_spectrum(outdir: Path, collider_label: str, delta_a_tau: float, sm: SampleData, bsm: SampleData, bins: int, mass_range: tuple[float, float]):
    centers, edges, sm_ds_dM, sm_err, sm_counts = histogram_ds_dx(sm.invariant_masses, sm.cross_section_pb, sm.events, bins=bins, value_range=mass_range)
    _, _, bsm_ds_dM, bsm_err, bsm_counts = histogram_ds_dx(bsm.invariant_masses, bsm.cross_section_pb, bsm.events, bins=bins, value_range=mass_range)
    fig, ax = plt.subplots(figsize=(8.2, 9.0))
    ax.step(edges[:-1], sm_ds_dM, where="post", linewidth=2.2, color="salmon", label=r"$\tau^+\tau^-\;({\rm SM})$")
    ax.step(edges[:-1], bsm_ds_dM, where="post", linewidth=2.2, linestyle="--", color="royalblue", label=rf"$\tau^+\tau^-\;(\delta a_\tau = {delta_a_tau:.3f})$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*mass_range)
    positive = np.concatenate([sm_ds_dM[sm_ds_dM > 0], bsm_ds_dM[bsm_ds_dM > 0]])
    ax.set_ylim(1e-5, 1.5 * np.max(positive) if positive.size else 1.0)
    ax.set_xlabel(r"$M_{\tau^+\tau^-}\;[\mathrm{GeV}]$")
    ax.set_ylabel(r"$d\sigma/dM_{\tau^+\tau^-}\;[\mathrm{pb/GeV}]$")
    ax.set_title(collider_label); ax.grid(True, which="both", linestyle="--", alpha=0.45); ax.legend(loc="upper right")
    fig.tight_layout(); fig.savefig(outdir / "Invariant_mass_tau_pair_SM_vs_BSM_delta_aTau.png"); plt.close(fig)
    return centers, edges, sm_ds_dM, sm_err, sm_counts, bsm_ds_dM, bsm_err, bsm_counts


def plot_rapidity_spectrum(outdir: Path, collider_label: str, delta_a_tau: float, sm: SampleData, bsm: SampleData, bins: int, rapidity_range: tuple[float, float]) -> None:
    _, edges, sm_ds_dY, _, _ = histogram_ds_dx(sm.pair_rapidities, sm.cross_section_pb, sm.events, bins=bins, value_range=rapidity_range)
    _, _, bsm_ds_dY, _, _ = histogram_ds_dx(bsm.pair_rapidities, bsm.cross_section_pb, bsm.events, bins=bins, value_range=rapidity_range)
    fig, ax = plt.subplots(figsize=(8.2, 9.0))
    ax.step(edges[:-1], sm_ds_dY, where="post", linewidth=2.2, color="salmon", label=r"$\tau^+\tau^-\;({\rm SM})$")
    ax.step(edges[:-1], bsm_ds_dY, where="post", linewidth=2.2, linestyle="--", color="royalblue", label=rf"$\tau^+\tau^-\;(\delta a_\tau = {delta_a_tau:.3f})$")
    ax.set_xlim(*rapidity_range)
    positive = np.concatenate([sm_ds_dY[sm_ds_dY > 0], bsm_ds_dY[bsm_ds_dY > 0]])
    if positive.size: ax.set_ylim(0, 1.25 * np.max(positive))
    ax.set_xlabel(r"$Y_{\tau^+\tau^-}$"); ax.set_ylabel(r"$d\sigma/dY_{\tau^+\tau^-}\;[\mathrm{pb}]$")
    ax.set_title(collider_label); ax.grid(True, linestyle="--", alpha=0.45); ax.legend(loc="best")
    fig.tight_layout(); fig.savefig(outdir / "Rapidity_tau_pair_SM_vs_BSM_delta_aTau.png"); plt.close(fig)


def plot_ratio_fit(outdir: Path, collider_label: str, delta_a_tau: float, lumi_fb: float, series: RatioSeries, fit: FitResult, filename: str, ymax: float = 3.0, scientific_lumi: bool = True) -> None:
    fig, ax = plt.subplots(figsize=(14.0, 5.4))
    ax.errorbar(series.centers[series.valid], series.ratio[series.valid], yerr=series.ratio_err[series.valid], fmt="o", color="red", markersize=8, label="Stat. uncertainty")
    fit_x = np.linspace(series.edges[0], series.edges[-1], 400)
    ax.plot(fit_x, linear_model(fit_x, fit.slope, fit.intercept), linestyle="--", color="magenta", linewidth=2.4,
            label=(rf"Fit: $y = ({fit.slope:.4f} \pm {fit.slope_err:.4f})x + "
                   rf"({fit.intercept:.4f} \pm {fit.intercept_err:.4f})$, ")) #                    rf"$\chi^2/\mathrm{{dof}} = {fit.chi2_dof:.2f}$"

    ax.axhline(1.0, color="green", linestyle="--", linewidth=2.0, label="y=1")
    ax.set_ylim(0, ymax); ax.set_xlabel(r"$M_{\tau^+\tau^-}\;[\mathrm{GeV}]$"); ax.set_ylabel(r"Ratio $(\mathrm{BSM}/\mathrm{SM})$")
    ax.set_title(_format_lumi_title(collider_label, delta_a_tau, lumi_fb, scientific=scientific_lumi))
    ax.grid(True, linestyle="--", alpha=0.45); ax.legend(loc="best")
    fig.tight_layout(); fig.savefig(outdir / filename); plt.close(fig)


def make_weighted_ratio_for_bins(sm: SampleData, bsm: SampleData, bins: int | Sequence[float], mass_range: tuple[float, float] | None, min_bin_events: int) -> RatioSeries:
    centers, edges, sm_ds, sm_err, sm_counts = histogram_ds_dx(sm.invariant_masses, sm.cross_section_pb, sm.events, bins=bins, value_range=mass_range)
    _, _, bsm_ds, bsm_err, bsm_counts = histogram_ds_dx(bsm.invariant_masses, bsm.cross_section_pb, bsm.events, bins=bins, value_range=mass_range)
    return build_ratio_series_from_weighted_spectra(bsm_ds, bsm_err, sm_ds, sm_err, centers, edges, bsm_counts, sm_counts, min_den_counts=min_bin_events, min_num_counts=min_bin_events)


def run_analysis(
    collider_label: str,
    sm_events_dir: Path,
    bsm_events_dir: Path,
    outdir: Path,
    delta_a_tau: float = 0.005,
    mass_plot_bins: int = 100,

    rapidity_plot_bins: int = 20,

    ratio_bins: int = 10,

    mass_range: tuple[float, float] = (10.0, 500.0),

    rapidity_range: tuple[float, float] = (-10.0, 10.0),

    custom_ratio_edges: Sequence[float] | None = None,
    display_lumi_fb: float | None = None,
    max_runs: int | None = None,
    min_bin_events: int = 10,
) -> None:
    configure_style(); outdir.mkdir(parents=True, exist_ok=True)
    sm = load_sample(sm_events_dir, "SM", max_runs=max_runs)
    bsm = load_sample(bsm_events_dir, "BSM", max_runs=max_runs)
    save_run_table(outdir, sm, bsm)

    effective_lumi_fb = average_effective_luminosity_fb(sm, bsm)
    shown_lumi_fb = effective_lumi_fb if display_lumi_fb is None else float(display_lumi_fb)

    centers, edges, sm_ds_dM, sm_err, sm_counts, bsm_ds_dM, bsm_err, bsm_counts = plot_mass_spectrum(outdir, collider_label, delta_a_tau, sm, bsm, mass_plot_bins, mass_range)
    plot_rapidity_spectrum(outdir, collider_label, delta_a_tau, sm, bsm, rapidity_plot_bins, rapidity_range)

    physical_series = build_ratio_series_from_weighted_spectra(bsm_ds_dM, bsm_err, sm_ds_dM, sm_err, centers, edges, bsm_counts, sm_counts, min_den_counts=min_bin_events, min_num_counts=min_bin_events)
    physical_fit = fit_ratio_direct(physical_series)
    plot_ratio_fit(outdir, collider_label, delta_a_tau, shown_lumi_fb, physical_series, physical_fit, "Ratio_BSM_over_SM_MtauTau_weighted_fine_fit.png", ymax=3.0, scientific_lumi=True)

    uniform_weighted = make_weighted_ratio_for_bins(sm, bsm, bins=ratio_bins, mass_range=mass_range, min_bin_events=min_bin_events)
    uniform_fit = fit_ratio_direct(uniform_weighted)
    plot_ratio_fit(outdir, collider_label, delta_a_tau, shown_lumi_fb, uniform_weighted, uniform_fit, "Ratio_Obs_Exp_with_Luminosity_Fit_uniform_WEIGHTED.png", ymax=3.0, scientific_lumi=True)

    custom_edges = np.asarray(custom_ratio_edges, dtype=float) if custom_ratio_edges is not None else default_custom_ratio_edges(mass_range)
    custom_weighted = make_weighted_ratio_for_bins(sm, bsm, bins=custom_edges, mass_range=None, min_bin_events=min_bin_events)
    custom_fit = fit_ratio_direct(custom_weighted)
    plot_ratio_fit(outdir, collider_label, delta_a_tau, shown_lumi_fb, custom_weighted, custom_fit, "Ratio_Obs_Exp_with_Luminosity_Fit_custom_bins_WEIGHTED.png", ymax=2.0, scientific_lumi=False)

    save_summary(outdir, collider_label, sm, bsm, delta_a_tau, effective_lumi_fb, physical_fit, uniform_fit, custom_fit)

    print("=" * 72)
    print(f"{collider_label} analysis completed")
    print(f"SM  : runs = {len(sm.file_paths)}, events = {sm.events}, sigma = {sm.cross_section_pb:.8f} pb, L_eff = {effective_luminosity_fb(sm):.6f} fb^-1")
    print(f"BSM : runs = {len(bsm.file_paths)}, events = {bsm.events}, sigma = {bsm.cross_section_pb:.8f} pb, L_eff = {effective_luminosity_fb(bsm):.6f} fb^-1")
    print(f"Average effective luminosity used for titles = {shown_lumi_fb:.6f} fb^-1")
    print(f"Weighted uniform fit: slope = {uniform_fit.slope:.6f} +/- {uniform_fit.slope_err:.6f}, intercept = {uniform_fit.intercept:.6f} +/- {uniform_fit.intercept_err:.6f}, chi2/dof = {uniform_fit.chi2_dof:.3f}")
    print(f"Weighted custom fit : slope = {custom_fit.slope:.6f} +/- {custom_fit.slope_err:.6f}, intercept = {custom_fit.intercept:.6f} +/- {custom_fit.intercept_err:.6f}, chi2/dof = {custom_fit.chi2_dof:.3f}")
    print(f"Output directory: {outdir}")
    print("=" * 72)


__all__ = ["run_analysis", "default_custom_ratio_edges"]
