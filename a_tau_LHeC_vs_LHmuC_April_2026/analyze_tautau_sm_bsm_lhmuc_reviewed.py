from __future__ import annotations

import argparse
from pathlib import Path

from tautau_sm_bsm_common_reviewed import run_analysis


def _parse_edges(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="SM vs BSM tau-pair analysis for the LHmuC sample set."
    )
    parser.add_argument(
        "--sm-events-dir",
        type=Path,
        default=Path(
            "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/"
            "aa_tautau_SM_NP_0_SMEFTsim_top_alphaScheme_UFO_LHmuC/Events"
        ),
        help="Path to the SM Events directory containing run_01, run_02, ...",
    )
    parser.add_argument(
        "--bsm-events-dir",
        type=Path,
        default=Path(
            "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/"
            "aa_tautau_BSM_NP_2_SMEFTsim_top_alphaScheme_UFO_LHmuC/Events"
        ),
        help="Path to the BSM Events directory containing run_01, run_02, ...",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("plots_LHmuC_SM_vs_BSM_delta_aTau_0p005_weighted"),
        help="Directory where output plots and summary files will be written.",
    )
    parser.add_argument(
        "--delta-a-tau",
        type=float,
        default=0.005,
        help="Label value used in legends/titles for the BSM benchmark.",
    )
    parser.add_argument("--mass-plot-bins", type=int, default=100)
    parser.add_argument("--rapidity-plot-bins", type=int, default=20)
    parser.add_argument(
        "--ratio-bins",
        type=int,
        default=10,
        help="Uniform bins used for the presentation-level weighted ratio fit.",
    )
    parser.add_argument(
        "--custom-ratio-edges",
        type=str,
        default=None,
        help="Comma-separated custom mass-bin edges, e.g. '10,35,60,100,150,250,500'.",
    )
    parser.add_argument(
        "--display-lumi-fb",
        type=float,
        default=None,
        help="Optional luminosity value to display in titles. Default: use average effective luminosity from the samples.",
    )
    parser.add_argument(
        "--max-runs",
        type=int,
        default=None,
        help="Optional debug option: read only the first N run_* folders. Leave unset to read all runs.",
    )
    parser.add_argument(
        "--min-bin-events",
        type=int,
        default=10,
        help="Require at least this many SM and BSM events in a bin before using it in the ratio/fit.",
    )
    parser.add_argument("--mass-min", type=float, default=10.0)
    parser.add_argument("--mass-max", type=float, default=500.0)
    parser.add_argument("--rapidity-min", type=float, default=-10.0)
    parser.add_argument("--rapidity-max", type=float, default=10.0)
    parser.add_argument(
        "--collider-label",
        type=str,
        default="LHmuC@5.3 TeV",
        help="Text label used in plot titles.",
    )
    args = parser.parse_args()

    if args.mass_min >= args.mass_max:
        raise ValueError("--mass-min must be smaller than --mass-max")
    if args.rapidity_min >= args.rapidity_max:
        raise ValueError("--rapidity-min must be smaller than --rapidity-max")

    custom_edges = _parse_edges(args.custom_ratio_edges) if args.custom_ratio_edges else None

    run_analysis(
        collider_label=args.collider_label,
        sm_events_dir=args.sm_events_dir,
        bsm_events_dir=args.bsm_events_dir,
        outdir=args.outdir,
        delta_a_tau=args.delta_a_tau,
        mass_plot_bins=args.mass_plot_bins,
        rapidity_plot_bins=args.rapidity_plot_bins,
        ratio_bins=args.ratio_bins,
        mass_range=(args.mass_min, args.mass_max),
        rapidity_range=(args.rapidity_min, args.rapidity_max),
        custom_ratio_edges=custom_edges,
        display_lumi_fb=args.display_lumi_fb,
        max_runs=args.max_runs,
        min_bin_events=args.min_bin_events,
    )


if __name__ == "__main__":
    main()
