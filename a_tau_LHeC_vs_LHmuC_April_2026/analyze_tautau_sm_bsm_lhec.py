from __future__ import annotations

import argparse
from pathlib import Path

from tautau_sm_bsm_common import default_custom_ratio_edges, run_analysis


def _parse_edges(text: str) -> list[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="SM vs BSM tau-pair analysis for the LHeC sample set."
    )
    parser.add_argument(
        "--sm-events-dir",
        type=Path,
        default=Path("/home/hamzeh-khanpour/MG5_aMC_v3_6_6/aa_tautau_SM_NP_0_SMEFTsim_top_alphaScheme_UFO_LHeC/Events"),
        help="Path to the SM Events directory.",
    )
    parser.add_argument(
        "--bsm-events-dir",
        type=Path,
        default=Path("/home/hamzeh-khanpour/MG5_aMC_v3_6_6/aa_tautau_BSM_NP_2_SMEFTsim_top_alphaScheme_UFO_LHeC/Events"),
        help="Path to the BSM Events directory.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("plots_LHeC_SM_vs_BSM_delta_aTau_0p005"),
        help="Directory where output plots will be written.",
    )
    parser.add_argument(
        "--delta-a-tau",
        type=float,
        default=0.005,
        help="Label value used in legends/titles for the BSM benchmark.",
    )
    parser.add_argument("--mass-plot-bins", type=int, default=100)
    parser.add_argument("--rapidity-plot-bins", type=int, default=20)
    parser.add_argument("--legacy-ratio-bins", type=int, default=10)
    parser.add_argument(
        "--custom-ratio-edges",
        type=str,
        default=None,
        help="Comma-separated custom mass-bin edges for the second legacy fit plot.",
    )
    parser.add_argument(
        "--display-lumi-fb",
        type=float,
        default=None,
        help="Optional luminosity value to display in the legacy-style titles. Default: use average effective luminosity from the samples.",
    )
    args = parser.parse_args()

    custom_edges = _parse_edges(args.custom_ratio_edges) if args.custom_ratio_edges else None

    run_analysis(
        collider_label="LHeC@1.2 TeV",
        sm_events_dir=args.sm_events_dir,
        bsm_events_dir=args.bsm_events_dir,
        outdir=args.outdir,
        delta_a_tau=args.delta_a_tau,
        mass_plot_bins=args.mass_plot_bins,
        rapidity_plot_bins=args.rapidity_plot_bins,
        legacy_ratio_bins=args.legacy_ratio_bins,
        custom_ratio_edges=custom_edges,
        display_lumi_fb=args.display_lumi_fb,
    )


if __name__ == "__main__":
    main()
