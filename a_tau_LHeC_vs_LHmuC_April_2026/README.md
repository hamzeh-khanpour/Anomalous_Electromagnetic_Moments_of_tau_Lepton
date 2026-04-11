# Anomalous Electromagnetic Moments of the Tau Lepton: LHeC vs LHmuC

This directory contains the April 2026 analysis setup for a direct comparison of photon-induced
$\gamma\gamma \to \tau^+\tau^-$ production at the **LHeC** and the **LHmuC**, focusing on the impact of a nonzero anomalous tau magnetic moment benchmark.

The current study compares:

- **LHeC @ 1.2 TeV**
- **LHmuC @ 3.7 TeV**
- **SM** vs **BSM** with **$\delta a_\tau = 0.008$**

using merged MadGraph event samples, differential distributions, and mass-dependent ratio fits.

---

## Directory contents

This subdirectory contains:

- `analyze_tautau_sm_bsm_lhec.py`  
  Driver script for the **LHeC** analysis.

- `analyze_tautau_sm_bsm_lhmuc.py`  
  Driver script for the **LHmuC** analysis.

- `tautau_sm_bsm_common.py`  
  Shared analysis utilities used by both collider-specific drivers.

- `plots_LHeC_SM_vs_BSM_delta_aTau_0p008/`  
  Output plots and summary text files for the LHeC benchmark.

- `plots_LHmuC_SM_vs_BSM_delta_aTau_0p008/`  
  Output plots and summary text files for the LHmuC benchmark.

---

## Physics goal

The purpose of this analysis is to compare how a benchmark anomalous tau magnetic moment modifies the
kinematics of photon-induced tau-pair production in two future lepton--hadron collider environments.

The main observables studied are:

- the invariant-mass spectrum
  $\frac{d\sigma}{dM_{\tau^+\tau^-}}$
- the pair-rapidity distribution
  $\frac{d\sigma}{dY_{\tau^+\tau^-}}$
- the ratio between BSM and SM predictions
  $\frac{\text{BSM}}{\text{SM}}$
- linear fits to the ratio as a function of $M_{\tau^+\tau^-}$

The physical motivation is that anomalous dipole-type couplings are expected to harden the invariant-mass spectrum, so the largest deviations from the Standard Model should appear in the high-mass tail.

---

## Event samples used in this study

For each collider and each theory hypothesis, the analysis is built from **five independent MadGraph runs**:

- `run_01`
- `run_02`
- `run_03`
- `run_04`
- `run_05`

Each run contains **1M unweighted events**, and the scripts automatically read all run directories found under the corresponding `Events/` folder, extract the LHE information, and merge the kinematic arrays into a single effective **5M-event sample**.

This is done separately for:

- **LHeC, SM**
- **LHeC, BSM**
- **LHmuC, SM**
- **LHmuC, BSM**

The integrated cross section is read directly from each LHE header, and the merged sample uses the event-count-weighted average cross section.

---

## Analysis strategy

The workflow implemented here is:

1. Read all `run_*` LHE or LHE.GZ files under the chosen `Events/` directory.
2. Parse final-state $\tau^+$ and $\tau^-$ four-momenta from each event.
3. Build the pair observables:
   - $M_{\tau^+\tau^-}$
   - $Y_{\tau^+\tau^-}$
4. Convert the binned event counts into weighted differential spectra using the total cross section read from the LHE headers.
5. Construct two types of ratio plots:
   - a **weighted physical ratio** based on differential spectra,
   - a **legacy-style count ratio** designed to reproduce the earlier plot style.
6. Fit the ratio with a simple linear form
   $y(M) = aM + b$
   to summarize whether the BSM enhancement grows with invariant mass.

The slope parameter is therefore used as a compact diagnostic of **tail hardening**.

---

## Generated outputs

For each collider, the scripts produce:

- `Invariant_mass_tau_pair_SM_vs_BSM_delta_aTau.png`
- `Rapidity_tau_pair_SM_vs_BSM_delta_aTau.png`
- `Ratio_BSM_over_SM_MtauTau_weighted.png`
- `Ratio_Obs_Exp_with_Luminosity.png`
- `Ratio_Obs_Exp_with_Luminosity_Fit_uniform.png`
- `Ratio_Obs_Exp_with_Luminosity_Fit_custom_bins.png`
- `analysis_summary.txt`

The `analysis_summary.txt` file reports:

- the number of runs successfully read,
- the total merged event count,
- the effective cross section of the merged sample,
- the effective luminosity implied by the sample,
- and the fit parameters for both weighted and legacy-style ratio analyses.

---

## Current benchmark results

For the currently committed benchmark $\delta a_\tau = 0.008$:

### LHeC

- SM cross section: `50.6950 pb`
- BSM cross section: `51.9226 pb`
- merged events per sample: `5,000,000`

### LHmuC

- SM cross section: `55.0944 pb`
- BSM cross section: `56.4608 pb`
- merged events per sample: `5,000,000`

In both machines, the BSM benchmark produces a ratio above unity in the hard region, indicating that the anomalous coupling makes the invariant-mass spectrum harder than the Standard Model prediction.

The LHmuC comparison shows a stronger rise of the ratio with invariant mass than the LHeC case, suggesting a more pronounced hard-tail enhancement for the same benchmark.

---

## Running the analysis

### LHeC

```bash
python3 analyze_tautau_sm_bsm_lhec.py
```

### LHmuC

```bash
python3 analyze_tautau_sm_bsm_lhmuc.py
```

Optional command-line arguments are available, for example:

```bash
python3 analyze_tautau_sm_bsm_lhec.py --display-lumi-fb 100
python3 analyze_tautau_sm_bsm_lhmuc.py --display-lumi-fb 100
```

The scripts currently assume local MG5 event directories of the form:

```text
/home/.../aa_tautau_SM_NP_0_SMEFTsim_top_alphaScheme_UFO_LHeC/Events
/home/.../aa_tautau_BSM_NP_2_SMEFTsim_top_alphaScheme_UFO_LHeC/Events
/home/.../aa_tautau_SM_NP_0_SMEFTsim_top_alphaScheme_UFO_LHmuC/Events
/home/.../aa_tautau_BSM_NP_2_SMEFTsim_top_alphaScheme_UFO_LHmuC/Events
```

If your files are stored elsewhere, edit the default paths in the two driver scripts or pass explicit command-line paths.

---

## Notes on the ratio plots

Two ratio styles are kept in this directory on purpose.

### 1. Weighted physical ratio

This uses

$\frac{(d\sigma/dM)_{\text{BSM}}}{(d\sigma/dM)_{\text{SM}}}$

and is the more directly physical comparison.

### 2. Legacy ratio plots

These are retained to reproduce the style of the earlier LHeC analysis plots and to facilitate visual continuity with previous versions of the study.

The linear fit

$y(M) = aM + b$

is not intended as a fundamental theory model. It is used as a compact summary of whether the BSM effect is approximately flat or grows toward high invariant mass.

---

## Requirements

The scripts use standard scientific Python packages:

```bash
pip install numpy matplotlib scipy
```

Python 3 is required.

---

## Reproducibility note

This directory is best viewed as a **research analysis workspace** rather than a polished standalone package.

The main reasons are:

- local MG5 paths are hard-coded by default,
- the analysis is tied to locally produced LHE event samples,
- the committed plot directories reflect a specific benchmark setup.

Still, the structure is now substantially cleaner than the earlier single-script exploratory workflow, since the collider-specific drivers and the common parsing/plotting utilities are separated.

---

## Related context inside the parent repository

This directory is part of the larger repository:

`Anomalous_Electromagnetic_Moments_of_tau_Lepton`

The parent repository also contains:

- earlier LHeC-only plotting scripts,
- SMEFT parameter utilities,
- UFO model archives,
- and CMS-related reference material used in connection with the anomalous tau moment study.

---

## Author

Hamzeh Khanpour (hamzeh.khanpour@agh.edu.pl)
