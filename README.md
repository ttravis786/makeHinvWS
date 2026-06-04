# makeHinvWS

Combine workspace and datacard creator for the VBF H→invisible analysis.

**Run3 scripts** (`*Run3.C`, `*_run3_*.py`, `run_pipeline.py`) are the active
versions. The non-`Run3` files are the working Run2/UL reference implementations
and should not be modified.

---

# Run3 workspace and datacard production

## Overview

Four scripts must be run **in order** to produce a Combine datacard from pyRAT
plotting output.  All input ROOT files are read from:
```
hinvisible_mtr_vtr/root_files/{year}/{region}_{cat}_VBF/{var}/VBF_shapes.root
```

| # | Script | What it produces | Environment |
|---|--------|-----------------|-------------|
| 1 | `derive_run3_nlo_ratios.py` | `run3_nlo_sf_{year}_{cat}_{var}.root` — per-bin NLO yield-variation histograms used by the parametric WS | pyRAT/LCG |
| 2 | `makeSignalAndMCBackgroundWSRun3.C` | `signal_mc_bkgs_ws_{cat}_{year}_VBF{var}.root` — signal + small MC backgrounds as `RooDataHist` | CMSSW |
| 3 | `makeWS_percategoryRun3.C` | `param_ws_{year}_{cat}_VBF{var}.root` — parametric W/Z transfer-factor workspace | CMSSW |
| 4 | `make_run3_datacard.py` | `datacards/{cat}_{year}/card_{cat}_{year}.txt` — Combine datacard | pyRAT/LCG |

Steps 1 & 4 need the pyRAT/LCG_104 environment (`setup_environment.sh`).
Steps 2, 3, 5 & 6 need CMSSW (`CMSSW_14_1_0_pre4`).
**These two environments are incompatible in the same shell** (conflicting
`ROOTSYS`, `LD_LIBRARY_PATH`, `PYTHONPATH`).  Use `make_ws_and_datacard.sh`
to handle the switching automatically.

The QCD data-driven (`_qcdDD.root`) and HF-noise (`_noiseDD.root`) workspaces are
produced automatically by the pyRAT plotting step (`hinvisible_mtr_vtr/region_analysis/main.py`)
and are copied into the datacard directory by `make_run3_datacard.py`.

---

## Quickstart — one command

`make_ws_and_datacard.sh` runs the full pipeline with automatic environment
switching.  **No pre-activation of CMSSW or pyRAT is required** — each step
spawns a fresh subprocess via `env -i` so the two incompatible environments
never coexist in the same shell.

```bash
cd /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT

# Full pipeline — VTR only, blind expected limit
bash makeHinvWS/make_ws_and_datacard.sh

# Both categories
CATS="VTR MTR" bash makeHinvWS/make_ws_and_datacard.sh

# Workspaces + datacard only (skip limits)
SKIP_LIMITS=1 bash makeHinvWS/make_ws_and_datacard.sh

# Observed limit
OBSERVED=1 bash makeHinvWS/make_ws_and_datacard.sh
```

Environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `CATS` | `VTR` | Space-separated list of categories |
| `VAR` | `Mjj` | Fit variable (`Mjj` or `SignalScore`) |
| `SKIP_LIMITS` | unset | Set to `1` to stop after Step 4 (datacard) |
| `OBSERVED` | unset | Set to `1` to run combine on observed data |

`run_pipeline.py` is an alternative Python wrapper that also orchestrates all
steps, but it still requires CMSSW to be sourced in the calling shell before
running.

Available `--steps` for `run_pipeline.py`: `derive_nlo`, `param_ws`,
`signal_ws`, `datacard`, `combine_dc`, `limits` (default: all).

---

## Step-by-step manual instructions

> **Environment note:** Steps 1 & 4 require the pyRAT/LCG environment;
> Steps 2, 3, 5 & 6 require CMSSW.  If running manually, open two terminals
> (one per environment) or use `make_ws_and_datacard.sh` to avoid the conflict.

### Step 1 — NLO ratio histograms (`derive_run3_nlo_ratios.py`)  [pyRAT/LCG]

Reads the SR `VBF_shapes.root` and extracts per-bin yield-variation histograms
for the QCD W and Z processes.  These replace the Run2 `fnlo_SF_*` histograms
embedded in the shapes files and must be produced **before** the parametric
workspace.

```bash
cd /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT
micromamba deactivate
source setup_environment.sh

python3 makeHinvWS/derive_run3_nlo_ratios.py \
    --year Run3Summer22_to_Run3Summer23BPix \
    --cat VTR \
    --var Mjj \
    --output-dir makeHinvWS/
```

Output: `makeHinvWS/run3_nlo_sf_Run3Summer22_to_Run3Summer23BPix_VTR_Mjj.root`

The four NLO histograms saved (one per process × four corrections):
- `{proc}_fnlo_SF_mufUp` — factorisation scale up (W-yield ratio, matching Run2)
- `{proc}_fnlo_SF_murUp` — renormalisation scale up (W-yield ratio)
- `{proc}_fnlo_SF_pdfUp` — PDF up (W/Z ratio variation; currently flat in pyRAT)
- `{proc}_fnlo_SF_EWK_corrUp` — EWK correction on QCD W/Z (W/Z ratio variation, per-bin)

---

### Step 2 — Signal + MC background workspace (`makeSignalAndMCBackgroundWSRun3.C`)  [CMSSW]

Packages signal (vbfH, ggH, WH, ZH) and small MC backgrounds (Top, Diboson,
QCD\_Zll\_NLONew, EWK\_Zll) as `RooDataHist` shapes with all systematic
variations.

```bash
cd /vols/cms/tt1020/Combine/CMSSW_14_1_0_pre4 && cmsenv
cd /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT/makeHinvWS
YEAR="Run3Summer22_to_Run3Summer23BPix"

root -l -b -q "makeSignalAndMCBackgroundWSRun3.C+(\"${YEAR}\",\"VTR\",false,true,true)"
root -l -b -q "makeSignalAndMCBackgroundWSRun3.C+(\"${YEAR}\",\"MTR\",false,true,true)"
```

Function signature:
```
makeSignalAndMCBackgroundWSRun3(year, cat, classifier, turn_off_jecs, turn_off_btag)
  classifier    false=Mjj  true=SignalScore
  turn_off_jecs true  = flat dummy JEC/JER shapes (use until Run3 JES file available)
                false = read shapes from vbf_shape_jes_uncs.root
  turn_off_btag true  = flat dummy b-tag shapes (Run3 b-tag not yet implemented)
```

---

### Step 3 — Parametric workspace (`makeWS_percategoryRun3.C`)

Builds per-bin `RooRealVar` Z yields, W/Z and CR/SR transfer-factor
`RooFormulaVar`s, and `RooParametricHist` PDFs for Z\_SR, W\_SR and each
control region.  Requires `RooParametricHist` from HiggsAnalysis/CombinedLimit.

**Reads** `run3_nlo_sf_*.root` from Step 1 — that file must exist in the same
directory before running this step.

```bash
# CMSSW required
cd /vols/cms/tt1020/Combine/CMSSW_14_1_0_pre4 && cmsenv
cd /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT/makeHinvWS
YEAR="Run3Summer22_to_Run3Summer23BPix"

root -l -b -q "makeWS_percategoryRun3.C+(\"${YEAR}\",\"VTR\",false,true,true)"
root -l -b -q "makeWS_percategoryRun3.C+(\"${YEAR}\",\"MTR\",false,true,true)"
```

Function signature:
```
makeWS_percategoryRun3(year, cat, classifier, turn_off_jecs, turn_off_btag)
  turn_off_jecs true  = omit JEC/JER terms from TF formula (use for preliminary fits)
                false = read TF corrections from vbf_jes_jer_tf_uncs.root
  turn_off_btag reserved for future use (no effect currently)
```

#### W/Z ratio transfer factor structure

The parametric workspace builds these nuisance terms on the W/Z SR ratio TF
(matching Run2 conventions):

| Term | Parameter type | Source |
|------|---------------|--------|
| muR scale (W-yield ratio) | Correlated across bins | `fnlo_SF_murUp` / W nominal |
| muF scale (W-yield ratio) | Correlated across bins | `fnlo_SF_mufUp` / W nominal |
| PDF (W/Z ratio variation) | Correlated across bins | `fnlo_SF_pdfUp` W/Z vs nominal |
| W/Z stat | Per-bin | MC stat errors |
| EWK correction on strong proc | Per-bin (like Run2) | `fnlo_SF_EWK_corrUp` W/Z ratio |
| JER / JES×11 | Correlated | `vbf_jes_jer_tf_uncs.root` or flat dummy |
| Tau veto, e-veto, μ-veto-id, μ-veto-iso | Correlated (shared with CR TF) | Histogram ratios; currently flat — activates when pyRAT veto weights land |

---

### Step 4 — Datacard (`make_run3_datacard.py`)  [pyRAT/LCG]

Assembles the Combine datacard, copies workspace files into the datacard
directory, and reads the number of bins directly from `VBF_shapes.root`.
Must run in the pyRAT/LCG environment — **do not run this under CMSSW** as the
LCG and CMSSW library paths are incompatible.

```bash
cd /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT
micromamba deactivate
source setup_environment.sh

python3 makeHinvWS/make_run3_datacard.py \
    --cat VTR \
    --campaigns Run3Summer22 Run3Summer22EE Run3Summer23 Run3Summer23BPix \
    --var Mjj
```

Output: `datacards/VTR_Run3Summer22_to_Run3Summer23BPix/card_VTR_....txt`

---

### Step 5 — Combine cards and run limits  [CMSSW]

```bash
cd /vols/cms/tt1020/Combine/CMSSW_14_1_0_pre4 && cmsenv
cd /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT/datacards
YEAR="Run3Summer22_to_Run3Summer23BPix"

# Combine VTR + MTR
combineCards.py \
    VTR_Run3=VTR_${YEAR}/card_VTR_${YEAR}.txt \
    MTR_Run3=MTR_${YEAR}/card_MTR_${YEAR}.txt \
    > run3_combined.txt

# Expected (blinded) limit
combine -M AsymptoticLimits -d run3_combined.txt -m 125 --run blind -t -1 -n Run3
```

---

## Run2 vs Run3 differences

### Workspace

| Aspect | Run2 | Run3 |
|--------|------|------|
| Input ROOT file layout | Per-region subdirectories | Flat file per region: `{proc}_calc_weight`, `{proc}_weight_{syst}_{up/down}` |
| NLO ratio source | `fnlo_SF_*` histograms embedded in shapes files | `derive_run3_nlo_ratios.py` → `run3_nlo_sf_*.root` |
| EWK correction on strong proc | Per-bin nuisance `wzratioEWK_on_strong[nT][nB]` | Same (fixed to match Run2; was wrongly correlated in original Run3) |
| muF/muR treatment | W-yield-only ratio (not W/Z ratio) | Same |
| Lepton veto in W/Z TF | Hardcoded `1/1.01` for e-veto + τ-veto | Histogram-derived (currently flat; activates when pyRAT weights land) |
| JEC/JER | Always from `vbf_jes_jer_tf_uncs.root` | `turn_off_jecs=true` flag (flat dummies for preliminary) |
| B-tagging | Full b-tag term | `turn_off_btag=true` (flat dummies) |
| Trigger SF | Not in workspace | `CMS_Trigger_{cat}_{year}` shape syst (per-event weight from pyRAT) |
| Prefiring | lnN nuisance | Not applicable in Run3 |

### Datacard

| Aspect | Run2 | Run3 |
|--------|------|------|
| Parametric WS params | `wzratioQCDcorrSyst_muR/muF/pdf` (correlated) + `wzratioEWK_on_strong_bin{i}` (per-bin) | Same naming scheme |
| Trigger | `trigger_SF_stat lnN 1.05` + `trigger_SF_syst lnN 1.01` (separate, flat) | Single `CMS_Trigger shape` from per-event weight variation (has mild Mjj shape ~1.5–2%) — no separate lnN needed |
| Photon CR | Included (external WS from another group) | Not yet included |
| Luminosity | 13 TeV values | 13.6 TeV, ~1.5% preliminary |
| Signal theory lnN | Updated Run2/YR4 values | YR4 values as placeholder — update for Run3 |

---

# Run2/UL instructions (legacy)

## Setting up fast_datacard (creating the env — only once)

```bash
source /afs/cern.ch/user/$U/$USER/miniconda2/etc/profile.d/conda.sh
conda create -n test_datacard_new
conda activate test_datacard_new

git clone ssh://git@gitlab.cern.ch:7999/fast-hep/public/fast-datacard.git
source /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-centos7-gcc8-opt/setup.sh
cd fast-datacard
make install-dev2

cd ..
export PATH=~/.local/bin:$PATH
```

## Running fast_datacard

```bash
git clone git@github.com:vukasinmilosevic/makeHinvWS.git
cd makeHinvWS/
# copy outputs of previous steps, normally located at IC under /vols/cms/VBFHinv/
. run_fast_datacard.sh
```

## Making the Run2 workspaces

From the `fast_datacard` output directory (requires CMSSW):

```bash
cd test_df_MTR_2017_2020v1/
root
root> .L ../makeWS_percategory.C++
root> makeWS_percategory("2017","MTR")
```

Limit steps and datacards: https://gitlab.cern.ch/cms-hcg/cadi/hig-20-003

---

# HIG-20-003

VBF H→invisible analysis datacards for HIG-20-003.

[![Cadi](https://amarini.web.cern.ch/amarini/ChargedHiggs/ci-status-images/img/icecave/regular/build-status/cadi-pas.png)](http://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-20-003)

Results produced with `CMSSW_10_2_13` and Combine `v8.1.0`.

## Combined datacards

Run `./mkCombinations.sh` or copy the snippet below:

```bash
c2017_MTR=MTR_2017/all_percategory.txt
c2017_VTR=VTR_2017/all_percategory.txt
p2017_MTR=photons/card_vbf_photons_2017.txt
c2018_MTR=MTR_2018/all_percategory.txt
c2018_VTR=VTR_2018/all_percategory.txt
p2018_MTR=photons/card_vbf_photons_2018.txt

combineCards.py MTR_2017=${c2017_MTR} VTR_2017=${c2017_VTR} photon_cr_2017=${p2017_MTR} > all_2017.txt
combineCards.py MTR_2018=${c2018_MTR} VTR_2018=${c2018_VTR} photon_cr_2018=${p2018_MTR} > all_2018.txt
combineCards.py MTR_2017=${c2017_MTR} MTR_2018=${c2018_MTR} VTR_2017=${c2017_VTR} VTR_2018=${c2018_VTR} \
    photon_cr_2017=${p2017_MTR} photon_cr_2018=${p2018_MTR} > alltime.txt
```
