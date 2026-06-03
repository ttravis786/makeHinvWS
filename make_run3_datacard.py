#!/usr/bin/env python3
"""
make_run3_datacard.py

Creates a CMS Combine datacard directory for the Run3 VBF H→invisible analysis.

Steps:
  1. Creates a labelled subdirectory under datacards/
  2. Copies workspace files from makeHinvWS/ and root_files/ into it
  3. Writes the .txt datacard

Usage (from /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT/):
    python3 makeHinvWS/make_run3_datacard.py \\
        --cat VTR \\
        --campaigns Run3Summer22 Run3Summer22EE Run3Summer23 Run3Summer23BPix \\
        --var Mjj

Known limitations / TODOs:
  - JEC and b-tag systematics are flat dummies (workspace produced with turn_off_jecs=True).
    Update the workspace macros and re-run once real JEC shapes are available.
  - Signal theory lnN values (QCDscale, pdf) are taken from Run2/YR4; update for Run3.
  - No photon CR included (Run2 used an external photon WS from another group).
  - Lumi uncertainty uses preliminary Run3 values.
"""

import os
import sys
import shutil
import argparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from utils.cms.year_run_utils import CampaignInfoMap  # noqa: E402

# ---------------------------------------------------------------------------
# Workspace object name helpers
# ---------------------------------------------------------------------------

PYRAT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def param_ws_name(cat, year, var):
    return f"param_ws_{year}_{cat}_VBF{var}.root"


def signal_ws_name(cat, year, var):
    return f"signal_mc_bkgs_ws_{cat}_{year}_VBF{var}.root"


def qcd_ws_name(cat, year):
    region = f"SR_{cat}_VBF"
    return f"out_{region}_{year}.root_qcdDD.root"


def noise_ws_name(cat, year):
    region = f"SR_{cat}_VBF"
    return f"{region}_{year}_noiseDD.root"


# Combine datacard process names → Run3 workspace object names
# (workspace objects from makeSignalAndMCBackgroundWSRun3 and makeWS_percategoryRun3)
SIGNAL_PROC_TO_WS = {
    "qqH_hinv": "vbfH",
    "ggH_hinv": "ggH",
    "WH_hinv":  "WH",
    "ZH_hinv":  "ZH",
}

MC_BKG_PROC_TO_WS = {
    "TOP":    "Top",
    "VV":     "Diboson",
    "DY":     "QCD_Zll_NLONew",
    "EWKZll": "EWK_Zll",
}

# Parametric PDF objects in param_ws — these have no $SYSTEMATIC in shapes lines
# (systematics are absorbed into the TF formula param nuisances)
PARAM_SR = {
    "ZJETS":    "QCDZ_SR",
    "EWKZNUNU": "EWKZ_SR",
    "WJETS":    "QCDW_SR",
    "EWKW":     "EWKW_SR",
}

# MC histogram objects for small backgrounds in Z CRs — come from param_ws
# (WJETS / EWKW contamination in Zee/Zmumu, with $SYSTEMATIC shapes)
PARAM_MC_PROCS = ["WJETS", "EWKW"]  # in Z CRs these use param_ws MC histo

# Data-driven processes — no $SYSTEMATIC shapes in datacard
DATA_DRIVEN = ["QCD", "HFNoise"]

# Region aliases matching workspace object names
REGION_ALIASES = {
    "SR":    "SR",
    "WENU":  "Wenu",
    "WMUNU": "Wmunu",
    "ZEE":   "Zee",
    "ZMUMU": "Zmumu",
}

# Process indices: ≤0 signal, >0 background
PROC_INDEX = {
    "qqH_hinv": -3, "ggH_hinv": -2, "WH_hinv": -1, "ZH_hinv": 0,
    "TOP": 1, "QCD": 2, "WJETS": 3, "EWKZll": 4, "VV": 5,
    "DY": 6, "EWKW": 7, "EWKZNUNU": 8, "HFNoise": 9, "ZJETS": 10,
}

# ---------------------------------------------------------------------------
# Systematic definitions
# ---------------------------------------------------------------------------

def build_syst_names(cat, year):
    """Returns the 14 weight-based shape systematic CMS names.

    CMS_Trigger replaces the Run2 lnN trigger_SF_stat + trigger_SF_syst pair:
    in Run3, per-event trigger weights from pyRAT capture the full uncertainty
    including mild Mjj shape (~1.5-2%), so no separate lnN is needed.
    """
    return {
        "CMS_Trigger":                   f"CMS_Trigger_{cat}_{year}",
        "CMS_eff_tauveto":               f"CMS_eff_tauveto_{year}",
        "CMS_eff_eTight_idiso":          f"CMS_eff_eTight_idiso_{year}",
        "CMS_eff_eTight_idiso_highpt":   f"CMS_eff_eTight_idiso_highpt_{year}",
        "CMS_eff_eVeto_idiso_veto":      f"CMS_eff_eVeto_idiso_veto_{year}",
        "CMS_eff_muTight_id":            "CMS_eff_muTight_id",
        "CMS_eff_muLoose_id_veto":       "CMS_eff_muLoose_id_veto",
        "CMS_eff_muLoose_iso_veto":      "CMS_eff_muLoose_iso_veto",
        "CMS_pdf":                       "CMS_pdf",
        "CMS_ewk_correction":            "CMS_ewk_correction",
        "CMS_qcd_correction":            "CMS_qcd_correction",
        "CMS_scale_mur":                 "CMS_scale_mur",
        "CMS_scale_muf":                 "CMS_scale_muf",
        "CMS_pileup":                    "CMS_pileup",
    }


def build_btag_syst_names(year):
    return {
        "CMS_fake_b": f"CMS_fake_b_{year}",  # most processes
        "CMS_eff_b":  f"CMS_eff_b_{year}",   # Top only
    }


def build_jec_syst_names(year):
    jec_systs = {}
    jec_systs["CMS_res_j"] = f"CMS_res_j_{year}"
    for name in ["jesAbsolute", f"jesAbsolute_{year}", "jesBBEC1", f"jesBBEC1_{year}",
                 "jesEC2", f"jesEC2_{year}", "jesFlavorQCD", "jesHF", f"jesHF_{year}",
                 "jesRelativeBal", f"jesRelativeSample_{year}"]:
        jec_systs[f"CMS_scale_j_{name}"] = f"CMS_scale_j_{name}"
    return jec_systs


def is_mc_shape(dc_proc, bin_alias):
    """Returns True if this (process, bin) combination uses MC histogram shapes.

    Parametric PDF processes (RooParametricHist) have no $SYSTEMATIC shapes —
    their systematics are absorbed into the TF formula param nuisances instead.
    The same Combine process name can be parametric in one bin and an MC histogram
    in another:
      - DY / EWKZll:  parametric in ZEE/ZMUMU, MC histogram in SR/WENU/WMUNU
      - WJETS / EWKW: always parametric (SR + W CRs)
    """
    if dc_proc in DATA_DRIVEN:
        return False
    if dc_proc in PARAM_SR:                                          # always parametric
        return False
    if dc_proc in ("DY", "EWKZll") and bin_alias in ("ZEE", "ZMUMU"):  # parametric in Z CRs
        return False
    return True


def has_weight_syst(dc_proc, syst_key, bin_alias):
    """Returns True if the workspace has up/down shape histograms for this combination."""
    if not is_mc_shape(dc_proc, bin_alias):
        return False
    # ewk_correction only for DY (QCD_Zll_NLONew)
    if syst_key == "CMS_ewk_correction" and dc_proc != "DY":
        return False
    # qcd_correction only for EWKZll
    if syst_key == "CMS_qcd_correction" and dc_proc != "EWKZll":
        return False
    # mur/muf skipped for signal processes
    if syst_key in ("CMS_scale_mur", "CMS_scale_muf") and dc_proc in SIGNAL_PROC_TO_WS:
        return False
    return True


def has_btag_syst(dc_proc, bin_alias):
    return is_mc_shape(dc_proc, bin_alias)


# ---------------------------------------------------------------------------
# Process / region layout
# ---------------------------------------------------------------------------

def get_bin_procs(bin_alias):
    """Returns ordered list of (datacard_process, ws_source) for a given bin alias."""
    if bin_alias == "SR":
        procs = (
            list(SIGNAL_PROC_TO_WS.keys()) +    # signals from signal_mc_bkgs_ws
            list(PARAM_SR.keys()) +             # parametric from param_ws
            list(MC_BKG_PROC_TO_WS.keys()) +    # MC from signal_mc_bkgs_ws
            DATA_DRIVEN
        )
    elif bin_alias in ("WENU", "WMUNU"):
        # parametric WJETS/EWKW + MC small bkgs
        procs = ["WJETS", "EWKW"] + list(MC_BKG_PROC_TO_WS.keys())
    elif bin_alias in ("ZEE", "ZMUMU"):
        # parametric DY/EWKZll + MC small bkgs (TOP, VV)
        # WJETS/EWKW contamination in Z CRs is negligible for VTR/MTR and
        # often has zero integral — exclude from preliminary fits to avoid
        # missing-object errors from the workspace.
        procs = ["DY", "EWKZll", "TOP", "VV"]
    else:
        raise ValueError(f"Unknown bin alias: {bin_alias}")
    return procs


def shapes_line(dc_proc, bin_label, ws_file_path, ws_obj, ws_syst_obj=None):
    if ws_syst_obj:
        return f"shapes {dc_proc:<16} {bin_label:<30} {ws_file_path} {ws_obj} {ws_syst_obj}"
    return f"shapes {dc_proc:<16} {bin_label:<30} {ws_file_path} {ws_obj}"


# ---------------------------------------------------------------------------
# Datacard writer
# ---------------------------------------------------------------------------

def write_datacard(args, label, ws_dir, lumi, energy, nB):
    cat, year, var = args.cat, args.year, args.var
    dc_path = os.path.join("datacards", label, f"card_{label}.txt")
    bins_short = ["ZEE", "ZMUMU", "WENU", "WMUNU", "SR"]
    bin_labels = [f"{label}_{b}" for b in bins_short]

    # Use bare filenames — the datacard and workspace files live in the same
    # directory, so paths must be relative to that directory (where combine runs).
    param_ws  = param_ws_name(cat, year, var)
    signal_ws = signal_ws_name(cat, year, var)
    qcd_ws    = qcd_ws_name(cat, year)
    noise_ws  = noise_ws_name(cat, year)

    syst_names  = build_syst_names(cat, year)
    btag_names  = build_btag_syst_names(year)
    jec_names   = build_jec_syst_names(year)

    # Build full (bin, process) table
    all_bin_proc = []  # list of (bin_short, bin_label, dc_proc)
    for b_short, b_label in zip(bins_short, bin_labels):
        for p in get_bin_procs(b_short):
            all_bin_proc.append((b_short, b_label, p))

    lines = []

    # ---- header -------------------------------------------------------------
    lines.append(f"# Run3 VBF H→invisible datacard — {cat}, {year}, {var}")
    lines.append(f"# Generated by makeHinvWS/make_run3_datacard.py")
    lines.append(f"imax {len(bins_short)}  number of channels")
    lines.append("jmax *  number of processes minus 1")
    lines.append("kmax *  number of nuisance parameters")
    lines.append("-" * 100)

    # ---- shapes lines -------------------------------------------------------
    # param_ws objects: labels from makeWS_percategoryRun3.C
    lcat  = f"{cat}_"
    lyear = f"{year}_"
    cr_pref = f"{lcat}{lyear}"

    for b_short, b_label in zip(bins_short, bin_labels):
        alias = REGION_ALIASES[b_short]
        for p in get_bin_procs(b_short):
            if p in SIGNAL_PROC_TO_WS:
                ws_obj = f"wspace_signal:{SIGNAL_PROC_TO_WS[p]}_hist_{alias}"
                ws_syst = f"wspace_signal:{SIGNAL_PROC_TO_WS[p]}_hist_{alias}_$SYSTEMATIC"
                lines.append(shapes_line(p, b_label, signal_ws, ws_obj, ws_syst))

            elif p in PARAM_SR and b_short == "SR":
                # Parametric PDF in SR — no $SYSTEMATIC
                lines.append(shapes_line(p, b_label, param_ws, f"wspace:{PARAM_SR[p]}"))

            elif p in ("WJETS", "EWKW") and b_short in ("WENU", "WMUNU"):
                # Parametric in W CRs
                ptype = "QCD" if p == "WJETS" else "EWK"
                obj = f"wspace:{cr_pref}{ptype}V_{alias}"
                lines.append(shapes_line(p, b_label, param_ws, obj))

            elif p in ("DY", "EWKZll") and b_short in ("ZEE", "ZMUMU"):
                # Parametric in Z CRs
                ptype = "QCD" if p == "DY" else "EWK"
                obj = f"wspace:{cr_pref}{ptype}V_{alias}"
                lines.append(shapes_line(p, b_label, param_ws, obj))

            elif p in MC_BKG_PROC_TO_WS and b_short != "SR":
                # MC small background in CRs (from signal_mc_bkgs_ws, only Top/Diboson/DY/EWKZll)
                if p in ("TOP", "VV"):
                    ws_n  = MC_BKG_PROC_TO_WS[p]
                    ws_obj  = f"wspace_signal:{ws_n}_hist_{alias}"
                    ws_syst = f"wspace_signal:{ws_n}_hist_{alias}_$SYSTEMATIC"
                    lines.append(shapes_line(p, b_label, signal_ws, ws_obj, ws_syst))
                elif p in ("DY", "EWKZll") and b_short in ("WENU", "WMUNU"):
                    ws_n  = MC_BKG_PROC_TO_WS[p]
                    ws_obj  = f"wspace_signal:{ws_n}_hist_{alias}"
                    ws_syst = f"wspace_signal:{ws_n}_hist_{alias}_$SYSTEMATIC"
                    lines.append(shapes_line(p, b_label, signal_ws, ws_obj, ws_syst))

            elif p in MC_BKG_PROC_TO_WS and b_short == "SR":
                ws_n  = MC_BKG_PROC_TO_WS[p]
                ws_obj  = f"wspace_signal:{ws_n}_hist_{alias}"
                ws_syst = f"wspace_signal:{ws_n}_hist_{alias}_$SYSTEMATIC"
                lines.append(shapes_line(p, b_label, signal_ws, ws_obj, ws_syst))

            elif p in ("WJETS", "EWKW") and b_short in ("ZEE", "ZMUMU"):
                # MC histogram contamination in Z CRs — from param_ws with $SYSTEMATIC
                ws_n  = "QCD_Wjet_NLONew" if p == "WJETS" else "EWK_W"
                ws_obj  = f"wspace:{ws_n}_hist_{alias}"
                ws_syst = f"wspace:{ws_n}_hist_{alias}_$SYSTEMATIC"
                lines.append(shapes_line(p, b_label, param_ws, ws_obj, ws_syst))

            elif p == "QCD":
                lines.append(shapes_line(p, b_label, qcd_ws, "qcd_wspace:QCD_DD"))

            elif p == "HFNoise":
                lines.append(shapes_line(p, b_label, noise_ws, "noise_wspace:QCD_noise"))

        # data_obs for each bin from param_ws
        lines.append(shapes_line("data_obs", b_label, param_ws, f"wspace:data_obs_{alias}"))

    lines.append("-" * 100)

    # ---- bin / observation --------------------------------------------------
    lines.append("bin          " + "  ".join(f"{b:<30}" for b in bin_labels))
    lines.append("observation  " + "  ".join("-1" + " " * 29 for _ in bin_labels))
    lines.append("-" * 100)

    # ---- process table ------------------------------------------------------
    bin_col    = "  ".join(f"{bp[1]:<30}" for bp in all_bin_proc)
    proc_col   = "  ".join(f"{bp[2]:<30}" for bp in all_bin_proc)
    idx_col    = "  ".join(f"{PROC_INDEX.get(bp[2], 99):<30}" for bp in all_bin_proc)
    rate_col   = "  ".join(f"{-1:<30}" for bp in all_bin_proc)
    lines.append("bin           " + bin_col)
    lines.append("process       " + proc_col)
    lines.append("process       " + idx_col)
    lines.append("rate          " + rate_col)
    lines.append("-" * 100)

    # ---- helper: build a systematics row ------------------------------------
    def syst_row(syst_label, syst_type, entries):
        """entries: list aligned to all_bin_proc, each '1.0' or '-'."""
        return f"{syst_label:<45} {syst_type}  " + "  ".join(f"{e:<30}" for e in entries)

    def entries_for(pred):
        return [("1.0" if pred(b_short, p) else "-") for b_short, _, p in all_bin_proc]

    # ---- shape systematics: Run3 weight-based -------------------------------
    for key, cms_name in syst_names.items():
        entries = entries_for(
            lambda b, p, k=key, cn=cms_name: has_weight_syst(p, cn, b)
        )
        if any(e != "-" for e in entries):
            lines.append(syst_row(cms_name, "shape", entries))

    # b-tag (dummy flat shapes)
    entries_fakeb = entries_for(lambda b, p: has_btag_syst(p, b) and p != "TOP")
    entries_effb  = entries_for(lambda b, p: p == "TOP" and is_mc_shape(p, b))
    if any(e != "-" for e in entries_fakeb):
        lines.append(syst_row(btag_names["CMS_fake_b"], "shape", entries_fakeb))
    if any(e != "-" for e in entries_effb):
        lines.append(syst_row(btag_names["CMS_eff_b"], "shape", entries_effb))

    # JEC (dummy flat shapes)
    for _, cms_jec in jec_names.items():
        entries = entries_for(lambda b, p: is_mc_shape(p, b))
        if any(e != "-" for e in entries):
            lines.append(syst_row(cms_jec, "shape", entries))

    # ---- lnN systematics ----------------------------------------------------
    # Luminosity (13.6 TeV Run3)
    lumi_unc = 1.015  # ~1.5% preliminary; update when official value available
    entries_all_mc = entries_for(lambda b, p: p not in DATA_DRIVEN)
    lines.append(syst_row(f"lumi_13p6TeV_{year}", "lnN",
                           [f"{lumi_unc}" if e == "1.0" else "-" for e in entries_all_mc]))

    # Signal theory (YR4 values at mH=125 GeV; use for preliminary)
    _sig_only = lambda b, p, sp: p == sp
    lines.append(syst_row("QCDscale_qqH", "lnN",
        ["0.997/1.004" if p == "qqH_hinv" else "-" for _, _, p in all_bin_proc]))
    lines.append(syst_row("QCDscale_ggH", "lnN",
        ["0.933/1.046" if p == "ggH_hinv" else "-" for _, _, p in all_bin_proc]))
    lines.append(syst_row("QCDscale_WH", "lnN",
        ["0.993/1.006" if p == "WH_hinv" else "-" for _, _, p in all_bin_proc]))
    lines.append(syst_row("QCDscale_ZH", "lnN",
        ["0.969/1.038" if p == "ZH_hinv" else "-" for _, _, p in all_bin_proc]))
    lines.append(syst_row("pdf_Higgs_qq", "lnN",
        ["1.021" if p in ("qqH_hinv", "WH_hinv", "ZH_hinv") else "-" for _, _, p in all_bin_proc]))
    lines.append(syst_row("pdf_Higgs_gg", "lnN",
        ["1.039" if p == "ggH_hinv" else "-" for _, _, p in all_bin_proc]))
    lines.append(syst_row("QCDscale_VV", "lnN",
        ["1.15" if p == "VV" else "-" for _, _, p in all_bin_proc]))
    lines.append(syst_row("QCDscale_ttbar", "lnN",
        ["1.1" if p == "TOP" else "-" for _, _, p in all_bin_proc]))

    # QCD multijet normalisation (data-driven uncertainty ~40%)
    lines.append(syst_row(f"VBF_MultiJetQCD_FitError_{label}", "lnN",
        ["1.4" if p == "QCD" else "-" for _, _, p in all_bin_proc]))

    # HF noise normalisation (~20%)
    lines.append(syst_row(f"VBF_NoiseHF_SysError_{label}", "lnN",
        ["1.2" if p == "HFNoise" else "-" for _, _, p in all_bin_proc]))

    # NOTE: Trigger uncertainty is covered by the CMS_Trigger_<cat>_<year> shape systematic
    # above (from per-event pyRAT weight variations with mild Mjj shape).
    # Run2 used separate lnN stat+syst because it had no per-event trigger weight — do
    # NOT add a redundant trigger_SF_stat lnN here.

    lines.append("-" * 100)

    # ---- param / flatParam lines -------------------------------------------
    # Z SR free bin yields
    for i in range(1, nB + 1):
        lines.append(f"{cr_pref}QCDZ_SR_bin{i}  flatParam")

    # WZ ratio theory nuisances correlated across bins (muR, muF, pdf)
    for t in ("QCD", "EWK"):
        for nuis in ("MURSyst", "MUFSyst", "PDFSyst"):
            lines.append(f"{lcat}{t}wzratio{nuis}  param  0.0  1")

    # WZ ratio EWK correction on strong proc — per-bin (matches Run2 wzratioEWK_on_strong)
    # Name: lCategory + lType + "wzratio_EWK_corr_on_Strong_bin" + iB  (no year in name)
    for t in ("QCD", "EWK"):
        for i in range(1, nB + 1):
            lines.append(f"{lcat}{t}wzratio_EWK_corr_on_Strong_bin{i}  param  0.0  1")

    # WZ ratio stat (per bin per type)
    for t in ("QCD", "EWK"):
        for i in range(1, nB + 1):
            lines.append(f"{cr_pref}{t}wzratio_stat_bin{i}  param  0.0  1")

    # EWK/QCD Z ratio stat (per bin)
    for i in range(1, nB + 1):
        lines.append(f"{cr_pref}ewkqcdratio_stat_bin{i}  param  0.0  1")

    # CR/SR TF stat nuisances (per type per CR per bin)
    for t in ("QCD", "EWK"):
        for cr_alias in ("Wenu", "Wmunu", "Zee", "Zmumu"):
            for i in range(1, nB + 1):
                lines.append(f"{cr_pref}{t}TF_{cr_alias}_stat_bin{i}  param  0.0  1")

    # JEC/JER params — always present in workspace (flat dummy 1.0 when turn_off_jecs=True)
    lines.append(f"CMS_res_j_{year}  param  0.0  1")
    for _, cms_jec in list(jec_names.items())[1:]:  # remaining JES (jer already above)
        lines.append(f"{cms_jec}  param  0.0  1")

    dc_text = "\n".join(lines) + "\n"
    with open(dc_path, "w") as f:
        f.write(dc_text)
    print(f"Wrote: {dc_path}")
    return dc_path


# ---------------------------------------------------------------------------
# File copying
# ---------------------------------------------------------------------------

def copy_workspace_files(args, label):
    cat, year, var = args.cat, args.year, args.var
    ws_dir = os.path.join("datacards", label)
    os.makedirs(ws_dir, exist_ok=True)

    issues = []

    def copy_file(src, dest_name):
        dest = os.path.join(ws_dir, dest_name)
        if not os.path.exists(src):
            issues.append(f"MISSING: {src}")
            return
        shutil.copy2(src, dest)
        print(f"  Copied: {os.path.basename(src)} → {ws_dir}/")

    # Workspace files from makeHinvWS/
    mkhws = os.path.join(PYRAT_ROOT, "makeHinvWS")
    copy_file(os.path.join(mkhws, param_ws_name(cat, year, var)),  param_ws_name(cat, year, var))
    copy_file(os.path.join(mkhws, signal_ws_name(cat, year, var)), signal_ws_name(cat, year, var))

    # QCD DD and HF noise workspaces from root_files/
    region = f"SR_{cat}_VBF"
    rf_dir = os.path.join(PYRAT_ROOT, "hinvisible_mtr_vtr", "root_files", year, region, var)
    copy_file(os.path.join(rf_dir, qcd_ws_name(cat, year)),   qcd_ws_name(cat, year))
    copy_file(os.path.join(rf_dir, noise_ws_name(cat, year)), noise_ws_name(cat, year))

    return ws_dir, issues


# ---------------------------------------------------------------------------
# Read number of bins from workspace
# ---------------------------------------------------------------------------

def get_nB(cat, year, var):
    """Read number of bins from the data_obs histogram in VBF_shapes.root.
    This avoids opening the workspace (which needs CombinedLimit loaded)."""
    try:
        import ROOT  # noqa: PLC0415
        ROOT.gROOT.SetBatch(True)
        shapes_path = os.path.join(
            PYRAT_ROOT, "hinvisible_mtr_vtr", "root_files",
            year, f"SR_{cat}_VBF", var, "VBF_shapes.root"
        )
        f = ROOT.TFile.Open(shapes_path)
        if not f or f.IsZombie():
            raise RuntimeError(f"Cannot open {shapes_path}")
        h = f.Get("data_obs")
        if not h:
            raise RuntimeError("data_obs not found")
        nB = h.GetNbinsX()
        f.Close()
        return nB
    except Exception as e:
        print(f"  WARNING: could not read nB from VBF_shapes.root ({e}), defaulting to 9")
        return 9


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a Run3 Combine datacard for VBF H→invisible."
    )
    parser.add_argument("--cat", required=True, choices=["VTR", "MTR"],
                        help="Analysis category")
    parser.add_argument("--campaigns", nargs="+", required=True,
                        help="One or more pyRAT campaign names, e.g. "
                             "Run3Summer22 Run3Summer22EE Run3Summer23 Run3Summer23BPix")
    parser.add_argument("--var", default="Mjj", choices=["Mjj", "SignalScore"],
                        help="Fit variable")
    parser.add_argument("--label", default=None,
                        help="Datacard directory label (default: {cat}_Run3)")
    args = parser.parse_args()

    # Derive combined year string (matches what the workspace macros use)
    if len(args.campaigns) == 1:
        args.year = args.campaigns[0]
    else:
        args.year = f"{args.campaigns[0]}_to_{args.campaigns[-1]}"

    label = args.label or f"{args.cat}_{args.year}"

    # Lumi from CampaignInfoMap
    cmap = CampaignInfoMap()
    total_lumi = sum(
        cmap.get_lumi(cmap.get_year(c)) for c in args.campaigns
    )
    total_lumi = round(total_lumi, 2)
    energy = cmap.get_energy(cmap.get_year(args.campaigns[0]))

    print(f"\n{'='*60}")
    print(f"  Category  : {args.cat}")
    print(f"  Year str  : {args.year}")
    print(f"  Variable  : {args.var}")
    print(f"  Label     : {label}")
    print(f"  Lumi      : {total_lumi} fb-1  ({energy} TeV)")
    print(f"{'='*60}\n")

    # Step 1: copy workspace files
    print("Copying workspace files...")
    ws_dir, issues = copy_workspace_files(args, label)

    if issues:
        print("\nISSUES — the following files were not found:")
        for msg in issues:
            print(f"  {msg}")
        print("\nMake sure you have run:")
        print(f"  makeSignalAndMCBackgroundWSRun3.C(\"{args.year}\", \"{args.cat}\", false, true, true)")
        print(f"  makeWS_percategoryRun3.C(\"{args.year}\", \"{args.cat}\", false, true, true)")
        print(f"  and that plotting has been run to produce the qcdDD/noiseDD workspaces\n")

    # Step 2: get number of bins from VBF_shapes.root
    nB = get_nB(args.cat, args.year, args.var)
    print(f"Number of bins: {nB}")

    # Step 3: write datacard
    print("\nWriting datacard...")
    dc_path = write_datacard(args, label, ws_dir, total_lumi, energy, nB)

    print(f"\nDone. Datacard: {dc_path}")
    print("\nTo test with Combine (from CMSSW):")
    print(f"  combine -M AsymptoticLimits -d {dc_path} -m 125 --run blind -t -1")


if __name__ == "__main__":
    main()
