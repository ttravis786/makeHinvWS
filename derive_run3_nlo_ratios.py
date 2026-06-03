#!/usr/bin/env python3
"""
derive_run3_nlo_ratios.py

Reads the Run3 SR VBF_shapes.root and extracts NLO yield-variation histograms for
QCD W and Z processes, saving them to a ROOT file for use by makeWS_percategoryRun3.C.

Run2 analogy (makeWS_percategory.C histosNLO indices):
  [0]  fnlo_SF_mufUp         -> QCD_Wjet_NLONew/QCD_Znunu_NLONew weight_muf_up
  [1]  fnlo_SF_murUp         -> weight_mur_up
  [2]  fnlo_SF_pdfUp         -> weight_pdf_up  (flat = nominal in current pyRAT)
  [3]  fnlo_SF_EWK_corrUp    -> weight_ewk_correction_up

The C++ reads these yield histograms and divides by nominal just like Run2:
  - muF, muR : WZratioSyst = W_syst_up / W_nominal  (W-yield ratio only, not W/Z)
  - pdf      : WZratioSyst = (W_nom/Z_nom) / (W_pdf/Z_pdf)  (W/Z ratio variation)
  - ewk_corr : WZratioSyst = (W_nom/Z_nom) / (W_ewk/Z_ewk)  (W/Z ratio variation)

Usage:
  python3 derive_run3_nlo_ratios.py --year Run3Summer22_to_Run3Summer23BPix --cat VTR --var Mjj
  python3 derive_run3_nlo_ratios.py --year Run3Summer22_to_Run3Summer23BPix --cat VTR --var SignalScore
"""

import argparse
import os
import sys
import numpy as np
import uproot
import boost_histogram as bh

ROOT_FILES_BASE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "hinvisible_mtr_vtr", "root_files"
)

# Processes for which NLO histograms are saved (Run2 histosNLO loop covers these four,
# but the TF formulas only use QCDW and QCDZnunu indices)
NLO_PROCS = [
    "QCD_Wjet_NLONew",   # PROCESS::QCDW    = 5
    "QCD_Znunu_NLONew",  # PROCESS::QCDZnunu = 3
]

# [index]: (Run3 syst name, output histogram suffix)
# Order must match the nNLO=4 array in makeWS_percategoryRun3.C
NLO_SYSTS = [
    ("muf",            "fnlo_SF_mufUp"),         # [0]
    ("mur",            "fnlo_SF_murUp"),         # [1]
    ("pdf",            "fnlo_SF_pdfUp"),         # [2] — currently flat in Run3
    ("ewk_correction", "fnlo_SF_EWK_corrUp"),   # [3]
]


def safe_ratio(num: np.ndarray, den: np.ndarray, fallback: float = 1.0) -> np.ndarray:
    return np.where(np.abs(den) > 1e-9, num / den, fallback)


def derive(year: str, cat: str, var: str, output_dir: str) -> None:
    sr_path = os.path.join(ROOT_FILES_BASE, year, f"SR_{cat}_VBF", var, "VBF_shapes.root")
    if not os.path.exists(sr_path):
        sys.exit(f"ERROR: SR file not found:\n  {sr_path}")

    out_path = os.path.join(output_dir, f"run3_nlo_sf_{year}_{cat}_{var}.root")
    print(f"Input :  {sr_path}")
    print(f"Output:  {out_path}\n")

    out_hists: dict[str, tuple[np.ndarray, np.ndarray]] = {}  # name -> (values, edges)

    with uproot.open(sr_path) as f:
        all_keys = {k.split(";")[0] for k in f.keys()}

        for proc in NLO_PROCS:
            nom_key = f"{proc}_calc_weight"
            if nom_key not in all_keys:
                sys.exit(f"ERROR: nominal histogram {nom_key!r} not found")

            nom_h = f[nom_key]
            nom_vals, edges = [np.array(x) for x in nom_h.to_numpy()]
            # to_numpy() returns (values, edges); values has nBins entries, edges has nBins+1

            print(f"  {proc}  nominal = {np.round(nom_vals, 2)}")

            for syst_name, out_suffix in NLO_SYSTS:
                up_key = f"{proc}_weight_{syst_name}_up"
                if up_key in all_keys:
                    up_vals: np.ndarray = np.array(f[up_key].to_numpy()[0])
                else:
                    print(f"    WARNING: {up_key!r} not found — using nominal")
                    up_vals = nom_vals.copy()

                ratio = safe_ratio(up_vals, nom_vals)
                flat = bool(np.allclose(ratio, 1.0, atol=1e-4))
                note = "  [flat — pdf not yet in pyRAT weights]" if flat else ""
                print(f"    {syst_name:22s} ratio = {np.round(ratio, 4)}{note}")

                out_name = f"{proc}_{out_suffix}"
                out_hists[out_name] = (up_vals, edges)  # yield histogram; C++ divides by nominal

            print()

    # Write ROOT file — pass (values, edges) tuple; uproot writes as TH1F
    with uproot.recreate(out_path) as fout:
        for name, (vals, edges) in out_hists.items():
            h = bh.Histogram(bh.axis.Variable(edges), storage=bh.storage.Weight())
            h.view().value[:] = vals
            h.view().variance[:] = 0.0
            fout[name] = h

    print(f"Written {len(out_hists)} histograms.")

    # Print summary of W/Z ratio effects for reference
    print("\n=== W/Z ratio systematic coefficients (per bin, for reference) ===")
    print(f"{'syst':25s}  {'W/Z ratio effect (syst/nom)':>40s}")

    with uproot.open(sr_path) as f:
        w_nom = np.array(f["QCD_Wjet_NLONew_calc_weight"].to_numpy()[0])
        z_nom = np.array(f["QCD_Znunu_NLONew_calc_weight"].to_numpy()[0])
        wz_nom = safe_ratio(w_nom, z_nom)
        all_keys_summary = {k.split(";")[0] for k in f.keys()}

        for syst_name, _ in NLO_SYSTS:
            w_up_key = f"QCD_Wjet_NLONew_weight_{syst_name}_up"
            z_up_key = f"QCD_Znunu_NLONew_weight_{syst_name}_up"
            w_up = np.array(f[w_up_key].to_numpy()[0]) if w_up_key in all_keys_summary else w_nom
            z_up = np.array(f[z_up_key].to_numpy()[0]) if z_up_key in all_keys_summary else z_nom
            wz_up = safe_ratio(w_up, z_up)
            # W/Z ratio change = wz_nom / wz_up  (matches Run2 ewk/pdf treatment)
            wz_ratio_effect = safe_ratio(wz_nom, wz_up)
            print(f"  {syst_name:22s}  {np.round(wz_ratio_effect, 4)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Derive Run3 NLO correction ratio histograms for makeWS_percategoryRun3.C"
    )
    parser.add_argument(
        "--year", default="Run3Summer22_to_Run3Summer23BPix",
        help="Campaign string (directory name under root_files/)"
    )
    parser.add_argument("--cat", default="VTR", choices=["VTR", "MTR"])
    parser.add_argument("--var", default="Mjj", choices=["Mjj", "SignalScore"])
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory in which to write run3_nlo_sf_*.root (default: current dir)"
    )
    args = parser.parse_args()
    derive(args.year, args.cat, args.var, args.output_dir)
