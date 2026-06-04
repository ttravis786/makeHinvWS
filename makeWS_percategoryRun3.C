#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "TSystem.h"
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "TH1F.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "../HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "RooAddition.h"
#include "RooMsgService.h"

// Run3 parametric workspace builder, adapted from the working Run2 makeWS_percategory.C.
//
// Key differences from Run2:
//   - Flat ROOT file layout: nominal key "{proc}_calc_weight" (or "data_obs" for data),
//     systematic key "{proc}_weight_{syst}_{up/down}".  No per-region subdirectories.
//   - Run3 process and region names.
//   - NLO corrections for WZ ratio TF come from weight_ewk_correction /
//     weight_qcd_correction / weight_pdf histograms already in the histos array —
//     no separate histosNLO array needed.
//   - turn_off_jecs / turn_off_btag flags: when true the corresponding terms are
//     omitted from the TF formula (preliminary fit mode).
//   - Mjj or SignalScore fit variable.

// ---- process indices — order must match lProcs[] below ----------------------
enum PROCESS {
    data      = 0,
    VBFH      = 1,   // vbfH
    GGH       = 2,   // ggH
    QCDZnunu  = 3,   // QCD_Znunu_NLONew
    EWKZnunu  = 4,   // EWK_Znunu
    QCDW      = 5,   // QCD_Wjet_NLONew
    EWKW      = 6,   // EWK_W
    QCDDYll   = 7,   // QCD_Zll_NLONew
    EWKZll    = 8,   // EWK_Zll
    WH        = 9,
    ZH        = 10,
};

// ---- small helpers ----------------------------------------------------------

double safeDiv(double num, double den, double fallback = 1.0) {
    return (std::fabs(den) > 1e-9) ? num / den : fallback;
}

double getBin(TH1F *h, int b, double fallback = 0.0) {
    return h ? h->GetBinContent(b) : fallback;
}

double getBinErr(TH1F *h, int b, double fallback = 0.0) {
    return h ? h->GetBinError(b) : fallback;
}

// ---- main function ----------------------------------------------------------

int makeWS_percategoryRun3(
    std::string year       = "Run3Summer22_to_Run3Summer23BPix",
    std::string cat        = "VTR",
    bool classifier        = false,   // false=Mjj  true=SignalScore
    bool turn_off_jecs     = true,    // omit JEC/JER terms from TF formula
    bool turn_off_btag     = true     // omit b-tag terms (not yet in Run3 files)
) {
    gSystem->Load("libHiggsAnalysisCombinedLimit.so");

    // Suppress RooFit INFO spam from RecycleConflictNodes imports
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // ---- variable and output naming -----------------------------------------
    const std::string lVarDir   = classifier ? "SignalScore" : "Mjj";
    const std::string lChannel  = "VBF";
    const std::string lCategory = cat + "_";
    const std::string lYear     = year + "_";
    const std::string lOutFileName = "param_ws_" + lYear + lCategory + lChannel + lVarDir + ".root";

    std::string varName, varTitle;
    double xmin, xmax;
    if (!classifier) {
        varName = "mjj_" + cat; varTitle = "M_{jj} (GeV)"; xmin = 200; xmax = 5000;
    } else {
        varName = "bdt_" + cat; varTitle = "Signal Score"; xmin = 0; xmax = 1;
    }
    RooRealVar lVarFit(varName.c_str(), varTitle.c_str(), xmin, xmax);

    // ---- regions ------------------------------------------------------------
    // lRegions: actual directory names used in the ROOT file paths
    // lRegionAlias: short labels used for workspace objects (datacard compatibility)
    const unsigned nR = 5;
    std::string lRegions[nR]     = {"SR", "singleElectron_CR", "singleMuon_CR", "diElectron_CR", "diMuon_CR"};
    std::string lRegionAlias[nR] = {"SR", "Wenu",              "Wmunu",          "Zee",           "Zmumu"};
    // iR < 3 are W-like CRs, iR >= 3 are Z-like CRs (same convention as Run2)

    // ---- processes ----------------------------------------------------------
    const unsigned nP = 11;
    // Run3 ROOT file keys (prefix before _calc_weight)
    std::string lProcs[nP] = {
        "data_obs",        // data  = 0
        "vbfH",            // VBFH  = 1
        "ggH",             // GGH   = 2
        "QCD_Znunu_NLONew",// QCDZnunu = 3
        "EWK_Znunu",       // EWKZnunu = 4
        "QCD_Wjet_NLONew", // QCDW  = 5
        "EWK_W",           // EWKW  = 6
        "QCD_Zll_NLONew",  // QCDDYll = 7
        "EWK_Zll",         // EWKZll  = 8
        "WH",              // WH    = 9
        "ZH",              // ZH    = 10
    };

    // ---- weight systematics -------------------------------------------------
    // Same 11 as in makeSignalAndMCBackgroundWSRun3.C (must be consistent).
    const unsigned nN = 14;
    std::string lNuisRun3[nN] = {
        "Trigger",
        "tau_veto",
        "electron_id",
        "electron_id_high_pt",
        "electron_veto",
        "muon_id",
        "muon_veto_id",
        "muon_veto_iso",
        "pdf",
        "ewk_correction",
        "qcd_correction",
        "mur",
        "muf",
        "pileup"
    };
    std::string lNuisCMS[nN] = {
        "CMS_Trigger",
        "CMS_eff_tauveto",
        "CMS_eff_eTight_idiso",
        "CMS_eff_eTight_idiso_highpt",
        "CMS_eff_eVeto_idiso_veto",
        "CMS_eff_muTight_id",
        "CMS_eff_muLoose_id_veto",
        "CMS_eff_muLoose_iso_veto",
        "CMS_pdf",
        "CMS_ewk_correction",
        "CMS_qcd_correction",
        "CMS_scale_mur",
        "CMS_scale_muf",
        "CMS_pileup"
    };
    const bool corrCat[nN]  = {false, true, true, true, true, true, true, true, true, true, true, true, true, true};
    const bool corrYear[nN] = {false, false, false, false, false, true, true, true, true, true, true, true, true, true};

    // Apply cat/year decorations to CMS names
    for (unsigned iN = 0; iN < nN; ++iN) {
        if (!corrCat[iN])  lNuisCMS[iN] += "_" + cat;
        if (!corrYear[iN]) lNuisCMS[iN] += "_" + year;
    }

    // Build syst suffix arrays: lSysts[0]="" (nominal), then up/down pairs
    // lSysts[2*iN+1] = nuisRun3[iN]+"_up", lSysts[2*iN+2] = nuisRun3[iN]+"_down"
    const unsigned nS = 2*nN + 1;
    std::string lSysts[nS];
    std::string lSystsCMS[nS];
    lSysts[0] = ""; lSystsCMS[0] = "";
    for (unsigned iN = 0; iN < nN; ++iN) {
        lSysts[2*iN+1]   = lNuisRun3[iN] + "_up";
        lSysts[2*iN+2]   = lNuisRun3[iN] + "_down";
        lSystsCMS[2*iN+1] = lNuisCMS[iN] + "Up";
        lSystsCMS[2*iN+2] = lNuisCMS[iN] + "Down";
    }

    // NLO correction histograms — yield variations read from derive_run3_nlo_ratios.py output.
    // Index order matches Run2 histosNLO: [0]=muF_up [1]=muR_up [2]=pdf_up [3]=ewk_correction_up.
    const unsigned nNLO = 4;
    const std::string lNLOname[nNLO] = {
        "fnlo_SF_mufUp",
        "fnlo_SF_murUp",
        "fnlo_SF_pdfUp",       // flat (=nominal) in current Run3 pyRAT weights
        "fnlo_SF_EWK_corrUp",
    };
    TH1F *histosNLO[nP][nNLO];
    for (unsigned iP = 0; iP < nP; ++iP)
        for (unsigned iS = 0; iS < nNLO; ++iS)
            histosNLO[iP][iS] = nullptr;

    // ---- vproc: which process index to use per type per region for TF -------
    // [QCD=0][iR]: process in region iR used for QCD W/Z ratio
    // [EWK=1][iR]: process in region iR used for EWK W/Z ratio
    // iR=0: SR Z (QCDZnunu / EWKZnunu)
    // iR=1,2: W CRs (QCDW / EWKW)
    // iR=3,4: Z CRs (QCDDYll / EWKZll)
    const unsigned nT = 2;
    std::string lType[nT]   = {"QCD", "EWK"};
    std::string lTypeLC[nT] = {"qcd", "ewk"};
    unsigned vproc[nT][nR] = {
        {PROCESS::QCDZnunu, PROCESS::QCDW,    PROCESS::QCDW,    PROCESS::QCDDYll, PROCESS::QCDDYll},
        {PROCESS::EWKZnunu, PROCESS::EWKW,    PROCESS::EWKW,    PROCESS::EWKZll,  PROCESS::EWKZll }
    };

    // ---- input file paths ---------------------------------------------------
    const std::string run3_path =
        "/vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT/hinvisible_mtr_vtr/root_files/";
    std::string lInFileName[nR];
    for (unsigned iR = 0; iR < nR; ++iR) {
        lInFileName[iR] = run3_path + year + "/" + lRegions[iR] + "_" + cat + "_VBF/" +
                          lVarDir + "/VBF_shapes.root";
        std::cout << "Input [" << iR << "]: " << lInFileName[iR] << std::endl;
    }

    // ---- JES/JER (same scheme as Run2, used only when turn_off_jecs=false) --
    const unsigned nJ = 11;
    std::string lJes[nJ] = {
        "jesAbsolute",
        Form("jesAbsolute_%s",  year.c_str()),
        "jesBBEC1",
        Form("jesBBEC1_%s",     year.c_str()),
        "jesEC2",
        Form("jesEC2_%s",       year.c_str()),
        "jesFlavorQCD",
        "jesHF",
        Form("jesHF_%s",        year.c_str()),
        "jesRelativeBal",
        Form("jesRelativeSample_%s", year.c_str())
    };

    // JEC/JER nuisance parameters are always created so the workspace
    // always contains them and the datacard can always reference them.
    // When turn_off_jecs=true the TF formula coefficients are set to 1.0
    // (flat dummy — no shape effect), matching makeSignalAndMCBackgroundWSRun3.C.
    TFile *finputJES = nullptr;
    if (!turn_off_jecs) {
        finputJES = TFile::Open("vbf_jes_jer_tf_uncs.root");
        if (!finputJES) {
            std::cout << "ERROR: vbf_jes_jer_tf_uncs.root not found. "
                      << "Re-run with turn_off_jecs=true to use flat dummies." << std::endl;
            return 1;
        }
    }
    RooRealVar *jer = new RooRealVar(("CMS_res_j_" + year).c_str(), "JER nuisance", 0);
    RooRealVar *jes[nJ];
    for (unsigned iJ = 0; iJ < nJ; ++iJ)
        jes[iJ] = new RooRealVar(Form("CMS_scale_j_%s", lJes[iJ].c_str()), "JES nuisance", 0);

    // ---- load all histograms ------------------------------------------------
    TH1F *histos[nR][nP][nS];
    for (unsigned iR = 0; iR < nR; ++iR)
        for (unsigned iP = 0; iP < nP; ++iP)
            for (unsigned iS = 0; iS < nS; ++iS)
                histos[iR][iP][iS] = nullptr;

    for (unsigned iR = 0; iR < nR; ++iR) {
        TFile *f = TFile::Open(lInFileName[iR].c_str());
        if (!f) { std::cout << "ERROR: could not open " << lInFileName[iR] << std::endl; return 1; }

        for (unsigned iP = 0; iP < nP; ++iP) {
            // Nominal histogram key
            std::string nomKey = (iP == PROCESS::data)
                ? "data_obs"
                : (lProcs[iP] + "_calc_weight");

            histos[iR][iP][0] = (TH1F *)f->Get(nomKey.c_str());
            if (!histos[iR][iP][0]) {
                if (iP == PROCESS::data) {
                    std::cout << "ERROR: data_obs not found in " << lInFileName[iR] << std::endl;
                    return 1;
                }
                // Process not present in this region (e.g. signal in CR) — leave null
                continue;
            }
            histos[iR][iP][0]->SetDirectory(0);

            // Systematic histograms
            for (unsigned iN = 0; iN < nN; ++iN) {
                for (unsigned dir = 0; dir < 2; ++dir) {
                    unsigned iS = 2*iN + 1 + dir;
                    std::string key = lProcs[iP] + "_weight_" + lNuisRun3[iN] +
                                      (dir == 0 ? "_up" : "_down");
                    TH1F *h = (TH1F *)f->Get(key.c_str());
                    if (h) { h->SetDirectory(0); histos[iR][iP][iS] = h; }
                    else   { histos[iR][iP][iS] = histos[iR][iP][0]; } // fall back to nominal
                }
            }
        }
        f->Close();

        // ---- diagnostic: print integral and negative-bin count per process ----
        std::cout << "\n[DIAG] Region " << lRegions[iR] << " (" << lRegionAlias[iR] << "):" << std::endl;
        for (unsigned iP = 0; iP < nP; ++iP) {
            TH1F *h = histos[iR][iP][0];
            if (!h) continue;
            double integral = h->Integral();
            int nNegBins = 0;
            double mostNeg = 0.0;
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                if (h->GetBinContent(b) < 0) {
                    ++nNegBins;
                    if (h->GetBinContent(b) < mostNeg) mostNeg = h->GetBinContent(b);
                }
            }
            std::string tag = (nNegBins > 0) ? " <<< HAS NEGATIVE BINS" : "";
            std::cout << "  [" << iP << "] " << (iP == PROCESS::data ? "data_obs" : lProcs[iP])
                      << "  integral=" << integral
                      << "  negBins=" << nNegBins
                      << "  mostNeg=" << mostNeg
                      << tag << std::endl;
        }
        std::cout << std::endl;
    }

    // ---- NLO correction file (generated by derive_run3_nlo_ratios.py) ----------
    // Contains yield variation histograms for QCD W and Z processes.
    // Must be regenerated if the shapes ROOT files change.
    {
        const std::string nloFileName = "run3_nlo_sf_" + year + "_" + cat + "_" + lVarDir + ".root";
        TFile *fNLO = TFile::Open(nloFileName.c_str());
        if (!fNLO) {
            std::cout << "ERROR: NLO file not found: " << nloFileName << "\n"
                      << "  Generate it with:\n"
                      << "    python3 derive_run3_nlo_ratios.py"
                      << " --year " << year << " --cat " << cat << " --var " << lVarDir << std::endl;
            return 1;
        }
        for (unsigned iP : {(unsigned)PROCESS::QCDZnunu, (unsigned)PROCESS::QCDW}) {
            for (unsigned iS = 0; iS < nNLO; ++iS) {
                std::string key = lProcs[iP] + "_" + lNLOname[iS];
                histosNLO[iP][iS] = (TH1F *)fNLO->Get(key.c_str());
                if (!histosNLO[iP][iS]) {
                    std::cout << "WARNING: NLO histo not found: " << key << " — using nominal" << std::endl;
                    histosNLO[iP][iS] = histos[0][iP][0];
                } else {
                    histosNLO[iP][iS]->SetDirectory(0);
                    std::cout << " NLO histo loaded: " << key
                              << "  integral=" << histosNLO[iP][iS]->Integral() << std::endl;
                }
            }
        }
        fNLO->Close();
    }

    // ---- binning from data_obs in SR ----------------------------------------
    TH1F *hData = histos[0][PROCESS::data][0];
    const unsigned nB = (unsigned)hData->GetNbinsX();
    double bins[nB+1];
    for (unsigned iB = 0; iB < nB; ++iB)
        bins[iB] = hData->GetXaxis()->GetBinLowEdge(iB+1);
    bins[nB] = hData->GetXaxis()->GetBinLowEdge(nB+1);
    TH1F dummyHist("dummyHist", "Dummy hist for binning", nB, bins);

    std::cout << nB << " bins: ";
    for (unsigned iB = 0; iB <= nB; ++iB) std::cout << bins[iB] << " ";
    std::cout << std::endl;

    // ---- output workspace ---------------------------------------------------
    TFile *fOut = new TFile(lOutFileName.c_str(), "RECREATE");
    RooWorkspace wspace("wspace", "wspace");
    RooArgList vars(lVarFit);

    // ---- import data_obs and histogram-based shapes -------------------------
    for (unsigned iR = 0; iR < nR; ++iR) {
        // data_obs
        RooDataHist *dh = new RooDataHist(("data_obs_" + lRegionAlias[iR]).c_str(),
                                          "Data observed", vars, histos[iR][PROCESS::data][0]);
        wspace.import(*dh);

        // All non-data processes: nominal + systematic shapes imported as RooDataHist.
        // These are used by the datacard for small MC backgrounds (Top, Diboson, etc.)
        // and also provide the histogram values used in TF formulas below.
        for (unsigned iP = 1; iP < nP; ++iP) {
            if (!histos[iR][iP][0]) continue;
            for (unsigned iS = 0; iS < nS; ++iS) {
                if (!histos[iR][iP][iS]) continue;
                std::ostringstream label;
                label << lProcs[iP] << "_hist_" << lRegionAlias[iR];
                if (iS > 0) label << "_" << lSystsCMS[iS];
                // Import regardless of integral so the datacard can always
                // reference these objects. Zero-yield processes are valid in
                // Combine (they simply contribute nothing to the likelihood).
                RooDataHist *bh = new RooDataHist(label.str().c_str(), "Background",
                                                   vars, histos[iR][iP][iS]);
                wspace.import(*bh);
            }
        }
    }

    // ---- nuisance parameters for TF formulas --------------------------------
    // Per-bin stat nuisances for EWK/QCD Z ratio in SR
    RooRealVar *ewkqcdratiostat[nB+1];
    // Per-bin W/Z ratio stat nuisances
    RooRealVar *wzratiostat[nT][nB+1];
    // Per-bin CR/SR stat nuisances
    RooRealVar *TFstat[nT][nB+1][nR-1];
    // Lepton SF nuisances for CR TF (one per systematic, shared across bins)
    RooRealVar *TFsysts[nN];

    // Initialise hardcoded nuisance overrides (< 0 means use histogram-derived value)
    double hardCodeNuisance[nR][nS];
    for (unsigned iR = 0; iR < nR; ++iR)
        for (unsigned iS = 0; iS < nS; ++iS)
            hardCodeNuisance[iR][iS] = -1.;

    // Create lepton SF nuisance RooRealVars
    for (unsigned iN = 0; iN < nN; ++iN) {
        std::ostringstream lname;
        lname << lNuisCMS[iN]; // already decorated with cat/year above
        TFsysts[iN] = new RooRealVar(lname.str().c_str(), "CR/SR syst nuisance", 0);
    }

    // NLO theory nuisances for W/Z ratio TF — same 4 as Run2:
    //   mur, muf: W-yield-only ratio, correlated across bins (not W/Z ratio, matching Run2)
    //   pdf:      W/Z ratio variation, correlated across bins (flat in current Run3)
    //   ewk_corr: W/Z ratio variation, PER-BIN (matches Run2's wzratioEWK_on_strong[nT][nB])
    RooRealVar *wzratioMURSyst[nT];
    RooRealVar *wzratioMUFSyst[nT];
    RooRealVar *wzratioPDFSyst[nT];
    RooRealVar *wzratioEWK_on_strong[nT][nB+1]; // per-bin, like Run2

    for (unsigned iT = 0; iT < nT; ++iT) {
        std::ostringstream n;
        n.str(""); n << lCategory << lType[iT] << "wzratioMURSyst";
        wzratioMURSyst[iT] = new RooRealVar(n.str().c_str(), "W/Z muR nuisance", 0);
        n.str(""); n << lCategory << lType[iT] << "wzratioMUFSyst";
        wzratioMUFSyst[iT] = new RooRealVar(n.str().c_str(), "W/Z muF nuisance", 0);
        n.str(""); n << lCategory << lType[iT] << "wzratioPDFSyst";
        wzratioPDFSyst[iT] = new RooRealVar(n.str().c_str(), "W/Z PDF nuisance", 0);
    }

    // ---- per-bin TF loop ----------------------------------------------------
    for (unsigned iB = 1; iB <= nB; ++iB) {
        std::cout << " -- bin " << iB << std::endl;
        std::ostringstream lname;

        RooFormulaVar *EWKQCDbin = nullptr;

        for (unsigned iT = 0; iT < nT; ++iT) {

            // ------ Z yield in SR (free parameter) ---------------------------
            unsigned zprocSR = (iT == 0) ? PROCESS::QCDZnunu : PROCESS::EWKZnunu;
            double zSRnom = getBin(histos[0][zprocSR][0], iB);

            lname.str(""); lname << lCategory << lType[iT] << "Z_SR_bin" << iB;
            RooRealVar binParZ(lname.str().c_str(),
                               (lType[iT] + " Z+jets yield in SR per bin").c_str(),
                               zSRnom, 0, 10.*zSRnom);
            wspace.import(binParZ, RooFit::RecycleConflictNodes());

            // ------ EWK/QCD Z link in SR (iT==0 only) -----------------------
            if (iT == 0) {
                double zEWKnom = std::max(0.0, getBin(histos[0][PROCESS::EWKZnunu][0], iB));
                double zQCDnom = std::max(0.0, getBin(histos[0][PROCESS::QCDZnunu][0], iB));
                double ewkqcdratio = safeDiv(zEWKnom, zQCDnom);
                double ewkqcdstat  = 1. + std::sqrt(
                    std::pow(safeDiv(getBinErr(histos[0][PROCESS::EWKZnunu][0], iB), zEWKnom), 2) +
                    std::pow(safeDiv(getBinErr(histos[0][PROCESS::QCDZnunu][0], iB), zQCDnom), 2));

                lname.str(""); lname << lCategory << "ewkqcdratio_stat_bin" << iB;
                ewkqcdratiostat[iB] = new RooRealVar(lname.str().c_str(), "EWK/QCD ratio stat", 0);

                lname.str(""); lname << lCategory << "TF_EWKQCDSR_bin" << iB;
                std::ostringstream fEWKQCD;
                fEWKQCD << ewkqcdratio << "*TMath::Power(" << ewkqcdstat << ",@0)";
                RooFormulaVar TFEWKQCD(lname.str().c_str(), "TF EWK/QCD Z",
                                       fEWKQCD.str().c_str(), RooArgList(*ewkqcdratiostat[iB]));
                wspace.import(TFEWKQCD, RooFit::RecycleConflictNodes());

                lname.str(""); lname << lCategory << "EWKQCD_SR_bin" << iB;
                EWKQCDbin = new RooFormulaVar(lname.str().c_str(),
                                              "EWK Z yield from QCD Z per bin",
                                              "@0*@1", RooArgList(TFEWKQCD, binParZ));
                wspace.import(*EWKQCDbin, RooFit::RecycleConflictNodes());
            }

            // ------ W/Z ratio in SR ------------------------------------------
            unsigned wprocSR = (iT == 0) ? PROCESS::QCDW : PROCESS::EWKW;

            // Clamp to zero: NLO samples can have negative total MC yield in
            // tight selections due to negative-weight events. A negative yield
            // is unphysical for the TF ratio — treat it as zero contribution.
            double wSRnom  = std::max(0.0, getBin(histos[0][wprocSR][0], iB));
            double zSRnom2 = std::max(0.0, getBin(histos[0][zprocSR][0], iB));
            double wzratio = safeDiv(wSRnom, zSRnom2);

            double wzratiostat_val = 1. + std::sqrt(
                std::pow(safeDiv(getBinErr(histos[0][wprocSR][0], iB), wSRnom), 2) +
                std::pow(safeDiv(getBinErr(histos[0][zprocSR][0], iB), zSRnom2), 2));

            lname.str(""); lname << lCategory << lType[iT] << "wzratio_stat_bin" << iB;
            wzratiostat[iT][iB] = new RooRealVar(lname.str().c_str(), "W/Z ratio stat", 0);

            // NLO scale and EWK corrections on the W/Z ratio — Run2 convention:
            //   muF, muR : W-yield ratio only  (histosNLO[QCDW][0/1] / W_nominal)
            //   pdf      : W/Z ratio variation  (W_nom/Z_nom) / (W_pdf/Z_pdf)
            //   ewk_corr : W/Z ratio variation  (W_nom/Z_nom) / (W_ewk/Z_ewk)
            // histosNLO always uses QCDW/QCDZnunu regardless of iT (same as Run2).
            double wNomNLO = getBin(histos[0][PROCESS::QCDW][0],    iB, 1.0);
            double zNomNLO = getBin(histos[0][PROCESS::QCDZnunu][0], iB, 1.0);

            double WZratioSyst_muf = safeDiv(
                getBin(histosNLO[PROCESS::QCDW][0], iB, wNomNLO), wNomNLO);
            double WZratioSyst_mur = safeDiv(
                getBin(histosNLO[PROCESS::QCDW][1], iB, wNomNLO), wNomNLO);
            double WZratioSyst_pdf = safeDiv(
                safeDiv(wNomNLO, zNomNLO),
                safeDiv(getBin(histosNLO[PROCESS::QCDW][2],    iB, wNomNLO),
                        getBin(histosNLO[PROCESS::QCDZnunu][2], iB, zNomNLO)));
            double WZratioSyst_ewk = safeDiv(
                safeDiv(wNomNLO, zNomNLO),
                safeDiv(getBin(histosNLO[PROCESS::QCDW][3],    iB, wNomNLO),
                        getBin(histosNLO[PROCESS::QCDZnunu][3], iB, zNomNLO)));

            // Per-bin EWK correction nuisance — matches Run2's wzratioEWK_on_strong[nT][nB]
            lname.str(""); lname << lCategory << lType[iT] << "wzratio_EWK_corr_on_Strong_bin" << iB;
            wzratioEWK_on_strong[iT][iB] = new RooRealVar(
                lname.str().c_str(), "W/Z EWK correction on strong proc per-bin nuisance", 0);

            std::cout << " [" << lType[iT] << " WZ SR bin " << iB << "]"
                      << " ratio=" << wzratio
                      << " muf=" << WZratioSyst_muf << " mur=" << WZratioSyst_mur
                      << " pdf=" << WZratioSyst_pdf << " ewk=" << WZratioSyst_ewk << std::endl;

            lname.str(""); lname << lCategory << lType[iT] << "TF_WZSR_bin" << iB;
            std::ostringstream fWZ;
            int atIdx = 0;
            fWZ << wzratio;
            fWZ << "*TMath::Power(" << WZratioSyst_mur << ",@" << atIdx++ << ")"; // muR  @0
            fWZ << "*TMath::Power(" << WZratioSyst_muf << ",@" << atIdx++ << ")"; // muF  @1
            fWZ << "*TMath::Power(" << WZratioSyst_pdf << ",@" << atIdx++ << ")"; // pdf  @2
            fWZ << "*TMath::Power(" << wzratiostat_val << ",@" << atIdx++ << ")"; // stat @3
            fWZ << "*TMath::Power(" << WZratioSyst_ewk << ",@" << atIdx++ << ")"; // ewk  @4 (per-bin)

            RooArgList wzVars(*wzratioMURSyst[iT], *wzratioMUFSyst[iT],
                              *wzratioPDFSyst[iT],
                              *wzratiostat[iT][iB],
                              *wzratioEWK_on_strong[iT][iB]);

            {
                // JER on W/Z ratio — flat dummy (1.0) when turn_off_jecs=true
                double jerWZ = 1.0;
                if (!turn_off_jecs && finputJES) {
                    TH1D *h = (TH1D *)finputJES->Get(
                        Form("znunu_over_wlnu_%s_%s_jerUp", lTypeLC[iT].c_str(), year.c_str()));
                    if (h) jerWZ = safeDiv(1., h->GetBinContent(1));
                }
                fWZ << "*TMath::Power(" << jerWZ << ",@" << atIdx++ << ")";
                wzVars.add(*jer);
                // JES on W/Z ratio — flat dummy (1.0) when turn_off_jecs=true
                for (unsigned iJ = 0; iJ < nJ; ++iJ) {
                    double jesWZ = 1.0;
                    if (!turn_off_jecs && finputJES) {
                        TH1D *h = (TH1D *)finputJES->Get(
                            Form("znunu_over_wlnu_%s_%s_%sUp",
                                 lTypeLC[iT].c_str(), year.c_str(), lJes[iJ].c_str()));
                        if (h) jesWZ = safeDiv(1., h->GetBinContent(1));
                    }
                    fWZ << "*TMath::Power(" << jesWZ << ",@" << atIdx++ << ")";
                    wzVars.add(*jes[iJ]);
                }
            }

            // Lepton and tau veto uncertainties on the W/Z ratio.
            // The SR applies tau, electron, and muon vetoes; these affect W and Z
            // differently since W+jets has a real lepton from the W decay.
            // lNuisRun3 indices: tau_veto=1, electron_veto=4, muon_veto_id=6, muon_veto_iso=7.
            // The same TFsysts[iN] parameters are shared with the CR TF (correlated).
            // Coefficient per bin = (W/Z nominal) / (W/Z with veto_syst_up).
            // Currently flat in Run3 ROOT files; skip if negligible until pyRAT weights land.
            {
                const unsigned vetoIdx[] = {1, 4, 6, 7};
                for (unsigned iN : vetoIdx) {
                    unsigned iSup = 2*iN + 1;
                    double wVeto = getBin(histos[0][PROCESS::QCDW][iSup],     iB, wSRnom);
                    double zVeto = getBin(histos[0][PROCESS::QCDZnunu][iSup], iB, zSRnom2);
                    double coeff = safeDiv(wzratio, safeDiv(wVeto, zVeto));
                    if (std::fabs(coeff - 1.) < 1e-4) continue;
                    fWZ << "*TMath::Power(" << coeff << ",@" << atIdx++ << ")";
                    wzVars.add(*TFsysts[iN]);
                }
            }

            RooFormulaVar *TFWZ = new RooFormulaVar(lname.str().c_str(), "TF W/Z SR",
                                                     fWZ.str().c_str(), wzVars);
            wspace.import(*TFWZ, RooFit::RecycleConflictNodes());

            lname.str(""); lname << lCategory << lType[iT] << "WZ_SR_bin" << iB;
            RooFormulaVar WZbin(lname.str().c_str(),
                                (lType[iT] + " W yield in SR per bin").c_str(),
                                "@0*@1",
                                (iT == 0) ? RooArgList(*TFWZ, binParZ)
                                          : RooArgList(*TFWZ, *EWKQCDbin));
            wspace.import(WZbin, RooFit::RecycleConflictNodes());

            // ------ CR transfer factors --------------------------------------
            for (unsigned iR = 1; iR < nR; ++iR) {

                unsigned crProc  = vproc[iT][iR];   // process in CR
                unsigned srProc  = (iR < 3) ? vproc[iT][iR]  // W CR: use same W proc in SR
                                             : vproc[iT][0];  // Z CR: use Z→νν in SR

                // Clamp negatives: same reason as W/Z ratio above.
                double crRaw = getBin(histos[iR][crProc][0], iB);
                double srRaw = getBin(histos[0][srProc][0], iB);
                double crNom = std::max(0.0, crRaw);
                double srNom = std::max(0.0, srRaw);
                double tfRatio = safeDiv(crNom, srNom);

                // ---- diagnostic ----
                std::cout << "[DIAG] " << lType[iT] << " CR=" << lRegionAlias[iR]
                          << " bin" << iB
                          << "  crRaw=" << crRaw << " srRaw=" << srRaw
                          << "  crNom=" << crNom << " srNom=" << srNom
                          << "  tfRatio=" << tfRatio << std::endl;

                double crErr = getBinErr(histos[iR][crProc][0], iB);
                double srErr = getBinErr(histos[0][srProc][0], iB);
                double tfStat = 1. + std::sqrt(
                    std::pow(safeDiv(crErr, crNom), 2) + std::pow(safeDiv(srErr, srNom), 2));

                lname.str(""); lname << lCategory << lType[iT] << "TF_" << lRegionAlias[iR] << "_stat_bin" << iB;
                TFstat[iT][iB][iR-1] = new RooRealVar(lname.str().c_str(), "CR/SR stat", 0);

                lname.str(""); lname << lCategory << lType[iT] << "TF_" << lRegionAlias[iR] << "_bin" << iB;
                std::ostringstream fCR;
                atIdx = 0;
                fCR << tfRatio;
                RooArgList crVars;

                {
                    // JER on CR/SR ratio — flat dummy (1.0) when turn_off_jecs=true
                    double jerCR = 1.0;
                    if (!turn_off_jecs && finputJES) {
                        const char *key = (iR < 3)
                            ? Form("wlnu_over_wmunu_%s_%s_jerUp", lTypeLC[iT].c_str(), year.c_str())
                            : Form("znunu_over_zmumu_%s_%s_jerUp", lTypeLC[iT].c_str(), year.c_str());
                        TH1D *h = (TH1D *)finputJES->Get(key);
                        if (h) jerCR = safeDiv(1., h->GetBinContent(1));
                    }
                    fCR << "*TMath::Power(" << jerCR << ",@" << atIdx++ << ")";
                    crVars.add(*jer);
                }

                // Stat uncertainty on CR/SR ratio
                fCR << "*TMath::Power(" << tfStat << ",@" << atIdx++ << ")";
                crVars.add(*TFstat[iT][iB][iR-1]);

                {
                    // JES on CR/SR ratio — flat dummy (1.0) when turn_off_jecs=true
                    for (unsigned iJ = 0; iJ < nJ; ++iJ) {
                        double jesCR = 1.0;
                        if (!turn_off_jecs && finputJES) {
                            const char *key = (iR < 3)
                                ? Form("wlnu_over_wmunu_%s_%s_%sUp",
                                       lTypeLC[iT].c_str(), year.c_str(), lJes[iJ].c_str())
                                : Form("znunu_over_zmumu_%s_%s_%sUp",
                                       lTypeLC[iT].c_str(), year.c_str(), lJes[iJ].c_str());
                            TH1D *h = (TH1D *)finputJES->Get(key);
                            if (h) jesCR = safeDiv(1., h->GetBinContent(1));
                        }
                        fCR << "*TMath::Power(" << jesCR << ",@" << atIdx++ << ")";
                        crVars.add(*jes[iJ]);
                    }
                }

                // Lepton SF terms: compute CR/SR ratio variation per systematic
                for (unsigned iN = 0; iN < nN; ++iN) {
                    unsigned iSup = 2*iN+1, iSdn = 2*iN+2;

                    TH1F *hCRup = histos[iR][crProc][iSup];
                    TH1F *hCRdn = histos[iR][crProc][iSdn];
                    TH1F *hSRup = histos[0][srProc][iSup];
                    TH1F *hSRdn = histos[0][srProc][iSdn];

                    double ratiovar_up = safeDiv(
                        safeDiv(getBin(hCRup, iB, crNom), getBin(hSRup, iB, srNom)),
                        tfRatio);
                    double ratiovar_dn = safeDiv(
                        safeDiv(getBin(hCRdn, iB, crNom), getBin(hSRdn, iB, srNom)),
                        tfRatio);

                    double crUpRaw = getBin(hCRup, iB, crNom);
                    double srUpRaw = getBin(hSRup, iB, srNom);
                    double crDnRaw = getBin(hCRdn, iB, crNom);
                    double srDnRaw = getBin(hSRdn, iB, srNom);

                    double rsUp = hardCodeNuisance[iR][iSup] >= 0
                        ? hardCodeNuisance[iR][iSup]
                        : 1. + (safeDiv(crUpRaw, srUpRaw) - tfRatio) / tfRatio;
                    double rsDn = hardCodeNuisance[iR][iSdn] >= 0
                        ? hardCodeNuisance[iR][iSdn]
                        : 1. + (safeDiv(crDnRaw, srDnRaw) - tfRatio) / tfRatio;

                    // ---- diagnostic: always print if unusual ----
                    bool unusual = (rsUp < 0 || rsDn < 0 || tfRatio == 0.0 ||
                                    crUpRaw < 0 || crDnRaw < 0);
                    if (unusual) {
                        std::cout << "[DIAG] LepSF syst=" << lNuisRun3[iN]
                                  << " " << lType[iT] << "/" << lRegionAlias[iR]
                                  << " bin" << iB
                                  << "  crUpRaw=" << crUpRaw << " srUpRaw=" << srUpRaw
                                  << "  crDnRaw=" << crDnRaw << " srDnRaw=" << srDnRaw
                                  << "  tfRatio=" << tfRatio
                                  << "  rsUp=" << rsUp << " rsDn=" << rsDn << std::endl;
                    }

                    if (rsUp < 0 || rsDn < 0) {
                        std::cout << " ERROR: negative ratiosyst for syst " << lNuisRun3[iN]
                                  << " bin " << iB << " region " << lRegionAlias[iR] << std::endl;
                        return 1;
                    }

                    // Skip if effect is negligible (< 0.1%)
                    if (std::fabs(rsUp - 1.) < 0.001 && std::fabs(1. - 1./rsDn) < 0.001) continue;

                    if (std::fabs(rsUp - 1./rsDn) > 0.001) {
                        fCR << "*((@" << atIdx << ">=0)*TMath::Power(" << rsUp << ",@" << atIdx
                            << ")+(@" << atIdx << "<0)*TMath::Power(" << safeDiv(1., rsDn) << ",@" << atIdx << "))";
                    } else {
                        fCR << "*TMath::Power(" << rsUp << ",@" << atIdx << ")";
                    }
                    crVars.add(*TFsysts[iN]);
                    atIdx++;
                }

                RooFormulaVar TF(lname.str().c_str(), "TF CR/SR", fCR.str().c_str(), crVars);
                wspace.import(TF, RooFit::RecycleConflictNodes());

                lname.str(""); lname << lCategory << lType[iT] << "V_" << lRegionAlias[iR] << "_bin" << iB;
                // W CRs scale off WZbin; Z CRs scale off Z SR yield
                RooFormulaVar CRbin(lname.str().c_str(),
                                    (lType[iT] + " V yield in CR per bin").c_str(),
                                    "@0*@1",
                                    (iR < 3) ? RooArgList(TF, WZbin)
                                    : (iT == 0) ? RooArgList(TF, binParZ)
                                               : RooArgList(TF, *EWKQCDbin));
                wspace.import(CRbin, RooFit::RecycleConflictNodes());

            } // CR loop
        } // type loop
    } // bin loop

    // ---- assemble RooParametricHist objects ---------------------------------
    std::cout << " - Assembling RooParametricHist objects." << std::endl;

    RooArgList Z_SR_bins[nT];
    RooArgList W_SR_bins[nT];
    RooArgList V_CR_bins[nT][nR-1];

    for (unsigned iT = 0; iT < nT; ++iT) {
        for (unsigned iB = 0; iB < nB; ++iB) {
            std::ostringstream lname;

            // Z SR bins
            if (iT == 0) {
                lname.str(""); lname << lCategory << lType[iT] << "Z_SR_bin" << iB+1;
                if (!wspace.var(lname.str().c_str())) {
                    std::cout << "ERROR: missing " << lname.str() << std::endl; return 1;
                }
                Z_SR_bins[iT].add(*wspace.var(lname.str().c_str()));
            } else {
                lname.str(""); lname << lCategory << "EWKQCD_SR_bin" << iB+1;
                if (!wspace.function(lname.str().c_str())) {
                    std::cout << "ERROR: missing " << lname.str() << std::endl; return 1;
                }
                Z_SR_bins[iT].add(*wspace.function(lname.str().c_str()));
            }

            // W SR bins
            lname.str(""); lname << lCategory << lType[iT] << "WZ_SR_bin" << iB+1;
            if (!wspace.function(lname.str().c_str())) {
                std::cout << "ERROR: missing " << lname.str() << std::endl; return 1;
            }
            W_SR_bins[iT].add(*wspace.function(lname.str().c_str()));

            // CR bins
            for (unsigned iR = 1; iR < nR; ++iR) {
                lname.str(""); lname << lCategory << lType[iT] << "V_" << lRegionAlias[iR] << "_bin" << iB+1;
                if (!wspace.function(lname.str().c_str())) {
                    std::cout << "ERROR: missing " << lname.str() << std::endl; return 1;
                }
                V_CR_bins[iT][iR-1].add(*wspace.function(lname.str().c_str()));
            }
        }
    }

    for (unsigned iT = 0; iT < nT; ++iT) {
        // SR parametric hists
        RooParametricHist p_Z((lType[iT] + "Z_SR").c_str(),
                              (lType[iT] + " Z+jets PDF in SR").c_str(),
                              lVarFit, Z_SR_bins[iT], dummyHist);
        RooAddition p_Z_norm((lType[iT] + "Z_SR_norm").c_str(),
                             ("N " + lType[iT] + " Z+jets SR").c_str(), Z_SR_bins[iT]);
        RooParametricHist p_W((lType[iT] + "W_SR").c_str(),
                              (lType[iT] + " W+jets PDF in SR").c_str(),
                              lVarFit, W_SR_bins[iT], dummyHist);
        RooAddition p_W_norm((lType[iT] + "W_SR_norm").c_str(),
                             ("N " + lType[iT] + " W+jets SR").c_str(), W_SR_bins[iT]);
        wspace.import(p_Z);
        wspace.import(p_Z_norm, RooFit::RecycleConflictNodes());
        wspace.import(p_W);
        wspace.import(p_W_norm, RooFit::RecycleConflictNodes());

        // CR parametric hists
        for (unsigned iR = 1; iR < nR; ++iR) {
            std::ostringstream lname;
            lname << lCategory << lType[iT] << "V_" << lRegionAlias[iR];
            RooParametricHist p_CR(lname.str().c_str(), "V+jets PDF in CR",
                                   lVarFit, V_CR_bins[iT][iR-1], dummyHist);
            lname << "_norm";
            RooAddition p_CR_norm(lname.str().c_str(), "N V+jets CR", V_CR_bins[iT][iR-1]);
            wspace.import(p_CR);
            wspace.import(p_CR_norm, RooFit::RecycleConflictNodes());
        }
    }

    // ---- write --------------------------------------------------------------
    fOut->cd();
    wspace.Write();
    fOut->Close();
    fOut->Delete();
    std::cout << "\nWritten: " << lOutFileName << std::endl;
    return 0;
}
