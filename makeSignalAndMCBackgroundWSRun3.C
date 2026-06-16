#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "TH1F.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooAddition.h"

// Run3 version of makeSignalAndMCBackgroundWS, built from the working Run2 version.
//
// Key differences from Run2:
//   - Flat ROOT file layout: nominal key is "{proc}_calc_weight",
//     systematic keys are "{proc}_weight_{syst}_{up/down}" — no per-region subdirectories.
//   - Input files live at: {run3_path}/{year}/{region}_{cat}_VBF/{var}/VBF_shapes.root
//   - Run3 process names and region names.
//   - Updated systematic list and CMS naming.
//   - Workspace object names use short region aliases (SR/Zee/Zmumu/Wenu/Wmunu) so
//     existing datacard templates remain compatible.
//   - jec_mode: "none" (Up=Down=nominal), "dummy" (flat Run2-derived, default),
//               "real" (per-bin shapes from vbf_shape_jes_uncs.root).
//   - turn_off_btag=true: flat dummy b-tag shapes; set false when real shapes available.

// Process indices — order must match lProcs[] below.
enum PROCESS {
    VBFH   = 0,  // vbfH
    GGH    = 1,  // ggH
    WH     = 2,  // WH
    ZH     = 3,  // ZH (combined qqZH + ggZH)
    TOP    = 4,  // Top
    VV     = 5,  // Diboson
    DY     = 6,  // QCD_Zll_NLONew
    EWKZll = 7,  // EWK_Zll
};

// ---- helpers ----------------------------------------------------------------

double findmax(TH1F *h) {
    double maxbv = h->GetBinContent(1);
    for (int b = 2; b <= h->GetNbinsX(); b++) {
        double v = h->GetBinContent(b);
        if (fabs(v) < 0.00001) continue;
        if (v > maxbv) maxbv = v;
    }
    return maxbv;
}

double findmin(TH1F *h) {
    double minbv = h->GetBinContent(1);
    for (int b = 2; b <= h->GetNbinsX(); b++) {
        double v = h->GetBinContent(b);
        if (fabs(v) < 0.00001) continue;
        if (v < minbv) minbv = v;
    }
    return minbv;
}

void makePlot(TDirectory *where, const std::string &name, const std::string &sys,
              TH1F *hC, TH1F *hU, TH1F *hD) {
    if (!hC || !hU || !hD) return;
    TCanvas *can = new TCanvas((name + "_" + sys).c_str(), "Syst plot", 700, 660);
    TPad p1("p1", "p1", 0.0, 0.3, 1, 0.98);
    TPad p2("p2", "p2", 0.0, 0.01, 1, 0.28);
    can->cd(); p1.Draw(); p2.Draw();
    TLegend leg(0.6, 0.78, 0.89, 0.89);
    leg.AddEntry(hC, "Nominal", "L");
    leg.AddEntry(hU, (sys + "Up").c_str(), "L");
    leg.AddEntry(hD, (sys + "Down").c_str(), "L");
    hC->SetTitle(name.c_str());
    hC->GetYaxis()->SetTitle("Events");
    hC->SetLineWidth(2); hC->SetLineColor(1);
    p1.cd(); hC->Draw("hist");
    hU->SetLineWidth(2); hU->SetLineColor(2); hU->Draw("histsame");
    hD->SetLineWidth(2); hD->SetLineColor(2); hD->SetLineStyle(2); hD->Draw("histsame");
    leg.Draw(); p1.SetLogy();
    TH1F *hUr = (TH1F *)hU->Clone(); hUr->Divide(hC);
    TH1F *hDr = (TH1F *)hD->Clone(); hDr->Divide(hC);
    p2.cd(); p2.SetGridy(); p2.SetGridx();
    hUr->SetTitle(""); hUr->GetYaxis()->SetTitle("Syst/Nominal");
    hUr->GetYaxis()->SetLabelSize(0.08); hUr->GetYaxis()->SetTitleSize(0.08);
    hUr->GetYaxis()->SetTitleOffset(0.6);
    double maxup = findmax(hUr), maxdn = findmax(hDr);
    double minup = findmin(hUr), mindn = findmin(hDr);
    hUr->SetMaximum(1.02 * (maxup > maxdn ? maxup : maxdn));
    hUr->SetMinimum(0.98 * (minup < mindn ? minup : mindn));
    hUr->Draw("hist"); hDr->Draw("histsame");
    where->WriteTObject(can);
}

// Helper: import a RooDataHist from histo into wspace, return true on success.
bool importHist(RooWorkspace &wspace, const std::string &name, RooArgList &vars, TH1F *h) {
    if (!h) { std::cout << " -- WARNING: null histogram for " << name << std::endl; return false; }
    RooDataHist *rdh = new RooDataHist(name.c_str(), name.c_str(), vars, h);
    wspace.import(*rdh);
    return true;
}

// ---- main function ----------------------------------------------------------

void makeSignalAndMCBackgroundWSRun3(
    std::string year       = "Run3Summer22_to_Run3Summer23BPix",
    std::string cat        = "VTR",
    bool classifier        = false,    // false=Mjj  true=SignalScore
    std::string jec_mode   = "dummy",  // "none" : Up=Down=nominal clone, no effect in fit
                                       // "dummy": flat Run2-derived shapes from vbf_shape_jes_uncs_run3.root
                                       // "real" : per-bin shapes from vbf_shape_jes_uncs.root
    bool turn_off_btag     = true      // true: flat dummy b-tag shapes; false: read from file (TODO)
) {
    // ---- variable and output naming -----------------------------------------
    std::string lVarDir    = classifier ? "SignalScore" : "Mjj";
    std::string lVarLabel  = lVarDir;  // used in RooRealVar name and output filename
    std::string lChannel   = "VBF";
    std::string lCategory  = cat + "_";
    std::string lYear      = year + "_";
    std::string lOutFileName = "signal_mc_bkgs_ws_" + lCategory + lYear + lChannel + lVarLabel + ".root";

    std::string varName, varTitle;
    double xmin, xmax;
    if (!classifier) {
        varName  = "mjj_" + cat;
        varTitle = "M_{jj} (GeV)";
        xmin = 200; xmax = 5000;
    } else {
        varName  = "bdt_" + cat;
        varTitle = "Signal Score";
        xmin = 0; xmax = 1;
    }
    RooRealVar lVarFit(varName.c_str(), varTitle.c_str(), xmin, xmax);
    RooArgList vars(lVarFit);

    TFile *fOut = new TFile(lOutFileName.c_str(), "RECREATE");
    RooWorkspace wspace("wspace_signal", "wspace_signal");
    TDirectory *outPlots = fOut->mkdir("Plots");

    // ---- input path ---------------------------------------------------------
    const std::string run3_path =
        "/vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT/hinvisible_mtr_vtr/root_files/";

    // ---- processes ----------------------------------------------------------
    // Run3 names as they appear in the ROOT file (prefix before _calc_weight)
    const unsigned nP = 8;
    std::string lProcs[nP] = {
        "vbfH",            // VBFH   = 0
        "ggH",             // GGH    = 1
        "WH",              // WH     = 2
        "ZH",              // ZH     = 3
        "Top",             // TOP    = 4
        "Diboson",         // VV     = 5
        "QCD_Zll_NLONew",  // DY     = 6
        "EWK_Zll"          // EWKZll = 7
    };

    // ---- regions ------------------------------------------------------------
    // lRegions: actual directory names in the ROOT file paths
    // lRegionAlias: short names used for workspace object labels (datacard compatibility)
    const unsigned nR = 5;
    std::string lRegions[nR]     = {"SR", "diElectron_CR", "diMuon_CR", "singleElectron_CR", "singleMuon_CR"};
    std::string lRegionAlias[nR] = {"SR", "Zee",           "Zmumu",     "Wenu",              "Wmunu"};

    // ---- weight systematics present in the Run3 ROOT files ------------------
    // Internal name (used to build histogram key: {proc}_weight_{name}_{up/down})
    // and corresponding CMS convention name (used in workspace object labels).
    // Photon systematics (photon_id, photon_id_high_pt) are deliberately excluded
    // since the processes in this workspace are not photon-related.
    const unsigned nN = 14;
    std::string lSystsRun3[nN] = {
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
    std::string lSystsCMS[nN] = {
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
    // corrCat: true = correlated across categories (no category suffix appended)
    const bool corrCat[nN]  = {false, true, true, true, true, true, true, true, true, true, true, true, true, true};
    // corrYear: true = correlated across years (no year suffix appended)
    const bool corrYear[nN] = {false, false, false, false, false, true, true, true, true, true, true, true, true, true};

    // Apply year/cat decorations to CMS names
    for (unsigned iN = 0; iN < nN; ++iN) {
        if (!corrCat[iN])  lSystsCMS[iN] += "_" + cat;
        if (!corrYear[iN]) lSystsCMS[iN] += "_" + year;
    }

    // Skip matrix: skipSyst[iP][iN] = true means don't add this systematic for this process.
    //
    // ewk_correction (index 9): NLO EWK correction applied to QCD-produced V+jets.
    //   Only meaningful for QCD_Zll_NLONew (DY=6). Skip for everything else.
    //
    // qcd_correction (index 10): NLO QCD correction applied to EWK-produced V+jets.
    //   Only meaningful for EWK_Zll (7). Skip for everything else.
    //
    // Higgs signal theory uncertainties (QCD scale, PDF) are handled as lnN in the
    // datacard rather than as shape systematics here, so ewk/qcd corrections are
    // skipped for vbfH, ggH, WH, ZH as well.
    bool skipSyst[nP][nN];
    for (unsigned iP = 0; iP < nP; ++iP)
        for (unsigned iN = 0; iN < nN; ++iN)
            skipSyst[iP][iN] = false;

    // ewk_correction (9): only apply to QCD_Zll_NLONew (DY=6)
    for (unsigned iP = 0; iP < nP; ++iP)
        if (iP != PROCESS::DY) skipSyst[iP][9] = true;

    // qcd_correction (10): only apply to EWK_Zll (EWKZll=7)
    for (unsigned iP = 0; iP < nP; ++iP)
        if (iP != PROCESS::EWKZll) skipSyst[iP][10] = true;

    // mur (11) and muf (12): QCD scale variations — skip for Higgs signals
    // (signal QCD scale uncertainties are lnN in the datacard, not shape)
    for (unsigned iP : {PROCESS::VBFH, PROCESS::GGH, PROCESS::WH, PROCESS::ZH}) {
        skipSyst[iP][11] = true;
        skipSyst[iP][12] = true;
    }
    // pileup (13): applies to all — no skips

    // ---- JEC/JER list (12 entries, same scheme as Run2) ---------------------
    const unsigned nJ = 12;
    std::string lJes[nJ] = {
        "jer",
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
    // JES process label: used for "dummy" and "real" modes to look up shapes in the JEC file
    std::string lJESLabel[nP] = {
        "VBF_HToInvisible_",
        "ZJetsToNuNu",
        "ZJetsToNuNu",
        "ZJetsToNuNu",
        "ZJetsToNuNu",
        "ZJetsToNuNu",
        "ZJetsToNuNu",
        "EWKZJetsToNuNu"
    };

    TFile *finputJES = nullptr;
    if (jec_mode == "dummy") {
        finputJES = TFile::Open("vbf_shape_jes_uncs_run3.root");
        if (!finputJES) {
            std::cout << "ERROR: vbf_shape_jes_uncs_run3.root not found. "
                      << "Run make_run3_jec_uncs.py first, or use jec_mode=\"none\"." << std::endl;
            return;
        }
    } else if (jec_mode == "real") {
        finputJES = TFile::Open("vbf_shape_jes_uncs.root");
        if (!finputJES) {
            std::cout << "ERROR: vbf_shape_jes_uncs.root not found. "
                      << "Use jec_mode=\"dummy\" or jec_mode=\"none\" instead." << std::endl;
            return;
        }
    }

    // ---- main loop over regions and processes --------------------------------
    for (unsigned iR = 0; iR < nR; ++iR) {

        std::string lInFileName = run3_path + year + "/" +
                                  lRegions[iR] + "_" + cat + "_VBF/" +
                                  lVarDir + "/VBF_shapes.root";
        std::cout << "\nOpening: " << lInFileName << std::endl;
        TFile *finput = TFile::Open(lInFileName.c_str());
        if (!finput) {
            std::cout << "ERROR: could not open " << lInFileName << std::endl;
            return;
        }

        for (unsigned iP = 0; iP < nP; ++iP) {

            // In control regions only include the small MC backgrounds
            // (Top, Diboson, QCD_Zll_NLONew, EWK_Zll) — not signals or VH.
            if (iR > 0 && iP < PROCESS::TOP) continue;

            // ---- nominal histogram ------------------------------------------
            std::string nomKey = lProcs[iP] + "_calc_weight";
            TH1F *Thist = (TH1F *)finput->Get(nomKey.c_str());
            if (!Thist) {
                std::cout << " ERROR: nominal histogram " << nomKey
                          << " not found in " << lInFileName << std::endl;
                return;
            }

            // Floor negative bins to zero: NLO samples with negative event weights
            // can produce negative bin content in tight kinematic regions (especially
            // the lowest Mjj bin). Combine FASTEXITs and Minuit cannot converge when
            // the total predicted yield goes negative, so we apply a hard zero floor.
            for (int iB = 1; iB <= Thist->GetNbinsX(); ++iB) {
                if (Thist->GetBinContent(iB) < 0.) {
                    std::cout << " -- WARNING: flooring negative bin " << iB
                              << " (" << Thist->GetBinContent(iB) << ") to 0 for "
                              << nomKey << " in " << lRegions[iR] << std::endl;
                    Thist->SetBinContent(iB, 0.);
                    Thist->SetBinError(iB, 0.);
                }
            }

            std::string wsObjBase = lProcs[iP] + "_hist_" + lRegionAlias[iR];
            importHist(wspace, wsObjBase, vars, Thist);
            std::cout << " Imported nominal: " << wsObjBase << std::endl;

            // ---- weight systematics -----------------------------------------
            for (unsigned iS = 0; iS < nN; ++iS) {
                if (skipSyst[iP][iS]) {
                    std::cout << " -- Skipping " << lSystsRun3[iS]
                              << " for " << lProcs[iP] << std::endl;
                    continue;
                }

                std::string upKey   = lProcs[iP] + "_weight_" + lSystsRun3[iS] + "_up";
                std::string downKey = lProcs[iP] + "_weight_" + lSystsRun3[iS] + "_down";
                TH1F *hUp   = (TH1F *)finput->Get(upKey.c_str());
                TH1F *hDown = (TH1F *)finput->Get(downKey.c_str());

                if (!hUp || !hDown) {
                    // Systematic not present for this sample (e.g. pileup not
                    // stored for WH/ZH). Add a flat dummy (up=down=nominal) so
                    // the workspace always contains these objects and the
                    // datacard can reference them without crashing Combine.
                    std::cout << " -- WARNING: systematic " << lSystsRun3[iS]
                              << " not found for " << lProcs[iP]
                              << " in " << lRegions[iR] << " — using flat dummy." << std::endl;
                    hUp   = (TH1F *)Thist->Clone(
                        Form("%s_flatdummy_%s_up",   lProcs[iP].c_str(), lSystsRun3[iS].c_str()));
                    hDown = (TH1F *)Thist->Clone(
                        Form("%s_flatdummy_%s_down", lProcs[iP].c_str(), lSystsRun3[iS].c_str()));
                }

                // Floor syst histograms too — a zero nominal can have non-zero syst
                for (int iB = 1; iB <= hUp->GetNbinsX(); ++iB) {
                    if (hUp->GetBinContent(iB)   < 0.) { hUp->SetBinContent(iB, 0.);   hUp->SetBinError(iB, 0.); }
                    if (hDown->GetBinContent(iB) < 0.) { hDown->SetBinContent(iB, 0.); hDown->SetBinError(iB, 0.); }
                }

                std::string cmsUpName   = wsObjBase + "_" + lSystsCMS[iS] + "Up";
                std::string cmsDownName = wsObjBase + "_" + lSystsCMS[iS] + "Down";
                importHist(wspace, cmsUpName,   vars, hUp);
                importHist(wspace, cmsDownName, vars, hDown);
                makePlot(outPlots, wsObjBase, lSystsRun3[iS], Thist, hUp, hDown);
                std::cout << " Imported syst: " << lSystsRun3[iS]
                          << " -> " << lSystsCMS[iS] << std::endl;
            }

            // ---- b-tagging systematic ---------------------------------------
            // b-tag efficiency for Top, fake-b for everything else.
            std::string btagCMSName = (iP == PROCESS::TOP)
                ? ("CMS_eff_b_" + year)
                : ("CMS_fake_b_" + year);

            if (turn_off_btag) {
                // Flat dummy: up = down = nominal clone (no shape effect).
                TH1F *hBUp   = (TH1F *)Thist->Clone(Form("%s_btag_dummy_up",   lProcs[iP].c_str()));
                TH1F *hBDown = (TH1F *)Thist->Clone(Form("%s_btag_dummy_down", lProcs[iP].c_str()));
                importHist(wspace, wsObjBase + "_" + btagCMSName + "Up",   vars, hBUp);
                importHist(wspace, wsObjBase + "_" + btagCMSName + "Down", vars, hBDown);
                std::cout << " b-tag dummy (flat) for " << btagCMSName << std::endl;
            } else {
                // TODO: read real b-tag shapes from dedicated input file.
                std::cout << " WARNING: real b-tag shapes not yet implemented for Run3."
                          << " Falling back to flat dummy for " << btagCMSName << std::endl;
                TH1F *hBUp   = (TH1F *)Thist->Clone(Form("%s_btag_dummy_up",   lProcs[iP].c_str()));
                TH1F *hBDown = (TH1F *)Thist->Clone(Form("%s_btag_dummy_down", lProcs[iP].c_str()));
                importHist(wspace, wsObjBase + "_" + btagCMSName + "Up",   vars, hBUp);
                importHist(wspace, wsObjBase + "_" + btagCMSName + "Down", vars, hBDown);
            }

            // ---- JEC / JER systematics --------------------------------------
            for (unsigned iJ = 0; iJ < nJ; ++iJ) {
                std::string jecCMSName = (iJ == 0)
                    ? ("CMS_res_j_" + year)
                    : ("CMS_scale_j_" + lJes[iJ]);

                TH1F *hJUp = nullptr, *hJDown = nullptr;

                if (jec_mode == "none") {
                    // Up = Down = nominal: zero effect in the fit.
                    hJUp   = (TH1F *)Thist->Clone(Form("%s_jecNone_%s_up",   lProcs[iP].c_str(), lJes[iJ].c_str()));
                    hJDown = (TH1F *)Thist->Clone(Form("%s_jecNone_%s_down", lProcs[iP].c_str(), lJes[iJ].c_str()));
                } else {
                    // Read ratio histograms from the JES input file and apply to nominal.
                    finputJES->cd();
                    TH1D *hRUp   = (TH1D *)finputJES->Get(Form("%s%s_%sUp",
                                       lJESLabel[iP].c_str(), year.c_str(), lJes[iJ].c_str()));
                    TH1D *hRDown = (TH1D *)finputJES->Get(Form("%s%s_%sDown",
                                       lJESLabel[iP].c_str(), year.c_str(), lJes[iJ].c_str()));
                    if (!hRUp || !hRDown) {
                        std::cout << " -- WARNING: JES ratio not found for "
                                  << lJes[iJ] << " process " << lProcs[iP]
                                  << " — using flat dummy." << std::endl;
                        hJUp   = (TH1F *)Thist->Clone();
                        hJDown = (TH1F *)Thist->Clone();
                    } else {
                        hJUp   = (TH1F *)Thist->Clone();
                        hJDown = (TH1F *)Thist->Clone();
                        for (int b = 1; b <= Thist->GetNbinsX(); ++b) {
                            double yv = hJUp->GetBinContent(b);
                            // Clamp FindBin result to [1, nBinsX] to avoid underflow/overflow
                            // bins (content=0) for non-Mjj fit variables (e.g. SignalScore [0,1]
                            // whose bin centres fall below the ratio histogram's Mjj range).
                            auto clampBin = [](TH1D *h, double x) {
                                int r = h->FindBin(x);
                                if (r < 1)            r = 1;
                                if (r > h->GetNbinsX()) r = h->GetNbinsX();
                                return r;
                            };
                            hJUp->SetBinContent(b,   yv * hRUp->GetBinContent(clampBin(hRUp,   hJUp->GetBinCenter(b))));
                            hJDown->SetBinContent(b, yv * hRDown->GetBinContent(clampBin(hRDown, hJDown->GetBinCenter(b))));
                        }
                    }
                    hJUp->SetName(Form("%s_%s_up",   wsObjBase.c_str(), lJes[iJ].c_str()));
                    hJDown->SetName(Form("%s_%s_down", wsObjBase.c_str(), lJes[iJ].c_str()));
                }

                importHist(wspace, wsObjBase + "_" + jecCMSName + "Up",   vars, hJUp);
                importHist(wspace, wsObjBase + "_" + jecCMSName + "Down", vars, hJDown);
                if (jec_mode != "none")
                    makePlot(outPlots, wsObjBase, lJes[iJ], Thist, hJUp, hJDown);
            }

        } // process loop

        finput->Close();

    } // region loop

    // ---- write output -------------------------------------------------------
    fOut->WriteTObject(&wspace);
    fOut->Close();
    std::cout << "\nWritten: " << lOutFileName << std::endl;
}
