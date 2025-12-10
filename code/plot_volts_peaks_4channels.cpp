// monitor_and_compare_pmts_mod.cpp

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

int main(int argc, char** argv) {
    // ------------ Configuration ------------
    std::vector<int> pulseWidths;
    // only 0 and 900 ps are processed
    pulseWidths.push_back(0);
    pulseWidths.push_back(900);

    std::vector<int> angleLaserBall;
    angleLaserBall.push_back(0);
    angleLaserBall.push_back(90);

    std::vector<std::string> PMT3_list;
    PMT3_list.push_back("-30");
    PMT3_list.push_back("H");

    std::vector<std::string> PMT4_list;
    PMT4_list.push_back("V");
    PMT4_list.push_back("+30");

    std::string treeNamePMT   = "PMT1data";  // monitor
    std::string treeNameTrig  = "PMT2data";  // trigger
    std::string treeNamePMT3  = "PMT3data";  // extra PMT
    std::string treeNamePMT4  = "PMT4data";  // extra PMT
    std::string branchName    = "volts";

    // ------------ Histogram settings (reduced bins) ------------
    int    nBins   = 200;   // was 400
    double xMin    = -0.05;
    double xMax    =  0.8;

    int    nBinsPeak   = 200;   // was 400
    double peakMin     = -0.05;
    double peakMax     =  0.8;

    int    nBinsSum    = 100;   // was 200
    double sumMin      = -0.05;
    double sumMax      =  1.0;

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    double sum  = 0.0;

    // Output names kept generic – you can adapt to your convention
    std::string outputRootName   = "../results/LaserBall_PMT1_Mon_PMT3_PMT4_0ps_900ps_spectrum.root";

    gStyle->SetOptStat(0);

    // ------------ Colors (fixed mapping by PMT, not by loop index) ------------
    // *** CHANGED: explicit mapping so colors are guaranteed different
    int colPMT1 = kRed + 1;
    int colPMT3 = kBlue + 1;
    int colPMT4 = kGreen + 2;

    // ------------ Separate canvases for 0 ps and 900 ps ------------
    // 0 ps
    TCanvas* c0_spec  = new TCanvas("c0_spec",  "Spectrum 0 ps",  1200, 1200);
    TCanvas* c0_peak  = new TCanvas("c0_peak",  "Peaks 0 ps",     1200, 800);
    TCanvas* c0_sum   = new TCanvas("c0_sum",   "Sum4 0 ps",      1200, 800);
    c0_spec->SetGrid(); c0_spec->SetLogy();
    c0_peak->SetGrid(); c0_sum->SetGrid();

    TLegend* l0_spec  = new TLegend(0.45, 0.50, 0.98, 0.98);
    TLegend* l0_peak  = new TLegend(0.45, 0.50, 0.98, 0.98);
    TLegend* l0_sum   = new TLegend(0.45, 0.50, 0.98, 0.98);
    l0_spec->SetBorderSize(0); l0_spec->SetFillStyle(0); l0_spec->SetTextSize(0.025);
    l0_peak->SetBorderSize(0); l0_peak->SetFillStyle(0); l0_peak->SetTextSize(0.025);
    l0_sum->SetBorderSize(0);  l0_sum->SetFillStyle(0);  l0_sum->SetTextSize(0.025);

    // 900 ps
    TCanvas* c900_spec  = new TCanvas("c900_spec",  "Spectrum 900 ps",  1200, 1200);
    TCanvas* c900_peak  = new TCanvas("c900_peak",  "Peaks 900 ps",     1200, 800);
    TCanvas* c900_sum   = new TCanvas("c900_sum",   "Sum4 900 ps",      1200, 800);
    c900_spec->SetGrid(); c900_spec->SetLogy();
    c900_peak->SetGrid(); c900_sum->SetGrid();

    TLegend* l900_spec  = new TLegend(0.45, 0.50, 0.98, 0.98);
    TLegend* l900_peak  = new TLegend(0.45, 0.50, 0.98, 0.98);
    TLegend* l900_sum   = new TLegend(0.45, 0.50, 0.98, 0.98);
    l900_spec->SetBorderSize(0); l900_spec->SetFillStyle(0); l900_spec->SetTextSize(0.025);
    l900_peak->SetBorderSize(0); l900_peak->SetFillStyle(0); l900_peak->SetTextSize(0.025);
    l900_sum->SetBorderSize(0);  l900_sum->SetFillStyle(0);  l900_sum->SetTextSize(0.025);

    // Internal flags to know if first histogram was drawn on a given canvas
    bool drawn0_spec  = false, drawn0_peak  = false, drawn0_sum  = false;
    bool drawn900_spec= false, drawn900_peak= false, drawn900_sum= false;

    std::vector<TH1D*> allHists; // for writing to ROOT later

    std::ofstream outTxt("../results/LaserBall_PMT_volts_summary.txt");
    if (!outTxt.is_open()) {
        std::cerr << "Error: could not open output summary for writing\n";
        return 1;
    }
    outTxt << "pw_ps angle_deg PMT3_cfg PMT4_cfg pw_sum1 min1 max1 mean1 sigma1 "
    << "meanSum3 sigmaSum3 meanSum4 sigmaSum4\n";

    // ======================= LOOP over configurations =======================
    for (size_t iPW = 0; iPW < pulseWidths.size(); ++iPW) {
        int pw = pulseWidths[iPW];

        for (size_t iAngle = 0; iAngle < angleLaserBall.size(); ++iAngle) {
            int ang = angleLaserBall[iAngle];

            for (size_t iPMT = 0; iPMT < PMT3_list.size(); ++iPMT) {
                std::string PMT3_name = PMT3_list[iPMT];
                std::string PMT4_name = PMT4_list[iPMT];

                std::string fileName =
                "../../data/LaserBall/hadded_PMT1_Mon_PMT2_Trig_PMT3_" + PMT3_name +
                "_PMT4_" + PMT4_name + "_results_" +
                std::to_string(pw) + "ps_" + std::to_string(ang) + "deg.root";

                std::cout << "Looking at file: " << fileName << std::endl;

                if (gSystem->AccessPathName(fileName.c_str())) {
                    std::cerr << "Warning: File not found: " << fileName << std::endl;
                    continue;
                }

                TFile* f = TFile::Open(fileName.c_str(), "READ");
                if (!f || f->IsZombie()) {
                    std::cerr << "Error: Could not open file: " << fileName << std::endl;
                    if (f) f->Close();
                    continue;
                }

                TTree* treePMT  = dynamic_cast<TTree*>(f->Get(treeNamePMT.c_str()));
                TTree* treeTrig = dynamic_cast<TTree*>(f->Get(treeNameTrig.c_str()));
                TTree* treePMT3 = dynamic_cast<TTree*>(f->Get(treeNamePMT3.c_str()));
                TTree* treePMT4 = dynamic_cast<TTree*>(f->Get(treeNamePMT4.c_str()));

                if (!treePMT || !treeTrig) {
                    std::cerr << "Error: Trees " << treeNamePMT << " or " << treeNameTrig
                    << " not found in file: " << fileName << std::endl;
                    f->Close();
                    continue;
                }

                bool hasPMT3 = (treePMT3 != nullptr);
                bool hasPMT4 = (treePMT4 != nullptr);
                if (!hasPMT3) std::cerr << "Note: " << treeNamePMT3 << " not found.\n";
                if (!hasPMT4) std::cerr << "Note: " << treeNamePMT4 << " not found.\n";

                Long64_t nEntriesTrig = treeTrig->GetEntries();
                Long64_t nEntriesPMT  = treePMT->GetEntries();
                Long64_t nEntriesPMT3 = hasPMT3 ? treePMT3->GetEntries() : nEntriesTrig;
                Long64_t nEntriesPMT4 = hasPMT4 ? treePMT4->GetEntries() : nEntriesTrig;

                Long64_t nEntries = std::min(std::min(nEntriesTrig, nEntriesPMT),
                                             std::min(nEntriesPMT3, nEntriesPMT4));
                if (nEntries == 0) {
                    std::cerr << "Warning: no entries in trees for file " << fileName << std::endl;
                    f->Close();
                    continue;
                }

                double voltsPMT = 0.0, voltsTrig = 0.0, voltsPMT3 = 0.0, voltsPMT4 = 0.0;
                treePMT->SetBranchAddress(branchName.c_str(), &voltsPMT);
                treeTrig->SetBranchAddress(branchName.c_str(), &voltsTrig);
                if (hasPMT3) treePMT3->SetBranchAddress(branchName.c_str(), &voltsPMT3);
                if (hasPMT4) treePMT4->SetBranchAddress(branchName.c_str(), &voltsPMT4);

                std::vector<double> pmt1Vals(nEntries), trigVals(nEntries), pmt3Vals, pmt4Vals;
                if (hasPMT3) pmt3Vals.resize(nEntries);
                if (hasPMT4) pmt4Vals.resize(nEntries);

                // ---------- Histograms for this configuration ----------
                std::string tag = Form("_%dps_%ddeg_P3%s_P4%s",
                                       pw, ang, PMT3_name.c_str(), PMT4_name.c_str());

                TH1D* h_spec1 = new TH1D(("h_Volts1" + tag).c_str(),
                                         ("PMT1 volts " + tag).c_str(),
                                         nBins, xMin, xMax);
                h_spec1->SetDirectory(nullptr);
                h_spec1->SetLineColor(colPMT1+iPMT+iAngle);
                h_spec1->SetLineWidth(2);

                TH1D* h_peak1 = new TH1D(("h_Peaks1" + tag).c_str(),
                                         ("PMT1 peaks " + tag).c_str(),
                                         nBinsPeak, peakMin, peakMax);
                h_peak1->SetDirectory(nullptr);
                h_peak1->SetLineColor(colPMT1+iPMT+iAngle);
                h_peak1->SetLineWidth(2);

                TH1D* h_sum1 = new TH1D(("h_sum4_PMT1" + tag).c_str(),
                                        ("PMT1 sum4 " + tag).c_str(),
                                        nBinsSum, sumMin, sumMax);
                h_sum1->SetDirectory(nullptr);
                h_sum1->SetLineColor(colPMT1+iPMT+iAngle);
                h_sum1->SetLineWidth(2);

                TH1D* h_peak3 = nullptr; TH1D* h_sum3 = nullptr;
                TH1D* h_peak4 = nullptr; TH1D* h_sum4 = nullptr;

                if (hasPMT3) {
                    h_peak3 = new TH1D(("h_Peaks3" + tag).c_str(),
                                       ("PMT3 peaks " + tag).c_str(),
                                       nBinsPeak, peakMin, peakMax);
                    h_peak3->SetDirectory(nullptr);
                    h_peak3->SetLineColor(colPMT3+iPMT+iAngle);
                    h_peak3->SetLineWidth(2);

                    h_sum3 = new TH1D(("h_sum4_PMT3" + tag).c_str(),
                                      ("PMT3 sum4 " + tag).c_str(),
                                      nBinsSum, sumMin, sumMax);
                    h_sum3->SetDirectory(nullptr);
                    h_sum3->SetLineColor(colPMT3+iPMT+iAngle);
                    h_sum3->SetLineWidth(2);
                }

                if (hasPMT4) {
                    h_peak4 = new TH1D(("h_Peaks4" + tag).c_str(),
                                       ("PMT4 peaks " + tag).c_str(),
                                       nBinsPeak, peakMin, peakMax);
                    h_peak4->SetDirectory(nullptr);
                    h_peak4->SetLineColor(colPMT4+iPMT+iAngle);
                    h_peak4->SetLineWidth(2);

                    h_sum4 = new TH1D(("h_sum4_PMT4" + tag).c_str(),
                                      ("PMT4 sum4 " + tag).c_str(),
                                      nBinsSum, sumMin, sumMax);
                    h_sum4->SetDirectory(nullptr);
                    h_sum4->SetLineColor(colPMT4+iPMT+iAngle);
                    h_sum4->SetLineWidth(2);
                }

                // Fill in-memory arrays + PMT1 spectrum
                for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
                    treePMT->GetEntry(iEntry);
                    treeTrig->GetEntry(iEntry);
                    pmt1Vals[iEntry] = voltsPMT;
                    trigVals[iEntry] = voltsTrig;
                    if (hasPMT3) { treePMT3->GetEntry(iEntry); pmt3Vals[iEntry] = voltsPMT3; }
                    if (hasPMT4) { treePMT4->GetEntry(iEntry); pmt4Vals[iEntry] = voltsPMT4; }

                    h_spec1->Fill(voltsPMT);
                    if (voltsPMT < vmin) vmin = voltsPMT;
                    if (voltsPMT > vmax) vmax = voltsPMT;
                    sum += voltsPMT;
                }

                // Trigger periods
                double thr = -0.6;
                std::vector<Long64_t> periodStarts;
                for (Long64_t iEntry = 1; iEntry < nEntries; ++iEntry) {
                    if (trigVals[iEntry-1] < thr && trigVals[iEntry] >= thr) {
                        periodStarts.push_back(iEntry);
                    }
                }

                auto computeSum4 = [&](const std::vector<double>& arr, Long64_t idx) {
                    Long64_t i0 = (idx > 0) ? (idx - 1) : idx;
                    Long64_t i1 = idx;
                    Long64_t i2 = (idx + 1 < nEntries) ? (idx + 1) : idx;
                    Long64_t i3 = (idx + 2 < nEntries) ? (idx + 2) : idx;
                    return arr[i0] + arr[i1] + arr[i2] + arr[i3];
                };

                for (size_t ip = 0; ip < periodStarts.size(); ++ip) {
                    Long64_t start = periodStarts[ip];
                    Long64_t end   = (ip + 1 < periodStarts.size()) ? periodStarts[ip+1] : nEntries;
                    if (end <= start) continue;

                    // PMT1 peak
                    double peak1 = -1e9; Long64_t idx1 = start;
                    for (Long64_t i = start; i < end; ++i) {
                        if (pmt1Vals[i] > peak1) { peak1 = pmt1Vals[i]; idx1 = i; }
                    }
                    h_peak1->Fill(peak1);
                    h_sum1->Fill(computeSum4(pmt1Vals, idx1));

                    // PMT3
                    if (hasPMT3) {
                        double peak3 = -1e9; Long64_t idx3 = start;
                        for (Long64_t i = start; i < end; ++i) {
                            if (pmt3Vals[i] > peak3) { peak3 = pmt3Vals[i]; idx3 = i; }
                        }
                        h_peak3->Fill(peak3);
                        h_sum3->Fill(computeSum4(pmt3Vals, idx3));
                    }

                    // PMT4
                    if (hasPMT4) {
                        double peak4 = -1e9; Long64_t idx4 = start;
                        for (Long64_t i = start; i < end; ++i) {
                            if (pmt4Vals[i] > peak4) { peak4 = pmt4Vals[i]; idx4 = i; }
                        }
                        h_peak4->Fill(peak4);
                        h_sum4->Fill(computeSum4(pmt4Vals, idx4));
                    }
                }

                // Fit PMT1 sum
                TF1* gaus1 = new TF1(("gaus1" + tag).c_str(), "gaus", sumMin, sumMax);
                gaus1->SetParameters(h_sum1->GetMaximum(), h_sum1->GetMean(), h_sum1->GetRMS());
                gaus1->SetLineColor(kBlack);
                gaus1->SetLineStyle(2);
                gaus1->SetLineWidth(2);
                h_sum1->Fit(gaus1, "Q0");
                double mean1  = gaus1->GetParameter(1);
                double sigma1 = gaus1->GetParameter(2);

                double meanSum3 = 0.0, sigmaSum3 = 0.0;
                double meanSum4 = 0.0, sigmaSum4 = 0.0;

                if (hasPMT3 && h_sum3->GetEntries() > 0) {
                    if (h_sum3->GetEntries() > 10) {
                        TF1* gaus3 = new TF1(("gaus3" + tag).c_str(), "gaus", sumMin, sumMax);
                        gaus3->SetParameters(h_sum3->GetMaximum(), h_sum3->GetMean(), h_sum3->GetRMS());
                        h_sum3->Fit(gaus3, "Q0");
                        meanSum3  = gaus3->GetParameter(1);
                        sigmaSum3 = gaus3->GetParameter(2);
                        delete gaus3;
                    } else {
                        meanSum3  = h_sum3->GetMean();
                        sigmaSum3 = h_sum3->GetRMS();
                    }
                }

                if (hasPMT4 && h_sum4->GetEntries() > 0) {
                    if (h_sum4->GetEntries() > 10) {
                        TF1* gaus4 = new TF1(("gaus4" + tag).c_str(), "gaus", sumMin, sumMax);
                        gaus4->SetParameters(h_sum4->GetMaximum(), h_sum4->GetMean(), h_sum4->GetRMS());
                        h_sum4->Fit(gaus4, "Q0");
                        meanSum4  = gaus4->GetParameter(1);
                        sigmaSum4 = gaus4->GetParameter(2);
                        delete gaus4;
                    } else {
                        meanSum4  = h_sum4->GetMean();
                        sigmaSum4 = h_sum4->GetRMS();
                    }
                }

                // ===================== DRAWING / LEGENDS =====================

                // Human-readable laser ball position string
                // (includes both angle and delay as you requested)
                std::string laserPosStr = Form("%d ps, %d deg", pw, ang);

                // Legend labels: include PMT name + laser ball position
                std::string labPMT1 = "PMT1, " + laserPosStr;
                std::string labPMT3 = "PMT3, " + laserPosStr +
                " (" + PMT3_name + ")";
                std::string labPMT4 = "PMT4, " + laserPosStr +
                " (" + PMT4_name + ")";

                // Choose which canvases to draw on based on pw
                bool is0ps   = (pw == 0);
                bool is900ps = (pw == 900);

                if (is0ps) {
                    // -------- 0 ps canvases --------
                    // Spectrum
                    c0_spec->cd();
                    if (!drawn0_spec) {
                        h_spec1->SetTitle("Spectrum at 0 ps;volts [V];entries");
                        h_spec1->Draw("HIST");
                        drawn0_spec = true;
                    } else {
                        h_spec1->Draw("HIST SAME");
                    }
                    l0_spec->AddEntry(h_spec1, labPMT1.c_str(), "l");

                    // Peaks
                    c0_peak->cd();
                    if (!drawn0_peak) {
                        h_peak1->SetTitle("Peaks at 0 ps;peak volts [V];entries");
                        h_peak1->Draw("HIST");
                        drawn0_peak = true;
                    } else {
                        h_peak1->Draw("HIST SAME");
                    }
                    l0_peak->AddEntry(h_peak1, labPMT1.c_str(), "l");
                    if (hasPMT3) { h_peak3->Draw("HIST SAME"); l0_peak->AddEntry(h_peak3, labPMT3.c_str(), "l"); }
                    if (hasPMT4) { h_peak4->Draw("HIST SAME"); l0_peak->AddEntry(h_peak4, labPMT4.c_str(), "l"); }

                    // Sums
                    c0_sum->cd();
                    if (!drawn0_sum) {
                        h_sum1->SetTitle("Sum4 at 0 ps;sum4 [V];entries");
                        h_sum1->Draw("HIST");
                        gaus1->Draw("SAME");
                        drawn0_sum = true;
                    } else {
                        h_sum1->Draw("HIST SAME");
                        gaus1->Draw("SAME");
                    }
                    l0_sum->AddEntry(h_sum1,
                                     Form("%s (mean=%.3f, #sigma=%.3f)", labPMT1.c_str(), mean1, sigma1),
                                     "l");
                    if (hasPMT3) {
                        h_sum3->Draw("HIST SAME");
                        l0_sum->AddEntry(h_sum3,
                                         Form("%s (mean=%.3f, #sigma=%.3f)", labPMT3.c_str(), meanSum3, sigmaSum3),
                                         "l");
                    }
                    if (hasPMT4) {
                        h_sum4->Draw("HIST SAME");
                        l0_sum->AddEntry(h_sum4,
                                         Form("%s (mean=%.3f, #sigma=%.3f)", labPMT4.c_str(), meanSum4, sigmaSum4),
                                         "l");
                    }
                }

                if (is900ps) {
                    // -------- 900 ps canvases --------
                    // Spectrum
                    c900_spec->cd();
                    if (!drawn900_spec) {
                        h_spec1->SetTitle("Spectrum at 900 ps;volts [V];entries");
                        h_spec1->Draw("HIST");
                        drawn900_spec = true;
                    } else {
                        h_spec1->Draw("HIST SAME");
                    }
                    l900_spec->AddEntry(h_spec1, labPMT1.c_str(), "l");

                    // Peaks
                    c900_peak->cd();
                    if (!drawn900_peak) {
                        h_peak1->SetTitle("Peaks at 900 ps;peak volts [V];entries");
                        h_peak1->Draw("HIST");
                        drawn900_peak = true;
                    } else {
                        h_peak1->Draw("HIST SAME");
                    }
                    l900_peak->AddEntry(h_peak1, labPMT1.c_str(), "l");
                    if (hasPMT3) { h_peak3->Draw("HIST SAME"); l900_peak->AddEntry(h_peak3, labPMT3.c_str(), "l"); }
                    if (hasPMT4) { h_peak4->Draw("HIST SAME"); l900_peak->AddEntry(h_peak4, labPMT4.c_str(), "l"); }

                    // Sums
                    c900_sum->cd();
                    if (!drawn900_sum) {
                        h_sum1->SetTitle("Sum4 at 900 ps;sum4 [V];entries");
                        h_sum1->Draw("HIST");
                        gaus1->Draw("SAME");
                        drawn900_sum = true;
                    } else {
                        h_sum1->Draw("HIST SAME");
                        gaus1->Draw("SAME");
                    }
                    l900_sum->AddEntry(h_sum1,
                                       Form("%s (mean=%.3f, #sigma=%.3f)", labPMT1.c_str(), mean1, sigma1),
                                       "l");
                    if (hasPMT3) {
                        h_sum3->Draw("HIST SAME");
                        l900_sum->AddEntry(h_sum3,
                                           Form("%s (mean=%.3f, #sigma=%.3f)", labPMT3.c_str(), meanSum3, sigmaSum3),
                                           "l");
                    }
                    if (hasPMT4) {
                        h_sum4->Draw("HIST SAME");
                        l900_sum->AddEntry(h_sum4,
                                           Form("%s (mean=%.3f, #sigma=%.3f)", labPMT4.c_str(), meanSum4, sigmaSum4),
                                           "l");
                    }
                }

                // keep pointers for writing later
                allHists.push_back(h_spec1);
                allHists.push_back(h_peak1);
                allHists.push_back(h_sum1);
                if (hasPMT3) { allHists.push_back(h_peak3); allHists.push_back(h_sum3); }
                if (hasPMT4) { allHists.push_back(h_peak4); allHists.push_back(h_sum4); }

                // write to text summary
                outTxt << pw << " " << ang << " " << PMT3_name << " " << PMT4_name << " "
                << sum << " " << vmin << " " << vmax << " "
                << mean1 << " " << sigma1 << " "
                << meanSum3 << " " << sigmaSum3 << " "
                << meanSum4 << " " << sigmaSum4 << "\n";

                f->Close();
                delete gaus1;

                std::cout << "Finished analysis for " << pw << " ps, " << ang
                << " deg, PMT3=" << PMT3_name
                << ", PMT4=" << PMT4_name << std::endl;
            } // iPMT
        } // iAngle
    } // iPW

    // ===================== SAVE CANVASES =====================

    // Draw legends & save PDFs for 0 ps
    c0_spec->cd();  l0_spec->Draw();
    c0_peak->cd();  l0_peak->Draw();
    c0_sum->cd();   l0_sum->Draw();

    c0_spec->SaveAs("../plots/summary/LaserBall_0ps_spectrum.pdf");
    c0_peak->SaveAs("../plots/summary/LaserBall_0ps_peaks.pdf");
    c0_sum->SaveAs("../plots/summary/LaserBall_0ps_Sum4Points.pdf");

    // Draw legends & save PDFs for 900 ps
    c900_spec->cd();  l900_spec->Draw();
    c900_peak->cd();  l900_peak->Draw();
    c900_sum->cd();   l900_sum->Draw();

    c900_spec->SaveAs("../plots/summary/LaserBall_900ps_spectrum.pdf");
    c900_peak->SaveAs("../plots/summary/LaserBall_900ps_peaks.pdf");
    c900_sum->SaveAs("../plots/summary/LaserBall_900ps_Sum4Points.pdf");

    // ===================== SAVE ROOT FILE =====================

    TFile* outFile = TFile::Open(outputRootName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: could not create output ROOT file: " << outputRootName << std::endl;
    } else {
        outFile->cd();
        c0_spec->Write("c0_spec");
        c0_peak->Write("c0_peak");
        c0_sum->Write("c0_sum");
        c900_spec->Write("c900_spec");
        c900_peak->Write("c900_peak");
        c900_sum->Write("c900_sum");

        for (auto* h : allHists) {
            if (h) h->Write();
        }
        outFile->Close();
        std::cout << "Wrote ROOT output to: " << outputRootName << std::endl;
    }

    outTxt.close();
    std::cout << "Done. Saved canvases and ROOT file." << std::endl;

    return 0;
}
