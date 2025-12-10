// monitor_and_compare_pmts.cpp
// Build with: g++ -O2 -std=c++17 `root-config --cflags --libs` monitor_and_compare_pmts.cpp -o monitor_and_compare_pmts

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
    for (int pw = 900; pw <= 900; pw += 10) pulseWidths.push_back(pw);
    pulseWidths.push_back(0);

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

    // Histogram settings
    int    nBins   = 400;
    double xMin    = -0.05;
    double xMax    =  0.8;

    int    nBinsPeak   = 400;
    double peakMin     = -0.05;
    double peakMax     =  0.8;

    int    nBinsSum    = 200;
    double sumMin      = -0.05;
    double sumMax      =  1;

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    double sum  = 0.0;

    std::string outputCanvasName = "volts_overlay";
    std::string outputPdfName    = "../plots/summary/LaserBall_PMT1_Mon_PMT3_-30_PMT4_V_spectrum.pdf";
    std::string outputPdfNamePeak    = "../plots/summary/LaserBall_PMT1_Mon_PMT3_-30_PMT4_V_peaks.pdf";
    std::string outputPdfNameSum    = "../plots/summary/LaserBall_PMT1_Mon_PMT3_-30_PMT4_V_Sum4Points.pdf";
    std::string outputRootName   = "../results/LaserBall_PMT1_Mon_PMT3_-30_PMT4_V_spectrum.root";

    gStyle->SetOptStat(0);

    // Canvas / legends
    TCanvas* c = new TCanvas("c", "Volts overlay", 1200, 1200);
    c->SetGrid();
    c->SetLogy();
    TLegend* leg = new TLegend(0.45, 0.50, 0.98, 0.98);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);

    TCanvas* c1 = new TCanvas("c1", "Sum4Points", 1200, 1200);
    c1->SetGrid();
    TLegend* leg1 = new TLegend(0.55, 0.60, 0.98, 0.98);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.025);

    TCanvas* c2 = new TCanvas("c2", "PeakValues", 1200, 800);
    c2->SetGrid();
    TLegend* leg2 = new TLegend(0.45, 0.50, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.025);

    std::vector<int> colors = { kRed+1, kBlue+1, kGreen+2, kMagenta+1, kCyan+1, kOrange+1, kViolet+1, kAzure+1, kPink+1, kTeal+1, kSpring+1, kGray+2, kBlack };

    std::vector<TH1D*> histos;        // full PMT1 distributions
    std::vector<TH1D*> histosPeaks;   // peak histograms (PMT1,PMT3,PMT4 as present)
    std::vector<TH1D*> histosSums;    // sum histograms (PMT1,PMT3,PMT4 as present)

    bool firstHist = true;

    std::ofstream outTxt("../results/LaserBall_PMT_volts_summary.txt");
    if (!outTxt.is_open()) {
        std::cerr << "Error: could not open output summary for writing\n";
        return 1;
    }
    outTxt << "pulse_width_ps pw_sum1 min1 max1 mean1 sigma1 meanSum1 sigmaSum1"
    << " meanSum3 sigmaSum3 meanSum4 sigmaSum4\n";

    for (size_t iPW = 0; iPW < pulseWidths.size(); ++iPW) {
        for (size_t iAngle = 0; iAngle < angleLaserBall.size(); ++iAngle){
            for (size_t iPMT = 0; iPMT < PMT3_list.size(); ++iPMT){
            int pw = pulseWidths[iPW];
            int ang = angleLaserBall[iAngle];
            std::string PMT3_name = PMT3_list[iPMT];
            std::string PMT4_name = PMT4_list[iPMT];


            std::string fileName = "../../data/LaserBall/hadded_PMT1_Mon_PMT2_Trig_PMT3_"+PMT3_name+"_PMT4_"+PMT4_name+"_results_"+ std::to_string(pw)+"ps_"+std::to_string(ang)+"deg.root";

            std::cout << "Looking at file: " << fileName << std::endl;

            //"../data/MonitorPMT/hadded_Monitor_PMT_and_OPM_results_" + std::to_string(pw) + "ps.root";

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

            TTree* treePMT = dynamic_cast<TTree*>(f->Get(treeNamePMT.c_str()));
            TTree* treeTrig = dynamic_cast<TTree*>(f->Get(treeNameTrig.c_str()));
            TTree* treePMT3 = dynamic_cast<TTree*>(f->Get(treeNamePMT3.c_str()));
            TTree* treePMT4 = dynamic_cast<TTree*>(f->Get(treeNamePMT4.c_str()));

            if (!treePMT || !treeTrig) {
                std::cerr << "Error: Trees " << treeNamePMT << " or " << treeNameTrig << " not found in file: " << fileName << std::endl;
                f->Close();
                continue;
            }

            bool hasPMT3 = (treePMT3 != nullptr);
            bool hasPMT4 = (treePMT4 != nullptr);
            if (!hasPMT3) std::cerr << "Note: " << treeNamePMT3 << " not found, skipping PMT3 sums for " << pw << " ps\n";
            if (!hasPMT4) std::cerr << "Note: " << treeNamePMT4 << " not found, skipping PMT4 sums for " << pw << " ps\n";

            // Prefer to use the minimum entries available among all present trees to be safe
            Long64_t nEntriesTrig = treeTrig->GetEntries();
            Long64_t nEntriesPMT  = treePMT->GetEntries();
            Long64_t nEntriesPMT3 = hasPMT3 ? treePMT3->GetEntries() : nEntriesTrig;
            Long64_t nEntriesPMT4 = hasPMT4 ? treePMT4->GetEntries() : nEntriesTrig;

            if (nEntriesPMT != nEntriesTrig) {
                std::cerr << "Warning: PMT1 and Trigger tree entry counts differ (" << nEntriesPMT << " vs " << nEntriesTrig << "). Using min.\n";
            }
            Long64_t nEntries = std::min(std::min(nEntriesTrig, nEntriesPMT), std::min(nEntriesPMT3, nEntriesPMT4));
            if (nEntries == 0) {
                std::cerr << "Warning: no entries in trees for file " << fileName << std::endl;
                f->Close();
                continue;
            }

            // Prepare branch addresses and in-memory vectors
            double voltsPMT = 0.0, voltsTrig = 0.0, voltsPMT3 = 0.0, voltsPMT4 = 0.0;
            treePMT->SetBranchAddress(branchName.c_str(), &voltsPMT);
            treeTrig->SetBranchAddress(branchName.c_str(), &voltsTrig);
            if (hasPMT3) treePMT3->SetBranchAddress(branchName.c_str(), &voltsPMT3);
            if (hasPMT4) treePMT4->SetBranchAddress(branchName.c_str(), &voltsPMT4);

            std::vector<double> pmt1Vals(nEntries), trigVals(nEntries), pmt3Vals, pmt4Vals;
            if (hasPMT3) pmt3Vals.resize(nEntries);
            if (hasPMT4) pmt4Vals.resize(nEntries);


            std::string histName = "h_Volts1_" + std::to_string(pw) + "ps";
            std::string histPeakName = "h_Peaks1_" + std::to_string(pw) + "ps";
            std::string histSumName = "h_sum4_pmt1_" + std::to_string(pw) + "ps";

            // Create histograms (PMT1 full, peaks and sums)
            TH1D* h = new TH1D(histName.c_str(),
                            ("Volts for " + std::to_string(pw) + " ps").c_str(),
                            nBins, xMin, xMax);
            h->SetDirectory(nullptr);
            h->SetLineColor(colors[iPW % colors.size()]);
            h->SetLineWidth(2);

            TH1D* hPeak = new TH1D(histPeakName.c_str(),
                                ("Peak volts for PMT1 " + std::to_string(pw) + " ps").c_str(),
                                nBinsPeak, peakMin, peakMax);
            hPeak->SetDirectory(nullptr);
            hPeak->SetLineColor(colors[iPW % colors.size()]);
            hPeak->SetLineWidth(2);

            TH1D* hSum = new TH1D(histSumName.c_str(),
                                ("4-sample sum around PMT1 peak for " + std::to_string(pw) + " ps").c_str(),
                                nBinsSum, sumMin, sumMax);
            hSum->SetDirectory(nullptr);
            hSum->SetLineColor(colors[iPW % colors.size()]);
            hSum->SetLineWidth(2);

            // histograms for PMT3 and PMT4 sums/peaks if present
            TH1D* hPeak3 = nullptr; TH1D* hSum3 = nullptr;
            TH1D* hPeak4 = nullptr; TH1D* hSum4 = nullptr;
            if (hasPMT3) {
                std::string histPeakName3 = "h_peak3_" + std::to_string(pw) + "ps";
                std::string histSumName3  = "h_sum4_pmt3_" + std::to_string(pw) + "ps";
                hPeak3 = new TH1D(histPeakName3.c_str(),
                                ("Peak volts for PMT3 " + std::to_string(pw) + " ps").c_str(),
                                nBinsPeak, peakMin, peakMax);
                hPeak3->SetDirectory(nullptr);
                hPeak3->SetLineColor(kGreen+2);
                hPeak3->SetLineWidth(2);
                hSum3 = new TH1D(histSumName3.c_str(),
                                ("4-sample sum around PMT3 peak for " + std::to_string(pw) + " ps").c_str(),
                                nBinsSum, sumMin, sumMax);
                hSum3->SetDirectory(nullptr);
                hSum3->SetLineColor(kGreen+2);
                hSum3->SetLineWidth(2);
            }
            if (hasPMT4) {
                std::string histPeakName4 = "h_peak4_" + std::to_string(pw) + "ps";
                std::string histSumName4  = "h_sum4_pmt4_" + std::to_string(pw) + "ps";
                hPeak4 = new TH1D(histPeakName4.c_str(),
                                ("Peak volts for PMT4 " + std::to_string(pw) + " ps").c_str(),
                                nBinsPeak, peakMin, peakMax);
                hPeak4->SetDirectory(nullptr);
                hPeak4->SetLineColor(kMagenta+1);
                hPeak4->SetLineWidth(2);
                hSum4 = new TH1D(histSumName4.c_str(),
                                ("4-sample sum around PMT4 peak for " + std::to_string(pw) + " ps").c_str(),
                                nBinsSum, sumMin, sumMax);
                hSum4->SetDirectory(nullptr);
                hSum4->SetLineColor(kMagenta+1);
                hSum4->SetLineWidth(2);
            }

            // Fill the in-memory arrays
            for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
                treePMT->GetEntry(iEntry);
                treeTrig->GetEntry(iEntry);
                pmt1Vals[iEntry] = voltsPMT;
                trigVals[iEntry] = voltsTrig;
                if (hasPMT3) {
                    treePMT3->GetEntry(iEntry);
                    pmt3Vals[iEntry] = voltsPMT3;
                }
                if (hasPMT4) {
                    treePMT4->GetEntry(iEntry);
                    pmt4Vals[iEntry] = voltsPMT4;
                }
                // Fill full PMT1 distribution
                h->Fill(voltsPMT);
                if (voltsPMT < vmin) vmin = voltsPMT;
                if (voltsPMT > vmax) vmax = voltsPMT;
                sum += voltsPMT;
            }

            // ---- Identify trigger periods via rising edges in Trigger::volts ----
            double thr = -0.6; // tune if necessary
            std::vector<Long64_t> periodStarts;
            for (Long64_t iEntry = 1; iEntry < nEntries; ++iEntry) {
                double prev = trigVals[iEntry - 1];
                double curr = trigVals[iEntry];
                if (prev < thr && curr >= thr) {
                    periodStarts.push_back(iEntry);
                }
            }

            // For optional comparison: compute PMT3/4 sums at PMT1 peak index
            bool compareOtherAtPMT1Peak = true;

            // Iterate over trigger periods and compute peaks & 4-sample sums
            for (size_t iPeriod = 0; iPeriod < periodStarts.size(); ++iPeriod) {
                Long64_t start = periodStarts[iPeriod];
                Long64_t end   = (iPeriod + 1 < periodStarts.size()) ? periodStarts[iPeriod + 1] : nEntries;
                if (end <= start) continue;

                // find peak for PMT1 in period
                double peakVal1 = -std::numeric_limits<double>::infinity();
                Long64_t peakIdx1 = start;
                for (Long64_t iEntry = start; iEntry < end; ++iEntry) {
                    if (pmt1Vals[iEntry] > peakVal1) {
                        peakVal1 = pmt1Vals[iEntry];
                        peakIdx1 = iEntry;
                    }
                }
                hPeak->Fill(peakVal1);

                // compute 4-sample sum around PMT1 peak
                auto computeSum4 = [&](const std::vector<double>& arr, Long64_t idx)->double {
                    Long64_t i0 = (idx > 0) ? (idx - 1) : idx;
                    Long64_t i1 = idx;
                    Long64_t i2 = (idx + 1 < nEntries) ? (idx + 1) : idx;
                    Long64_t i3 = (idx + 2 < nEntries) ? (idx + 2) : idx;
                    return arr[i0] + arr[i1] + arr[i2] + arr[i3];
                };

                double sum4_1 = computeSum4(pmt1Vals, peakIdx1);
                hSum->Fill(sum4_1);

                // PMT3: find its own peak and sum, and optionally sum at PMT1 peak
                if (hasPMT3) {
                    double peakVal3 = -std::numeric_limits<double>::infinity();
                    Long64_t peakIdx3 = start;
                    for (Long64_t iEntry = start; iEntry < end; ++iEntry) {
                        if (pmt3Vals[iEntry] > peakVal3) {
                            peakVal3 = pmt3Vals[iEntry];
                            peakIdx3 = iEntry;
                        }
                    }
                    hPeak3->Fill(peakVal3);
                    double sum4_3 = computeSum4(pmt3Vals, peakIdx3);
                    hSum3->Fill(sum4_3);
                    if (compareOtherAtPMT1Peak) {
                        double sum4_3_at1 = computeSum4(pmt3Vals, peakIdx1);
                        // optionally store into another histogram if needed
                        // For now, also fill hSum3 (mixing might be confusing) -- better to create a separate hist if desired.
                    }
                }

                // PMT4
                if (hasPMT4) {
                    double peakVal4 = -std::numeric_limits<double>::infinity();
                    Long64_t peakIdx4 = start;
                    for (Long64_t iEntry = start; iEntry < end; ++iEntry) {
                        if (pmt4Vals[iEntry] > peakVal4) {
                            peakVal4 = pmt4Vals[iEntry];
                            peakIdx4 = iEntry;
                        }
                    }
                    hPeak4->Fill(peakVal4);
                    double sum4_4 = computeSum4(pmt4Vals, peakIdx4);
                    hSum4->Fill(sum4_4);
                    if (compareOtherAtPMT1Peak) {
                        double sum4_4_at1 = computeSum4(pmt4Vals, peakIdx1);
                        // as above, optionally store separately
                    }
                }
            } // end periods

            // Fit PMT1 sum histogram with gaussian to extract mean/sigma (quiet)
            TF1* gausFit1 = new TF1(("gaus1_" + std::to_string(pw)).c_str(), "gaus", sumMin, sumMax);
            gausFit1->SetLineColor(kBlack);
            gausFit1->SetLineStyle(2);
            gausFit1->SetLineWidth(2);
            gausFit1->SetParameters(hSum->GetMaximum(), hSum->GetMean(), hSum->GetRMS());
            hSum->Fit(gausFit1, "Q0");
            double mean1  = gausFit1->GetParameter(1);
            double sigma1 = gausFit1->GetParameter(2);

            double meanSum3 = 0.0, sigmaSum3 = 0.0, meanSum4 = 0.0, sigmaSum4 = 0.0;
            if (hasPMT3) {
                // try fit but only if entries > a small number
                if (hSum3->GetEntries() > 10) {
                    TF1* gaus3 = new TF1(("gaus3_" + std::to_string(pw)).c_str(), "gaus", sumMin, sumMax);
                    gaus3->SetParameters(hSum3->GetMaximum(), hSum3->GetMean(), hSum3->GetRMS());
                    hSum3->Fit(gaus3, "Q0");
                    meanSum3  = gaus3->GetParameter(1);
                    sigmaSum3 = gaus3->GetParameter(2);
                    delete gaus3;
                } else {
                    meanSum3 = hSum3->GetMean();
                    sigmaSum3 = hSum3->GetRMS();
                }
            }
            if (hasPMT4) {
                if (hSum4->GetEntries() > 10) {
                    TF1* gaus4 = new TF1(("gaus4_" + std::to_string(pw)).c_str(), "gaus", sumMin, sumMax);
                    gaus4->SetParameters(hSum4->GetMaximum(), hSum4->GetMean(), hSum4->GetRMS());
                    hSum4->Fit(gaus4, "Q0");
                    meanSum4  = gaus4->GetParameter(1);
                    sigmaSum4 = gaus4->GetParameter(2);
                    delete gaus4;
                } else {
                    meanSum4 = hSum4->GetMean();
                    sigmaSum4 = hSum4->GetRMS();
                }
            }

            // Draw histograms (first vs same)
            if (firstHist) {
                c->cd();
                h->SetTitle("Monitor PMT volts distribution;volts [V];entries");
                h->Draw("HIST");
                firstHist = false;

                c1->cd();
                hSum3->SetTitle("4-sample sums (PMT1/PMT3/PMT4);volts[V];entries");
                // hSum->Draw("HIST");
                if (hasPMT3) { hSum3->Draw("HIST SAME"); }
                if (hasPMT4) { hSum4->Draw("HIST SAME"); }

                c2->cd();
                hPeak->SetTitle("Peak voltages;peak volts[V];entries");
                hPeak->Draw("HIST");
                if (hasPMT3) hPeak3->Draw("HIST SAME");
                if (hasPMT4) hPeak4->Draw("HIST SAME");
            } else {
                c->cd(); h->Draw("HIST SAME");
                c1->cd(); if (hasPMT3) hSum3->Draw("HIST SAME"); if (hasPMT4) hSum4->Draw("HIST SAME");
                c2->cd(); hPeak->Draw("HIST SAME"); if (hasPMT3) hPeak3->Draw("HIST SAME"); if (hasPMT4) hPeak4->Draw("HIST SAME");
            }

            // Legends
            c->cd();
            std::string legEntryText = Form("%d ps: sum=%.3g, max=%.3g", pw, sum, vmax);
            leg->AddEntry(h, legEntryText.c_str(), "l");

            c1->cd();
            leg1->AddEntry(hSum, Form("%d ps: PMT1 mean=%.3f sigma=%.3f", pw, mean1, sigma1), "l");
            if (hasPMT3) leg1->AddEntry(hSum3, Form("%d ps: PMT3 mean=%.3f sigma=%.3f", pw, meanSum3, sigmaSum3), "l");
            if (hasPMT4) leg1->AddEntry(hSum4, Form("%d ps: PMT4 mean=%.3f sigma=%.3f", pw, meanSum4, sigmaSum4), "l");

            c2->cd();
            leg2->AddEntry(hPeak, Form("%d ps", pw), "l");

            // Store hist pointers for writing later
            histos.push_back(h);
            histosPeaks.push_back(hPeak);
            histosSums.push_back(hSum);
            if (hasPMT3) { histosPeaks.push_back(hPeak3); histosSums.push_back(hSum3); }
            if (hasPMT4) { histosPeaks.push_back(hPeak4); histosSums.push_back(hSum4); }

            // write out the information line
            outTxt << pw << " " << sum << " " << vmin << " " << vmax << " "
            << mean1 << " " << sigma1 << " "
            << meanSum3 << " " << sigmaSum3 << " "
            << meanSum4 << " " << sigmaSum4 << "\n";

            f->Close();
            delete gausFit1;

            std::cout << "Finished the analysis for run at " << pw << " ps" << std::endl;
            } // end of PMT config loop
        }//end laser ball position loop
    } // end pulse width loop

    // Draw and save legends and canvases
    // Draw and save legends and canvases
    c->cd(); leg->Draw();
    c1->cd(); leg1->Draw();
    c2->cd(); leg2->Draw();

    // Save canvases to PDF files
    c->SaveAs(outputPdfName.c_str());
    c1->SaveAs(outputPdfNamePeak.c_str());
    c2->SaveAs(outputPdfNameSum.c_str());

    // Optionally, save histograms and canvases to a ROOT file
    TFile* outFile = TFile::Open(outputRootName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: could not create output ROOT file: " << outputRootName << std::endl;
    } else {
        outFile->cd();
        c->Write(outputCanvasName.c_str());
        c1->Write("Sum4Points_canvas");
        c2->Write("PeakValues_canvas");

        // write all histograms we collected
        for (auto* h : histos) {
            if (h) h->Write();
        }
        for (auto* hP : histosPeaks) {
            if (hP) hP->Write();
        }
        for (auto* hS : histosSums) {
            if (hS) hS->Write();
        }

        outFile->Close();
        std::cout << "Wrote ROOT output to: " << outputRootName << std::endl;
    }

    outTxt.close();

    std::cout << "Done. Saved canvases to PDFs and histograms to ROOT file." << std::endl;

    // Clean up: delete canvases and legends (optional; OS will reclaim on exit)
    delete c;
    delete c1;
    delete c2;
    delete leg;
    delete leg1;
    delete leg2;

    // Note: we intentionally did not delete the histograms written to file because
    // ROOT may have taken ownership when writing. If you prefer explicit deletion:
    // for (auto* h : histos) delete h;
    // for (auto* hP : histosPeaks) delete hP;
    // for (auto* hS : histosSums) delete hS;

    return 0;
}

