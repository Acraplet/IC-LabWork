#include <iostream>
#include <vector>
#include <string>
#include <limits>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include <fstream>

int main(int argc, char** argv) {
    // ------------ Configuration ------------
    // Pulse widths (ps)
    std::vector<int> pulseWidths;
    for (int pw = 760; pw <= 900; pw += 10) {
        pulseWidths.push_back(pw);
    }
    // pulseWidths.push_back(900);

    // Tree & branch names
    std::string treeNamePMT   = "PMT1data";
    std::string treeNameTrig  = "Trigger";
    std::string branchName    = "volts";

    // Histogram settings
    int    nBins   = 200;
    double xMin    = -0.2;
    double xMax    =  0.8;

    // Histograms for peaks and 4-sample sums
    int    nBinsPeak   = 200;
    double peakMin     = -0.2;   // adjust if needed
    double peakMax     =  0.8;   // adjust if needed

    int    nBinsSum    = 200;
    double sumMin      = -0.2;   // roughly 4 * single-sample range if needed
    double sumMax      =  2.5;   // adjust to your expected integral range

    // Output
    std::string outputCanvasName = "volts_overlay";
    std::string outputPdfName    = "../plots/summary/Monitor_PMT_volts_overlay_spectrum.pdf";
    std::string outputPdfNamePeak    = "../plots/summary/Monitor_PMT_volts_overlay_peaks.pdf";
    std::string outputPdfNameSum    = "../plots/summary/Monitor_PMT_volts_overlay_Sum4Points.pdf";
    std::string outputRootName   = "../results/Monitor_PMT_volts_overlay_spectrum.root";

    // ------------ ROOT Style ------------
    gStyle->SetOptStat(0);

    // ------------ Prepare canvas & legend ------------
    TCanvas* c = new TCanvas("c", "Volts overlay", 1200, 1200);
    c->SetGrid();

    TLegend* leg = new TLegend(0.55, 0.60, 0.98, 0.98);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);
    c->SetLogy();

    // ------------ Prepare canvas & legend ------------
    TCanvas* c1 = new TCanvas("c1", "Sum4Points", 1200, 1200);
    c1->SetGrid();
    c1->cd();

    TLegend* leg1 = new TLegend(0.55, 0.60, 0.98, 0.98);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.025);

    // ------------ Prepare canvas & legend ------------
    TCanvas* c2 = new TCanvas("c2", "PeakValues", 1200, 800);
    c2->SetGrid();
    c2->cd();

    TLegend* leg2 = new TLegend(0.45, 0.50, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.025);

    // ------------ Colors for histograms ------------
    std::vector<int> colors = {
        kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
        kCyan + 1, kOrange + 1, kViolet + 1, kAzure + 1,
        kPink + 1, kTeal + 1, kSpring + 1, kGray + 2,
        kBlack
    };

    // Store histograms to keep them alive
    std::vector<TH1D*> histos;
    std::vector<TH1D*> histosPeaks;
    std::vector<TH1D*> histosSums;

    bool firstHist = true;

    std::ofstream outTxt("../results/Monitor_PMT_volts_summary.txt");
    if (!outTxt.is_open()) {
        std::cerr << "Error: could not open volts_summary.txt for writing\n";
        return 1;
    }
    outTxt << "pulse_width_ps sum min max mean4points sigma4points\n";

    // ------------ Loop over pulse widths / files ------------
    for (size_t iPW = 0; iPW < pulseWidths.size(); ++iPW) {
        int pw = pulseWidths[iPW];

        // Build file name
        // std::string fileName =
        //     "/home/ac4317/Laptops/Year1/WCTE/LabWork/2025_characterisations/data/"
        //     "MonitorPMT/hadded_Monitor_PMT_results_" + std::to_string(pw) + "ps.root";

        std::string fileName =
            "../data/MonitorPMT/hadded_Monitor_PMT_and_OPM_results_" + std::to_string(pw) + "ps.root";

        // Check file existence
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

        // Get trees
        TTree* treePMT = dynamic_cast<TTree*>(f->Get(treeNamePMT.c_str()));
        TTree* treeTrig = dynamic_cast<TTree*>(f->Get(treeNameTrig.c_str()));
        if (!treePMT || !treeTrig) {
            std::cerr << "Error: Trees " << treeNamePMT << " or " << treeNameTrig
                      << " not found in file: " << fileName << std::endl;
            f->Close();
            continue;
        }

        // Check entry counts match
        Long64_t nEntriesPMT  = treePMT->GetEntries();
        Long64_t nEntriesTrig = treeTrig->GetEntries();
        if (nEntriesPMT != nEntriesTrig) {
            std::cerr << "Error: PMT1data and Trigger trees have different entry counts in "
                      << fileName << " (" << nEntriesPMT << " vs " << nEntriesTrig << ")\n";
            f->Close();
            continue;
        }
        Long64_t nEntries = nEntriesPMT;
        if (nEntries == 0) {
            std::cerr << "Warning: Trees have zero entries in " << fileName << std::endl;
            f->Close();
            continue;
        }

        // Set up branch addresses
        double voltsPMT = 0.0;
        double voltsTrig = 0.0;
        treePMT->SetBranchAddress(branchName.c_str(), &voltsPMT);
        treeTrig->SetBranchAddress(branchName.c_str(), &voltsTrig);

        // Create histogram for full PMT distribution
        std::string histName = "h_volts_" + std::to_string(pw) + "ps";
        TH1D* h = new TH1D(histName.c_str(),
                           ("Volts for " + std::to_string(pw) + " ps").c_str(),
                           nBins, xMin, xMax);
        h->SetDirectory(nullptr);
        int color = colors[iPW % colors.size()];
        h->SetLineColor(color);
        h->SetLineWidth(2);

        // Histograms for peak values and 4-sample sums
        std::string histPeakName = "h_peak_" + std::to_string(pw) + "ps";
        std::string histSumName  = "h_sum4_" + std::to_string(pw) + "ps";

        TH1D* hPeak = new TH1D(histPeakName.c_str(),
                               ("Peak volts for " + std::to_string(pw) + " ps").c_str(),
                               nBinsPeak, peakMin, peakMax);
        hPeak->SetDirectory(nullptr);
        hPeak->SetLineColor(color);
        hPeak->SetLineWidth(2);

        TH1D* hSum = new TH1D(histSumName.c_str(),
                              ("4-sample sum around peak for " + std::to_string(pw) + " ps").c_str(),
                              nBinsSum, sumMin, sumMax);
        hSum->SetDirectory(nullptr);
        hSum->SetLineColor(color);

        hSum->SetLineWidth(2);

        // Track min and max of PMT volts, and sum
        double vmin = std::numeric_limits<double>::infinity();
        double vmax = -std::numeric_limits<double>::infinity();
        double sum  = 0.0;

        // First pass: fill global PMT histogram and copy Trigger values into memory
        std::vector<double> trigVals(nEntries);
        for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
            treePMT->GetEntry(iEntry);
            treeTrig->GetEntry(iEntry);

            h->Fill(voltsPMT);
            sum += voltsPMT;
            if (voltsPMT < vmin) vmin = voltsPMT;
            if (voltsPMT > vmax) vmax = voltsPMT;

            trigVals[iEntry] = voltsTrig;
        }

        // ---- Identify trigger periods via rising edges in Trigger::volts ----
        // Simple threshold-based rising-edge detection; tune this as needed.
        // Example: assume low ~0, high > 0.5
        double thr = -0.6;
        std::vector<Long64_t> periodStarts;
        for (Long64_t iEntry = 1; iEntry < nEntries; ++iEntry) {
            double prev = trigVals[iEntry - 1];
            double curr = trigVals[iEntry];
            if (prev < thr && curr >= thr) {
                periodStarts.push_back(iEntry);
            }
        }

        // If we found N periods, define period ranges as [start[k], start[k+1]) and
        // last period [start[N-1], nEntries).
        for (size_t iPeriod = 0; iPeriod < periodStarts.size(); ++iPeriod) {
            Long64_t start = periodStarts[iPeriod];
            Long64_t end   = (iPeriod + 1 < periodStarts.size())
                           ? periodStarts[iPeriod + 1]
                           : nEntries;

            if (end <= start) continue;

            // Find peak PMT value in this period
            double peakVal = -std::numeric_limits<double>::infinity();
            Long64_t peakIdx = start;

            for (Long64_t iEntry = start; iEntry < end; ++iEntry) {
                treePMT->GetEntry(iEntry); // read PMT for this entry
                if (voltsPMT > peakVal) {
                    peakVal = voltsPMT;
                    peakIdx = iEntry;
                }
            }

            // Fill peak histogram
            hPeak->Fill(peakVal);

            // Compute sum: sample before, peak, and two after: (peakIdx-1, peakIdx, peakIdx+1, peakIdx+2)
            double sum4 = 0.0;
            // We need to be careful at boundaries
            Long64_t i0 = (peakIdx > 0) ? (peakIdx - 1) : peakIdx;
            Long64_t i1 = peakIdx;
            Long64_t i2 = (peakIdx + 1 < nEntries) ? (peakIdx + 1) : peakIdx;
            Long64_t i3 = (peakIdx + 2 < nEntries) ? (peakIdx + 2) : peakIdx;

            // Read and sum these indices
            double v0, v1, v2, v3;
            treePMT->GetEntry(i0); v0 = voltsPMT;
            treePMT->GetEntry(i1); v1 = voltsPMT;
            treePMT->GetEntry(i2); v2 = voltsPMT;
            treePMT->GetEntry(i3); v3 = voltsPMT;

            sum4 = v0 + v1 + v2 + v3;
            hSum->Fill(sum4);
        }

        TF1* gausFit = new TF1("fit", "gaus", sumMin, sumMax);
        gausFit->SetLineColor(kBlack);
        gausFit->SetLineStyle(2);  // dashed
        gausFit->SetLineWidth(2);

        // Optional: reasonable initial guesses (amplitude, mean, sigma)
        gausFit->SetParameters(hSum->GetMaximum(), hSum->GetMean(), hSum->GetRMS());

        // Do the fit quietly ("Q") and not draw automatically ("0")
        hSum->Fit(gausFit, "Q0");
        double amp     = gausFit->GetParameter(0);
        double mean    = gausFit->GetParameter(1);
        double sigma   = gausFit->GetParameter(2);


        // Draw the full PMT histogram on the main canvas
        if (firstHist) {
            c->cd();
            h->SetTitle("Monitor PMT volts distribution;volts [V];entries");
            h->Draw("HIST");
            firstHist = false;
            c1->cd();
            hPeak->SetTitle("Monitor PMT peak voltage;peak volts[V];entries");
            hPeak->Draw("HIST");
            c2->cd();
            hSum->SetTitle("Monitor PMT sum 4 points around peak;volts[V];entries");
            hSum->Draw("HIST");
            gausFit->Draw("SAME");
        } else {
            c->cd();
            h->Draw("HIST SAME");
            c1->cd();
            hPeak->Draw("HIST SAME");
            c2->cd();
            hSum->Draw("HIST SAME");
            gausFit->Draw("SAME");
        }

        // Legend entry (you can change what you report)
        c->cd();
        std::string legEntryText = Form("%d ps: sum=%.3g, max=%.3g",
                                        pw, sum, vmax);
        leg->AddEntry(h, legEntryText.c_str(), "l");


        c1->cd();
        // Legend entry (you can change what you report)
        std::string legEntryText1 = Form("%d ps",
                                        pw);
        leg1->AddEntry(hPeak, legEntryText1.c_str(), "l");


        c2->cd();
        std::string legEntryText2 = Form("%d ps: mean = %.3f V, sigma = %.3f V",
                                         pw, mean, sigma);
        leg2->AddEntry(hSum, legEntryText2.c_str(), "l");

        histos.push_back(h);
        histosPeaks.push_back(hPeak);
        histosSums.push_back(hSum);

        // write out the information
        outTxt << pw << " " << sum << " " << vmin << " " << vmax << " " << mean << " " << sigma << "\n";

        f->Close();

        std::cout << "Finished the analysis for run at " << pw << " ps" << std::endl;
    }

    // Draw legend
    c->cd();
    leg->Draw();
    c1->cd();
    leg1->Draw();
    c2->cd();
    leg2->Draw();

    // Save canvas (full PMT distributions)
    c->SaveAs(outputPdfName.c_str());
    c1->SaveAs(outputPdfNamePeak.c_str());
    c2->SaveAs(outputPdfNameSum.c_str());

    // Optionally, save histograms and canvas to a ROOT file
    TFile* outFile = TFile::Open(outputRootName.c_str(), "RECREATE");
    c->Write(outputCanvasName.c_str());
    for (auto* h : histos)      h->Write();
    for (auto* hP : histosPeaks) hP->Write();
    for (auto* hS : histosSums)  hS->Write();
    outFile->Close();

    outTxt.close();

    std::cout << "Done. Saved canvas to " << outputPdfName
              << " and objects to " << outputRootName << std::endl;

    return 0;
}
