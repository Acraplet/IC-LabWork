// monitor_and_compare_pmts_mod.cpp
//This is a simple code to extract the gain of each of the PMTs
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
    pulseWidths.push_back(0);

    std::vector<int> angleLaserBall;
    angleLaserBall.push_back(0);

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
    //only look at the integrated sum 
    int    nBinsSum    = 50;   
    double sumMin      = -0.05;
    double sumMax      =  0.1;

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    double sum  = 0.0;

    // Output names kept generic – you can adapt to your convention
    std::string outputRootName   = "../results/LaserBall_2025_gain_calibration.root";

    gStyle->SetOptStat(0);

    // ------------ Colors (fixed mapping by PMT, not by loop index) ------------
    // *** CHANGED: explicit mapping so colors are guaranteed different
    int colPMT1 = kRed + 1;

    std::vector<int> PMT3_cols;
    PMT3_cols.push_back(kBlue+1);
    PMT3_cols.push_back(kRed+2);

    std::vector<int> PMT4_cols;
    PMT4_cols.push_back(kGreen+1);
    PMT4_cols.push_back(kRed+4);

    //The order here is PMT3 (-30) PMT3(H) PMT4(V) PMT4(+30)
    std::vector<double> min_mean = {0.03, 0.025, 3.31470e-02,  0.045};
    std::vector<double> max_mean = {0.05, 0.045, 3.51470e-02, 0.055};
    std::vector<double> min_std  = {9.7e-03, 9.1e-03,  6.8706e-03, 1.20936e-02};
    std::vector<double> max_std  = {9.8e-03,  9.3e-03, 8.06706e-03, 1.5936e-02};

    // ------------ Separate canvases for 0 ps and 900 ps ------------
    // 0 ps
    TCanvas* c0_sum   = new TCanvas("c1_sum",   "Sum4 0 ps",      1200, 800);
    c0_sum->SetGrid(); c0_sum->SetLogy();

    TCanvas* c1_sum   = new TCanvas("c2_sum",   "Sum4 0 ps",      1200, 800);
    c1_sum->SetGrid(); c1_sum->SetLogy();
    

    std::vector<TCanvas*> canvases_0;
    canvases_0.push_back(c0_sum);
    canvases_0.push_back(c1_sum);

    TCanvas* c2_sum   = new TCanvas("c3_sum",   "Sum4 0 ps",      1200, 800);
    c2_sum->SetGrid(); c2_sum->SetLogy();

    TCanvas* c3_sum   = new TCanvas("c4_sum",   "Sum4 0 ps",      1200, 800);
    c3_sum->SetGrid(); c3_sum->SetLogy();
    
    std::vector<TCanvas*> canvases_1;
    canvases_1.push_back(c2_sum);
    canvases_1.push_back(c3_sum);

    TLegend* l0_sum   = new TLegend(0.25, 0.50, 0.98, 0.98);
    l0_sum->SetBorderSize(0);  l0_sum->SetFillStyle(0);  l0_sum->SetTextSize(0.025);


    TLegend* l1_sum   = new TLegend(0.25, 0.50, 0.98, 0.98);
    l1_sum->SetBorderSize(0);  l1_sum->SetFillStyle(0);  l1_sum->SetTextSize(0.025);

    TLegend* l2_sum   = new TLegend(0.25, 0.50, 0.98, 0.98);
    l2_sum->SetBorderSize(0);  l2_sum->SetFillStyle(0);  l2_sum->SetTextSize(0.025);

    TLegend* l3_sum   = new TLegend(0.25, 0.50, 0.98, 0.98);
    l3_sum->SetBorderSize(0);  l3_sum->SetFillStyle(0);  l3_sum->SetTextSize(0.025);

    std::vector<TLegend*> legends_0;
    std::vector<TLegend*> legends_1;

    legends_0.push_back(l0_sum);
    legends_0.push_back(l1_sum);
    legends_1.push_back(l2_sum);
    legends_1.push_back(l3_sum);



    // Internal flags to know if first histogram was drawn on a given canvas
    bool drawn0_sum  = false;

    std::vector<TH1D*> allHists; // for writing to ROOT later
    std::vector<TF1*>  allFunc; // for writing to ROOT later

    std::ofstream outTxt("../results/LaserBall_PMT_gain_summary.txt");
    if (!outTxt.is_open()) {
        std::cerr << "Error: could not open output summary for writing\n";
        return 1;
    }
    outTxt << "pw_ps angle_deg PMT1_pos PMT3_pos PMT4_pos meanSum1 sigmaSum1 "
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
                "/vols/hyperk/users/jnugent/540_data/LaserBall/LaserBall/hadded_PMT1_Mon_PMT2_Trig_PMT3_" + PMT3_name +
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


                Long64_t nEntriesTrig = treeTrig->GetEntries();
                Long64_t nEntriesPMT  = treePMT->GetEntries();
                Long64_t nEntriesPMT3 = treePMT3->GetEntries();
                Long64_t nEntriesPMT4 = treePMT4->GetEntries();

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
                treePMT3->SetBranchAddress(branchName.c_str(), &voltsPMT3);
                treePMT4->SetBranchAddress(branchName.c_str(), &voltsPMT4);

                std::vector<double> pmt1Vals(nEntries), trigVals(nEntries), pmt3Vals, pmt4Vals;
                pmt3Vals.resize(nEntries);
                pmt4Vals.resize(nEntries);

                // ---------- Histograms for this configuration ----------
                std::string tag = Form("_%dps_%ddeg_P3%s_P4%s",
                                       pw, ang, PMT3_name.c_str(), PMT4_name.c_str());

                TH1D* h_sum1 = new TH1D(("h_sum4_PMT1" + tag).c_str(),
                                        ("PMT1 sum4 " + tag).c_str(),
                                        nBinsSum, sumMin, sumMax);

                h_sum1->SetDirectory(nullptr);
                h_sum1->SetLineColor(colPMT1+iPMT);
                h_sum1->SetLineWidth(2);

                TH1D* h_sum3 = nullptr;
                TH1D* h_sum4 = nullptr;

                h_sum3 = new TH1D(("h_sum4_PMT3" + tag).c_str(),
                                  ("PMT3 sum4 " + tag).c_str(),
                                  nBinsSum, sumMin, sumMax);
                h_sum3->SetDirectory(nullptr);
                h_sum3->SetLineColor(PMT3_cols[iPMT]);
                h_sum3->SetLineWidth(2);

                h_sum4 = new TH1D(("h_sum4_PMT4" + tag).c_str(),
                                  ("PMT4 sum4 " + tag).c_str(),
                                  nBinsSum, sumMin, sumMax);
                h_sum4->SetDirectory(nullptr);
                h_sum4->SetLineColor(PMT4_cols[iPMT]);
                h_sum4->SetLineWidth(2);

                // Fill in-memory arrays + PMT1 spectrum
                for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
                    treePMT->GetEntry(iEntry);
                    treeTrig->GetEntry(iEntry);
                    treePMT3->GetEntry(iEntry); 
                    treePMT4->GetEntry(iEntry); 
                    
		    pmt1Vals[iEntry] = voltsPMT;
                    trigVals[iEntry] = voltsTrig;
		    pmt3Vals[iEntry] = voltsPMT3;
		    pmt4Vals[iEntry] = voltsPMT4;

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
                    //Simple function that returns the sum of the four points around the peak,
                    //can modify to only return the peak value or do a fit and give that value
                    Long64_t i0 = (idx > 0) ? (idx - 1) : idx;
                    Long64_t i1 = idx;
                    Long64_t i2 = (idx + 1 < nEntries) ? (idx + 1) : idx;
                    Long64_t i3 = (idx + 2 < nEntries) ? (idx + 2) : idx;
                    return arr[i0] + arr[i1] + arr[i2] + arr[i3];
                };

                for (size_t ip = 0; ip < periodStarts.size(); ++ip) {
                    //run over all periods of the laser, finding the peaks and integrating them
                    Long64_t start = periodStarts[ip];
                    Long64_t end   = (ip + 1 < periodStarts.size()) ? periodStarts[ip+1] : nEntries;
                    if (end <= start) continue;

                    double peak1 = -1e9; double peak3 = -1e9; double peak4 = -1e9; 
                    Long64_t idx1 = start; Long64_t idx3 = start; Long64_t idx4 = start;
                    for (Long64_t i = start; i < end; ++i) {
                        if (pmt1Vals[i] > peak1) { peak1 = pmt1Vals[i]; idx1 = i; }
                        if (pmt3Vals[i] > peak3) { peak3 = pmt3Vals[i]; idx3 = i; }
                        if (pmt4Vals[i] > peak4) { peak4 = pmt4Vals[i]; idx4 = i; }
                    }
                    h_sum1->Fill(computeSum4(pmt1Vals, idx1));
                    h_sum3->Fill(computeSum4(pmt3Vals, idx3));
                    h_sum4->Fill(computeSum4(pmt4Vals, idx4));
                }

                // Fit PMT1 sum fitting the sum with a gaussian might not be optimal
                // check what i was doing in the previous years.
                TF1* gaus1 = new TF1(("gaus1" + tag).c_str(), "gaus", sumMin, sumMax);


                TF1* doublegaus3 = new TF1("doubleGaus", "gaus(0)+gaus(3)", sumMin, sumMax);
                doublegaus3->SetLineColor(PMT3_cols[iPMT]);
                doublegaus3->SetLineStyle(2);
                doublegaus3->SetLineWidth(2);


                TF1* doublegaus4 = new TF1("doubleGaus4", "gaus(0)+gaus(3)", sumMin, sumMax);
                doublegaus4->SetLineColor(PMT4_cols[iPMT]);
                doublegaus4->SetLineStyle(2);
                doublegaus4->SetLineWidth(2);


                //TF1* gaus3 = new TF1(("gaus3" + tag).c_str(), "gaus", sumMin, sumMax);
                //TF1* gaus4 = new TF1(("gaus4" + tag).c_str(), "gaus", sumMin, sumMax);
                gaus1->SetParameters(h_sum1->GetMaximum(), h_sum1->GetMean(), h_sum1->GetRMS());
                gaus1->SetLineColor(kBlack);
                gaus1->SetLineStyle(2);
                gaus1->SetLineWidth(2);
                h_sum1->Fit(gaus1, "Q0");
                
		        std::cout << "Doing the fits" << std::endl;
		        double mean1  = gaus1->GetParameter(1);
                double sigma1 = gaus1->GetParameter(2);

                double ped_mean3 = 0.005;
                double onepe_mean3 = (min_mean[iPMT] + max_mean[iPMT])/2;
                double onepe_std3 = (min_std[iPMT] + max_std[iPMT])/2;

                doublegaus3->SetParameters(h_sum3->GetMaximum(), ped_mean3, h_sum3->GetRMS(), h_sum3->GetMaximum()/100., onepe_mean3, onepe_std3);

                doublegaus3->SetParLimits(4, min_mean[iPMT], max_mean[iPMT]);
                doublegaus3->SetParLimits(5, min_std[iPMT], max_std[iPMT]);
                h_sum3->Fit(doublegaus3, "Q0");
                double meanSum3     = doublegaus3->GetParameter(1);
                double sigmaSum3    = doublegaus3->GetParameter(2);
                double mean1peSum3  = doublegaus3->GetParameter(4);
                double sigma1peSum3  = doublegaus3->GetParameter(5);
                
                double ped_mean4 = 0.013;
                double onepe_mean4 = (min_mean[iPMT+2] + max_mean[iPMT+2])/2;
                double onepe_std4 = (min_std[iPMT+2] + max_std[iPMT+2])/2;
                doublegaus4->SetParameters(h_sum4->GetMaximum(), ped_mean4, h_sum4->GetRMS(), h_sum4->GetMaximum()/100., onepe_mean4, onepe_std4);

//                gaus4->SetParameters(h_sum4->GetMaximum(), h_sum4->GetMean(), h_sum4->GetRMS());
                doublegaus4->SetParLimits(4, min_mean[iPMT+2], max_mean[iPMT+2]);
                doublegaus4->SetParLimits(5, min_std[iPMT+2], max_std[iPMT+2]);
                h_sum4->Fit(doublegaus4, "Q0");
                double meanSum4     = doublegaus3->GetParameter(1);
                double sigmaSum4    = doublegaus4->GetParameter(2);
                double mean1peSum4  = doublegaus4->GetParameter(4);
                double sigma1peSum4  = doublegaus4->GetParameter(5);

                // ===================== DRAWING / LEGENDS =====================

                // Human-readable laser ball position string
                // (includes both angle and delay as you requested)
                std::string laserPosStr = Form("Dark rate: ");

                // Legend labels: include PMT name + laser ball position
                std::string labPMT1 = "Monitor PMT " + laserPosStr;
                std::string labPMT3 = laserPosStr + PMT3_name + " PMT";
                std::string labPMT4 = laserPosStr + PMT4_name + " PMT";

		        std::cout << "Finished the fits" << std::endl;

                // Sums
                canvases_0[iPMT]->cd();
                
                if (!drawn0_sum) {
                    h_sum3->SetTitle("Sum4 at 0 ps;sum4 [V];entries");
                    h_sum3->Draw("HIST");
                    doublegaus3->Draw("SAME");
                    drawn0_sum = true;
                } else {
                    h_sum3->Draw("HIST SAME");
                    doublegaus3->Draw("SAME");
                }
                legends_0[iPMT]->AddEntry(h_sum3,
                                 Form("%s (ped mean=%.3f, ped #sigma=%.3f, 1PE mean=%.3f, 1PE #sigma=%.3f)", labPMT3.c_str(), meanSum3, sigmaSum3, mean1peSum3, sigma1peSum3),
                                 "l");

                canvases_1[iPMT]->cd();
                h_sum4->Draw("HIST SAME");
                legends_1[iPMT]->AddEntry(h_sum4,
                                 Form("%s (ped mean=%.3f, ped #sigma=%.3f, 1PE mean=%.3f, 1PE #sigma=%.3f)", labPMT4.c_str(), meanSum4, sigmaSum4, mean1peSum4, sigma1peSum4),
                                 "l");
                doublegaus4->Draw("SAME");


                allHists.push_back(h_sum1);
                allHists.push_back(h_sum3); 
                allHists.push_back(h_sum4); 

                // write to text summary
                outTxt << pw << " " << ang << " Mon " << PMT3_name << " " << PMT4_name << " "
                << mean1 << " " << sigma1 << " "
                << meanSum3 << " " << sigmaSum3 << " "
                << meanSum4 << " " << sigmaSum4 << "\n";

                f->Close();
//                delete gaus1;
//                delete doublegaus3;
//                delete gaus4;

                std::cout << "Finished analysis for " << pw << " ps, " << ang
                << " deg, PMT3=" << PMT3_name
                << ", PMT4=" << PMT4_name << std::endl;
            } // iPMT
        } // iAngle
    } // iPW

    // ===================== SAVE CANVASES =====================

    for (int i=0; i<2; i++){
        canvases_0[i]->cd();   legends_0[i]->Draw();
        canvases_0[i]->SaveAs(Form("../plots/summary/LaserBall_darkrate_0ps_Sum4Points_%c.pdf", PMT3_list[i]));
        canvases_1[i]->cd();   legends_1[i]->Draw();
        canvases_1[i]->SaveAs(Form("../plots/summary/LaserBall_darkrate_0ps_Sum4Points_%c.pdf", PMT4_list[i]));
    }


    // ===================== SAVE ROOT FILE =====================

    TFile* outFile = TFile::Open(outputRootName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: could not create output ROOT file: " << outputRootName << std::endl;
    } else {
        outFile->cd();

        // Draw legends & save PDFs for 0 ps



        canvases_0[0]->Write(Form("c1_sum"));
        canvases_0[1]->Write(Form("c2_sum"));
        canvases_1[0]->Write(Form("c3_sum"));
        canvases_1[1]->Write(Form("c4_sum"));


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
