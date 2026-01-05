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
    // pulseWidths.push_back(0);
    pulseWidths.push_back(0);

    std::vector<int> angleLaserBall;
    angleLaserBall.push_back(0);
    // angleLaserBall.push_back(90);

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
    int    nBinsSum    = 55;   // was 200
    double sumMin      = -0.5;
    double sumMax      =  3;

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    double sum  = 0.0;

    // Output names kept generic – you can adapt to your convention
    std::string outputRootName   = Form("../results/LaserBall_PMT1_Mon_PMT3_PMT4_0ps_%ips_%ideg_spectrum.root", pulseWidths[0], angleLaserBall[0]);

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

    //Saving the PMT gain for each of the PMTs so we can make the conversion, we have integrated the peak which should be ok
    //-30, H, V, +30
    std::vector<double> gains_per_PMT = {3.14e-2, 2.60e-2, 1.99e-2, 3.37e-2};
    std::vector<double> pedestal_pos_per_PMT = {7.375e-3, 6.887e-3, 1.320e-2, 1.291e-2};

    // std::vector<double> gains_per_PMT = {2.60e-2, 3.14e-2, 2.58e-2, 3.97e-2};
    // std::vector<double> pedestal_pos_per_PMT = {6.887e-3, 7.375e-3, 7.375e-3, 6.887e-3};

    // ------------ Separate canvases for 0 ps and 900 ps ------------
    




    TCanvas* c0_sum   = new TCanvas("c0_sum",   Form("Sum 4 points %ips", pulseWidths[0]),      1200, 800);
    TCanvas* c1_sum   = new TCanvas("c1_sum",   Form("Sum 4 points %ips", pulseWidths[0]),      1200, 800);
    TCanvas* c2_sum   = new TCanvas("c2_sum",   Form("Sum 4 points %ips", pulseWidths[0]),      1200, 800);
    TCanvas* c3_sum   = new TCanvas("c3_sum",   Form("Sum 4 points %ips", pulseWidths[0]),      1200, 800);

    c0_sum->SetGrid();
    c1_sum->SetGrid();
    c2_sum->SetGrid();
    c3_sum->SetGrid();

    std::vector<TCanvas*> canvases_0;
    canvases_0.push_back(c0_sum);
    canvases_0.push_back(c1_sum);

     std::vector<TCanvas*> canvases_1;
    canvases_1.push_back(c2_sum);
    canvases_1.push_back(c3_sum);


    TLegend* l0_sum   = new TLegend(0.15, 0.50, 0.98, 0.98);
    l0_sum->SetBorderSize(0);  l0_sum->SetFillStyle(0);  l0_sum->SetTextSize(0.025);


    TLegend* l1_sum   = new TLegend(0.15, 0.50, 0.98, 0.98);
    l1_sum->SetBorderSize(0);  l1_sum->SetFillStyle(0);  l1_sum->SetTextSize(0.025);

    TLegend* l2_sum   = new TLegend(0.15, 0.50, 0.98, 0.98);
    l2_sum->SetBorderSize(0);  l2_sum->SetFillStyle(0);  l2_sum->SetTextSize(0.025);

    TLegend* l3_sum   = new TLegend(0.15, 0.50, 0.98, 0.98);
    l3_sum->SetBorderSize(0);  l3_sum->SetFillStyle(0);  l3_sum->SetTextSize(0.025);

    std::vector<TLegend*> legends_0;
    std::vector<TLegend*> legends_1;

    legends_0.push_back(l0_sum);
    legends_0.push_back(l1_sum);
    legends_1.push_back(l2_sum);
    legends_1.push_back(l3_sum);


    std::vector<TH1D*> allHists; // for writing to ROOT later
    std::vector<TF1*>  allFunc;

    std::ofstream outTxt(Form("../results/LaserBall_values_%i.txt", pulseWidths[0]));
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
                int nWindows = 0;
                std::string PMT3_name = PMT3_list[iPMT];
                std::string PMT4_name = PMT4_list[iPMT];

                std::string fileName =
                "/vols/hyperk/users/jnugent/540_data/LaserBall/LaserBall/hadded_PMT1_Mon_PMT2_Trig_PMT3_" + PMT3_name +
                "_PMT4_" + PMT4_name + "_results_" +
                std::to_string(pw) + "ps_" + std::to_string(ang) + "deg.root";

                if (pw == 900){fileName =
                "/vols/hyperk/users/jnugent/540_data/LaserBall/91225_data/hadded_PMT1_Mon_PMT2_Trig_PMT3_" + PMT3_name +
                "_PMT4_" + PMT4_name + "_results_" +
                std::to_string(pw) + "ps_" + std::to_string(ang) + "deg.root";}
                

                

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
                h_sum1->SetLineColor(colPMT1+iPMT+iAngle);
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
                        nWindows += 1;
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
                    double peak1 = -1e9; double peak3 = -1e9; double peak4 = -1e9; 
                    Long64_t idx1 = start; Long64_t idx3 = start; Long64_t idx4 = start;
                    for (Long64_t i = start; i < end; ++i) {
                        if (pmt1Vals[i] > peak1) { peak1 = pmt1Vals[i]; idx1 = i; }
                        if (pmt3Vals[i] > peak3) { peak3 = pmt3Vals[i]; idx3 = i; }
                        if (pmt4Vals[i] > peak4) { peak4 = pmt4Vals[i]; idx4 = i; }
                    }

                    //convert the sum to units of PE

                    
                    double integrated_pulse_PE_PMT3 = (computeSum4(pmt3Vals, idx3)- pedestal_pos_per_PMT[iPMT])/gains_per_PMT[iPMT];
                    double integrated_pulse_PE_PMT4 = (computeSum4(pmt4Vals, idx4)- pedestal_pos_per_PMT[iPMT+2])/gains_per_PMT[iPMT+2];

                    h_sum1->Fill(computeSum4(pmt1Vals, idx1));
                    h_sum3->Fill(integrated_pulse_PE_PMT3);
                    h_sum4->Fill(integrated_pulse_PE_PMT4);
                }

                // Fit PMT1 sum
                TF1* gaus1 = new TF1(("gaus1" + tag).c_str(), "gaus", sumMin, sumMax);

                // TF1* doublegaus3 = new TF1("doubleGaus", "gaus(0)+gaus(3)", sumMin, sumMax);
                // doublegaus3->SetLineColor(PMT3_cols[iPMT]);
                // doublegaus3->SetLineStyle(2);
                // doublegaus3->SetLineWidth(2);


                // TF1* doublegaus4 = new TF1("doubleGaus4", "gaus(0)+gaus(3)", sumMin, sumMax);
                // doublegaus4->SetLineColor(PMT4_cols[iPMT]);
                // doublegaus4->SetLineStyle(2);
                // doublegaus4->SetLineWidth(2);


                // TF1* gaus3 = new TF1(("gaus3" + tag).c_str(), "gaus", sumMin, sumMax);
                // TF1* gaus4 = new TF1(("gaus4" + tag).c_str(), "gaus", sumMin, sumMax);
                // gaus1->SetParameters(h_sum1->GetMaximum(), h_sum1->GetMean(), h_sum1->GetRMS());
                // gaus1->SetLineColor(kBlack);
                // gaus1->SetLineStyle(2);
                // gaus1->SetLineWidth(2);
                // h_sum1->Fit(gaus1, "Q0");
                // double mean1  = gaus1->GetParameter(1);
                // double sigma1 = gaus1->GetParameter(2);

                // double meanSum3 = 0.0, sigmaSum3 = 0.0;
                // double meanSum4 = 0.0, sigmaSum4 = 0.0;

                // if (hasPMT3 && h_sum3->GetEntries() > 0) {
                //     if (h_sum3->GetEntries() > 10) {

                //         gaus3->SetParameters(h_sum3->GetMaximum(), h_sum3->GetMean(), h_sum3->GetRMS());
                //         h_sum3->Fit(gaus3, "Q0");
                //         meanSum3  = gaus3->GetParameter(1);
                //         sigmaSum3 = gaus3->GetParameter(2);
                //         // delete gaus3;
                //     } else {
                //         meanSum3  = h_sum3->GetMean();
                //         sigmaSum3 = h_sum3->GetRMS();
                //     }
                // }

                // if (hasPMT4 && h_sum4->GetEntries() > 0) {
                //     if (h_sum4->GetEntries() > 10) {

                //         gaus4->SetParameters(h_sum4->GetMaximum(), h_sum4->GetMean(), h_sum4->GetRMS());
                //         h_sum4->Fit(gaus4, "Q0");
                //         meanSum4  = gaus4->GetParameter(1);
                //         sigmaSum4 = gaus4->GetParameter(2);
                //         // delete gaus4;
                //     } else {
                //         meanSum4  = h_sum4->GetMean();
                //         sigmaSum4 = h_sum4->GetRMS();
                //     }
                // }

                // ===================== DRAWING / LEGENDS =====================

                // Human-readable laser ball position string
                // (includes both angle and delay as you requested)
                std::string laserPosStr = Form("%d ps, %d deg", pw, ang);

                // Legend labels: include PMT name + laser ball position
                std::string labPMT1 = "Monitor PMT " + laserPosStr;
                std::string labPMT3 = laserPosStr + PMT3_name + " PMT";
                std::string labPMT4 = laserPosStr + PMT4_name + " PMT";

                canvases_0[iPMT]->cd();
                
                
                h_sum3->SetTitle("Sum4 at 0 ps;sum4 [PE];entries");
                h_sum3->Draw("HIST");
                // doublegaus3->Draw("SAME");

                legends_0[iPMT]->AddEntry(h_sum3,
                                Form("%s; %i ps; %i deg", labPMT3.c_str(), pulseWidths[0], angleLaserBall[0]),
                                "l");
                
                // legends_0[iPMT]->AddEntry(h_sum3,
                                //  Form("%s (ped mean=%.3e, ped #sigma=%.3e, 1PE mean=%.3e, 1PE #sigma=%.3e)", labPMT3.c_str(), meanSum3, sigmaSum3, mean1peSum3, sigma1peSum3),
                                //  "l");
                // legends_0[iPMT]->AddEntry((TObject*)0, Form("Gain = %.2e #pm %.2e", gain_PMT3, err_gain_PMT3), "");

                canvases_1[iPMT]->cd();

        
                h_sum4->SetTitle("Sum4 at 0 ps;sum4 [PE];entries");
                h_sum4->Draw("HIST");
                // doublegaus4->Draw("SAME");

                // legends_1[iPMT]->AddEntry(h_sum4,
                                //  Form("%s (ped mean=%.3e, ped #sigma=%.3e, 1PE mean=%.3e, 1PE #sigma=%.3e)", labPMT4.c_str(), meanSum4, sigmaSum4, mean1peSum4, sigma1peSum4),
                                //  "l");

                legends_1[iPMT]->AddEntry(h_sum4,
                                Form("%s; %i ps; %i deg", labPMT4.c_str(), pulseWidths[0], angleLaserBall[0]),
                                "l");
                // legends_1[iPMT]->AddEntry((TObject*)0, Form("Gain = %.2e #pm %.2e", gain_PMT4, err_gain_PMT4), "");

                


                allHists.push_back(h_sum1);
                allHists.push_back(h_sum3); 
                allHists.push_back(h_sum4); 

                // write to text summary
                // outTxt << pw << " " << ang << " Mon " << PMT3_name << " " << PMT4_name << " "
                // << mean1 << " " << sigma1 << " "
                // << meanSum3 << " " << sigmaSum3 << " "
                // << meanSum4 << " " << sigmaSum4 << "\n";

                f->Close();
//                delete gaus1;
//                delete doublegaus3;
//                delete gaus4;

                std::cout << "Finished analysis for " << pw << " ps, " << ang
                << " deg, PMT3=" << PMT3_name
                << ", PMT4=" << PMT4_name << "total number of windows: " << nWindows << std::endl;
            } // iPMT
        } // iAngle
    } // iPW


    //             if (is0ps) {
    //                 // -------- 0 ps canvases --------
    //                 // Spectrum
    //                 c0_spec->cd();
    //                 if (!drawn0_spec) {
    //                     h_spec1->SetTitle("Spectrum at 0 ps;volts [V];entries");
    //                     h_spec1->Draw("HIST");
    //                     drawn0_spec = true;
    //                 } else {
    //                     h_spec1->Draw("HIST SAME");
    //                 }
    //                 l0_spec->AddEntry(h_spec1, labPMT1.c_str(), "l");

    //                 // Peaks
    //                 c0_peak->cd();
    //                 if (!drawn0_peak) {
    //                     h_peak1->SetTitle("Peaks at 0 ps;peak volts [V];entries");
    //                     h_peak1->Draw("HIST");
    //                     drawn0_peak = true;
    //                 } else {
    //                     h_peak1->Draw("HIST SAME");
    //                 }
    //                 l0_peak->AddEntry(h_peak1, labPMT1.c_str(), "l");
    //                 if (hasPMT3) { h_peak3->Draw("HIST SAME"); l0_peak->AddEntry(h_peak3, labPMT3.c_str(), "l"); }
    //                 if (hasPMT4) { h_peak4->Draw("HIST SAME"); l0_peak->AddEntry(h_peak4, labPMT4.c_str(), "l"); }

    //                 // Sums
    //                 c0_sum->cd();
    //                 if (!drawn0_sum) {
    //                     h_sum3->SetTitle("Sum4 at 0 ps;sum4 [V];entries");
    //                     h_sum3->Draw("HIST");
    //                     gaus3->Draw("SAME");
    //                     drawn0_sum = true;
    //                 } else {
    //                     h_sum3->Draw("HIST SAME");
    //                     gaus3->Draw("SAME");
    //                 }
    //                 l0_sum->AddEntry(h_sum3,
    //                                  Form("%s (mean=%.3f, #sigma=%.3f)", labPMT3.c_str(), meanSum3, sigmaSum3),
    //                                  "l");

    //                 if (hasPMT4) {
    //                     h_sum4->Draw("HIST SAME");
    //                     l0_sum->AddEntry(h_sum4,
    //                                      Form("%s (mean=%.3f, #sigma=%.3f)", labPMT4.c_str(), meanSum4, sigmaSum4),
    //                                      "l");
    //                     gaus4->Draw("SAME");

    //                 }
    //             }

    //             if (is900ps) {
    //                 // -------- 900 ps canvases --------
    //                 // Spectrum
    //                 c900_spec->cd();
    //                 if (!drawn900_spec) {
    //                     h_spec1->SetTitle("Spectrum at 900 ps;volts [V];entries");
    //                     h_spec1->Draw("HIST");
    //                     drawn900_spec = true;
    //                 } else {
    //                     h_spec1->Draw("HIST SAME");
    //                 }
    //                 l900_spec->AddEntry(h_spec1, labPMT1.c_str(), "l");

    //                 // Peaks
    //                 c900_peak->cd();
    //                 if (!drawn900_peak) {
    //                     h_peak1->SetTitle("Peaks at 900 ps;peak volts [V];entries");
    //                     h_peak1->Draw("HIST");
    //                     drawn900_peak = true;
    //                 } else {
    //                     h_peak1->Draw("HIST SAME");
    //                 }
    //                 l900_peak->AddEntry(h_peak1, labPMT1.c_str(), "l");
    //                 if (hasPMT3) { h_peak3->Draw("HIST SAME"); l900_peak->AddEntry(h_peak3, labPMT3.c_str(), "l"); }
    //                 if (hasPMT4) { h_peak4->Draw("HIST SAME"); l900_peak->AddEntry(h_peak4, labPMT4.c_str(), "l"); }

    //                 // Sums
    //                 c900_sum->cd();
    //                 if (!drawn900_sum) {
    //                     h_sum3->SetTitle("Sum4 at 900 ps;sum4 [V];entries");
    //                     h_sum3->Draw("HIST");
    //                     gaus3->Draw("SAME");
    //                     drawn900_sum = true;
    //                 } else {
    //                     h_sum3->Draw("HIST SAME");
    //                     gaus3->Draw("SAME");
    //                 }
    //                 l900_sum->AddEntry(h_sum3,
    //                                    Form("%s (mean=%.3f, #sigma=%.3f)", labPMT3.c_str(), meanSum3, sigmaSum3),
    //                                    "l");

    //                 if (hasPMT4) {
    //                     h_sum4->Draw("HIST SAME");
    //                     l900_sum->AddEntry(h_sum4,
    //                                        Form("%s (mean=%.3f, #sigma=%.3f)", labPMT4.c_str(), meanSum4, sigmaSum4),
    //                                        "l");
    //                 }
    //             }

    //             // keep pointers for writing later
    //             allHists.push_back(h_spec1);
    //             allHists.push_back(h_peak1);
    //             allHists.push_back(h_sum1);
    //             if (hasPMT3) { allHists.push_back(h_peak3); allHists.push_back(h_sum3); }
    //             if (hasPMT4) { allHists.push_back(h_peak4); allHists.push_back(h_sum4); }

    //             // write to text summary
    //             outTxt << pw << " " << ang << " " << PMT3_name << " " << PMT4_name << " "
    //             << sum << " " << vmin << " " << vmax << " "
    //             << mean1 << " " << sigma1 << " "
    //             << meanSum3 << " " << sigmaSum3 << " "
    //             << meanSum4 << " " << sigmaSum4 << "\n";

    //             f->Close();
    //             delete gaus1;
    //             delete gaus3;
    //             delete gaus4;

    //             std::cout << "Finished analysis for " << pw << " ps, " << ang
    //             << " deg, PMT3=" << PMT3_name
    //             << ", PMT4=" << PMT4_name << std::endl;
    //         } // iPMT
    //     } // iAngle
    // } // iPW

    // ===================== SAVE CANVASES =====================

    
    for (int i=0; i<2; i++){
        canvases_0[i]->cd();   legends_0[i]->Draw();

        std::string filename =
        "../plots/summary/Characterisation_"+
        std::to_string(angleLaserBall[0]) +"deg_0ps_Sum4Points_" +
        PMT3_list[i] +
        ".pdf";

        canvases_0[i]->SaveAs(filename.c_str());
        canvases_1[i]->cd();   legends_1[i]->Draw();
        std::string filename2 =
        "../plots/summary/Characterisation_"+
        std::to_string(angleLaserBall[0]) +"deg_0ps_Sum4Points_" +
        PMT4_list[i] +
        ".pdf";
        canvases_1[i]->SaveAs(filename2.c_str());
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
