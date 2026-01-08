#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"

void plotLaserBallMap_withTimingCuts() {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // ================= CONFIGURATION =================

    std::string filename =
    "../results/LaserBall_calibration_allAngles_810ps_060126_testWithTimes.root";

    std::vector<int> angles = {0, 90, 180, -90};

    std::vector<std::string> pmtNames = {
        "PMT3_m30",
        "PMT3_H",
        "PMT4_p30",
        "PMT4_V"
    };

    std::vector<std::string> treeNames = {
        "pulseTree_0",
        "pulseTree_1",
        "pulseTree_1",
        "pulseTree_0"
    };

    std::vector<std::string> monitorNames = {
        "PMT1_Monitor_0",
        "PMT1_Monitor_1",
        "PMT1_Monitor_1",
        "PMT1_Monitor_0"
    };

    // Initial fit guesses
    std::vector<double> ped_guess = {-0.002, -0.004, -0.002, -0.001};
    std::vector<double> ped_sigma_guess = {0.001, 0.001, 0.003, 0.004};
    std::vector<double> onepe_guess = {0.032, 0.025, 0.042, 0.024};
    std::vector<double> onepe_sigma_guess = {0.009, 0.008, 0.013, 0.006};


    //In case it's needed, we can also apply a fixed pedestal and gain value, deciding it from
    //a combination of fits. If the fits and gain are stable this should not modify much, if not
    //then it is better to understand the reason behind the PMT instability...
    std::vector<double> ped_mean = {-0.00242, -0.00219, -0.00184, 0.0044};
    std::vector<double> gain_mean = {0.03353, 0.02891, 0.04337, 0.02271};
    // ================= OUTPUT FILE =================

    std::ofstream gainFile("PMT_gain_summary.txt");
    gainFile << "# angle  PMT  gain[V/PE]  gain_err[V/PE] pedestal[V] pedestal_sigma[V]\n";

    // ================= OPEN ROOT FILE =================

    TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open input file\n";
        return;
    }

    // ================= LASER BALL MAP =================

    TH2D* hMap = new TH2D(
        "LaserBallMap",
        "Sum(#SignalPE/MonitorPMT)/nTriggersTot;PMT;Angle",
        pmtNames.size(), 0, pmtNames.size(),
                          angles.size(), 0, angles.size()
    );

    for (size_t i = 0; i < pmtNames.size(); ++i)
        hMap->GetXaxis()->SetBinLabel(i+1, pmtNames[i].c_str());

    for (size_t i = 0; i < angles.size(); ++i)
        hMap->GetYaxis()->SetBinLabel(i+1, Form("%d deg", angles[i]));

    // ================= MAIN LOOP =================

    for (size_t ia = 0; ia < angles.size(); ++ia) {

        std::string angleDirName = Form("angle_%ideg", angles[ia]);
        TDirectory* angleDir = (TDirectory*)f->Get(angleDirName.c_str());
        if (!angleDir) continue;

        for (size_t ip = 0; ip < pmtNames.size(); ++ip) {

            TTree* T = (TTree*)angleDir->Get(treeNames[ip].c_str());
            if (!T) continue;

            double pulse, monitor;
            int hit_index;

            T->SetBranchAddress(pmtNames[ip].c_str(), &pulse);
            T->SetBranchAddress(monitorNames[ip].c_str(), &monitor);

            std::string hitBranch = pmtNames[ip] + "hit_index";
            T->SetBranchAddress(hitBranch.c_str(), &hit_index);

            Long64_t nEntries = T->GetEntries();

            // ---------- Find laser arrival time ----------

            TH1I* hHit = new TH1I("hHit", "", 2000, 4000, 6000);

            for (Long64_t i = 0; i < nEntries; ++i) {
                T->GetEntry(i);
                if (hit_index >= 4000 && hit_index <= 6000)
                    hHit->Fill(hit_index);
            }

            int mean_arrival_time =
            hHit->GetBinCenter(hHit->GetMaximumBin());

            // ---------- Separate samples ----------

            std::vector<double> dark_pulses;
            std::vector<double> laser_pulses;
            std::vector<double> laser_monitor;

            for (Long64_t i = 0; i < nEntries; ++i) {
                T->GetEntry(i);
                if (std::abs(hit_index - mean_arrival_time) <= 5) {
                    laser_pulses.push_back(pulse);
                    laser_monitor.push_back(monitor);
                } else {
                    dark_pulses.push_back(pulse);
                }
            }

            // ---------- Gain fit (dark-rate only) ----------

            TH1D* hDark = new TH1D(
                "hDark",
                Form("Dark rate: %s at %d deg; sum4 [V]; Entries",
                     pmtNames[ip].c_str(), angles[ia]),
                                   50, -0.03, 0.1
            );

            for (double v : dark_pulses)
                hDark->Fill(v);

            // Double Gaussian: pedestal + 1 PE
            TF1* fG = new TF1("fG", "gaus(0)+gaus(3)", -0.03, 0.1);
            fG->SetParameters(
                hDark->GetMaximum(), ped_guess[ip], ped_sigma_guess[ip],
                              hDark->GetMaximum()*0.05, onepe_guess[ip], onepe_sigma_guess[ip]
            );

            // fG->SetParLimits(0, 0.9*p0, 1.1*p0);
            fG->SetParLimits(1, -0.005, 0.005);
            fG->SetParLimits(2, 0.8* ped_sigma_guess[ip], 1.2* ped_sigma_guess[ip]);
            // fG->SetParLimits(3, 0.9*p3, 1.1*p3);
            fG->SetParLimits(4, 0.8*onepe_guess[ip], 1.2*onepe_guess[ip]);
            fG->SetParLimits(5, 0.8*onepe_sigma_guess[ip], 1.2*onepe_sigma_guess[ip]);


            TFitResultPtr res = hDark->Fit(fG, "SRQ");

            // ---------- Extract fit parameters ----------

            double ped     = fG->GetParameter(1);
            double ped_sig = fG->GetParameter(2);
            double onepe   = fG->GetParameter(4);
            double onepe_sig = fG->GetParameter(5);

            double gain = onepe - ped;

            double ped_err   = fG->GetParError(1);
            double onepe_err = fG->GetParError(4);
            double cov = res->GetCovarianceMatrix()(1,4);

            double gain_err = std::sqrt(
                onepe_err*onepe_err +
                ped_err*ped_err -
                2.0*cov
            );

            // ---------- Save gain to text file ----------

            gainFile << angles[ia] << " "
            << pmtNames[ip] << " "
            << gain << " "
            << gain_err << " " << ped << " " << ped_sig << "\n";

            // ---------- Plot ----------

            TCanvas* cDark = new TCanvas(
                Form("cDark_%s_%d", pmtNames[ip].c_str(), angles[ia]),
                                         "", 800, 600
            );

            cDark ->SetLogy();
            hDark->SetLineWidth(2);
            hDark->Draw("HIST");

            fG->SetLineColor(kRed);
            fG->SetLineWidth(2);
            fG->Draw("SAME");

            // ---------- Legend with fit results ----------

            TLegend* leg = new TLegend(0.35, 0.60, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.032);

            leg->AddEntry(hDark, "Dark-rate spectrum", "l");
            leg->AddEntry(fG, "2-Gaussian fit", "l");

            leg->AddEntry((TObject*)0,
                          Form("Pedestal = %.4f #pm %.4f V", ped, ped_sig),
                          "");
            leg->AddEntry((TObject*)0,
                          Form("1 PE = %.4f #pm %.4f V", onepe, onepe_sig),
                          "");
            leg->AddEntry((TObject*)0,
                          Form("Gain = %.4f #pm %.4f V/PE", gain, gain_err),
                          "");

            leg->Draw();

            // ---------- Save ----------

            cDark->SaveAs(Form("../plots/060126/DarkRate/DarkRateGain_%s_%ddeg.png",
                               pmtNames[ip].c_str(), angles[ia]));

            // ---------- PE overlap plot ----------

            TH1D* hDarkPE  = new TH1D(
                "hDarkPE",
                Form("%s at %d deg;PE;Entries", pmtNames[ip].c_str(), angles[ia]),
                                      100, -0.1, 3
            );

            TH1D* hLaserPE = new TH1D(
                "hLaserPE",
                Form("%s at %d deg;PE;Entries", pmtNames[ip].c_str(), angles[ia]),
                                      100, -0.1, 3
            );

            for (double v : dark_pulses)
                hDarkPE->Fill((v - ped) / gain);

            for (double v : laser_pulses)
                hLaserPE->Fill((v - ped) / gain);

            hDarkPE->SetLineColor(kBlue);
            hDarkPE->SetLineWidth(2);

            hLaserPE->SetLineColor(kRed);
            hLaserPE->SetLineWidth(2);

            TCanvas* cPE = new TCanvas(
                Form("cPE_%s_%d", pmtNames[ip].c_str(), angles[ia]),
                                       "", 800, 600
            );

            // Draw
            cPE->SetLogy();
            hDarkPE->Draw("HIST");
            hLaserPE->Draw("HIST SAME");

            // ---------- Legend ----------
            TLegend* leg1 = new TLegend(0.60, 0.70, 0.88, 0.88);
            leg1->SetBorderSize(0);
            leg1->SetFillStyle(0);
            leg1->SetTextSize(0.035);

            leg1->AddEntry(hDarkPE,  "Dark-rate events",  "l");
            leg1->AddEntry(hLaserPE, "Laser signal",      "l");

            leg1->Draw();

            // Save
            cPE->SaveAs(Form("../plots/060126/Signal/PE_overlap_%s_%ddeg.png",
                             pmtNames[ip].c_str(), angles[ia]));

            // ---------- Laser ball observable ----------

            double sum_norm = 0.0;

            for (size_t i = 0; i < laser_pulses.size(); ++i) {
                //here can use the ped and gain calculated for this specific entry or the mean vals
                double pe = (laser_pulses[i] - ped) / gain;
                sum_norm += pe / laser_monitor[i];
            }

            double laserBallValue =
            sum_norm / laser_pulses.size();

            hMap->SetBinContent(ip+1, ia+1, laserBallValue);
        }
    }

    // ================= FINAL MAP =================

    TCanvas* cMap = new TCanvas("cLaserBallMap", "", 900, 700);
    hMap->Draw("COLZ TEXT");
    cMap->SaveAs("LaserBallMap_clean.png");

    gainFile.close();
}
