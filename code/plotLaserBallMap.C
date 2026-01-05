#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

void plotLaserBallMap() {

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    // ---------- Configuration ----------

    std::string filename =
    "../results/LaserBall_calibration_allAngles_810ps_withTTree.root";

    std::vector<int> angles = {0, 90, 180, -90};

    // Order matters: defines the map layout
    std::vector<std::string> pmtDirs = {
        "PMT3_m30_norm",
        "PMT3_H_norm",
        "PMT4_p30_norm",
        "PMT4_V_norm",

    };

    std::vector<std::string> pmtDirs_raw = {
        "PMT3_m30",
        "PMT3_H",
        "PMT4_p30",
        "PMT4_V",

    };

    //actually we should read the TTree with the entries, fit it with gaussians and then extract the gain. Does it have to be stable ? not really ? Maybe print out the values, just to ensure it's ok 


    std::vector<double> pmtGains = {
        0.0333144,
        0.029094,
        0.042336,
        0.0222327
        ,

    };

    std::vector<double> pmtPedestal= {
        0.00371834,
        0.00430222,
        0.0043416775,
        0.004803675,

    };

    std::vector<double> pmtIntLims = {
        0.0021896,
        0.0031439,
        0.0055,
        0.00794,

    };

    std::vector<std::string> tree_name = {
        "pulseTree_0",
        "pulseTree_1",
        "pulseTree_1",
        "pulseTree_0",

    };


    std::vector<std::string> monitor_raw = {
        "PMT1_Monitor_0",
        "PMT1_Monitor_1",
        "PMT1_Monitor_1",
        "PMT1_Monitor_0",

    };





    // ---------- Open file ----------

    TFile* f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open file\n";
        return;
    }

    // ---------- Create map ----------

    TH2D* hMap = new TH2D(
        "LaserBallMap",
        "Average number of PE per Trigger \n normalised to Monitor (3sig. ped. cut);PMT position;Angle",
                          pmtDirs.size(), 0, pmtDirs.size(),
                          angles.size(), 0, angles.size()
    );

    // Axis labels
    for (size_t i = 0; i < pmtDirs.size(); ++i)
        hMap->GetXaxis()->SetBinLabel(i+1, pmtDirs[i].c_str());

    for (size_t i = 0; i < angles.size(); ++i)
        hMap->GetYaxis()->SetBinLabel(i+1,
                                      Form("%d deg", angles[i]));

        // ---------- Fill map ----------

        for (size_t ia = 0; ia < angles.size(); ++ia) {

            std::string angleDirName =
            Form("angle_%ideg", angles[ia]);

            TDirectory* angleDir =
            (TDirectory*)f->Get(angleDirName.c_str());

            if (!angleDir) {
                std::cerr << "Missing directory "
                << angleDirName << std::endl;
                continue;
            }

            



            for (size_t ip = 0; ip < pmtDirs.size(); ++ip) {

                TDirectory* d =
                (TDirectory*)angleDir->Get(pmtDirs_raw[ip].c_str());

                if (!d) {
                    std::cerr << "Missing directory "
                    << pmtDirs[ip] << std::endl;
                    continue;
                }

                // Get the only histogram in the directory
                TH1D* h = nullptr;
                d->GetObject(d->GetListOfKeys()->At(0)->GetName(), h);

                
                //some quick plotting to validate
                TCanvas* c = new TCanvas(
                    Form("c_%s_%s", pmtDirs_raw[ip].c_str(), angleDirName.c_str()),
                    pmtDirs[ip].c_str(),
                    800, 600
                );

                h->SetLineWidth(2);
                h->Draw("HIST");  

                

                if (!h) {
                    std::cerr << "No histogram in "
                    << pmtDirs[ip] << std::endl;
                    continue;
                }
                //Fit the histogram with two gaussians
                TF1* doublegaus = new TF1("doubleGaus", "gaus(0)+gaus(3)", 0.0, 0.1);

                // Inititalise the fit
                doublegaus->SetParameters(
                    h->GetMaximum(), 0.0037, 0.0001,
                                          h->GetMaximum()*0.01, 0.04, 0.01
                );

                //V
                if (ip == 3) doublegaus->SetParameters(
                    h->GetMaximum(), 0.005, 0.002, h->GetMaximum()*0.01, 0.026, 0.007
                );

                //p30
                if (ip == 2) doublegaus->SetParameters(
                    h->GetMaximum(), 0.0047, 0.002, h->GetMaximum()*0.01, 0.047, 0.0136
                );

                //H
                if (ip == 1) doublegaus->SetParameters(
                    h->GetMaximum(), 0.0046, 0.0012, h->GetMaximum()*0.01, 0.0332, 0.010
                );


                TFitResultPtr res = h->Fit(doublegaus, "SRQ0");
                h->SetTitle(Form("c_%s_%s; sum4 [V]; entries", pmtDirs_raw[ip].c_str(), angleDirName.c_str()));

                double mean_ped  = doublegaus->GetParameter(1);
                double sigma_ped = doublegaus->GetParameter(2);
                double mean_1pe  = doublegaus->GetParameter(4);
                double sigma_1pe = doublegaus->GetParameter(5);

                double gain_PMT = mean_1pe - mean_ped;

                double ped_err   = doublegaus->GetParError(1);
                double onepe_err = doublegaus->GetParError(4);
                double cov = res->GetCovarianceMatrix()(1,4);

                double err_gain_PMT3 = std::sqrt(
                    onepe_err*onepe_err +
                    ped_err*ped_err -
                    2.0*cov
                );

                doublegaus->SetLineColor(kRed);
                doublegaus->SetLineWidth(2);
                doublegaus->Draw("SAME");     // overlay fit



                c->Update();

                std::cout << "\n\n  "<< angleDirName << ":\nFor histogram " << pmtDirs_raw[ip] << " the fit parameters are \nmean_ped = " << mean_ped << "\nsigma_ped = " << sigma_ped << "\n mean_1pe = " << mean_1pe << "\nsigma_1pe = " << sigma_1pe << "\nGAIN = " << gain_PMT << " +/- "<< err_gain_PMT3 << " V/PE" <<  std::endl;


                // double value = h->Integral();   // or GetEntries()
                // double value = h->Integral(h->FindBin(pmtIntLims[ip]), h->GetNbinsX()+1);
                //
                // if (angles[ia] == 180 & ip == 3){
                //     value = 0;
                // }


                TTree* T  = dynamic_cast<TTree*>(angleDir->Get(tree_name[ip].c_str()));
                if (!T) {
                    std::cerr << "ERROR: Tree \"" << tree_name[ip] << "\" not found in file!" << std::endl;
                    continue;
                }

                double value_PMT;
                double value_Monitor;
                std::cout << pmtDirs_raw[ip].c_str() << " " << tree_name[ip] << std::endl;


                const char* brName = pmtDirs_raw[ip].c_str();

                TBranch* br = T->GetBranch(brName);
                TBranch* brMonitor = T->GetBranch(monitor_raw[ip].c_str());
                if (!br) {
                    std::cerr << "ERROR: Branch \"" << brName
                    << "\" not found in tree!" << std::endl;
                    continue;
                }

                T->SetBranchAddress(brName, &value_PMT);
                T->SetBranchAddress(monitor_raw[ip].c_str(), &value_Monitor);

                Long64_t nEntries = T->GetEntries();
                double sum_entries = 0 ;
                double sum_entries_normalised = 0 ;

                for (Long64_t i = 0; i < nEntries; ++i) {
                    T->GetEntry(i);
                    // if (value_PMT > pmtIntLims[ip] && ip == 2) std::cout << "Entry i: " << i << " has value:  " << value_PMT << " i.e. " <<  (value_PMT - pmtPedestal[ip])/pmtGains[ip] << " PE." << std::endl;
                    double entryPE = (value_PMT - pmtPedestal[ip])/pmtGains[ip];
                    if (value_PMT > (pmtIntLims[ip]+pmtPedestal[ip])){
                        sum_entries += entryPE;
                        sum_entries_normalised += entryPE/value_Monitor;
                    }
                }



                hMap->SetBinContent(ip+1, ia+1, sum_entries_normalised/nEntries);
            }

            //here we are reading the entries, we'll have to normalise to 1pe before dividing 
         /*
            TTree* treePMT_0  = dynamic_cast<TTree*>(f->Get("pulseTree_0"));
            Long64_t nTriggers_0 = treePMT_0->GetEntries();*/
        }

        // ---------- Draw ----------

        TCanvas* c = new TCanvas("cLaserBallMap",
                                 "LaserBall map", 900, 700);

        hMap->Draw("COLZ TEXT");

        c->SaveAs("LaserBallMap.png");
}
