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

    //actually we should read the TTree with the entries, fit it with gaussians and then extract the gain. Does it have to be stable ? not really ? Maybe print out the values, just to ensure it's ok 


    std::vector<double> pmtIntLims = {
        0.03,
        0.035,
        0.04,
        0.04,

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
        "LaserBall intensity map (normalised);PMT position;Angle",
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
                (TDirectory*)angleDir->Get(pmtDirs[ip].c_str());

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
                    Form("c_%s", pmtDirs[ip].c_str()),
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
                    h->GetMaximum(), h->GetMean()*0.2, h->GetRMS()*0.5,
                    h->GetMaximum()*0.1, h->GetMean(), h->GetRMS()
                );
                h->Fit(doublegaus, "R");

                double mean_ped  = doublegaus->GetParameter(1);
                double sigma_ped = doublegaus->GetParameter(2);
                double mean_1pe  = doublegaus->GetParameter(4);
                double sigma_1pe = doublegaus->GetParameter(5);

                doublegaus->SetLineColor(kRed);
                doublegaus->SetLineWidth(2);
                doublegaus->Draw("SAME");     // overlay fit

                c->Update();

                std::cout << "For histogram " << pmtDirs[ip] << " the fit parameters are \nmean_ped = " << mean_ped << "\nsigma_ped = " << sigma_ped << "\n mean_1pe = " << mean_1pe << "\nsigma_1pe = " << sigma_1pe << std::endl;


                // double value = h->Integral();   // or GetEntries()
                double value = h->Integral(h->FindBin(pmtIntLims[ip]), h->GetNbinsX()+1);

                if (angles[ia] == 180 & ip == 3){
                    value = 0;
                }

                hMap->SetBinContent(ip+1, ia+1, value);
            }

            //here we are reading the entries, we'll have to normalise to 1pe before dividing 
            
            TTree* treePMT_0  = dynamic_cast<TTree*>(f->Get("pulseTree_0"));
            Long64_t nTriggers_0 = treePMT_0->GetEntries();
        }

        // ---------- Draw ----------

        TCanvas* c = new TCanvas("cLaserBallMap",
                                 "LaserBall map", 900, 700);

        hMap->Draw("COLZ TEXT");

        c->SaveAs("LaserBallMap.png");
}
