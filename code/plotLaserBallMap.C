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
    "../results/LaserBall_calibration_allAngles_810ps_full.root";

    std::vector<int> angles = {0, 90, 180, -90};

    // Order matters: defines the map layout
    std::vector<std::string> pmtDirs = {
        "PMT3_-30_norm",
        "PMT3_H_norm",
        "PMT4_+30_norm",
        "PMT4_V_norm",

    };


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

                if (!h) {
                    std::cerr << "No histogram in "
                    << pmtDirs[ip] << std::endl;
                    continue;
                }

                // double value = h->Integral();   // or GetEntries()
                double value = h->Integral(h->FindBin(pmtIntLims[ip]), h->GetNbinsX()+1);

                if (angles[ia] == 180 & ip == 3){
                    value = 0;
                }

                hMap->SetBinContent(ip+1, ia+1, value);
            }
        }

        // ---------- Draw ----------

        TCanvas* c = new TCanvas("cLaserBallMap",
                                 "LaserBall map", 900, 700);

        hMap->Draw("COLZ TEXT");

        c->SaveAs("LaserBallMap.png");
}
