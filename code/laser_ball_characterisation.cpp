// LaserBall characterisation
// Saves histograms by ANGLE → PMT directories with legends

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TDirectory.h"

int main(int argc, char** argv) {

    // ================= CONFIGURATION =================

    std::vector<int> pulseWidths = {810};
    std::vector<int> angleLaserBall = {0, 90, 180, -90};   // ADD MORE ANGLES HERE

    std::vector<std::string> PMT3_list = {"-30", "H"};
    std::vector<std::string> PMT4_list = {"V", "+30"};

    std::string treeNamePMT   = "PMT1data";
    std::string treeNameTrig  = "PMT2data";
    std::string treeNamePMT3  = "PMT3data";
    std::string treeNamePMT4  = "PMT4data";
    std::string branchName    = "volts";

    int    nBinsSum = 400;
    double sumMin   = -0.02;
    double sumMax   =  0.20;
    double monitorMax = 0.5;


    std::string outputRootName =
        Form("../results/LaserBall_calibration_allAngles_%dps_full.root", pulseWidths[0]);

    gStyle->SetOptStat(1111111);

    // ================= OUTPUT ROOT FILE =================

    TFile* outFile = TFile::Open(outputRootName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "ERROR: Cannot create output ROOT file\n";
        return 1;
    }

    // ================= LOOP =================

    for (int pw : pulseWidths) {

        for (int ang : angleLaserBall) {

            // -------- Angle directory --------
            outFile->cd();
            std::string angleDirName = Form("angle_%ideg", ang);
            TDirectory* angleDir = outFile->mkdir(angleDirName.c_str());
            angleDir->cd();

            for (size_t iPMT = 0; iPMT < PMT3_list.size(); ++iPMT) {

                std::string PMT3_name = PMT3_list[iPMT];
                std::string PMT4_name = PMT4_list[iPMT];

                std::string fileName =
                    "/vols/hyperk/540PC/171225/hadded_PMT1_Mon_PMT2_Trig_PMT3_" +
                    PMT3_name + "_PMT4_" + PMT4_name +
                    "_results_" + std::to_string(pw) +
                    "ps_" + std::to_string(ang) + "deg.root";

                if (gSystem->AccessPathName(fileName.c_str())) {
                    std::cerr << "Missing file: " << fileName << std::endl;
                    continue;
                }

                TFile* f = TFile::Open(fileName.c_str(), "READ");
                if (!f || f->IsZombie()) continue;

                TTree* tPMT1 = (TTree*)f->Get(treeNamePMT.c_str());
                TTree* tTrig = (TTree*)f->Get(treeNameTrig.c_str());
                TTree* tPMT3 = (TTree*)f->Get(treeNamePMT3.c_str());
                TTree* tPMT4 = (TTree*)f->Get(treeNamePMT4.c_str());

                if (!tPMT1 || !tTrig || !tPMT3 || !tPMT4) {
                    f->Close();
                    continue;
                }

                Long64_t nEntries = std::min(
                    std::min(tPMT1->GetEntries(), tTrig->GetEntries()),
                    std::min(tPMT3->GetEntries(), tPMT4->GetEntries())
                );

                double v1, vTrig, v3, v4;
                tPMT1->SetBranchAddress(branchName.c_str(), &v1);
                tTrig->SetBranchAddress(branchName.c_str(), &vTrig);
                tPMT3->SetBranchAddress(branchName.c_str(), &v3);
                tPMT4->SetBranchAddress(branchName.c_str(), &v4);

                std::vector<double> p1(nEntries), tr(nEntries),
                                    p3(nEntries), p4(nEntries);

                for (Long64_t i=0;i<nEntries;i++) {
                    tPMT1->GetEntry(i);
                    tTrig->GetEntry(i);
                    tPMT3->GetEntry(i);
                    tPMT4->GetEntry(i);
                    p1[i]=v1; tr[i]=vTrig; p3[i]=v3; p4[i]=v4;
                }

                // -------- Histograms --------
                TH1D* h1 = new TH1D("h1_sum4", "PMT1 sum4",
                                    nBinsSum, sumMin, monitorMax);
                TH1D* h3 = new TH1D("h3_sum4", Form("PMT %s sum4", PMT3_name.c_str()),
                                    nBinsSum, sumMin, sumMax);
                TH1D* h4 = new TH1D("h4_sum4", Form("PMT %s sum4", PMT4_name.c_str()),
                                    nBinsSum, sumMin, sumMax);
                TH1D* h3_norm = new TH1D("h3_sum4_norm", Form("PMT %s sum4 div. mon. PMT", PMT3_name.c_str()),
                                    nBinsSum, sumMin, sumMax/monitorMax);
                TH1D* h4_norm = new TH1D("h4_sum4_norm", Form("PMT %s sum4 div. mon. PMT", PMT4_name.c_str()),
                                    nBinsSum, sumMin, sumMax/monitorMax);
                h1->SetDirectory(nullptr);
                h3->SetDirectory(nullptr);
                h4->SetDirectory(nullptr);
                h3_norm->SetDirectory(nullptr);
                h4_norm->SetDirectory(nullptr);


                auto sum4 = [&](const std::vector<double>& a, Long64_t i){
                    //returns the sum of the four points around the max that was found
                    return a[i] +
                           (i>0 ? a[i-1]:0) +
                           (i+1<nEntries ? a[i+1]:0) +
                           (i+2<nEntries ? a[i+2]:0);
                };

                double thr = -0.6;
                std::vector<Long64_t> starts;
                //finding the laser periods - threshold crossing
                for (Long64_t i=1;i<nEntries;i++)
                    if (tr[i-1]<thr && tr[i]>=thr)
                        starts.push_back(i);

                for (size_t ip=0; ip<starts.size(); ip++) {
                    //look at each individual period
                    Long64_t s = starts[ip];
                    Long64_t e = (ip+1<starts.size()) ? starts[ip+1] : nEntries;

                    double m1=-1e9, m3=-1e9, m4=-1e9;
                    Long64_t i1=s, i3=s, i4=s;
                    //find the maximums 
                    for (Long64_t i=s;i<e;i++) {
                        if (p1[i]>m1){m1=p1[i]; i1=i;}
                        if (p3[i]>m3){m3=p3[i]; i3=i;}
                        if (p4[i]>m4){m4=p4[i]; i4=i;}
                    }
                    double s1 = std::max(sum4(p1,i1), 1e-9);
                    double s3 = sum4(p3,i3);
                    double s4 = sum4(p4,i4);
                    h1->Fill(s1);
                    h3->Fill(s3);
                    h4->Fill(s4);
                    h3_norm->Fill(std::min(s3/s1, 999.));
                    h4_norm->Fill(std::min(s4/s1, 999.));

                }

                // -------- PMT directories --------
                angleDir->cd();
                TDirectory* d1 = angleDir->mkdir(("PMT1_"+PMT3_name).c_str());
                TDirectory* d3 = angleDir->mkdir(("PMT3_"+PMT3_name).c_str());
                TDirectory* d4 = angleDir->mkdir(("PMT4_"+PMT4_name).c_str());
                TDirectory* d3_norm = angleDir->mkdir(("PMT3_"+PMT3_name+"_norm").c_str());
                TDirectory* d4_norm = angleDir->mkdir(("PMT4_"+PMT4_name+"_norm").c_str());


                // -------- Write PMT1 --------
                d1->cd();
                h1->Write(Form("Monitor PMT %i: %d deg",iPMT, ang));

                TLegend* l1 = new TLegend(0.15,0.65,0.9,0.9);
                l1->SetBorderSize(0);
                l1->SetFillStyle(0);
                l1->AddEntry(h1, Form("Monitor PMT %i: %d deg",iPMT, ang), "l");
                l1->AddEntry((TObject*)0,
                    Form("Integral = %.0f", h1->Integral()), "");
                l1->Write("legend");

                // -------- Write PMT3 --------
                d3->cd();
                h3->Write(Form("PMT %s: %d deg",PMT3_name.c_str(), ang));

                TLegend* l3 = new TLegend(0.15,0.65,0.9,0.9);
                l3->SetBorderSize(0);
                l3->SetFillStyle(0);
                l3->AddEntry(h3, Form("PMT %s: %d deg",PMT3_name.c_str(), ang), "l");
                l3->AddEntry((TObject*)0,
                    Form("Integral = %.0f", h3->Integral()), "");
                l3->Write("legend");

                // -------- Write PMT4 --------
                d4->cd();
                h4->Write(Form("PMT %s: %d deg",PMT4_name.c_str(), ang));

                TLegend* l4 = new TLegend(0.15,0.65,0.9,0.9);
                l4->SetBorderSize(0);
                l4->SetFillStyle(0);
                l4->AddEntry(h4, Form("PMT %s: %d deg",PMT4_name.c_str(), ang), "l");
                l4->AddEntry((TObject*)0,
                    Form("Integral = %.0f", h4->Integral()), "");
                l4->Write("legend");

                // -------- Write PMT3 divided by monitor PMT --------
                d3_norm->cd();
                h3_norm->Write(Form("PMT %s norm. to mon.: %d deg",PMT3_name.c_str(), ang));

                TLegend* l3_norm = new TLegend(0.15,0.65,0.9,0.9);
                l3_norm->SetBorderSize(0);
                l3_norm->SetFillStyle(0);
                l3_norm->AddEntry(h3_norm, Form("PMT %s norm. to mon.: %d deg",PMT3_name.c_str(), ang), "l");
                l3_norm->AddEntry((TObject*)0,
                    Form("Integral = %.0f", h3_norm->Integral()), "");
                l3_norm->Write("legend");

                // -------- Write PMT4 divided by monitor PMT --------
                d4_norm->cd();
                h4_norm->Write(Form("PMT %s norm. to mon.: %d deg",PMT4_name.c_str(),ang));

                TLegend* l4_norm = new TLegend(0.15,0.65,0.9,0.9);
                l4_norm->SetBorderSize(0);
                l4_norm->SetFillStyle(0);
                l4_norm->AddEntry(h4_norm, Form("PMT %s norm. to mon.: %d deg",PMT4_name.c_str(), ang), "l");
                l4_norm->AddEntry((TObject*)0,
                    Form("Integral = %.0f", h4_norm->Integral()), "");
                l4_norm->Write("legend");

                f->Close();
            }
        }
    }

    outFile->Write();
    outFile->Close();

    std::cout << "DONE. Output written to "
              << outputRootName << std::endl;

    return 0;
}
