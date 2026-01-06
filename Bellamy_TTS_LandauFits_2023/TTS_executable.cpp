#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TText.h"
#include "TStyle.h"

//This is a code that takes in laser data and then gets the TTS out
bool passVeto(double y[], double vetoHigh, int lenWindow){
	for (int i = 0; i<lenWindow; i++){
		if (y[i]>vetoHigh){
			std::cout << "passed the " << vetoHigh << " veto" << std::endl;
		       	return true;
		}	
	}
	return false;
	
}

void HelpMessage(){
	std::cout << "This is the code to optain the TTS for the LASER BALL \n"
		  << "USAGE: \n"
		  << "-f the PMTform (Losange, Square, Circle, Triangle) \n"
		  << "-p the laser ball orientation and pulse current (e.g. 180deg154mA) \n"
		  << "-r the number of 1ms windows that we have (default 500) \n"
		  << "-t TTS or TTSLB -> if you want to look at LB data or fibre (default LB) \n"
		  << "careful, this relies on correct location of your files and naming... \n";
}

void CalculateTTS(TFile *f, char *name, TH1D *hTTS, TH2F *hcol1, int channelTrigger, int rep=500){
	//variable definition
	double pe1; //this is the number of pe of the laser file
	double volts_trigger; //this is the adc of the trigger from the file
	int lenWindow = 2000; //this is one single scope window = one trigger -> CAREFULL, it depends
	//on your sampling rate and window size (on the scope)
	//int rep = 10; //the number of 1ms windows that we have i
	int numberWindCheck = 1000 * rep;
	double x[lenWindow], y[lenWindow], trigger[lenWindow]; //this is the x and y 
	//that we want to zoom in to have an idea
	double threshold = -0.003;
// 	double distToLocalTrig; //for each pulse
	//int n_largeWindows = 0; // the number of large windows we have circled through
	double laserT0 = -9999.;
	double pmtT0 = -9999.;
// 	double laserRate = 1000; //interval in ns between two laser pulses

	//double TTS[numberWindCheck];
       	//double pe_max[numberWindCheck];
	double max_pe_in_window = 10;
	int offsetWind = 1000; //sometimes the window starts just after the laser trigger -> that causes a mess
	//adding an offset is a simple fix
	bool laserTriggerFound = false;
	bool pmtTriggerFound = false;
	double laser_thresh = -0.5; //the drop after which we consider the laser has lased

	
	TTree *t1 = (TTree*)f->Get(name);
	t1->SetBranchAddress("volts",&pe1);
    Int_t nentries = (Int_t)t1->GetEntries();
    
    TTree *T = (TTree*)f->Get(Form("PMT%idata", channelTrigger));
	T->SetBranchAddress("volts",&volts_trigger);
    
    std::cout << "Read the data hopefully correctly " << std::endl;
	
	//Where we store the raw data
	TH1D *h = new TH1D("hLASnormal","pe distribution normal",100,0,10);

	for (Int_t k=0; k<numberWindCheck; k++){
		//re-initialise the points at which we cross thresholds 
		laserTriggerFound = false;
		pmtTriggerFound = false;
		max_pe_in_window = 10.;
		for (Int_t i=0; i<lenWindow+1; i++) {
	                t1->GetEntry(i + k*lenWindow + offsetWind);
	                T->GetEntry(i + k*lenWindow + offsetWind);
	                h->Fill(pe1);
			int j = i%lenWindow; 
			x[j] = j; //look at the raw output of the file
            		y[j] = -pe1;
			if (-pe1 < max_pe_in_window) max_pe_in_window = -pe1;
			trigger[j] = volts_trigger;
			//note: here we assume that we stat at 0 and then drop to -1 during the period:
			//really might not always be the case
//             std::cout << /*volts_trigger*/ << std::endl;
			if (volts_trigger > laser_thresh && j>0 && trigger[j-1] < laser_thresh && !laserTriggerFound){
                 		//std::cout << "laser trigger found here: " << j << " window " << i << std::endl; 
				laserTriggerFound = true;
			        //laserT0 = j; //-> without interpolation
				//linear interpolation part
				double a_interpolation = (trigger[j]-trigger[j-1]);
				double b_interpolation = trigger[j] - a_interpolation * j;
				laserT0 = (laser_thresh - b_interpolation) / a_interpolation;

			}
			//std::cout << pe1 << std::endl; 
			if (pe1 > threshold && j>0 && y[j-1] < threshold && !pmtTriggerFound){
				pmtTriggerFound = true;
				double a_interpolation = -(y[j]-y[j-1]);
				double b_interpolation = y[j] - a_interpolation * j;
				pmtT0 = (threshold - b_interpolation) / a_interpolation;
			}
		}
		if (laserTriggerFound and pmtTriggerFound){
			double TTS = -9999.;
			
            
            		if ((pmtT0-laserT0)> 0) TTS = (pmtT0-laserT0) * 0.5;

			hTTS->Fill(TTS);
            hcol1->Fill(TTS,max_pe_in_window);
			}
	}
	delete h;
}


int main(int argc, char **argv){
	char * PMTform = NULL;
	char * face_position = NULL;
	char option;
	int rep = 10; //the number of 1ms windows that we have i
	int channelTrigger = 0;
	const char * fibre_or_LB;
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.5);	
    
    	char * path_to_input_file = NULL;
    	char * path_to_output_file = "../data/output_TTS.root";
	
	//can change this for a more user-chosen approach, these values will work for our
	//set-up @ Imperial
	double minHist = 540;
	double maxHist = 650;
	
	while ((option = getopt(argc, argv, "f:p:r:t:b:i:o:e:h")) != -1 )
	{
		switch(option)
		{
			case 'f':
				PMTform = optarg;
				break;
			case 'p':
                face_position = optarg;
                break;
			case 'r':
				rep = std::stoi(optarg);
				break;
			case 't':
				channelTrigger = std::stoi(optarg);
				break;
			case 'b':
				minHist = std::stod(optarg);
				break;
			case 'e':
				maxHist = std::stod(optarg);
				break;
			case 'h':
				HelpMessage();
				break;
            case 'i':
				path_to_input_file = optarg;
				break;
			case 'o':
				path_to_output_file = optarg;
				break;
			default:
				return 0;
		}
	}

	if (path_to_input_file == NULL || channelTrigger == 0) HelpMessage();

   	auto hcol1 = new TH2F("QvsT_PMT1","Collected charge vs hit time",100,400,minHist,maxHist,-0.1,0.03);
	TH1D *hTTS = new TH1D("hTTS_PMT1","hTTS_PMT1", int((maxHist-minHist)/0.5),minHist, maxHist);
   	auto h2col1 = new TH2F("QvsT_PMT2","Collected charge vs hit time",100,400,minHist, maxHist,-0.1,0.03);
	TH1D *h2TTS = new TH1D("hTTS_PMT2","hTTS_PMT2", int((maxHist-minHist)/0.5),minHist, maxHist);
   	auto h3col1 = new TH2F("QvsT_PMT3","Collected charge vs hit time",100,400,minHist, maxHist,-0.1,0.03);
	TH1D *h3TTS = new TH1D("hTTS_PMT3","hTTS_PMT3", int((maxHist-minHist)/0.5),minHist, maxHist);
    auto h4col1 = new TH2F("QvsT_PMT4","Collected charge vs hit time",100,400,minHist, maxHist,-0.1,0.03);
	TH1D *h4TTS = new TH1D("hTTS_PMT4","hTTS_PMT4", int((maxHist-minHist)/0.5),minHist, maxHist);
    
    std::vector<TH2F*> vector_TH2F;
    vector_TH2F.push_back(hcol1);
    vector_TH2F.push_back(h2col1);
    vector_TH2F.push_back(h3col1);
    vector_TH2F.push_back(h4col1);
    
    std::vector<TH1D*> vector_TH1D;
    vector_TH1D.push_back(hTTS);
    vector_TH1D.push_back(h2TTS);
    vector_TH1D.push_back(h3TTS);
    vector_TH1D.push_back(h4TTS);
        
    TFile *f = new TFile(Form("%s", path_to_input_file));
    
    for (int channel=1; channel<=4; channel++){
      if (channel!=channelTrigger) {
        	CalculateTTS(f, Form("PMT%idata", channel), vector_TH1D[channel-1], vector_TH2F[channel-1], channelTrigger, rep);
         }
    }
	
//set up the plots in a cleaner way
	for (int channel=1; channel<=4; channel++){
       		 auto h = vector_TH2F[channel-1];
       		 h->Draw("COLZ");
       		 h->GetXaxis()->SetTitle("Delay between laser and pmt threshold crossing (ns)");
       		 h->GetYaxis()->SetTitle("Max depth of PMT pulse (arbitrary.)");
       		 auto h1 = vector_TH1D[channel-1];
       		 h1->Draw();
       		 h1->GetXaxis()->SetTitle("Delay between laser and pmt threshold crossing (ns)");
       		 h1->GetYaxis()->SetTitle(Form("Occurencies/0.5ns"));
       		 h1->Fit("gaus");
    }
    TFile* o = new TFile(Form("%s",path_to_output_file), "RECREATE");
    //h_LAS1sPD->Write();
    
    for (int channel=1; channel<=4; channel++){
        auto h = vector_TH2F[channel-1];
        auto h1 = vector_TH1D[channel-1];
        h->Write();
        h1->Write();
        
    }
	
	f->Close();
	o->Close();
	}


