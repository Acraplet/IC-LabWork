/// PMT Response Function
/// Based on:
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
/// and "Fitting Single Photo-electron peak," Qing He, Princeton, 2010

#include "pmt_response_function.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TText.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include <fstream>

const double pi = acos(-1.0);



/// factorial using TMath::Factorial
double factorial( unsigned n ){
  return TMath::Factorial( n );
}


/// Poisson term
double poisson_term( double mu, unsigned n ){
  return std::pow( mu, n ) * std::exp( -mu ) / factorial(n);
}

/// Gaussian part of signal for n'th photoelectron
/// This gaussian is defined only for x>0 and normalized accordingly
/// such that
///  integra( 0, infty, gaussian ) = 1
double gaussian( double x, double Qn, double sigman ){
  //double norm = 1.0 / ( sigman * std::sqrt( 2*pi ) );
  double norm = sqrt(2/pi) / ( sigman * (1 + std::erf( Qn / ( sqrt(2)*sigman ) ) ) );

  return norm * std::exp( -std::pow( (x-Qn)/(sqrt(2)*sigman), 2 ) );
}

double sign( double x ){
  if (x > 0.) return 1.;
  if (x < 0.) return -1.;
  return 0.;
}

/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// *** Exclude Pedestal term. ***
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q1, charge of single photoelectron
/// p[2] is sigma1, width of single photoelectron gaussian
/// p[3] is mu, the poisson mean number of photoelectrons
/// p[4] is w, the fraction of the background signal that is exponential
/// p[5] is alpha, the exponential constant
double pmtresponse( double * xx, double * p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];

  /// how many terms to include?
  /// go to 3 sigma above mean
  double approx_sigma = std::sqrt( mu );
  unsigned nmax = std::max(4, int(approx_sigma*3+mu) );

  double signal = 0;
  for (unsigned n=1; n<=nmax; ++n){
    double Qn = n * Q1;
    double sigman = std::sqrt(n) * s1;
    double poisson = poisson_term( mu, n );
    double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) );
    double gaus_term = gaussian( x, Qn, sigman );
    double erf1 = std::erf( fabs(n*Q1 + sigman*sigman*alpha) / ( sigman*sqrt(2.) ) );
    double erfarg = x - Qn - sigman*sigman*alpha;
    double erf2 = sign(erfarg)*std::erf( fabs(erfarg) / ( sigman*sqrt(2.) ) );
    double IGn_term = exp_term*( erf1 + erf2 );
    signal += poisson*( (1-w)*gaus_term + w*IGn_term );
  }
  return N*signal;
}

///  Background response
/// Just the background portionof the function
double pmtbackgroundresponse( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];

  /// how many terms to include?
  /// go to 3 sigma above mean
  double approx_sigma = std::sqrt( mu );
  //unsigned nmax = std::max(2, int(approx_sigma*3+mu) );
  unsigned nmax = 4;

  double signal = 0;
  for (unsigned n=1; n<=nmax; n++){
    double Qn = n * Q1;
    double sigman = std::sqrt(n) * s1;
    double poisson = poisson_term( mu, n );
    double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) );
    double gaus_term = gaussian( x, Qn, sigman );
    double erf1 = std::erf( fabs(n*Q1 + sigman*sigman*alpha) / ( sigman*sqrt(2.) ) );
    double erfarg = x - Qn - sigman*sigman*alpha;
    double erf2 = sign(erfarg)*std::erf( fabs(erfarg) / ( sigman*sqrt(2.) ) );
    double IGn_term = exp_term*( erf1 + erf2 );
    signal += poisson* w * IGn_term ;
  }
  return N*signal;
}

/// N'th photoelectron ideal response
/// make seventh parameter N (the photoelectron number)
double pmtG_n( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q1    = p[1];
  double s1    = p[2];
  double mu    = p[3];
  double w     = p[4];
  double alpha = p[5];
  double n     = p[6];

  double signal = 0;

  double Qn = n * Q1;
  double sigman = std::sqrt(n) * s1;
  double poisson = poisson_term( mu, n );
  double exp_term = alpha / 2 * std::exp( -alpha * ( x - Qn - alpha*sigman*sigman/2 ) );
  double gaus_term = gaussian( x, Qn, sigman );
  signal += poisson*( (1-w)*gaus_term );

  return N*signal;
}

/// Return a vector of functions that sums to be approxametely the Bellamy function
/// First term is the background, then 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components( double * p ){
  double pars[ 7 ];
  for ( unsigned i = 0; i < 6; ++i){
    pars[i] = p[i];
  }
  pars[6]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", pmtbackgroundresponse, 0., 5000., 6 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  fcns.push_back( bg );

  double mu=p[3];
  double approx_sigma = std::sqrt( mu );
  unsigned nmax = std::max(5, int(approx_sigma*3+mu) );

  for (unsigned n=1; n<=nmax; ++n ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), pmtG_n, 0., 5000., 7 );
    pars[6] = n;
    ftmp->SetNpx(1000);
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    fcns.push_back( ftmp );
  }
  return fcns;
}



/////****************************************************************
///// Repeat all of above functions but including pedestal
/////****************************************************************
///
/// PMT response Sreal(x) from
/// E.H. Bellamy, et. al., NIM A 339 (1994) 468.
///
/// *** INclude Pedestal term. ***
/// x[0] is charge
/// p[0] is N, overall normalization
/// p[1] is Q0, charge of pedestal
/// p[2] is sigma0, width of pedestal
/// p[3] is Q1, charge of single photoelectron
/// p[4] is sigma1, width of single photoelectron gaussian
/// p[5] is mu, the poisson mean number of photoelectrons
/// p[6] is w, the fraction of the background signal that is exponential
/// p[7] is alpha, the exponential constant
double pmtresponseped( double * xx, double * p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];
  double A2    = p[8];
  double A3    = p[9];

  /// how many terms to include?
  /// go to 3 sigma above mean
  double approx_sigma = std::sqrt( mu );
  //unsigned nmax = std::max(5, int(approx_sigma*3+mu) ); //number of pe? //change 1 for 3
  unsigned nmax = 4;

  double signal = 0;
  //A2 = 1.;
  //A3 = 1.;

  // Use approx equation from eq. 10 of Bellamy
  double exp0term = w * alpha * std::exp( -alpha * ( x - Q0 ) );
  if ( x < Q0 ) exp0term = 0.;
  double gaus0term = (1-w) * gaussian( x, Q0, s0 );
  signal = poisson_term( mu, 0 ) * ( gaus0term + exp0term );


  for (unsigned n=1; n<=nmax; n++){
    double gaus_term = gaussian( x, Q0 + n*Q1 + w/alpha, std::sqrt(n) * s1 );
    if (n==2) signal += A2 * poisson_term(mu,n) * gaus_term ;
    if (n==3) signal += TMath::Abs(A3)* poisson_term(mu,n) * gaus_term ;

    else signal += 1 * poisson_term(mu,n) * gaus_term ;
    //signal += poisson_term(mu,n) * gaus_term ;
  }
  return N*signal;
}

///  Background response *** with pedesal ***
/// Just the background portion of the function, excludes pedestal term
double pmtbackgroundresponseped( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];

  double signal = 0;

  // Use approx equation from eq. 10 of Bellamy
  double exp0term = TMath::Abs(w) * alpha * std::exp( -alpha * ( x - Q0 ) );
  if ( x < Q0 ) exp0term = 0.;
  double gaus0term = (1-w) * gaussian( x, Q0, s0 );
  signal = poisson_term( mu, 0 ) * ( gaus0term + exp0term );
  //std::cout <<"The aplitude " <<  N  << "the vallue of the bg at x = " << x << " is " << N*signal << std::endl;
  return N*signal;
}


/// N'th photoelectron ideal response
/// make ninth parameter n (the photoelectron number)
/// 0 photoelectron term is 1pe
double pmtG_n_ped( double * xx, double *p ){
  double x     = xx[0];
  double N     = p[0];
  double Q0    = p[1];
  double s0    = p[2];
  double Q1    = p[3];
  double s1    = p[4];
  double mu    = p[5];
  double w     = p[6];
  double alpha = p[7];
  unsigned n   = unsigned(p[10]);
  double A2    = p[8];
  double A3    = p[9];  

  double signal = 0;

  //A2 = 1; 
  //A3 = 1;
  double gaus_term = gaussian( x, Q0 + n*Q1 + w/alpha, std::sqrt(n) * s1 );
  //signal += poisson_term(mu,n) * gaus_term ;
  if (n==1) signal += poisson_term(mu,n) * gaus_term ;
  else if (n==3) {
	  signal += TMath::Abs(A3) * poisson_term(mu,n) * gaus_term ;
	  //std::cout << "Position of the 3 pe peak :" << Q0 + n*Q1 + w/alpha << std::endl;
  }
  else if (n==2){
	  signal += A2 * poisson_term(mu,n) * gaus_term ;
	  //std::cout << "Position of the " << n << " pe peak " << Q0 + n*Q1 + w/alpha << " n: " << n << std::endl;
  }
  else signal += poisson_term(mu,n) * gaus_term ;

  return N*signal;
}

/// Return a vector of functions that sums to be approxametely the Bellamy function
/// First term is the background, then 0 pe, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_pmt_response_components_ped( double * p, double hist_binning ){
  double pars[ 11 ];
  double peToVolts = 0.01;
  for ( unsigned i = 0; i < 10; ++i){
    pars[i] = p[i];
  }
  pars[10]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", pmtbackgroundresponseped, -0.1 * peToVolts, 6 * peToVolts, 10 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  std::cout << "Integral BG " << " " << bg->Integral(-0.1, 1)/hist_binning << std::endl;
  //std::cout << "Integral BG " << " " << bg->Integral(-0.1 * peToVolts, 6 * peToVolts) * pars[0] << std::endl;
  fcns.push_back( bg );

  double mu=p[5];
  double approx_sigma = std::sqrt( mu );
  //unsigned nmax = std::max(5, int(approx_sigma*3+mu) );
  unsigned nmax = 4;
  for (unsigned n=1; n<=nmax; n++ ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), pmtG_n_ped, -0.1 * peToVolts, 6 * peToVolts, 11 );
    std::cout << " n " << n << " mu " << mu << " approx sigma " << approx_sigma << std::endl; 
    pars[10] = n;
    ftmp->SetNpx(1000);
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    //std::cout << ftmp->Eval(0.013) << std::endl;
    std::cout << "Integral " << n << " " << ftmp->Integral(-0.1 * peToVolts, 6 * peToVolts)/hist_binning << std::endl;
    fcns.push_back( ftmp );
  }
  return fcns;
}


PMTResponsePed * PMTResponsePed::fInstance = nullptr;
TH1D* PMTResponsePed::fPDF = nullptr;
double PMTResponsePed::fWid = 0;
TFile* PMTResponsePed::fin = nullptr;

void PMTResponsePed::set_pedestal( const std::string& fname, const std::string& histname ){
  if ( fInstance == nullptr ) fInstance = new PMTResponsePed( 0.0 );
  fWid = 0.0;
  if ( fin == nullptr ) fin = new TFile( fname.c_str(), "read" );
  if ( !fin ) {
    std::cout<<"Could not find background root file: "<<fname<<std::endl;
    exit(0);
  }
  fPDF = (TH1D*) fin->Get( histname.c_str() );
  if ( !fPDF ){
    std::cout<<"Could not find histogram: "<<histname<<" in root file: "<<fname<<std::endl;
    exit(0);
  }
}


/*
Model 1:
x[0] = charge
p[0] = normalization (count)
p[1] = q1 single pe mean charge
p[2] = s1 single pe charge width
p[3] = mu mean number of pe
p[4] = w weight of exponential bg (0-1)
p[5] = alpha exponential constant bg
Change background to be part of pedestal signal!
 */
double model1( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];

  double retval = 0.0;

  retval = w * alpha * exp( -alpha*x[0] );

  retval += (1-w)*poisson_term( mu, 0 ) * PMTResponsePed::get_prob_density( x[0] );

  for ( unsigned npe =1; npe<5; ++npe ){
    retval += (1-w)*poisson_term( mu, npe ) * gaussian( x[0], npe*q1, sqrt( npe )*s1 )
      ;
  }

  return retval * N;
}


double model1bg( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];

  double retval = 0.0;

  //double binwid = PMTResponsePed_BinWid::GetBinWid();
  //double expnorm = alpha / (1.0 - exp( -alpha*binwid ) );
  //if ( x[0] < binwid ) retval = (1-w)*poisson_term( mu, 0 )/binwid;


  retval = w * alpha * exp( -alpha*x[0] );
  retval += (1-w)* poisson_term( mu, 0 ) * PMTResponsePed::get_prob_density( x[0] );


  return retval * N;
}


double model1npe( double * x, double *p ){
  double N = p[0];
  double q1 = p[1];
  double s1 = p[2];
  double mu = p[3];
  double w = p[4];
  double alpha = p[5];

  unsigned npe = p[6];

  double retval = 0.0;

  retval += (1-w) * poisson_term( mu, npe ) * gaussian( x[0], npe*q1, sqrt( npe ) * s1 );

  return retval * N;
}



/// Return a vector of functions that sums to be approxametely model1
/// First term is 0 pe + bg, 1 pe, 2 pe, ... etc
std::vector< TF1* > get_model1_components( double * p ){
  double pars[ 7 ];
  for ( unsigned i = 0; i < 7; ++i){
    pars[i] = p[i];
  }
  pars[6]=0.;

  std::vector< TF1* > fcns;
  TF1* bg = new TF1("Background", model1bg, 0., 0.065, 6 );
  bg->SetParameters( pars );
  //bg->SetLineColor( 52 );
  fcns.push_back( bg );

  for (unsigned n=1; n<6; ++n ){
    std::ostringstream fname;
    fname <<"fnpe_"<<n;
    TF1* ftmp = new TF1( fname.str().c_str(), model1npe, 0., 0.065, 7 );
    pars[6] = n;
    ftmp->SetLineWidth(3);
    ftmp->SetParameters( pars );
    ftmp->SetLineColor( 52 + 5*n );
    fcns.push_back( ftmp );
  }
  for ( TF1 * fcn : fcns ) fcn->SetNpx(1000);
  return fcns;
}




/// Test case code to test the fitting ***with pedestal***
//int main(void){
//void pmt_response_function_tweakedCircle_volts(){
int main(int n_argv, const char** args){
    
    const char * PMTform = "None";
    const char * face_position = "None";
    const char * TTS =  "None";
    
  //TFile* fin = new TFile( Form("/home/a/WCTE540/PMT_data/%s/%sPMT/combined_datasets/hadded_all_%s_TTS_%s_1pe_1ms_LandauCorrectedVolts.root", TTS, PMTform, PMTform, face_position), "read" );

  //Get the inputs
  std::string runNumber = args[1]; 
  std::string channelTrigger = args[2];
  std::string laserPower = args[3];
  std::cout << laserPower << std::endl;
  std::string PMT1name = args[4];
  std::string PMT1theta = args[5];
  std::string PMT2name = args[6];
  std::string PMT2theta = args[7];
  std::string PMT3name = args[8];
  std::string PMT3theta = args[9];
  std::string PMT4name = args[10];
  std::string PMT4theta = args[11];
  std::string LBphi = args[12];
  int PMT1number = std::stoi(args[13]);
  int PMT2number = std::stoi(args[14]);
  int PMT3number = std::stoi(args[15]);
  int PMT4number = std::stoi(args[16]);
    
  
  std::vector<int> PMTnumbers;
  PMTnumbers.push_back(PMT1number);
  PMTnumbers.push_back(PMT2number);
  PMTnumbers.push_back(PMT3number);
  PMTnumbers.push_back(PMT4number);

  std::vector<std::string> PMTnames;
  PMTnames.push_back(PMT1name);
  PMTnames.push_back(PMT2name);
  PMTnames.push_back(PMT3name);
  PMTnames.push_back(PMT4name);

  std::vector<std::string> PMTthetas;
  PMTthetas.push_back(PMT1theta);
  PMTthetas.push_back(PMT2theta);
  PMTthetas.push_back(PMT3theta);
  PMTthetas.push_back(PMT4theta);

  TFile* fin = new TFile(Form("../data/landau_run%s_triggerChannel%s_laserPower%s_Ch1%s%s_Ch2%s%s_Ch3%s%s_Ch4%s%s_phi%s.root", runNumber.c_str(), channelTrigger.c_str(), laserPower.c_str(), PMT1name.c_str(), PMT1theta.c_str(), PMT2name.c_str(), PMT2theta.c_str(), PMT3name.c_str(), PMT3theta.c_str(), PMT4name.c_str(), PMT4theta.c_str(), LBphi.c_str()));
  std::cout << "The data that is read in is "<< Form("../data/landau_run%s_triggerChannel%s_laserPower%s_Ch1%s%s_Ch2%s%s_Ch3%s%s_Ch4%s%s_phi%s.root", runNumber.c_str(), channelTrigger.c_str(), laserPower.c_str(), PMT1name.c_str(), PMT1theta.c_str(), PMT2name.c_str(), PMT2theta.c_str(), PMT3name.c_str(), PMT3theta.c_str(), PMT4name.c_str(), PMT4theta.c_str(), LBphi.c_str()) << std::endl;
  std::string runNb = runNumber;

  //We'll have to read the characteristics of the run from reading the txt file 
  //here a loop to read all of the PMTs (that aren't the trigger)
  std::vector<double> list_gains;
  std::vector<double> list_gain_errors;
  std::vector<double> list_chi2perPt;
  std::vector<double> list_Q1;
  std::vector<double> list_1peNumber;
  std::vector<double> list_2peNumber;
  int numberOfPointsPerWindow = 200;
  double samplingFrequency = 2;//in Gsample/s
  double durationWindow = numberOfPointsPerWindow/(samplingFrequency * 1e9);
  double totalDataDuration; // = numberOfLandauFittedPoints*durationWindow; 
  for (int pmt_nb=1; pmt_nb<5; pmt_nb++) {
	  if (pmt_nb != std::stoi(channelTrigger)){
		//get that PMT data 
	  	std::cout << "\n\n\n\n\n\n "<< Form("PMT%d_Corrected", pmt_nb) << " is being processed. " << std::endl;
	  	TH1D * h = (TH1D*)fin->Get(Form("PMT%d_Corrected", pmt_nb));
	 	std::cout << "File is read in" << std::endl;
	  	double gain = -9999;
	  	double peToVolts = 0.01; //an arbitrary scaling that could be cleaned up later on, used for debugging
		
		//Draw the histogram 
	  	TCanvas *c=new TCanvas();
	  	c->cd();
	  	h->SetLineWidth(3);
	  	h->Draw("E1");
          	h->SetTitle(Form("Run%s triggerCh%s laserPow%s Ch1%s%s Ch2%s%s Ch3%s%s Ch4%s%s phi%s", runNumber.c_str(), channelTrigger.c_str(), laserPower.c_str(), PMT1name.c_str(), PMT1theta.c_str(), PMT2name.c_str(), PMT2theta.c_str(), PMT3name.c_str(), PMT3theta.c_str(), PMT4name.c_str(), PMT4theta.c_str(), LBphi.c_str()));
          	h->GetXaxis()->SetTitle("Volts");
          	h->GetYaxis()->SetTitle("Events");
	  	h->SetMarkerStyle(20);
	  	h->SetMarkerSize(0.5);
	
	  	double Nfix = h->Integral( 0, h->GetNbinsX()+1, "width" );
		std::cout << "The bin width is "<< h->GetBinWidth(4) << std::endl;
		double hist_binning = h->GetBinWidth(4);
	  	int numberOfLandauFittedPoints = h->GetEntries();
  	  	totalDataDuration = numberOfLandauFittedPoints*durationWindow; //in units of seconds, needed for calculating frequencies 
	
	 	TF1* ff = new TF1( "pmt_response", pmtresponseped, -0.2 * peToVolts, 6 * peToVolts, 10 );
	
	  	double N0 = h->Integral(0,4);
	 	double Ntot = h->Integral(0,1000);
	  	double Nrest = h->Integral(4,100); // beyond bin 15 dominated by exp background
	  	double height_onepe = h->Integral(5,1000); //4th bin is just above DR
	  	double mufix = Nrest / N0;
	  	mufix = log( mufix + 1 ); // correction factor to go from estimated mu, to true
	  	//double mufix = 5;
	  	std::cout << "Nfix = "<<Nfix<<" Nrest="<<Nrest<<" / N0 = "<<N0<<" = "<<mufix<<std::endl;
	  	std::cout << "mufix = " << mufix << std::endl;
	
	  	ff->SetNpx(1000);
	  	//Best values are square, circle, losange, triangle
	  	double best_Q0[5]= {5.07777e-4, 4.94150e-4, 4.94368e-4, 5.04369e-4, 5.04369e-4};
	  	//double best_s0[4] = {2.72169e-4, 2.7621e-4, 2.72831e-4, 2.77949e-4};
	  	//double best_Q1[4]={1.57682e-2, 9.9694e-3, 1.10866e-2, 1.22630e-2};
	  	double best_Q1[5]={1.9482e-2, 9.094e-3, 1.0866e-2, 1.10866e-2, 5.10866e-3};
	  	double best_s1[5] = {3.43079e-3, 2.82959e-3, 2.81860e-3, 3.5586e-3, 1.5586e-3};
	  	//double best_w[4] = {3.76437e-3, 3.42151e-3, -1.27453e-5, 2.75938e-3};
	  	double best_w[5] = {50, 0.8, 50, 0.5, 5000};
	  	//double best_a[4] = {1.07855e-1, 3.12277e-1, 1.16969e-1, 1.34141e-2};
	  	double best_a[5] = {5e-1, 5e-1, 4e-1, 5e-1, 5e-1};
	  	double best_s0[5] = {0.00018, 0.00046, 0.00039, 0.00055, 0.0004};
	
	  	int PMTnumber = PMTnumbers[pmt_nb-1]-1; //-1 free fit, 0 Square fixed, 1 Circle, 2 Losange. 3 Tirangle
	
	  	//dec ent-ish guess for square
	  	ff->SetParNames("N", "Q_{0}", "#sigma_{0}", "Q_{1}", "#sigma_{1}", "#mu", "w", "#alpha", "A2", "A3" );
	  	ff->SetParameters( Nfix, 0.00005, best_s0[PMTnumber], best_Q1[PMTnumber], best_s1[PMTnumber], mufix, best_w[PMTnumber], best_a[PMTnumber], 1, 1  );
		std::cout << ff->GetParameter(2) << " is s0 for PMT " << pmt_nb << std::endl;
	  	ff->FixParameter( 0, Nfix );
	  	//ff->FixParameter( 1, ff->GetParameter(1) );
	  	//ff->FixParameter( 2, ff->GetParameter(2));
	  	//ff->FixParameter( 3, ff->GetParameter(3) );
	  	//ff->FixParameter( 4, ff->GetParameter(4) );
	  	ff->FixParameter( 8, ff->GetParameter(8) ) ;
	  	ff->FixParameter( 9, ff->GetParameter(9) );
	  	//ff->FixParameter( 6, ff->GetParameter(6) ) ;
	  	//ff->FixParameter( 7, ff->GetParameter(7) ) ;
	  	ff->SetParLimits( 6, 0, 1 );
	  	ff->SetParLimits( 3, 0, 0.02 );
	  	ff->SetParLimits( 5, 0, 10 );
	  	ff->SetParLimits( 7, 0, 100 );
	
	  	std::cout << "\n The First fit I think " << std::endl;	
	  	h->Fit(ff, "", "", 0.1 * peToVolts, 6 * peToVolts ); //background at the end?
	
	  	for ( unsigned ipar=0 ; ipar < 10; ++ipar ) ff->ReleaseParameter(ipar);
	  	ff->FixParameter( 0, ff->GetParameter(0) );
	  	//ff->FixParameter( 1, ff->GetParameter(1) );
	  	//ff->FixParameter( 3, ff->GetParameter(3) ) ;
	  	ff->FixParameter( 2, ff->GetParameter(2) ) ;
	  	ff->FixParameter( 4, ff->GetParameter(4) ) ;
	  	ff->SetParLimits( 6, 1e-7, 1e-4 );
	  	ff->SetParLimits( 1, -0.004, 0.004 );
	  	ff->SetParLimits( 5, 0, 10 );
	  	ff->SetParLimits( 7, 0, 50 );
	  	ff->SetParLimits( 8, 0, 1000 );
	  	//ff->SetParLimits( 9, 0, 50 );
	  	ff->SetParLimits( 3, ff->GetParameter(3) * 0.8, ff->GetParameter(3)*1.2 );
	
	  	ff->SetLineWidth(3);
	  	ff->SetLineColor(kRed+2);
	
	  	std::cout << "\n The whole fit I think " << std::endl;	
	  	TFitResultPtr r = h->Fit(ff, "S", "",-0.1 * peToVolts, 6 * peToVolts ); //whole fit? here in  pe
	//	These are all of the n-pe responses
	  	std::vector< TF1* > fcmp = get_pmt_response_components_ped( ff->GetParameters(), hist_binning );
	  	//std::vector< TF1* > fcmp_err = ff->GetErrors();
	  	TF1* bg = fcmp[0];
	  	int count = 0;
	  	//on changing this.    
	  	c->cd();
	  	double totalNumberPE = 0;
	  	double totalNumberTwoPE = 0;
	  	double totalNumberPedestal = bg->Integral(-0.001,1)/hist_binning;
	  	std::cout << "Number of counts in the pedestal "<< totalNumberPedestal<< std::endl;
	  	for ( TF1* f : fcmp ){
	    		f->Draw("same");
	    		std::cout << "Integral of this function " << count << " is "<< f->Integral(-0.015, 10)/hist_binning << std::endl;
	    	  if (count == 1){
		  	totalNumberPE = f->Integral(-0.001, 10)/hist_binning;
	    	  }
	    	  if (count == 2){
		  	totalNumberTwoPE = f->Integral(-0.001, 10)/hist_binning;
	    	  }
	    	  //totalNumberPedestal = f->Integral(0, 1000)/hist_binning;	
	    	  f->GetXaxis()->SetTitle("volts");
	    	  count += 1;
	  	}
	  
	  	c->cd();
	
	  	double Q0 = ff->GetParameter(1);
	  	double Q0_err = r->ParError(1);
	  	double Q1 = ff->GetParameter(3);
	  	double Q1_err = r->ParError(3);
	  	double chi2 = r->Chi2();
	  	double Q_err = TMath::Sqrt(Q0_err*Q0_err + Q1_err*Q1_err);
	  	double R = 50.;
	  	double Tr = 0.8e-9;
	  	double e = 1.6e-19;
	  	double gain_guess = (((Q1-Q0) * Tr)/R)/e ;
	  	double gain_error = ((Q_err * Tr)/R)/e;
	
	  	std::cout << PMTform << " " << face_position  << " Q0 = " << Q0 << " +/- " << Q0_err << " volts Q1 = " << Q1 << " +/- " << Q1_err << " volts gain = " << gain_guess << " +/- " << gain_error << std::endl;
	  	
	  	TPad *newpad=new TPad("newpad","a transparent pad",0.04,0.22,0.87,0.90);
	  	newpad->SetFillStyle(4000);
	  	newpad->Draw();
	  	newpad->cd();
	  	TText *t = new TText(.5,.5,Form("Run %s  PMT %s theta %s phi %s laser current %smA", runNumber.c_str(), PMTnames[pmt_nb-1].c_str(), PMTthetas[pmt_nb-1].c_str(), LBphi.c_str(), laserPower.c_str()));
	  	t->SetX(0.14107383);
	  	t->SetY(0.950383);
	  	t->SetTextColor(1);
	  	t->SetTextFont(22);
	  	t->SetTextSize(0.05);
	  	t->SetTextAngle(0);
	  	t->Draw();
	  	TText *t2 = new TText(.5,.5,Form("Counts in pedestal: %.3e, 1pe peak: %.0f, 2pe peak: %.0f",totalNumberPedestal, totalNumberPE, totalNumberTwoPE));
	  	t2->SetX(0.13107383);
	  	t2->SetY(0.900383);
	  	t2->SetTextColor(1);
	  	t2->SetTextFont(22);
	  	t2->SetTextSize(0.05);
	  	t2->SetTextAngle(0);
	  	t2->Draw();
	  	//in landau code we take one maximum per window of 200points
	  	std::cout << "totalDataDuration " << totalDataDuration << std::endl;

	  	TText *t3 = new TText(.5,.5,Form("1PE event frequency  %.2f kHz, 2PE event frequency %.2f kHz",totalNumberPE/totalDataDuration/1000, totalNumberTwoPE/totalDataDuration/1000));
	  	t3->SetX(0.14107383);
	  	t3->SetY(0.850383);
	  	t3->SetTextColor(1);
	  	t3->SetTextFont(22);
	  	t3->SetTextSize(0.05);
	  	t3->SetTextAngle(0);
	  	t3->Draw();
	  	      
	  	 TText *t4 = new TText(.5,.5,Form("Gain %.3e +/- %.3e \n chi2 per pt = %.4f", PMTform, face_position, gain_guess, gain_error, chi2/Ntot ));
	  	t4->SetX(0.14107383);
	  	t4->SetY(0.800383);
	  	t4->SetTextColor(1);
	  	t4->SetTextFont(22);
	  	t4->SetTextSize(0.05);
	  	t4->SetTextAngle(0);
	  	t4->Draw();
	
	  	h->GetXaxis()->SetRangeUser(-0.1 * peToVolts, 6 * peToVolts);
	  	h->GetYaxis()->SetRangeUser(1, N0*1.5);
	  	c->SetLogy();
	  	c->SaveAs(Form("./Bellamy_plots/BellamyFit_Run%s_PMT%d.pdf", runNb.c_str(), pmt_nb));
	  	std::cout << Form("%s %s - gain %.3e +/- %.3e \n chi2 per pt = %.4f", PMTform, face_position, gain_guess, gain_error, chi2/Ntot ) << std::endl;

 		list_gains.push_back(gain_guess);
		list_gain_errors.push_back(gain_error);
        	list_chi2perPt.push_back(chi2/Ntot);
		list_Q1.push_back(Q1);
		list_1peNumber.push_back(totalNumberPE);
		list_2peNumber.push_back(totalNumberTwoPE);
	} //for everything that is not the trigger
	else { //for the trigger
		list_gains.push_back(-9999);
        	list_gain_errors.push_back(-9999);
        	list_chi2perPt.push_back(-9999);
        	list_Q1.push_back(-9999);
        	list_1peNumber.push_back(-9999);
        	list_2peNumber.push_back(-9999);
		std::cout << "Not processing the trigger channel which is channel " << channelTrigger << std::endl;
	}
	
	}//we have done all of the PMTs in that file (except the trigger)
	

	//read in the TTS from the other files 
	std::vector<double> TTSresults;
	std::vector<double> list_FWHM;
	std::vector<double> list_FWHM_err;
	std::vector<double> list_nHits_TTS;
	std::vector<double> list_nHits_TTS_err;
	std::vector<double> list_nHits_DR;
	std::vector<double> list_nHits_DR_err;
	if (std::stoi(channelTrigger) != 0){
		TFile* finTTS = new TFile(Form("../data/TTS_run%s_triggerChannel%s_laserPower%s_Ch1%s%s_Ch2%s%s_Ch3%s%s_Ch4%s%s_phi%s.root", runNumber.c_str(), channelTrigger.c_str(), laserPower.c_str(), PMT1name.c_str(), PMT1theta.c_str(), PMT2name.c_str(), PMT2theta.c_str(), PMT3name.c_str(), PMT3theta.c_str(), PMT4name.c_str(), PMT4theta.c_str(), LBphi.c_str()));	
		for (int channel = 1; channel < 5; channel ++){
			if (channel!=channelTrigger){ // break;
			TH1D * hTTS = (TH1D*)finTTS->Get(Form("hTTS_PMT%d", channel));
			double TTS_binning = hTTS->GetBinWidth(3);
			double binMin = hTTS->GetXaxis()->GetBinLowEdge(0);
			double binMax = hTTS->GetXaxis()->GetBinUpEdge(hTTS->GetNbinsX());
			TF1 *fitFunction = new TF1("fitFunction", "landau + [3]", binMin, binMax);
			fitFunction->SetNpx(10000);
			fitFunction->SetParameters(hTTS->GetBinContent(hTTS->GetMaximumBin()), hTTS->GetBinCenter(hTTS->GetMaximumBin()), 0.4, 0.55);
			fitFunction->SetParLimits(2, 0, 5);
			//TF1 *fitFunction = new TF1("fitFunction", "[0]/2 * [2] *exp(2*[1] + [2] * [3]*[3]- 2 * x) *ROOT::Math::erfc((([1]+[2] * [3]*[3]-x)/(sqrt(2.0)*[3]))/ROOT::Math::erfc(x)) + [4]", binMin, binMax);
			//fitFunction->SetParameters(hTTS->GetBinContent(hTTS->GetMaximumBin()), hTTS->GetBinCenter(hTTS->GetMaximumBin()), 0.4, 5, 4, 0.55);
			double amplitude = hTTS->GetBinContent(hTTS->GetMaximumBin());
			double mean_guess = hTTS->GetBinCenter(hTTS->GetMaximumBin());
			std::cout << mean_guess << std::endl; 
			double sigma_guess = 1.2;                    
			double decay_rate_exp = 0.05;      
			//fitFunction->SetParameters(2*amplitude, mean_guess, sigma_guess, decay_rate_exp, 0);
			//fitFunction->SetParLimits(1, mean_guess*0.9, mean_guess*1.1);
			fitFunction->SetParLimits(0, amplitude*0.9, amplitude*50);
			fitFunction->SetParLimits(4, 0, 2);
			//fitFunction->SetParLimits(2, 0.3, 1.2);

			hTTS->Fit("fitFunction", "B");
			hTTS->GetYaxis()->SetRangeUser(0, amplitude * 1.7);
			double mean = fitFunction->GetParameter(1);
			double sigma = fitFunction->GetParameter(2);
			double mean_err = fitFunction->GetParError(1);
			double sigma_err = fitFunction->GetParError(2);
			TCanvas *cTTS=new TCanvas();
			cTTS->SetMargin(0.12, 0.1, 0.1, 0.1);
          		cTTS->cd();
          		hTTS->Draw("E1");
			fitFunction-> Draw("same");
			fitFunction->SetLineColor(kRed);
			TPad *newpad=new TPad("newpad","a transparent pad",0.04,0.22,0.87,0.90);
          		newpad->SetFillStyle(4000);
          		newpad->Draw();
          		newpad->cd();
          		TText *tTTS = new TText(.5,.5,Form("Run %s  PMT %s theta %s phi %s laser current %smA", runNumber.c_str(), PMTnames[channel-1].c_str(), PMTthetas[channel-1].c_str(), LBphi.c_str(), laserPower.c_str()));
        	  	tTTS->SetX(0.14107383);
        	   	tTTS->SetY(0.950383);
        	    	tTTS->SetTextColor(1);
        	    	tTTS->SetTextFont(22);
        	    	tTTS->SetTextSize(0.05);
        	    	tTTS->SetTextAngle(0);
        	    	tTTS->Draw();
          		TText *tTTS2 = new TText(.5,.5,Form("Mean delay: %.1f+/-%.1fns, TTS (FWHM): %.2f +/- %.2fns", mean, mean_err, sigma*2, sigma_err*2));
			//need to multiply the landau sigma by two to get the FWHM 
        	  	tTTS2->SetX(0.14107383);
        	   	tTTS2->SetY(0.900383);
        	    	tTTS2->SetTextColor(1);
        	    	tTTS2->SetTextFont(22);
        	    	tTTS2->SetTextSize(0.05);
        	    	tTTS2->SetTextAngle(0);
        	    	tTTS2->Draw();
			//double nHits_TTS = fitFunction->Integral(mean-7*sigma, mean+7*sigma)/TTS_binning;

			//need to get the number of hits in the TTS and then in the DR
			TF1 *landauTTS = new TF1("fitFunction", "landau", binMin, binMax);
			landauTTS->SetParameters(fitFunction->GetParameter(0), fitFunction->GetParameter(1), fitFunction->GetParameter(2));
			double nHits_TTS = landauTTS->Integral(binMin, binMax)/TTS_binning;
			double nHits_DR = hTTS->Integral(0, (binMax-binMin)/TTS_binning) - nHits_TTS;
			double err_TTS = TMath::Sqrt(nHits_TTS);
			double err_DR = TMath::Sqrt(nHits_TTS + hTTS->Integral(0, (binMin-binMax)/TTS_binning));
			if (mean > 10 and sigma > 0.1){
				std::cout << "YEAYAAYAYAY " << std::endl;
				//If we actually have a resonable signal we integrate from teh 10th bin to 10sigma before the mean signal time to 
				//make the average number of DR hits per bin (better than integrating) and then we extrapolate to the end, need to have the Bin width is the delay between two triggers
				//double meanDRhitsPerBins = hTTS->Integral(10, (binMin + mean - 10 * sigma)/TTS_binning)/((binMin + mean - 10 * sigma)/TTS_binning - 10);
				int nBinsDRaverage = (mean - 15 - binMin)/TTS_binning;
				double totalNumberOfWindows = totalDataDuration/((binMax-binMin) * 10e-9);
				double meanDRhitsPerBins = hTTS->Integral(10, nBinsDRaverage)/(nBinsDRaverage);
				std::cout << "nBinsDRaverage " << nBinsDRaverage << " and the position of the cut is " << mean - 15  - binMin  << " bin width is " << TTS_binning << "\n The average number of hits per bin is " << meanDRhitsPerBins << " which implies a DR of " << meanDRhitsPerBins * (binMax-binMin)/(totalDataDuration * 1000 * TTS_binning)  << "Hz. " << " Total number of bins " << (binMax-binMin)/TTS_binning  << " and the total nmber of hits in first half of DR: " <<  hTTS->Integral(10, nBinsDRaverage/2) << " vs in the second half " << hTTS->Integral(nBinsDRaverage/2, nBinsDRaverage-10)<< std::endl;
				nHits_DR = meanDRhitsPerBins * (binMax-binMin)/TTS_binning;
				nHits_TTS = hTTS->Integral(0, (binMax-binMin)/TTS_binning) - nHits_DR;
				err_DR = TMath::Sqrt(meanDRhitsPerBins * nBinsDRaverage) * (binMax-binMin)/(mean - 15  - binMin);
				err_TTS = TMath::Sqrt(err_DR*err_DR + hTTS->Integral(0, (binMax-binMin)/TTS_binning));
			}

          		TText *tTTS3 = new TText(.5,.5,Form("#hits in TTS peak: %.0f+/-%.1f (%.2f+/-%.2f kHz)", nHits_TTS, err_TTS, nHits_TTS/(totalDataDuration*1000), err_TTS/(totalDataDuration*1000)));
        	  	tTTS3->SetX(0.14107383);
        	   	tTTS3->SetY(0.850383);
        	    	tTTS3->SetTextColor(1);
        	    	tTTS3->SetTextFont(22);
        	    	tTTS3->SetTextSize(0.05);
        	    	tTTS3->SetTextAngle(0);
        	    	tTTS3->Draw();
          		TText *tTTS4 = new TText(.5,.5,Form("#hits in DR: %.0f +/-%.1f (%.2f +/- %.2fkHz)", nHits_DR, err_DR,  nHits_DR/(totalDataDuration*1000), err_DR/(totalDataDuration*1000)));
        	  	tTTS4->SetX(0.14107383);
        	   	tTTS4->SetY(0.800383);
        	    	tTTS4->SetTextColor(1);
        	    	tTTS4->SetTextFont(22);
        	    	tTTS4->SetTextSize(0.05);
        	    	tTTS4->SetTextAngle(0);
        	    	tTTS4->Draw();
			cTTS->Draw();
			cTTS->SaveAs(Form("./TTS_plots/TTS_Run%s_PMT%d.pdf", runNb.c_str(), channel));	
			hTTS->GetXaxis()->SetRangeUser(mean-35*sigma, mean+ 70 * sigma);
			cTTS->SaveAs(Form("./TTS_plots/TTS_Run%s_PMT%d_peak.pdf", runNb.c_str(), channel));	
			list_nHits_TTS.push_back(nHits_TTS);
			list_nHits_TTS_err.push_back(err_TTS);
			list_nHits_DR.push_back(nHits_DR);
			list_nHits_DR_err.push_back(err_DR);
			list_FWHM.push_back(sigma*2);
			// list_FWHM.push_back(sigma*2);
			list_FWHM_err.push_back(sigma_err*2);
			}

			else{//the trigger channel if we do have the TTS
			list_nHits_TTS.push_back(-9999);
			list_nHits_TTS_err.push_back(-9999);
			list_nHits_DR.push_back(-9999);
			list_nHits_DR_err.push_back(-9999);
			list_FWHM.push_back(-9999);
			list_FWHM_err.push_back(-9999);
			}
		}
		}
		
		else { //if we do not have TTS at all 	
		for (int channel = 1; channel < 5; channel ++){
			list_nHits_TTS.push_back(-9999);
			list_nHits_TTS_err.push_back(-9999);
			list_nHits_DR.push_back(-9999);
			list_nHits_DR_err.push_back(-9999);
                	list_FWHM.push_back(-9999);
                	list_FWHM_err.push_back(-9999);
		}
		}
	//output
	std::ofstream outputFile;
	outputFile.open(Form("output_run%s.txt", runNumber.c_str()), std::ios::app);
	if (outputFile.is_open()) {
	outputFile << runNumber << " " << channelTrigger << " " << laserPower << " " << PMT1name << " " << PMT1number << " " << PMT1theta << " " << PMT2name << " " << PMT2number << " " << PMT2theta << " "<< PMT3name << " " << PMT3number << " " << PMT3theta << " " << PMT4name << " " << PMT4number << " " << PMT4theta << " "<< LBphi << " " << list_gains[0] << " " << list_gains[1] << " " << list_gains[2] << " " << list_gains[3] << " " << list_gain_errors[0] << " " << list_gain_errors[1] << " " << list_gain_errors[2] << " " << list_gain_errors[3] << " " << list_Q1[0] << " " << list_Q1[1] << " " << list_Q1[2] << " " << list_Q1[3] << " " << list_1peNumber[0] << " " << list_1peNumber[1] <<" " << list_1peNumber[2] <<" " << list_1peNumber[3] << " " << list_2peNumber[0] << " " << list_2peNumber[1] <<" " << list_2peNumber[2] <<" " << list_2peNumber[3] << " " << list_chi2perPt[0] << " " << list_chi2perPt[1] << " " << list_chi2perPt[2] << " " << list_chi2perPt[3] << " " << totalDataDuration << " " << list_nHits_TTS[0] << " " << list_nHits_TTS[1] <<" " << list_nHits_TTS[2] <<" " << list_nHits_TTS[3] << " " << list_FWHM[0] <<  " " << list_FWHM[1] << " " << list_FWHM[2] << " " << list_FWHM[3] << " " << list_FWHM_err[0] <<" " << list_FWHM_err[1] << " " << list_FWHM_err[2] << " " << list_FWHM_err[3] <<  std::endl;
	outputFile.close();
	} 
	std::ofstream outputFile2;
	outputFile2.open(Form("general_output_new.txt"), std::ios::app);
	if (outputFile2.is_open()) {
	outputFile2 << runNumber << " " << channelTrigger << " " << laserPower << " " << PMT1name << " " << PMT1number << " " << PMT1theta << " " << PMT2name << " " << PMT2number << " " << PMT2theta << " "<< PMT3name << " " << PMT3number << " " << PMT3theta << " " << PMT4name << " " << PMT4number << " " << PMT4theta << " "<< LBphi << " " << list_gains[0] << " " << list_gains[1] << " " << list_gains[2] << " " << list_gains[3] << " " << list_gain_errors[0] << " " << list_gain_errors[1] << " " << list_gain_errors[2] << " " << list_gain_errors[3] << " " << list_Q1[0] << " " << list_Q1[1] << " " << list_Q1[2] << " " << list_Q1[3] << " " << list_1peNumber[0] << " " << list_1peNumber[1] <<" " << list_1peNumber[2] <<" " << list_1peNumber[3] << " " << list_2peNumber[0] << " " << list_2peNumber[1] <<" " << list_2peNumber[2] <<" " << list_2peNumber[3] << " " << list_chi2perPt[0] << " " << list_chi2perPt[1] << " " << list_chi2perPt[2] << " " << list_chi2perPt[3] << " " << totalDataDuration << " " << list_nHits_TTS[0] << " " << list_nHits_TTS[1] <<" " << list_nHits_TTS[2] <<" " << list_nHits_TTS[3] << " " << list_FWHM[0] <<  " " << list_FWHM[1] << " " << list_FWHM[2] << " " << list_FWHM[3] << " " << list_FWHM_err[0] <<" " << list_FWHM_err[1] << " " << list_FWHM_err[2] << " " << list_FWHM_err[3] << std::endl;
	outputFile2.close();

	//here do the acrylic rod scan
	if (std::stoi(runNumber.c_str()) >= 76 and std::stoi(runNumber.c_str()) <= 93){
	std::ofstream outputFile3;
	outputFile3.open(Form("acrylicRodScan.txt"), std::ios::app);
	if (outputFile3.is_open()) {
	outputFile3 << runNumber << " " << channelTrigger << " " << laserPower << " " << PMT1number << " " << PMT1theta << " " <<  PMT2number << " " << PMT2theta << " "<< PMT3number << " " << PMT3theta << " " << PMT4name << " " <<  PMT4theta << " "<< LBphi << " " << list_gains[0] << " " << list_gains[1] << " " << list_gains[2] << " " << list_gains[3] << " " << list_gain_errors[0] << " " << list_gain_errors[1] << " " << list_gain_errors[2] << " " << list_gain_errors[3] << " " << list_Q1[0] << " " << list_Q1[1] << " " << list_Q1[2] << " " << list_Q1[3] << " " << list_1peNumber[0] << " " << list_1peNumber[1] <<" " << list_1peNumber[2] <<" " << list_1peNumber[3] << " " << list_2peNumber[0] << " " << list_2peNumber[1] <<" " << list_2peNumber[2] <<" " << list_2peNumber[3] << " " << list_chi2perPt[0] << " " << list_chi2perPt[1] << " " << list_chi2perPt[2] << " " << list_chi2perPt[3] << " " << totalDataDuration << " " << list_nHits_TTS[0] << " " << list_nHits_TTS[1] <<" " << list_nHits_TTS[2] <<" " << list_nHits_TTS[3] << " " << list_nHits_TTS_err[0] << " " << list_nHits_TTS_err[1] <<" " << list_nHits_TTS_err[2] <<" " << list_nHits_TTS_err[3] << " " << list_nHits_DR[0] << " " << list_nHits_DR[1] <<" " << list_nHits_DR[2] <<" " << list_nHits_DR[3] << " " << list_nHits_DR_err[0] << " " << list_nHits_DR_err[1] <<" " << list_nHits_DR_err[2] <<" " << list_nHits_DR_err[3] << " "<< list_FWHM[0] <<  " " << list_FWHM[1] << " " << list_FWHM[2] << " " << list_FWHM[3] << " " << list_FWHM_err[0] <<" " << list_FWHM_err[1] << " " << list_FWHM_err[2] << " " << list_FWHM_err[3] << std::endl;
        outputFile3.close();
	}
	}
	if (std::stoi(runNumber.c_str()) >= 94 ){
	std::ofstream outputFile4;
	outputFile4.open(Form("LaserBallPrototype1Scan.txt"), std::ios::app);
	if (outputFile4.is_open()) {
	outputFile4 << runNumber << " " << channelTrigger << " " << laserPower << " " << PMT1number << " " << PMT1theta << " " <<  PMT2number << " " << PMT2theta << " "<< PMT3number << " " << PMT3theta << " " << PMT4name << " " <<  PMT4theta << " "<< LBphi << " " << list_gains[0] << " " << list_gains[1] << " " << list_gains[2] << " " << list_gains[3] << " " << list_gain_errors[0] << " " << list_gain_errors[1] << " " << list_gain_errors[2] << " " << list_gain_errors[3] << " " << list_Q1[0] << " " << list_Q1[1] << " " << list_Q1[2] << " " << list_Q1[3] << " " << list_1peNumber[0] << " " << list_1peNumber[1] <<" " << list_1peNumber[2] <<" " << list_1peNumber[3] << " " << list_2peNumber[0] << " " << list_2peNumber[1] <<" " << list_2peNumber[2] <<" " << list_2peNumber[3] << " " << list_chi2perPt[0] << " " << list_chi2perPt[1] << " " << list_chi2perPt[2] << " " << list_chi2perPt[3] << " " << totalDataDuration << " " << list_nHits_TTS[0] << " " << list_nHits_TTS[1] <<" " << list_nHits_TTS[2] <<" " << list_nHits_TTS[3] << " " << list_nHits_TTS_err[0] << " " << list_nHits_TTS_err[1] <<" " << list_nHits_TTS_err[2] <<" " << list_nHits_TTS_err[3] << " " << list_nHits_DR[0] << " " << list_nHits_DR[1] <<" " << list_nHits_DR[2] <<" " << list_nHits_DR[3] << " " << list_nHits_DR_err[0] << " " << list_nHits_DR_err[1] <<" " << list_nHits_DR_err[2] <<" " << list_nHits_DR_err[3] << " "<< list_FWHM[0] <<  " " << list_FWHM[1] << " " << list_FWHM[2] << " " << list_FWHM[3] << " " << list_FWHM_err[0] <<" " << list_FWHM_err[1] << " " << list_FWHM_err[2] << " " << list_FWHM_err[3] << std::endl;
        outputFile4.close();
	}
	}


 
	std::cout << "Data successfully written to the file." << std::endl;
	}
	else { std::cerr << "Error opening the file." << std::endl;}
	
  return 0;
}
