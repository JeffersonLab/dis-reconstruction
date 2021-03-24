R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

void kinematic_resolution_ycut(){

  int energy_set(0);

  cout << "Choose Which Beam Energies for ep" << endl;
  cout << "1) e = 5 GeV, p = 41 GeV" << endl;
  cout << "2) e = 18 GeV, p = 275 GeV"<< endl;
  cin >> energy_set;

  //Constants                                                                                                        
  const double alpha(1./137.036);
  double s_cm(0);
  const double fbgev(1./(0.3894E12));
  const double Mp(0.9383);

  if(energy_set == 1) s_cm = 4.*5*41;
  if(energy_set == 2) s_cm = 4.*18*275;

  //We will apply cuts on the true y in this version
  double y_max = 0.95;double y_min = 1e-2;

  //Histograms

  //Set Style
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  //gStyle->SetLabelOffset(0.01,"X");
  //gStyle->SetLabelOffset(0.01,"Y");
  gStyle->SetTitleXSize(0.04);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYSize(0.04);
  gStyle->SetTitleYOffset(0.9);

  //x Binning                                                                                                                       
  double x_min = 1E-4;
  double x_max = 1;
  const int nbins_x = 100;
  double log_bw_x = (log10(x_max) - log10(x_min))/(nbins_x); //Determine bin width                                                             
  double log_x_div;
  double x_bins[nbins_x+1];
  for(int i=0;i<nbins_x+1;i++){
    log_x_div = log10(x_min) + (i*log_bw_x);
    x_bins[i] = pow(10,log_x_div);
  }

  //Electron Method using ecal energy
  TH2 *h1a_1 = new TH2D("h1a_1","Q^{2} Resolution vs. Q^{2}",200,0,200,200,-10,10);
  h1a_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h1a_1->GetXaxis()->CenterTitle();
  h1a_1->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h1a_1->GetYaxis()->CenterTitle();
  TH2 *h1a_2 = new TH2D("h1a_2","Q^{2} Resolution vs. y",200,0,1,200,-10,10);
  h1a_2->GetXaxis()->SetTitle("True y");h1a_2->GetXaxis()->CenterTitle();
  h1a_2->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h1a_2->GetYaxis()->CenterTitle();
  TH2 *h1a_3 = new TH2D("h1a_3","Q^{2} Resolution vs. x",100,x_bins,200,-10,10);
  h1a_3->GetXaxis()->SetTitle("True x");h1a_3->GetXaxis()->CenterTitle();
  h1a_3->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h1a_3->GetYaxis()->CenterTitle();

  TH2 *h1b_1 = new TH2D("h1b_1","y Resolution vs. Q^{2}",200,0,200,200,-20,20);
  h1b_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h1b_1->GetXaxis()->CenterTitle();
  h1b_1->GetYaxis()->SetTitle("Difference from true y [%]");h1b_1->GetYaxis()->CenterTitle();
  TH2 *h1b_2 = new TH2D("h1b_2","y Resolution vs. y",200,0,1,200,-20,20);
  h1b_2->GetXaxis()->SetTitle("True y");h1b_2->GetXaxis()->CenterTitle();
  h1b_2->GetYaxis()->SetTitle("Difference from true y [%]");h1b_2->GetYaxis()->CenterTitle();
  TH2 *h1b_3 = new TH2D("h1b_3","y Resolution vs. x",100,x_bins,200,-20,20);
  h1b_3->GetXaxis()->SetTitle("True x");h1b_3->GetXaxis()->CenterTitle();
  h1b_3->GetYaxis()->SetTitle("Difference from true y [%]");h1b_3->GetYaxis()->CenterTitle();

  TH2 *h1c_1 = new TH2D("h1c_1","x Resolution vs. Q^{2}",200,0,200,200,-20,20);
  h1c_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h1c_1->GetXaxis()->CenterTitle();
  h1c_1->GetYaxis()->SetTitle("Difference from true x [%]");h1c_1->GetYaxis()->CenterTitle();
  TH2 *h1c_2 = new TH2D("h1c_2","x Resolution vs. y",200,0,1,200,-20,20);
  h1c_2->GetXaxis()->SetTitle("True y");h1c_2->GetXaxis()->CenterTitle();
  h1c_2->GetYaxis()->SetTitle("Difference from true x [%]");h1c_2->GetYaxis()->CenterTitle();
  TH2 *h1c_3 = new TH2D("h1c_3","x Resolution vs. x",100,x_bins,200,-20,20);
  h1c_3->GetXaxis()->SetTitle("True x");h1c_3->GetXaxis()->CenterTitle();
  h1c_3->GetYaxis()->SetTitle("Difference from true x [%]");h1c_3->GetYaxis()->CenterTitle();

  TH2 *h1d_1 = new TH2D("h1d_1","Reconstructed vs. True Q^{2}",200,0,200,200,0,200);
  h1d_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h1d_1->GetXaxis()->CenterTitle();
  h1d_1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");h1d_1->GetYaxis()->CenterTitle();
  TH2 *h1d_2 = new TH2D("h1d_2","Reconstructed vs. True y",200,0,1,200,0,1);
  h1d_2->GetXaxis()->SetTitle("True y");h1d_2->GetXaxis()->CenterTitle();
  h1d_2->GetYaxis()->SetTitle("Reconstructed y");h1d_2->GetYaxis()->CenterTitle();
  TH2 *h1d_3 = new TH2D("h1d_3","Reconstructed vs. True x",100,x_bins,100,x_bins);
  h1d_3->GetXaxis()->SetTitle("True x");h1d_3->GetXaxis()->CenterTitle();
  h1d_3->GetYaxis()->SetTitle("Reconstructed x");h1d_3->GetYaxis()->CenterTitle();

  //Electron Method using momentum as energy
  TH2 *h2a_1 = new TH2D("h2a_1","Q^{2} Resolution vs. Q^{2}",200,0,200,200,-10,10);
  h2a_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h2a_1->GetXaxis()->CenterTitle();
  h2a_1->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h2a_1->GetYaxis()->CenterTitle();
  TH2 *h2a_2 = new TH2D("h2a_2","Q^{2} Resolution vs. y",200,0,1,200,-10,10);
  h2a_2->GetXaxis()->SetTitle("True y");h2a_2->GetXaxis()->CenterTitle();
  h2a_2->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h2a_2->GetYaxis()->CenterTitle();
  TH2 *h2a_3 = new TH2D("h2a_3","Q^{2} Resolution vs. x",100,x_bins,200,-10,10);
  h2a_3->GetXaxis()->SetTitle("True x");h2a_3->GetXaxis()->CenterTitle();
  h2a_3->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h2a_3->GetYaxis()->CenterTitle();

  TH2 *h2b_1 = new TH2D("h2b_1","y Resolution vs. Q^{2}",200,0,200,200,-20,20);
  h2b_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h2b_1->GetXaxis()->CenterTitle();
  h2b_1->GetYaxis()->SetTitle("Difference from true y [%]");h2b_1->GetYaxis()->CenterTitle();
  TH2 *h2b_2 = new TH2D("h2b_2","y Resolution vs. y",200,0,1,200,-20,20);
  h2b_2->GetXaxis()->SetTitle("True y");h2b_2->GetXaxis()->CenterTitle();
  h2b_2->GetYaxis()->SetTitle("Difference from true y [%]");h2b_2->GetYaxis()->CenterTitle();
  TH2 *h2b_3 = new TH2D("h2b_3","y Resolution vs. x",100,x_bins,200,-20,20);
  h2b_3->GetXaxis()->SetTitle("True x");h2b_3->GetXaxis()->CenterTitle();
  h2b_3->GetYaxis()->SetTitle("Difference from true y [%]");h2b_3->GetYaxis()->CenterTitle();

  TH2 *h2c_1 = new TH2D("h2c_1","x Resolution vs. Q^{2}",200,0,200,200,-20,20);
  h2c_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h2c_1->GetXaxis()->CenterTitle();
  h2c_1->GetYaxis()->SetTitle("Difference from true x [%]");h2c_1->GetYaxis()->CenterTitle();
  TH2 *h2c_2 = new TH2D("h2c_2","x Resolution vs. y",200,0,1,200,-20,20);
  h2c_2->GetXaxis()->SetTitle("True y");h2c_2->GetXaxis()->CenterTitle();
  h2c_2->GetYaxis()->SetTitle("Difference from true x [%]");h2c_2->GetYaxis()->CenterTitle();
  TH2 *h2c_3 = new TH2D("h2c_3","x Resolution vs. x",100,x_bins,200,-20,20);
  h2c_3->GetXaxis()->SetTitle("True x");h2c_3->GetXaxis()->CenterTitle();
  h2c_3->GetYaxis()->SetTitle("Difference from true x [%]");h2c_3->GetYaxis()->CenterTitle();

  TH2 *h2d_1 = new TH2D("h2d_1","Reconstructed vs. True Q^{2}",200,0,200,200,0,200);
  h2d_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h2d_1->GetXaxis()->CenterTitle();
  h2d_1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");h2d_1->GetYaxis()->CenterTitle();
  TH2 *h2d_2 = new TH2D("h2d_2","Reconstructed vs. True y",200,0,1,200,0,1);
  h2d_2->GetXaxis()->SetTitle("True y");h2d_2->GetXaxis()->CenterTitle();
  h2d_2->GetYaxis()->SetTitle("Reconstructed y");h2d_2->GetYaxis()->CenterTitle();
  TH2 *h2d_3 = new TH2D("h2d_3","Reconstructed vs. True x",100,x_bins,100,x_bins);
  h2d_3->GetXaxis()->SetTitle("True x");h2d_3->GetXaxis()->CenterTitle();
  h2d_3->GetYaxis()->SetTitle("Reconstructed x");h2d_3->GetYaxis()->CenterTitle();

  //JB Method -- using jet
  TH2 *h3a_1 = new TH2D("h3a_1","Q^{2} Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h3a_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h3a_1->GetXaxis()->CenterTitle();
  h3a_1->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h3a_1->GetYaxis()->CenterTitle();
  TH2 *h3a_2 = new TH2D("h3a_2","Q^{2} Resolution vs. y",200,0,1,200,-120,120);
  h3a_2->GetXaxis()->SetTitle("True y");h3a_2->GetXaxis()->CenterTitle();
  h3a_2->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h3a_2->GetYaxis()->CenterTitle();
  TH2 *h3a_3 = new TH2D("h3a_3","Q^{2} Resolution vs. x",100,x_bins,200,-120,120);
  h3a_3->GetXaxis()->SetTitle("True x");h3a_3->GetXaxis()->CenterTitle();
  h3a_3->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h3a_3->GetYaxis()->CenterTitle();

  TH2 *h3b_1 = new TH2D("h3b_1","y Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h3b_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h3b_1->GetXaxis()->CenterTitle();
  h3b_1->GetYaxis()->SetTitle("Difference from true y [%]");h3b_1->GetYaxis()->CenterTitle();
  TH2 *h3b_2 = new TH2D("h3b_2","y Resolution vs. y",200,0,1,200,-120,120);
  h3b_2->GetXaxis()->SetTitle("True y");h3b_2->GetXaxis()->CenterTitle();
  h3b_2->GetYaxis()->SetTitle("Difference from true y [%]");h3b_2->GetYaxis()->CenterTitle();
  TH2 *h3b_3 = new TH2D("h3b_3","y Resolution vs. x",100,x_bins,200,-120,120);
  h3b_3->GetXaxis()->SetTitle("True x");h3b_3->GetXaxis()->CenterTitle();
  h3b_3->GetYaxis()->SetTitle("Difference from true y [%]");h3b_3->GetYaxis()->CenterTitle();

  TH2 *h3c_1 = new TH2D("h3c_1","x Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h3c_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h3c_1->GetXaxis()->CenterTitle();
  h3c_1->GetYaxis()->SetTitle("Difference from true x [%]");h3c_1->GetYaxis()->CenterTitle();
  TH2 *h3c_2 = new TH2D("h3c_2","x Resolution vs. y",200,0,1,200,-120,120);
  h3c_2->GetXaxis()->SetTitle("True y");h3c_2->GetXaxis()->CenterTitle();
  h3c_2->GetYaxis()->SetTitle("Difference from true x [%]");h3c_2->GetYaxis()->CenterTitle();
  TH2 *h3c_3 = new TH2D("h3c_3","x Resolution vs. x",100,x_bins,200,-120,120);
  h3c_3->GetXaxis()->SetTitle("True x");h3c_3->GetXaxis()->CenterTitle();
  h3c_3->GetYaxis()->SetTitle("Difference from true x [%]");h3c_3->GetYaxis()->CenterTitle();

  TH2 *h3d_1 = new TH2D("h3d_1","Reconstructed vs. True Q^{2}",200,0,200,200,0,200);
  h3d_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h3d_1->GetXaxis()->CenterTitle();
  h3d_1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");h3d_1->GetYaxis()->CenterTitle();
  TH2 *h3d_2 = new TH2D("h3d_2","Reconstructed vs. True y",200,0,1,200,0,1);
  h3d_2->GetXaxis()->SetTitle("True y");h3d_2->GetXaxis()->CenterTitle();
  h3d_2->GetYaxis()->SetTitle("Reconstructed y");h3d_2->GetYaxis()->CenterTitle();
  TH2 *h3d_3 = new TH2D("h3d_3","Reconstructed vs. True x",100,x_bins,100,x_bins);
  h3d_3->GetXaxis()->SetTitle("True x");h3d_3->GetXaxis()->CenterTitle();
  h3d_3->GetYaxis()->SetTitle("Reconstructed x");h3d_3->GetYaxis()->CenterTitle();
  
  //JB Method -- Summing over all particles
  TH2 *h4a_1 = new TH2D("h4a_1","Q^{2} Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h4a_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h4a_1->GetXaxis()->CenterTitle();
  h4a_1->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h4a_1->GetYaxis()->CenterTitle();
  TH2 *h4a_2 = new TH2D("h4a_2","Q^{2} Resolution vs. y",200,0,1,200,-120,120);
  h4a_2->GetXaxis()->SetTitle("True y");h4a_2->GetXaxis()->CenterTitle();
  h4a_2->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h4a_2->GetYaxis()->CenterTitle();
  TH2 *h4a_3 = new TH2D("h4a_3","Q^{2} Resolution vs. x",100,x_bins,200,-120,120);
  h4a_3->GetXaxis()->SetTitle("True x");h4a_3->GetXaxis()->CenterTitle();
  h4a_3->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h4a_3->GetYaxis()->CenterTitle();

  TH2 *h4b_1 = new TH2D("h4b_1","y Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h4b_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h4b_1->GetXaxis()->CenterTitle();
  h4b_1->GetYaxis()->SetTitle("Difference from true y [%]");h4b_1->GetYaxis()->CenterTitle();
  TH2 *h4b_2 = new TH2D("h4b_2","y Resolution vs. y",200,0,1,200,-120,120);
  h4b_2->GetXaxis()->SetTitle("True y");h4b_2->GetXaxis()->CenterTitle();
  h4b_2->GetYaxis()->SetTitle("Difference from true y [%]");h4b_2->GetYaxis()->CenterTitle();
  TH2 *h4b_3 = new TH2D("h4b_3","y Resolution vs. x",100,x_bins,200,-120,120);
  h4b_3->GetXaxis()->SetTitle("True x");h4b_3->GetXaxis()->CenterTitle();
  h4b_3->GetYaxis()->SetTitle("Difference from true y [%]");h4b_3->GetYaxis()->CenterTitle();

  TH2 *h4c_1 = new TH2D("h4c_1","x Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h4c_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h4c_1->GetXaxis()->CenterTitle();
  h4c_1->GetYaxis()->SetTitle("Difference from true x [%]");h4c_1->GetYaxis()->CenterTitle();
  TH2 *h4c_2 = new TH2D("h4c_2","x Resolution vs. y",200,0,1,200,-120,120);
  h4c_2->GetXaxis()->SetTitle("True y");h4c_2->GetXaxis()->CenterTitle();
  h4c_2->GetYaxis()->SetTitle("Difference from true x [%]");h4c_2->GetYaxis()->CenterTitle();
  TH2 *h4c_3 = new TH2D("h4c_3","x Resolution vs. x",100,x_bins,200,-120,120);
  h4c_3->GetXaxis()->SetTitle("True x");h4c_3->GetXaxis()->CenterTitle();
  h4c_3->GetYaxis()->SetTitle("Difference from true x [%]");h4c_3->GetYaxis()->CenterTitle();

  TH2 *h4d_1 = new TH2D("h4d_1","Reconstructed vs. True Q^{2}",200,0,200,200,0,200);
  h4d_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h4d_1->GetXaxis()->CenterTitle();
  h4d_1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");h4d_1->GetYaxis()->CenterTitle();
  TH2 *h4d_2 = new TH2D("h4d_2","Reconstructed vs. True y",200,0,1,200,0,1);
  h4d_2->GetXaxis()->SetTitle("True y");h4d_2->GetXaxis()->CenterTitle();
  h4d_2->GetYaxis()->SetTitle("Reconstructed y");h4d_2->GetYaxis()->CenterTitle();
  TH2 *h4d_3 = new TH2D("h4d_3","Reconstructed vs. True x",100,x_bins,100,x_bins);
  h4d_3->GetXaxis()->SetTitle("True x");h4d_3->GetXaxis()->CenterTitle();
  h4d_3->GetYaxis()->SetTitle("Reconstructed x");h4d_3->GetYaxis()->CenterTitle();

  //JB 4-Vector Method
  TH2 *h5a_1 = new TH2D("h5a_1","Q^{2} Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h5a_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h5a_1->GetXaxis()->CenterTitle();
  h5a_1->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h5a_1->GetYaxis()->CenterTitle();
  TH2 *h5a_2 = new TH2D("h5a_2","Q^{2} Resolution vs. y",200,0,1,200,-120,120);
  h5a_2->GetXaxis()->SetTitle("True y");h5a_2->GetXaxis()->CenterTitle();
  h5a_2->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h5a_2->GetYaxis()->CenterTitle();
  TH2 *h5a_3 = new TH2D("h5a_3","Q^{2} Resolution vs. x",100,x_bins,200,-120,120);
  h5a_3->GetXaxis()->SetTitle("True x");h5a_3->GetXaxis()->CenterTitle();
  h5a_3->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h5a_3->GetYaxis()->CenterTitle();

  TH2 *h5b_1 = new TH2D("h5b_1","y Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h5b_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h5b_1->GetXaxis()->CenterTitle();
  h5b_1->GetYaxis()->SetTitle("Difference from true y [%]");h5b_1->GetYaxis()->CenterTitle();
  TH2 *h5b_2 = new TH2D("h5b_2","y Resolution vs. y",200,0,1,200,-120,120);
  h5b_2->GetXaxis()->SetTitle("True y");h5b_2->GetXaxis()->CenterTitle();
  h5b_2->GetYaxis()->SetTitle("Difference from true y [%]");h5b_2->GetYaxis()->CenterTitle();
  TH2 *h5b_3 = new TH2D("h5b_3","y Resolution vs. x",100,x_bins,200,-120,120);
  h5b_3->GetXaxis()->SetTitle("True x");h5b_3->GetXaxis()->CenterTitle();
  h5b_3->GetYaxis()->SetTitle("Difference from true y [%]");h5b_3->GetYaxis()->CenterTitle();

  TH2 *h5c_1 = new TH2D("h5c_1","x Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h5c_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h5c_1->GetXaxis()->CenterTitle();
  h5c_1->GetYaxis()->SetTitle("Difference from true x [%]");h5c_1->GetYaxis()->CenterTitle();
  TH2 *h5c_2 = new TH2D("h5c_2","x Resolution vs. y",200,0,1,200,-120,120);
  h5c_2->GetXaxis()->SetTitle("True y");h5c_2->GetXaxis()->CenterTitle();
  h5c_2->GetYaxis()->SetTitle("Difference from true x [%]");h5c_2->GetYaxis()->CenterTitle();
  TH2 *h5c_3 = new TH2D("h5c_3","x Resolution vs. x",100,x_bins,200,-120,120);
  h5c_3->GetXaxis()->SetTitle("True x");h5c_3->GetXaxis()->CenterTitle();
  h5c_3->GetYaxis()->SetTitle("Difference from true x [%]");h5c_3->GetYaxis()->CenterTitle();

  TH2 *h5d_1 = new TH2D("h5d_1","Reconstructed vs. True Q^{2}",200,0,200,200,0,200);
  h5d_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h5d_1->GetXaxis()->CenterTitle();
  h5d_1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");h5d_1->GetYaxis()->CenterTitle();
  TH2 *h5d_2 = new TH2D("h5d_2","Reconstructed vs. True y",200,0,1,200,0,1);
  h5d_2->GetXaxis()->SetTitle("True y");h5d_2->GetXaxis()->CenterTitle();
  h5d_2->GetYaxis()->SetTitle("Reconstructed y");h5d_2->GetYaxis()->CenterTitle();
  TH2 *h5d_3 = new TH2D("h5d_3","Reconstructed vs. True x",100,x_bins,100,x_bins);
  h5d_3->GetXaxis()->SetTitle("True x");h5d_3->GetXaxis()->CenterTitle();
  h5d_3->GetYaxis()->SetTitle("Reconstructed x");h5d_3->GetYaxis()->CenterTitle();
  
  //DA Method -- using jet
  TH2 *h6a_1 = new TH2D("h6a_1","Q^{2} Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h6a_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h6a_1->GetXaxis()->CenterTitle();
  h6a_1->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h6a_1->GetYaxis()->CenterTitle();
  TH2 *h6a_2 = new TH2D("h6a_2","Q^{2} Resolution vs. y",200,0,1,200,-120,120);
  h6a_2->GetXaxis()->SetTitle("True y");h6a_2->GetXaxis()->CenterTitle();
  h6a_2->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h6a_2->GetYaxis()->CenterTitle();
  TH2 *h6a_3 = new TH2D("h6a_3","Q^{2} Resolution vs. x",100,x_bins,200,-120,120);
  h6a_3->GetXaxis()->SetTitle("True x");h6a_3->GetXaxis()->CenterTitle();
  h6a_3->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h6a_3->GetYaxis()->CenterTitle();

  TH2 *h6b_1 = new TH2D("h6b_1","y Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h6b_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h6b_1->GetXaxis()->CenterTitle();
  h6b_1->GetYaxis()->SetTitle("Difference from true y [%]");h6b_1->GetYaxis()->CenterTitle();
  TH2 *h6b_2 = new TH2D("h6b_2","y Resolution vs. y",200,0,1,200,-120,120);
  h6b_2->GetXaxis()->SetTitle("True y");h6b_2->GetXaxis()->CenterTitle();
  h6b_2->GetYaxis()->SetTitle("Difference from true y [%]");h6b_2->GetYaxis()->CenterTitle();
  TH2 *h6b_3 = new TH2D("h6b_3","y Resolution vs. x",100,x_bins,200,-120,120);
  h6b_3->GetXaxis()->SetTitle("True x");h6b_3->GetXaxis()->CenterTitle();
  h6b_3->GetYaxis()->SetTitle("Difference from true y [%]");h6b_3->GetYaxis()->CenterTitle();

  TH2 *h6c_1 = new TH2D("h6c_1","x Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h6c_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h6c_1->GetXaxis()->CenterTitle();
  h6c_1->GetYaxis()->SetTitle("Difference from true x [%]");h6c_1->GetYaxis()->CenterTitle();
  TH2 *h6c_2 = new TH2D("h6c_2","x Resolution vs. y",200,0,1,200,-120,120);
  h6c_2->GetXaxis()->SetTitle("True y");h6c_2->GetXaxis()->CenterTitle();
  h6c_2->GetYaxis()->SetTitle("Difference from true x [%]");h6c_2->GetYaxis()->CenterTitle();
  TH2 *h6c_3 = new TH2D("h6c_3","x Resolution vs. x",100,x_bins,200,-120,120);
  h6c_3->GetXaxis()->SetTitle("True x");h6c_3->GetXaxis()->CenterTitle();
  h6c_3->GetYaxis()->SetTitle("Difference from true x [%]");h6c_3->GetYaxis()->CenterTitle();

  TH2 *h6d_1 = new TH2D("h6d_1","Reconstructed vs. True Q^{2}",200,0,200,200,0,200);
  h6d_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h6d_1->GetXaxis()->CenterTitle();
  h6d_1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");h6d_1->GetYaxis()->CenterTitle();
  TH2 *h6d_2 = new TH2D("h6d_2","Reconstructed vs. True y",200,0,1,200,0,1);
  h6d_2->GetXaxis()->SetTitle("True y");h6d_2->GetXaxis()->CenterTitle();
  h6d_2->GetYaxis()->SetTitle("Reconstructed y");h6d_2->GetYaxis()->CenterTitle();
  TH2 *h6d_3 = new TH2D("h6d_3","Reconstructed vs. True x",100,x_bins,100,x_bins);
  h6d_3->GetXaxis()->SetTitle("True x");h6d_3->GetXaxis()->CenterTitle();
  h6d_3->GetYaxis()->SetTitle("Reconstructed x");h6d_3->GetYaxis()->CenterTitle();

  //DA Method -- Summing over all particles
  TH2 *h7a_1 = new TH2D("h7a_1","Q^{2} Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h7a_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h7a_1->GetXaxis()->CenterTitle();
  h7a_1->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h7a_1->GetYaxis()->CenterTitle();
  TH2 *h7a_2 = new TH2D("h7a_2","Q^{2} Resolution vs. y",200,0,1,200,-120,120);
  h7a_2->GetXaxis()->SetTitle("True y");h7a_2->GetXaxis()->CenterTitle();
  h7a_2->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h7a_2->GetYaxis()->CenterTitle();
  TH2 *h7a_3 = new TH2D("h7a_3","Q^{2} Resolution vs. x",100,x_bins,200,-120,120);
  h7a_3->GetXaxis()->SetTitle("True x");h7a_3->GetXaxis()->CenterTitle();
  h7a_3->GetYaxis()->SetTitle("Difference from true Q^{2} [%]");h7a_3->GetYaxis()->CenterTitle();

  TH2 *h7b_1 = new TH2D("h7b_1","y Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h7b_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h7b_1->GetXaxis()->CenterTitle();
  h7b_1->GetYaxis()->SetTitle("Difference from true y [%]");h7b_1->GetYaxis()->CenterTitle();
  TH2 *h7b_2 = new TH2D("h7b_2","y Resolution vs. y",200,0,1,200,-120,120);
  h7b_2->GetXaxis()->SetTitle("True y");h7b_2->GetXaxis()->CenterTitle();
  h7b_2->GetYaxis()->SetTitle("Difference from true y [%]");h7b_2->GetYaxis()->CenterTitle();
  TH2 *h7b_3 = new TH2D("h7b_3","y Resolution vs. x",100,x_bins,200,-120,120);
  h7b_3->GetXaxis()->SetTitle("True x");h7b_3->GetXaxis()->CenterTitle();
  h7b_3->GetYaxis()->SetTitle("Difference from true y [%]");h7b_3->GetYaxis()->CenterTitle();

  TH2 *h7c_1 = new TH2D("h7c_1","x Resolution vs. Q^{2}",200,0,200,200,-120,120);
  h7c_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h7c_1->GetXaxis()->CenterTitle();
  h7c_1->GetYaxis()->SetTitle("Difference from true x [%]");h7c_1->GetYaxis()->CenterTitle();
  TH2 *h7c_2 = new TH2D("h7c_2","x Resolution vs. y",200,0,1,200,-120,120);
  h7c_2->GetXaxis()->SetTitle("True y");h7c_2->GetXaxis()->CenterTitle();
  h7c_2->GetYaxis()->SetTitle("Difference from true x [%]");h7c_2->GetYaxis()->CenterTitle();
  TH2 *h7c_3 = new TH2D("h7c_3","x Resolution vs. x",100,x_bins,200,-120,120);
  h7c_3->GetXaxis()->SetTitle("True x");h7c_3->GetXaxis()->CenterTitle();
  h7c_3->GetYaxis()->SetTitle("Difference from true x [%]");h7c_3->GetYaxis()->CenterTitle();

  TH2 *h7d_1 = new TH2D("h7d_1","Reconstructed vs. True Q^{2}",200,0,200,200,0,200);
  h7d_1->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");h7d_1->GetXaxis()->CenterTitle();
  h7d_1->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");h7d_1->GetYaxis()->CenterTitle();
  TH2 *h7d_2 = new TH2D("h7d_2","Reconstructed vs. True y",200,0,1,200,0,1);
  h7d_2->GetXaxis()->SetTitle("True y");h7d_2->GetXaxis()->CenterTitle();
  h7d_2->GetYaxis()->SetTitle("Reconstructed y");h7d_2->GetYaxis()->CenterTitle();
  TH2 *h7d_3 = new TH2D("h7d_3","Reconstructed vs. True x",100,x_bins,100,x_bins);
  h7d_3->GetXaxis()->SetTitle("True x");h7d_3->GetXaxis()->CenterTitle();
  h7d_3->GetYaxis()->SetTitle("Reconstructed x");h7d_3->GetYaxis()->CenterTitle();

  //Line for reconstructed vs. true plots
  TF1 *func_1 = new TF1("func_1","x",0,200);
  func_1->SetLineColor(kRed);func_1->SetLineWidth(2);func_1->SetLineStyle(2);

  //--------------------------------//
  //   Analyse PYTHIA Simulation    //
  //--------------------------------//

  //Load ROOT Files
  erhic::EventPythia *event(NULL);    //Event Class
  erhic::ParticleMC *particle(NULL);  //Particle Class

  Smear::Event *event_s(NULL);           //Smeared Event Class
  Smear::ParticleMCS *particle_s(NULL);  //Smeared Particle Class

  TChain *tree = new TChain("EICTree");
  TChain *tree_s = new TChain("Smeared");
  
  if(energy_set == 1){
    for(int i=0;i<15;i++){
        tree->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/yellow_report/5_41/ep_5_41_newtune_%d.root",i));
        tree_s->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/yellow_report/5_41/ep_5_41_newtune_%d_matrixsmear.root",i));
    }
  }

  if(energy_set == 2){
    for(int i=0;i<15;i++){
        tree->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/yellow_report/18_275/ep_18_275_newtune_%d.root",i));
        tree_s->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/yellow_report/18_275/ep_18_275_newtune_%d_matrixsmear.root",i));
    }
  }

  tree->AddFriend(tree_s,"Smeared"); //Smeared Tree
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("eventS",&event_s);

  //MC Variables
  Int_t nParticles(0);
  Int_t Status[500];
  Int_t id[500];
  Int_t orig[500];
  Double_t Ef_e(0), pxf_e(0), pyf_e(0), pzf_e(0); //Scattered Electron
  TLorentzVector ef;
  Double_t Ei_e(0),pxi_e(0),pyi_e(0),pzi_e(0); //Initial Electron
  TLorentzVector ei;
  Double_t E_pn(0),px_pn(0),py_pn(0),pz_pn(0),Mpn(0); //Initial Proton/Neutron
  TLorentzVector pni;
  TLorentzVector q_e; //q 4-vector using scattered electron
  Double_t Q2_e(0),y_e(0),x_e(0),W_e(0); //Scattered electron (no assumptions)
  Double_t s_nm(0); //CM Energy Squared (massless e,p)
  Int_t electronIndex(0);

  //Smeared Variables
  Int_t nParticles_s(0);
  Int_t Status_s[500];
  Int_t id_s[500];
  Double_t Ef_e_s(0), pxf_e_s(0), pyf_e_s(0), pzf_e_s(0),theta_e_s(0); //Scattered Electron
  Double_t Ef_e_p_s(0); //Scattered electron, using total momentum as final energy
  Double_t pztot_jet_s(0),pttot_jet_s(0),Etot_jet_s(0); // Total added hadronic variables via. FastJet
  Double_t pxtot_sumh_s(0),pytot_sumh_s(0),pztot_sumh_s(0),pttot_sumh_s(0),Etot_sumh_s(0); // Total added hadronic variables by summing over hadrons
  Double_t holdpx_h(0),holdpy_h(0),holdpz_h(0),holdE_h(0); //Hold variables for summing over hadrons
  TLorentzVector ef_s;
  TLorentzVector q_e_s; //q 4-vector using scattered electron
  Double_t Q2_e_nm_s(0),y_e_nm_s(0),x_e_nm_s(0); //Scattered electron (massless e,e',p)
  Double_t Q2_e_nm_p_s(0),y_e_nm_p_s(0),x_e_nm_p_s(0); //Scattered electron using total momentum as final energy(massless e,e',p)
  Double_t Theta_jb_jet_s(0),Q2_jb_jet_s(0),y_jb_jet_s(0),x_jb_jet_s(0); // JB Method (Using Jets)
  Double_t Theta_jb_sumh_s(0),Q2_jb_sumh_s(0),y_jb_sumh_s(0),x_jb_sumh_s(0); // JB Method (Summing Over Hadrons)
  Double_t Theta_h_nm_s(0), Q2_da_s(0), y_da_s(0), x_da_s(0); //DA Method
  Double_t mass[500];
  bool detected_elec(false);
  TLorentzVector eh_s_tot;
  TLorentzVector q_h_s;

  //FastJet Variables
  vector<fastjet::PseudoJet> holdParticles;
  valarray<double> fourvec(4);
  TLorentzVector virtualPhoton;
  TLorentzVector jet;
  Double_t R = 1.0;     //For now, going to use the paramters and cuts from help-script
  fastjet::JetDefinition jetDef(fastjet::kt_algorithm,R,fastjet::E_scheme,fastjet::Best); 
  Double_t etaMax = 4.0;
  Double_t etaMin = -etaMax;
  
  //Not going to apply pt cuts to jet. Don't see a reason to yet
  //auto selectJetEta = fastjet::SelectorEtaRange(etaMin+R,etaMax+R);

  Int_t nevents = 1e6;//tree->GetEntries();
  
  // Loop over all events
  for(int i=0;i<nevents;i++){
    if(i%10000==0) cout<<"Events Analysed = "<<i<<"!"<<endl;

    //Reset FastJet Variables
    holdParticles.clear();

    //Reset Summing Over Hadrons Variables
    Etot_sumh_s = 0;
    pxtot_sumh_s = 0;
    pytot_sumh_s = 0;
    pztot_sumh_s = 0;
    pttot_sumh_s = 0;

    tree->GetEntry(i);
    
    //MC Tree First
    nParticles = event->GetNTracks();

    //First loop to find scattered electron and beam particles
    for(int j=0;j<nParticles;j++){
      
      particle = event->GetTrack(j);
      Status[j] = (Int_t) particle->GetStatus();
      id[j] = (Int_t) particle->Id();
      mass[j] = (Double_t) particle->GetM();
      orig[j] = (Int_t) particle->GetParentIndex();

      if(Status[j]==1 && id[j]==11 && orig[j]==3){
	        Ef_e = particle->GetE();
	        pxf_e = particle->GetPx();
	        pyf_e = particle->GetPy();
	        pzf_e = particle->GetPz();
	        ef.SetPxPyPzE(pxf_e,pyf_e,pzf_e,Ef_e);
	        electronIndex = j;
      }

      // Get initial state election
      if(j==0 && (id[j]==11 && Status[j]==21)){
	      Ei_e = particle->GetE();
	      pxi_e = particle->GetPx();
	      pyi_e = particle->GetPy();
	      pzi_e = particle->GetPz();
	      ei.SetPxPyPzE(pxi_e,pyi_e,pzi_e,Ei_e);
      }
	 
      // Get initial state proton
      if(j==1 && (id[j]==2212 && Status[j]==21)){
	      E_pn = particle->GetE();
	      px_pn = particle->GetPx();
	      py_pn = particle->GetPy();
	      pz_pn = particle->GetPz();
	      Mpn = particle->GetM();
	      pni.SetPxPyPzE(px_pn,py_pn,pz_pn,E_pn);
      }
    }
    
    //Second loop over all particles to get smeared information
    for(int j=0;j<nParticles;j++){
      
      particle = event->GetTrack(j);
      Status[j] = (Int_t) particle->GetStatus();
      id[j] = (Int_t) particle->Id();
      mass[j] = (Double_t) particle->GetM();
      orig[j] = (Int_t) particle->GetParentIndex();

      particle_s = event_s->GetTrack(j);
      
      if(particle_s){ //make sure not null pointer
        Status_s[j] = (Int_t) particle_s->GetStatus();
	      id_s[j] = (Int_t) particle_s->Id();
 
	      //Get Smeared-Scattered Electron
	      if(Status_s[j]==1 && Status[j]==1 && id[j]==11 && orig[j]==3){
	        Ef_e_s = particle_s->GetE();
	        pxf_e_s = particle_s->GetPx();
	        pyf_e_s = particle_s->GetPy();
	        pzf_e_s = particle_s->GetPz();
          theta_e_s = particle_s->GetTheta();
	        ef_s.SetPxPyPzE(pxf_e_s,pyf_e_s,pzf_e_s,Ef_e_s);
          detected_elec = true;
	      }

	      //Put hadrons into PseudoJet object
	      if(j!=electronIndex && Status_s[j]==1){
          
          //Summing over all particles
          holdE_h = particle_s->GetE();

          //handle neutral particles
          if(id[j]==22 || id[j]==130 || id[j]==2112){

            double neutral_mom = 0;

            switch(id[j]){
              case 22:
                neutral_mom = holdE_h;
              case 130:
                neutral_mom = sqrt(holdE_h*holdE_h - 0.4976*0.4976);
              case 2112:
                neutral_mom = sqrt(holdE_h*holdE_h - 0.9396*0.9396);
            }

            holdpx_h = neutral_mom * sin( particle_s->GetTheta() ) * cos( particle_s->GetPhi() );
            holdpy_h = neutral_mom * sin( particle_s->GetTheta() ) * sin( particle_s->GetPhi() );
            holdpz_h = neutral_mom * cos( particle_s->GetTheta() );

          }else{
            holdpx_h = particle_s->GetPx();
            holdpy_h = particle_s->GetPy();
            holdpz_h = particle_s->GetPz();
          }

          pxtot_sumh_s = pxtot_sumh_s + holdpx_h;
          pytot_sumh_s = pytot_sumh_s + holdpy_h;
          pztot_sumh_s = pztot_sumh_s + holdpz_h;
          Etot_sumh_s =  Etot_sumh_s + holdE_h;

          //Fast Jet
	        fourvec[0] = particle_s->GetPx();
	        fourvec[1] = particle_s->GetPy();
	        fourvec[2] = particle_s->GetPz();
	        fourvec[3] = particle_s->GetE();

	        fastjet::PseudoJet finalStateParticle(fourvec);
	        finalStateParticle.set_user_index(j);
	        holdParticles.push_back(finalStateParticle);
	      }
      }  
    }//Finished Loop Over Particles

    //-------------Calculate invariants using 4-vectors------------------//
    //Since we are comparing to smeared tree here, only calculate
    //variables using scattered electron method.
    //We could have instead used the functions:
    //----
    //Q2_event = event->GetQ2();
    //x_event = event->GetX();
    //y_event = event->GetY();
    //W2_event = event->GetW2();
    //----
    //These calculate the kinematic quantities using the scattered electron ignoring
    //the electron mass, which is fine for the Q2 range being considered.
    //If we need to consider electron mass or QED radiative effects are ON, may want to use:
    //----
    //Q2_true = event->GetTrueQ2();
    //x_true = event->GetTrueX();
    //y_true = event->GetTrueY();
    //W2_true = event->GetTrueW2();
    //----
    q_e = (ei - ef);
    Q2_e = -1.0*(q_e*q_e);
    y_e = (pni*q_e)/(pni*ei);
    x_e = Q2_e/(2*pni*q_e); 
    s_nm = 4*Ei_e*E_pn;

    //-------------Calculate *Smeared* invariants using 4-vectors------------------//
  
    if(detected_elec && y_e>y_min && y_e<y_max){
      //1.1) Smeared scattered election using Ecal energy the final energy
      Q2_e_nm_s = 4.*Ei_e*Ef_e_s*TMath::Cos(theta_e_s/2.)*TMath::Cos(theta_e_s/2.);
      y_e_nm_s = 1. - ( (Ef_e_s/(2.*Ei_e))*(1. - TMath::Cos(theta_e_s)) );
      x_e_nm_s = Q2_e_nm_s/(s_nm*y_e_nm_s);

      h1a_1->Fill(Q2_e, 100.*(Q2_e-Q2_e_nm_s)/Q2_e );
      h1a_2->Fill(y_e, 100.*(Q2_e-Q2_e_nm_s)/Q2_e );
      h1a_3->Fill(x_e, 100.*(Q2_e-Q2_e_nm_s)/Q2_e );

      h1b_1->Fill(Q2_e, 100.*(y_e-y_e_nm_s)/y_e );
      h1b_2->Fill(y_e, 100.*(y_e-y_e_nm_s)/y_e );
      h1b_3->Fill(x_e, 100.*(y_e-y_e_nm_s)/y_e );

      h1c_1->Fill(Q2_e, 100.*(x_e-x_e_nm_s)/x_e );
      h1c_2->Fill(y_e, 100.*(x_e-x_e_nm_s)/x_e );
      h1c_3->Fill(x_e, 100.*(x_e-x_e_nm_s)/x_e );

      h1d_1->Fill(Q2_e,Q2_e_nm_s);
      h1d_2->Fill(y_e,y_e_nm_s);
      h1d_3->Fill(x_e,x_e_nm_s);

      //1.2) Using Smeared scattered election using track momentum as the total momentum as final energy.  
      Ef_e_p_s = TMath::Sqrt( TMath::Power(pxf_e_s,2)+TMath::Power(pyf_e_s,2)+TMath::Power(pzf_e_s,2) + (0.511e-3)*(0.511e-3) ); 
      Q2_e_nm_p_s = 4.*Ei_e*Ef_e_p_s*TMath::Cos(theta_e_s/2.)*TMath::Cos(theta_e_s/2.);
      y_e_nm_p_s = 1. - ( (Ef_e_p_s/(2.*Ei_e))*(1. - TMath::Cos(theta_e_s)) );
      x_e_nm_p_s = Q2_e_nm_p_s/(s_nm*y_e_nm_p_s);

      h2a_1->Fill(Q2_e, 100.*(Q2_e-Q2_e_nm_p_s)/Q2_e );
      h2a_2->Fill(y_e, 100.*(Q2_e-Q2_e_nm_p_s)/Q2_e );
      h2a_3->Fill(x_e, 100.*(Q2_e-Q2_e_nm_p_s)/Q2_e );

      h2b_1->Fill(Q2_e, 100.*(y_e-y_e_nm_p_s)/y_e );
      h2b_2->Fill(y_e, 100.*(y_e-y_e_nm_p_s)/y_e );
      h2b_3->Fill(x_e, 100.*(y_e-y_e_nm_p_s)/y_e );

      h2c_1->Fill(Q2_e, 100.*(x_e-x_e_nm_p_s)/x_e );
      h2c_2->Fill(y_e, 100.*(x_e-x_e_nm_p_s)/x_e );
      h2c_3->Fill(x_e, 100.*(x_e-x_e_nm_p_s)/x_e );

      h2d_1->Fill(Q2_e,Q2_e_nm_p_s);
      h2d_2->Fill(y_e,y_e_nm_p_s);
      h2d_3->Fill(x_e,x_e_nm_p_s);

    } //detected electron

    if(!(holdParticles.empty()) && y_e>y_min && y_e<y_max){
      pttot_sumh_s = TMath::Sqrt( TMath::Power(pxtot_sumh_s,2)+TMath::Power(pytot_sumh_s,2) );

      fastjet::ClusterSequence cs(holdParticles, jetDef);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

      if(!(jets.empty())){
        //For now, I am going to calculate the DIS Kinematics via. JB Method by using the highest pt jet
        pztot_jet_s = jets[0].pz();
        pttot_jet_s = jets[0].pt();
        Etot_jet_s = jets[0].E();

        //2.1) Using Smeared JB Method                                                                                                           
        y_jb_jet_s = (Etot_jet_s - pztot_jet_s)/(2*Ei_e);
        Q2_jb_jet_s = (pttot_jet_s*pttot_jet_s)/(1. - y_jb_jet_s);
        x_jb_jet_s = Q2_jb_jet_s/(s_nm*y_jb_jet_s);

        h3a_1->Fill(Q2_e, 100.*(Q2_e-Q2_jb_jet_s)/Q2_e );
        h3a_2->Fill(y_e, 100.*(Q2_e-Q2_jb_jet_s)/Q2_e );
        h3a_3->Fill(x_e, 100.*(Q2_e-Q2_jb_jet_s)/Q2_e );

        h3b_1->Fill(Q2_e, 100.*(y_e-y_jb_jet_s)/y_e );
        h3b_2->Fill(y_e, 100.*(y_e-y_jb_jet_s)/y_e );
        h3b_3->Fill(x_e, 100.*(y_e-y_jb_jet_s)/y_e );

        h3c_1->Fill(Q2_e, 100.*(x_e-x_jb_jet_s)/x_e );
        h3c_2->Fill(y_e, 100.*(x_e-x_jb_jet_s)/x_e );
        h3c_3->Fill(x_e, 100.*(x_e-x_jb_jet_s)/x_e );

        h3d_1->Fill(Q2_e,Q2_jb_jet_s);
        h3d_2->Fill(y_e,y_jb_jet_s);
        h3d_3->Fill(x_e,x_jb_jet_s);
        
        //3.1) Using Smeared DA Method
        if(detected_elec){
          Theta_h_nm_s = 2.* TMath::ATan( (Etot_jet_s - pztot_jet_s)/pttot_jet_s);    
          Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(theta_e_s/2.) )*
                   ( 1./(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
          y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
          x_da_s = Q2_da_s/(s_nm*y_da_s);  

          h6a_1->Fill(Q2_e, 100.*(Q2_e-Q2_da_s)/Q2_e );
          h6a_2->Fill(y_e, 100.*(Q2_e-Q2_da_s)/Q2_e );
          h6a_3->Fill(x_e, 100.*(Q2_e-Q2_da_s)/Q2_e );

          h6b_1->Fill(Q2_e, 100.*(y_e-y_da_s)/y_e );
          h6b_2->Fill(y_e, 100.*(y_e-y_da_s)/y_e );
          h6b_3->Fill(x_e, 100.*(y_e-y_da_s)/y_e );

          h6c_1->Fill(Q2_e, 100.*(x_e-x_da_s)/x_e );
          h6c_2->Fill(y_e, 100.*(x_e-x_da_s)/x_e );
          h6c_3->Fill(x_e, 100.*(x_e-x_da_s)/x_e );

          h6d_1->Fill(Q2_e,Q2_da_s);
          h6d_2->Fill(y_e,y_da_s);
          h6d_3->Fill(x_e,x_da_s);
        }

      }//jets.empty()

      //2.2) Using Smeared JB Method (Summing Over all particles)
      y_jb_sumh_s = (Etot_sumh_s - pztot_sumh_s)/(2*Ei_e);
      Q2_jb_sumh_s = (pttot_sumh_s*pttot_sumh_s)/(1. - y_jb_sumh_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s);

      h4a_1->Fill(Q2_e, 100.*(Q2_e-Q2_jb_sumh_s)/Q2_e );
      h4a_2->Fill(y_e, 100.*(Q2_e-Q2_jb_sumh_s)/Q2_e );
      h4a_3->Fill(x_e, 100.*(Q2_e-Q2_jb_sumh_s)/Q2_e );

      h4b_1->Fill(Q2_e, 100.*(y_e-y_jb_sumh_s)/y_e );
      h4b_2->Fill(y_e, 100.*(y_e-y_jb_sumh_s)/y_e );
      h4b_3->Fill(x_e, 100.*(y_e-y_jb_sumh_s)/y_e );

      h4c_1->Fill(Q2_e, 100.*(x_e-x_jb_sumh_s)/x_e );
      h4c_2->Fill(y_e, 100.*(x_e-x_jb_sumh_s)/x_e );
      h4c_3->Fill(x_e, 100.*(x_e-x_jb_sumh_s)/x_e );

      h4d_1->Fill(Q2_e,Q2_jb_sumh_s);
      h4d_2->Fill(y_e,y_jb_sumh_s);
      h4d_3->Fill(x_e,x_jb_sumh_s);

      //2.3) Using 4-Vector of outgoing X
      eh_s_tot.SetPxPyPzE(pxtot_sumh_s,pytot_sumh_s,pztot_sumh_s,Etot_sumh_s);

      q_h_s= (eh_s_tot - pni);
      Q2_jb_sumh_s = -1.0*(q_h_s*q_h_s);
      y_jb_sumh_s = (pni*q_h_s)/(pni*ei);
      //x_jb_sumh_s = Q2_jb_sumh_s/(2*pni*q_h_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s); //Using Q2 = xsy, instead of 4-vector product Q2_jb_sumh_s/(2*pni*q_h_s); 

      h5a_1->Fill(Q2_e, 100.*(Q2_e-Q2_jb_sumh_s)/Q2_e );
      h5a_2->Fill(y_e, 100.*(Q2_e-Q2_jb_sumh_s)/Q2_e );
      h5a_3->Fill(x_e, 100.*(Q2_e-Q2_jb_sumh_s)/Q2_e );

      h5b_1->Fill(Q2_e, 100.*(y_e-y_jb_sumh_s)/y_e );
      h5b_2->Fill(y_e, 100.*(y_e-y_jb_sumh_s)/y_e );
      h5b_3->Fill(x_e, 100.*(y_e-y_jb_sumh_s)/y_e );

      h5c_1->Fill(Q2_e, 100.*(x_e-x_jb_sumh_s)/x_e );
      h5c_2->Fill(y_e, 100.*(x_e-x_jb_sumh_s)/x_e );
      h5c_3->Fill(x_e, 100.*(x_e-x_jb_sumh_s)/x_e );

      h5d_1->Fill(Q2_e,Q2_jb_sumh_s);
      h5d_2->Fill(y_e,y_jb_sumh_s);
      h5d_3->Fill(x_e,x_jb_sumh_s);

      //3.2) Using Smeared DA Method (Summing Over Hadrons)
      if(detected_elec){
        Theta_h_nm_s = 2.* TMath::ATan( (Etot_sumh_s - pztot_sumh_s)/pttot_sumh_s);
        Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(theta_e_s/2.) )*
          ( 1./(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
        y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
        x_da_s = Q2_da_s/(s_nm*y_da_s);

        h7a_1->Fill(Q2_e, 100.*(Q2_e-Q2_da_s)/Q2_e );
        h7a_2->Fill(y_e, 100.*(Q2_e-Q2_da_s)/Q2_e );
        h7a_3->Fill(x_e, 100.*(Q2_e-Q2_da_s)/Q2_e );

        h7b_1->Fill(Q2_e, 100.*(y_e-y_da_s)/y_e );
        h7b_2->Fill(y_e, 100.*(y_e-y_da_s)/y_e );
        h7b_3->Fill(x_e, 100.*(y_e-y_da_s)/y_e );

        h7c_1->Fill(Q2_e, 100.*(x_e-x_da_s)/x_e );
        h7c_2->Fill(y_e, 100.*(x_e-x_da_s)/x_e );
        h7c_3->Fill(x_e, 100.*(x_e-x_da_s)/x_e );

        h7d_1->Fill(Q2_e,Q2_da_s);
        h7d_2->Fill(y_e,y_da_s);
        h7d_3->Fill(x_e,x_da_s);

      }

    }//holdParticles.empty()

    //reset variables
    detected_elec = false;

  }//Finished Event Loop

  //Make Latex
  TPaveText* tex_energy = new TPaveText(0.1,0.7,0.9,0.9,"NDCNB");

  if(energy_set == 1){
    tex_energy->AddText("5 GeV e^{-} on 41 GeV p, #sqrt{s}=28.6 GeV");
	  tex_energy->SetFillStyle(4000);tex_energy->SetTextFont(63);tex_energy->SetTextSize(15);
  }

  if(energy_set == 2){
    tex_energy->AddText("18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
	  tex_energy->SetFillStyle(4000);tex_energy->SetTextFont(63);tex_energy->SetTextSize(15);
  }

  TPaveText *texy = new TPaveText(0.1,0.3,0.9,0.5,"NDCNB");
  texy->AddText(Form("%.2f < y_{true} < %.2f",y_min,y_max));
  texy->SetFillStyle(4000);texy->SetTextFont(63);texy->SetTextSize(10);

  TPaveText *tex1_1 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex1_1->AddText("#frac{True - Electron Method (using ECal Energy)}{True} vs. True");
  tex1_1->SetFillStyle(4000);tex1_1->SetTextFont(63);tex1_1->SetTextSize(10);

  TPaveText *tex1_2 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex1_2->AddText("Electron Method (using ECal Energy) vs. True");
  tex1_2->SetFillStyle(4000);tex1_2->SetTextFont(63);tex1_2->SetTextSize(10);

  TPaveText *tex2_1 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex2_1->AddText("#frac{True - Electron Method (using track momentum)}{True} vs. True");
  tex2_1->SetFillStyle(4000);tex2_1->SetTextFont(63);tex2_1->SetTextSize(10);

  TPaveText *tex2_2 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex2_2->AddText("Electron Method (using track momentum) vs. True");
  tex2_2->SetFillStyle(4000);tex2_2->SetTextFont(63);tex2_2->SetTextSize(10);

  TPaveText *tex3_1 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex3_1->AddText("#frac{True - J.B. Method (using jet)}{True} vs. True");
  tex3_1->SetFillStyle(4000);tex3_1->SetTextFont(63);tex3_1->SetTextSize(10);

  TPaveText *tex3_2 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex3_2->AddText("J.B. Method (using jet) vs. True");
  tex3_2->SetFillStyle(4000);tex3_2->SetTextFont(63);tex3_2->SetTextSize(10);

  TPaveText *tex4_1 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex4_1->AddText("#frac{True - J.B. Method (summing all particles)}{True} vs. True");
  tex4_1->SetFillStyle(4000);tex4_1->SetTextFont(63);tex4_1->SetTextSize(10);

  TPaveText *tex4_2 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex4_2->AddText("J.B. Method (summing all particles) vs. True");
  tex4_2->SetFillStyle(4000);tex4_2->SetTextFont(63);tex4_2->SetTextSize(10);

  TPaveText *tex5_1 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex5_1->AddText("#frac{True - J.B. 4-Vector Method}{True} vs. True");
  tex5_1->SetFillStyle(4000);tex5_1->SetTextFont(63);tex5_1->SetTextSize(10);

  TPaveText *tex5_2 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex5_2->AddText("J.B. 4-Vector Method vs. True");
  tex5_2->SetFillStyle(4000);tex5_2->SetTextFont(63);tex5_2->SetTextSize(10);

  TPaveText *tex6_1 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex6_1->AddText("#frac{True - D.A. Method (using jet)}{True} vs. True");
  tex6_1->SetFillStyle(4000);tex6_1->SetTextFont(63);tex6_1->SetTextSize(10);

  TPaveText *tex6_2 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex6_2->AddText("D.A. Method (using jet) vs. True");
  tex6_2->SetFillStyle(4000);tex6_2->SetTextFont(63);tex6_2->SetTextSize(10);

  TPaveText *tex7_1 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex7_1->AddText("#frac{True - D.A. Method (summing all particles)}{True} vs. True");
  tex7_1->SetFillStyle(4000);tex7_1->SetTextFont(63);tex7_1->SetTextSize(10);

  TPaveText *tex7_2 = new TPaveText(0.1,0.5,0.9,0.7,"NDCNB");
  tex7_2->AddText("D.A. Method (summing all particles) vs. True");
  tex7_2->SetFillStyle(4000);tex7_2->SetTextFont(63);tex7_2->SetTextSize(10);

  //Make Plots

  //----------------------------------
  //Electron Method using ecal energy
  //----------------------------------
  TCanvas *c1a = new TCanvas("c1a");
  c1a->Divide(2,2);
  c1a->cd(1);
  h1a_1->Draw("colz");
  c1a->cd(2);
  h1a_2->Draw("colz");
  c1a->cd(3);
  gPad->SetLogx();h1a_3->Draw("colz");
  c1a->cd(4);tex_energy->Draw();tex1_1->Draw();texy->Draw();

  TCanvas *c1b = new TCanvas("c1b");
  c1b->Divide(2,2);
  c1b->cd(1);
  h1b_1->Draw("colz");
  c1b->cd(2);
  h1b_2->Draw("colz");
  c1b->cd(3);
  gPad->SetLogx();h1b_3->Draw("colz");
  c1b->cd(4);tex_energy->Draw();tex1_1->Draw();texy->Draw();

  TCanvas *c1c = new TCanvas("c1c");
  c1c->Divide(2,2);
  c1c->cd(1);
  h1c_1->Draw("colz");
  c1c->cd(2);
  h1c_2->Draw("colz");
  c1c->cd(3);
  gPad->SetLogx();h1c_3->Draw("colz");
  c1c->cd(4);tex_energy->Draw();tex1_1->Draw();texy->Draw();

  TCanvas *c1d = new TCanvas("c1d");
  c1d->Divide(2,2);
  c1d->cd(1);h1d_1->Draw("colz");func_1->Draw("same");
  c1d->cd(2);h1d_2->Draw("colz");func_1->Draw("same");
  c1d->cd(3);gPad->SetLogx();gPad->SetLogy();h1d_3->Draw("colz");func_1->Draw("same");
  c1d->cd(4);tex_energy->Draw();tex1_2->Draw();texy->Draw();

  //----------------------------------
  //Electron Method using momentum as energy
  //----------------------------------
  TCanvas *c2a = new TCanvas("c2a");
  c2a->Divide(2,2);
  c2a->cd(1);
  h2a_1->Draw("colz");
  c2a->cd(2);
  h2a_2->Draw("colz");
  c2a->cd(3);
  gPad->SetLogx();h2a_3->Draw("colz");
  c2a->cd(4);tex_energy->Draw();tex2_1->Draw();texy->Draw();

  TCanvas *c2b = new TCanvas("c2b");
  c2b->Divide(2,2);
  c2b->cd(1);
  h2b_1->Draw("colz");
  c2b->cd(2);
  h2b_2->Draw("colz");
  c2b->cd(3);
  gPad->SetLogx();h2b_3->Draw("colz");
  c2b->cd(4);tex_energy->Draw();tex2_1->Draw();texy->Draw();

  TCanvas *c2c = new TCanvas("c2c");
  c2c->Divide(2,2);
  c2c->cd(1);
  h2c_1->Draw("colz");
  c2c->cd(2);
  h2c_2->Draw("colz");
  c2c->cd(3);
  gPad->SetLogx();h2c_3->Draw("colz");
  c2c->cd(4);tex_energy->Draw();tex2_1->Draw();texy->Draw();

  TCanvas *c2d = new TCanvas("c2d");
  c2d->Divide(2,2);
  c2d->cd(1);h2d_1->Draw("colz");func_1->Draw("same");
  c2d->cd(2);h2d_2->Draw("colz");func_1->Draw("same");
  c2d->cd(3);gPad->SetLogx();gPad->SetLogy();h2d_3->Draw("colz");func_1->Draw("same");
  c2d->cd(4);tex_energy->Draw();tex2_2->Draw();texy->Draw();

  //----------------------------------
  //JB Method -- using jet
  //----------------------------------
  TCanvas *c3a = new TCanvas("c3a");
  c3a->Divide(2,2);
  c3a->cd(1);
  h3a_1->Draw("colz");
  c3a->cd(2);
  h3a_2->Draw("colz");
  c3a->cd(3);
  gPad->SetLogx();h3a_3->Draw("colz");
  c3a->cd(4);tex_energy->Draw();tex3_1->Draw();texy->Draw();

  TCanvas *c3b = new TCanvas("c3b");
  c3b->Divide(2,2);
  c3b->cd(1);
  h3b_1->Draw("colz");
  c3b->cd(2);
  h3b_2->Draw("colz");
  c3b->cd(3);
  gPad->SetLogx();h3b_3->Draw("colz");
  c3b->cd(4);tex_energy->Draw();tex3_1->Draw();texy->Draw();

  TCanvas *c3c = new TCanvas("c3c");
  c3c->Divide(2,2);
  c3c->cd(1);
  h3c_1->Draw("colz");
  c3c->cd(2);
  h3c_2->Draw("colz");
  c3c->cd(3);
  gPad->SetLogx();h3c_3->Draw("colz");
  c3c->cd(4);tex_energy->Draw();tex3_1->Draw();texy->Draw();

  TCanvas *c3d = new TCanvas("c3d");
  c3d->Divide(2,2);
  c3d->cd(1);h3d_1->Draw("colz");func_1->Draw("same");
  c3d->cd(2);h3d_2->Draw("colz");func_1->Draw("same");
  c3d->cd(3);gPad->SetLogx();gPad->SetLogy();h3d_3->Draw("colz");func_1->Draw("same");
  c3d->cd(4);tex_energy->Draw();tex3_2->Draw();texy->Draw();

  //----------------------------------
  //JB Method -- summing all particles
  //----------------------------------
  TCanvas *c4a = new TCanvas("c4a");
  c4a->Divide(2,2);
  c4a->cd(1);
  h4a_1->Draw("colz");
  c4a->cd(2);
  h4a_2->Draw("colz");
  c4a->cd(3);
  gPad->SetLogx();h4a_3->Draw("colz");
  c4a->cd(4);tex_energy->Draw();tex4_1->Draw();texy->Draw();

  TCanvas *c4b = new TCanvas("c4b");
  c4b->Divide(2,2);
  c4b->cd(1);
  h4b_1->Draw("colz");
  c4b->cd(2);
  h4b_2->Draw("colz");
  c4b->cd(3);
  gPad->SetLogx();h4b_3->Draw("colz");
  c4b->cd(4);tex_energy->Draw();tex4_1->Draw();texy->Draw();

  TCanvas *c4c = new TCanvas("c4c");
  c4c->Divide(2,2);
  c4c->cd(1);
  h4c_1->Draw("colz");
  c4c->cd(2);
  h4c_2->Draw("colz");
  c4c->cd(3);
  gPad->SetLogx();h4c_3->Draw("colz");
  c4c->cd(4);tex_energy->Draw();tex4_1->Draw();texy->Draw();

  TCanvas *c4d = new TCanvas("c4d");
  c4d->Divide(2,2);
  c4d->cd(1);h4d_1->Draw("colz");func_1->Draw("same");
  c4d->cd(2);h4d_2->Draw("colz");func_1->Draw("same");
  c4d->cd(3);gPad->SetLogx();gPad->SetLogy();h4d_3->Draw("colz");func_1->Draw("same");
  c4d->cd(4);tex_energy->Draw();tex4_2->Draw();texy->Draw();

  //----------------------------------
  //JB 4-Vector Method
  //----------------------------------
  TCanvas *c5a = new TCanvas("c5a");
  c5a->Divide(2,2);
  c5a->cd(1);
  h5a_1->Draw("colz");
  c5a->cd(2);
  h5a_2->Draw("colz");
  c5a->cd(3);
  gPad->SetLogx();h5a_3->Draw("colz");
  c5a->cd(4);tex_energy->Draw();tex5_1->Draw();texy->Draw();

  TCanvas *c5b = new TCanvas("c5b");
  c5b->Divide(2,2);
  c5b->cd(1);
  h5b_1->Draw("colz");
  c5b->cd(2);
  h5b_2->Draw("colz");
  c5b->cd(3);
  gPad->SetLogx();h5b_3->Draw("colz");
  c5b->cd(4);tex_energy->Draw();tex5_1->Draw();texy->Draw();

  TCanvas *c5c = new TCanvas("c5c");
  c5c->Divide(2,2);
  c5c->cd(1);
  h5c_1->Draw("colz");
  c5c->cd(2);
  h5c_2->Draw("colz");
  c5c->cd(3);
  gPad->SetLogx();h5c_3->Draw("colz");
  c5c->cd(4);tex_energy->Draw();tex5_1->Draw();texy->Draw();

  TCanvas *c5d = new TCanvas("c5d");
  c5d->Divide(2,2);
  c5d->cd(1);h5d_1->Draw("colz");func_1->Draw("same");
  c5d->cd(2);h5d_2->Draw("colz");func_1->Draw("same");
  c5d->cd(3);gPad->SetLogx();gPad->SetLogy();h5d_3->Draw("colz");func_1->Draw("same");
  c5d->cd(4);tex_energy->Draw();tex5_2->Draw();texy->Draw();

  //----------------------------------
  //DA Method -- using jet
  //----------------------------------
  TCanvas *c6a = new TCanvas("c6a");
  c6a->Divide(2,2);
  c6a->cd(1);
  h6a_1->Draw("colz");
  c6a->cd(2);
  h6a_2->Draw("colz");
  c6a->cd(3);
  gPad->SetLogx();h6a_3->Draw("colz");
  c6a->cd(4);tex_energy->Draw();tex6_1->Draw();texy->Draw();

  TCanvas *c6b = new TCanvas("c6b");
  c6b->Divide(2,2);
  c6b->cd(1);
  h6b_1->Draw("colz");
  c6b->cd(2);
  h6b_2->Draw("colz");
  c6b->cd(3);
  gPad->SetLogx();h6b_3->Draw("colz");
  c6b->cd(4);tex_energy->Draw();tex6_1->Draw();texy->Draw();

  TCanvas *c6c = new TCanvas("c6c");
  c6c->Divide(2,2);
  c6c->cd(1);
  h6c_1->Draw("colz");
  c6c->cd(2);
  h6c_2->Draw("colz");
  c6c->cd(3);
  gPad->SetLogx();h6c_3->Draw("colz");
  c6c->cd(4);tex_energy->Draw();tex6_1->Draw();texy->Draw();

  TCanvas *c6d = new TCanvas("c6d");
  c6d->Divide(2,2);
  c6d->cd(1);h6d_1->Draw("colz");func_1->Draw("same");
  c6d->cd(2);h6d_2->Draw("colz");func_1->Draw("same");
  c6d->cd(3);gPad->SetLogx();gPad->SetLogy();h6d_3->Draw("colz");func_1->Draw("same");
  c6d->cd(4);tex_energy->Draw();tex6_2->Draw();texy->Draw();

  //----------------------------------
  //DA Method -- summing all particles
  //----------------------------------
  TCanvas *c7a = new TCanvas("c7a");
  c7a->Divide(2,2);
  c7a->cd(1);
  h7a_1->Draw("colz");
  c7a->cd(2);
  h7a_2->Draw("colz");
  c7a->cd(3);
  gPad->SetLogx();h7a_3->Draw("colz");
  c7a->cd(4);tex_energy->Draw();tex7_1->Draw();texy->Draw();

  TCanvas *c7b = new TCanvas("c7b");
  c7b->Divide(2,2);
  c7b->cd(1);
  h7b_1->Draw("colz");
  c7b->cd(2);
  h7b_2->Draw("colz");
  c7b->cd(3);
  gPad->SetLogx();h7b_3->Draw("colz");
  c7b->cd(4);tex_energy->Draw();tex7_1->Draw();texy->Draw();

  TCanvas *c7c = new TCanvas("c7c");
  c7c->Divide(2,2);
  c7c->cd(1);
  h7c_1->Draw("colz");
  c7c->cd(2);
  h7c_2->Draw("colz");
  c7c->cd(3);
  gPad->SetLogx();h7c_3->Draw("colz");
  c7c->cd(4);tex_energy->Draw();tex7_1->Draw();texy->Draw();

  TCanvas *c7d = new TCanvas("c7d");
  c7d->Divide(2,2);
  c7d->cd(1);h7d_1->Draw("colz");func_1->Draw("same");
  c7d->cd(2);h7d_2->Draw("colz");func_1->Draw("same");
  c7d->cd(3);gPad->SetLogx();gPad->SetLogy();h7d_3->Draw("colz");func_1->Draw("same");
  c7d->cd(4);tex_energy->Draw();tex7_2->Draw();texy->Draw();

  //Print to File
  if(energy_set == 1){
    c1a->Print("./plots/kinematic_resolution_ycut_5_41.pdf[");
    c1a->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c1b->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c1c->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c1d->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c2a->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c2b->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c2c->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c2d->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c3a->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c3b->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c3c->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c3d->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c4a->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c4b->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c4c->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c4d->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c5a->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c5b->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c5c->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c5d->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c6a->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c6b->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c6c->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c6d->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c7a->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c7b->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c7c->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c7d->Print("./plots/kinematic_resolution_ycut_5_41.pdf");
    c7d->Print("./plots/kinematic_resolution_ycut_5_41.pdf]");
  }
  if(energy_set == 2){
    c1a->Print("./plots/kinematic_resolution_ycut_18_275.pdf[");
    c1a->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c1b->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c1c->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c1d->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c2a->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c2b->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c2c->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c2d->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c3a->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c3b->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c3c->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c3d->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c4a->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c4b->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c4c->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c4d->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c5a->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c5b->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c5c->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c5d->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c6a->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c6b->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c6c->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c6d->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c7a->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c7b->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c7c->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c7d->Print("./plots/kinematic_resolution_ycut_18_275.pdf");
    c7d->Print("./plots/kinematic_resolution_ycut_18_275.pdf]");
  }
}