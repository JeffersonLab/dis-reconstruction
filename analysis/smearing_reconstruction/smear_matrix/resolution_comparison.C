R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "includes/energy_resolution_optimization.h"
#include "includes/PadMxN.h"

void resolution_comparison(){

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

  //Profile Histograms
  //We want to resolution, which we get by doing GetBinError() with the "s" option.
  //Otherwise it is the error on the mean.
  //All events (Q2>0.5)
  TProfile *prof1 = new TProfile("prof1","",30,0,1,-300,300,"s"); //Electron % deviation vs. y_true
  TProfile *prof2 = new TProfile("prof2","",30,0,1,-300,300,"s"); //JB % deviation vs. y_true
  TProfile *prof3 = new TProfile("prof3","",30,0,1,-300,300,"s"); //DA % deviation vs. y_true

  TProfile *prof4 = new TProfile("prof4","",30,0,180,-300,300,"s"); //Electron % deviation vs. theta_e
  TProfile *prof5 = new TProfile("prof5","",30,0,180,-300,300,"s"); //JB % deviation vs. theta_e
  TProfile *prof6 = new TProfile("prof6","",30,0,180,-300,300,"s"); //DA % deviation vs. theta_e

  TProfile *prof7 = new TProfile("prof7","",30,0,180,-300,300,"s"); //Electron % deviation vs. theta_h
  TProfile *prof8 = new TProfile("prof8","",30,0,180,-300,300,"s"); //JB % deviation vs. theta_h
  TProfile *prof9 = new TProfile("prof9","",30,0,180,-300,300,"s"); //DA % deviation vs. theta_h

  //Q2>10
  TProfile *prof1a = new TProfile("prof1a","",30,0,1,-300,300,"s"); //Electron % deviation vs. y_true
  TProfile *prof2a = new TProfile("prof2a","",30,0,1,-300,300,"s"); //JB % deviation vs. y_true
  TProfile *prof3a = new TProfile("prof3a","",30,0,1,-300,300,"s"); //DA % deviation vs. y_true

  TProfile *prof4a = new TProfile("prof4a","",30,0,180,-300,300,"s"); //Electron % deviation vs. theta_e
  TProfile *prof5a = new TProfile("prof5a","",30,0,180,-300,300,"s"); //JB % deviation vs. theta_e
  TProfile *prof6a = new TProfile("prof6a","",30,0,180,-300,300,"s"); //DA % deviation vs. theta_e

  TProfile *prof7a = new TProfile("prof7a","",30,0,180,-300,300,"s"); //Electron % deviation vs. theta_h
  TProfile *prof8a = new TProfile("prof8a","",30,0,180,-300,300,"s"); //JB % deviation vs. theta_h
  TProfile *prof9a = new TProfile("prof9a","",30,0,180,-300,300,"s"); //DA % deviation vs. theta_h

  //Q2>100
  TProfile *prof1b = new TProfile("prof1b","",30,0,1,-300,300,"s"); //Electron % deviation vs. y_true
  TProfile *prof2b = new TProfile("prof2b","",30,0,1,-300,300,"s"); //JB % deviation vs. y_true
  TProfile *prof3b = new TProfile("prof3b","",30,0,1,-300,300,"s"); //DA % deviation vs. y_true

  TProfile *prof4b = new TProfile("prof4b","",30,0,180,-300,300,"s"); //Electron % deviation vs. theta_e
  TProfile *prof5b = new TProfile("prof5b","",30,0,180,-300,300,"s"); //JB % deviation vs. theta_e
  TProfile *prof6b = new TProfile("prof6b","",30,0,180,-300,300,"s"); //DA % deviation vs. theta_e

  TProfile *prof7b = new TProfile("prof7b","",30,0,180,-300,300,"s"); //Electron % deviation vs. theta_h
  TProfile *prof8b = new TProfile("prof8b","",30,0,180,-300,300,"s"); //JB % deviation vs. theta_h
  TProfile *prof9b = new TProfile("prof9b","",30,0,180,-300,300,"s"); //DA % deviation vs. theta_h

  //Two Dimensional profiles and histograms

  //Q2 Binning                                                                                                                                        
  double Q2_min = 1E-1;  
  double Q2_max = 1E4;
  const int nbins_Q2 = 25;
  double log_bw_Q2 = (log10(Q2_max) - log10(Q2_min))/(nbins_Q2);   //Determine bin width                                                              
  double log_Q2_div;
  double Q2_bins[nbins_Q2+1];
  for(int i=0;i<nbins_Q2+1;i++){
    log_Q2_div = log10(Q2_min) + (i*log_bw_Q2);
    Q2_bins[i] = pow(10,log_Q2_div);
  }

  //x Binning                                                                                                                                        
  double x_min = 1E-5;
  double x_max = 1;
  const int nbins_x = 25;
  double log_bw_x = (log10(x_max) - log10(x_min))/(nbins_x);   //Determine bin width                                                                  
  double log_x_div;
  double x_bins[nbins_x+1];
  for(int i=0;i<nbins_x+1;i++){
    log_x_div = log10(x_min) + (i*log_bw_x);
    x_bins[i] = pow(10,log_x_div);
  }

  TProfile2D *prof10 = new TProfile2D("prof10","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //Electron Method
  TProfile2D *prof11 = new TProfile2D("prof11","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //JB Method
  TProfile2D *prof12 = new TProfile2D("prof12","",nbins_x,x_bins,nbins_Q2,Q2_bins,"s"); //DA Method

  TH2 *hh1 = new TH2D("hh1","",nbins_x,x_bins,nbins_Q2,Q2_bins); //Electron Method
  hh1->GetXaxis()->SetTitle("x");hh1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh1->GetXaxis()->SetLabelFont(63);hh1->GetYaxis()->SetLabelFont(63);
  hh1->GetXaxis()->SetLabelSize(25);hh1->GetYaxis()->SetLabelSize(25);
  hh1->GetXaxis()->SetLabelOffset(0.01);hh1->GetYaxis()->SetLabelOffset(0.01);
  hh1->GetXaxis()->CenterTitle(1);hh1->GetYaxis()->CenterTitle(1);
  hh1->GetXaxis()->SetTitleSize(40);hh1->GetXaxis()->SetTitleOffset(2.5); 
  hh1->GetYaxis()->SetTitleSize(40);hh1->GetYaxis()->SetTitleOffset(3.0);

  TH2 *hh2 = new TH2D("hh2","",nbins_x,x_bins,nbins_Q2,Q2_bins); //JB Method
  hh2->GetXaxis()->SetTitle("x");hh2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh2->GetXaxis()->SetLabelFont(63);hh2->GetYaxis()->SetLabelFont(63);
  hh2->GetXaxis()->SetLabelSize(25);hh2->GetYaxis()->SetLabelSize(25);
  hh2->GetXaxis()->SetLabelOffset(0.01);hh2->GetYaxis()->SetLabelOffset(0.01);
  hh2->GetXaxis()->CenterTitle(1);hh2->GetYaxis()->CenterTitle(1);
  hh2->GetXaxis()->SetTitleSize(40);hh2->GetXaxis()->SetTitleOffset(2.5); 
  hh2->GetYaxis()->SetTitleSize(40);hh2->GetYaxis()->SetTitleOffset(3.0);
  
  TH2 *hh3 = new TH2D("hh3","",nbins_x,x_bins,nbins_Q2,Q2_bins); //DA Method
  hh3->GetXaxis()->SetTitle("x");hh3->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	hh3->GetXaxis()->SetLabelFont(63);hh3->GetYaxis()->SetLabelFont(63);
  hh3->GetXaxis()->SetLabelSize(25);hh3->GetYaxis()->SetLabelSize(25);
  hh3->GetXaxis()->SetLabelOffset(0.01);hh3->GetYaxis()->SetLabelOffset(0.01);
  hh3->GetXaxis()->CenterTitle(1);hh3->GetYaxis()->CenterTitle(1);
  hh3->GetXaxis()->SetTitleSize(40);hh3->GetXaxis()->SetTitleOffset(2.5); 
  hh3->GetYaxis()->SetTitleSize(40);hh3->GetYaxis()->SetTitleOffset(3.0);

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
  Double_t Theta_h_calc(0);
  Double_t Theta_e_true(0);

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
  Double_t px_sum_noh(0),py_sum_noh(0),pz_sum_noh(0),pt_sum_noh(0),E_sum_noh(0); //Sum without barrel HCAL
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

  Int_t nevents = tree->GetEntries();
  
  // Loop over all events
  for(int i=0;i<nevents;i++){
    if(i%10000==0) cout<<"Events Analysed = "<<i<<"!"<<endl;

    //Reset FastJet Variables
    holdParticles.clear();

    //Reset Summing Over Hadrons Variables
    Etot_sumh_s = 0; pxtot_sumh_s = 0; pytot_sumh_s = 0; pztot_sumh_s = 0; pttot_sumh_s = 0;
    E_sum_noh = 0; px_sum_noh = 0; py_sum_noh = 0; pz_sum_noh = 0; pt_sum_noh = 0;

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
          Theta_e_true = particle->GetTheta()*TMath::RadToDeg();
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

        int sigma_best;
 
	      //Get Smeared-Scattered Electron
	      if(Status_s[j]==1 && Status[j]==1 && id[j]==11 && orig[j]==3){
	        	        
          //If return 0, we use calorimeter. If return 1, we use tracker.
          sigma_best = energy_resolution_optimization(particle->GetE(), particle->GetP(), particle->GetEta(), id[j], particle_s->IsESmeared(), particle_s->IsPSmeared());

          if(sigma_best == 0){
            if(particle_s->GetE() > mass[j]){
              Ef_e_s = particle_s->GetE();
              pxf_e_s = sqrt(Ef_e_s*Ef_e_s - mass[j]*mass[j]) * sin( particle_s->GetTheta() ) * cos( particle_s->GetPhi() );
              pyf_e_s = sqrt(Ef_e_s*Ef_e_s - mass[j]*mass[j]) * sin( particle_s->GetTheta() ) * sin( particle_s->GetPhi() );
              pzf_e_s = sqrt(Ef_e_s*Ef_e_s - mass[j]*mass[j]) * cos( particle_s->GetTheta() );
            } else{
              Ef_e_s = particle_s->GetE();
              pxf_e_s = 0;
              pyf_e_s = 0;
              pzf_e_s = 0;
            }

            theta_e_s = particle_s->GetTheta();
	          ef_s.SetPxPyPzE(pxf_e_s,pyf_e_s,pzf_e_s,Ef_e_s);
            detected_elec = true;

          } else if(sigma_best == 1){
            pxf_e_s = particle_s->GetPx();
            pyf_e_s = particle_s->GetPy();
            pzf_e_s = particle_s->GetPz();
            Ef_e_s = sqrt( pxf_e_s*pxf_e_s + pyf_e_s* pyf_e_s + pzf_e_s*pzf_e_s + mass[j]*mass[j] );

            theta_e_s = particle_s->GetTheta();
	          ef_s.SetPxPyPzE(pxf_e_s,pyf_e_s,pzf_e_s,Ef_e_s);
            detected_elec = true;

          }else{
            detected_elec = false;
          }
	      }

	      //Put hadrons into PseudoJet object
	      if(j!=electronIndex && Status_s[j]==1){
                   
          //If return 0, we use calorimeter. If return 1, we use tracker.
          sigma_best = energy_resolution_optimization(particle->GetE(), particle->GetP(), particle->GetEta(), id[j], particle_s->IsESmeared(), particle_s->IsPSmeared());

          if(sigma_best == 0){
            if(particle_s->GetE() > mass[j]){
              holdE_h = particle_s->GetE();
              holdpx_h = sqrt(holdE_h*holdE_h - mass[j]*mass[j]) * sin( particle_s->GetTheta() ) * cos( particle_s->GetPhi() );
              holdpy_h = sqrt(holdE_h*holdE_h - mass[j]*mass[j]) * sin( particle_s->GetTheta() ) * sin( particle_s->GetPhi() );
              holdpz_h = sqrt(holdE_h*holdE_h - mass[j]*mass[j]) * cos( particle_s->GetTheta() );
            } else{
              holdE_h = particle_s->GetE();
              holdpx_h = 0;
              holdpy_h = 0;
              holdpz_h = 0;
            }
          } else if(sigma_best == 1){
            holdpx_h = particle_s->GetPx();
            holdpy_h = particle_s->GetPy();
            holdpz_h = particle_s->GetPz();
            holdE_h = sqrt( holdpx_h*holdpx_h + holdpy_h* holdpy_h + holdpz_h*holdpz_h + mass[j]*mass[j] );
          } else{
            continue;
          }

          pxtot_sumh_s = pxtot_sumh_s + holdpx_h;
          pytot_sumh_s = pytot_sumh_s + holdpy_h;
          pztot_sumh_s = pztot_sumh_s + holdpz_h;
          Etot_sumh_s =  Etot_sumh_s + holdE_h;

          //Fast Jet
	        fourvec[0] = (double) holdpx_h;
	        fourvec[1] = (double) holdpy_h;
	        fourvec[2] = (double) holdpz_h;
	        fourvec[3] = (double) holdE_h;

	        fastjet::PseudoJet finalStateParticle(fourvec);
	        finalStateParticle.set_user_index(j);
	        holdParticles.push_back(finalStateParticle);

          //Calculation if barrel HCAL (eta=[-1,1]) is removed
          //Do not include neutrons (id=2112) or K0_L (id=130)
          if( !(id[j]==2112 || id[j]==130) || !(particle->GetEta()>-1.0 && particle->GetEta()<1.0 ) ){
            px_sum_noh = px_sum_noh + holdpx_h;
            py_sum_noh = py_sum_noh + holdpy_h;
            pz_sum_noh = pz_sum_noh + holdpz_h;
            E_sum_noh =  E_sum_noh + holdE_h;
          }

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

    Theta_h_calc = TMath::RadToDeg()*TMath::ACos( ( -y_e*Ei_e +(1.-y_e)*x_e*E_pn ) / (y_e*Ei_e + (1.-y_e)*x_e*E_pn) );

    //-------------Calculate *Smeared* invariants using 4-vectors------------------//
  
    if(detected_elec && y_e>y_min && y_e<y_max){
      //1.1) Smeared scattered election using optimization
      Q2_e_nm_s = 4.*Ei_e*Ef_e_s*TMath::Cos(theta_e_s/2.)*TMath::Cos(theta_e_s/2.);
      y_e_nm_s = 1. - ( (Ef_e_s/(2.*Ei_e))*(1. - TMath::Cos(theta_e_s)) );
      x_e_nm_s = Q2_e_nm_s/(s_nm*y_e_nm_s);

      prof1->Fill(y_e, 100.*(y_e-y_e_nm_s)/y_e);
      prof4->Fill(Theta_e_true,100.*(y_e-y_e_nm_s)/y_e);
      prof7->Fill(Theta_h_calc,100.*(y_e-y_e_nm_s)/y_e);

      if(Q2_e>10){
        prof1a->Fill(y_e, 100.*(y_e-y_e_nm_s)/y_e);
        prof4a->Fill(Theta_e_true,100.*(y_e-y_e_nm_s)/y_e);
        prof7a->Fill(Theta_h_calc,100.*(y_e-y_e_nm_s)/y_e);
      }
      if(Q2_e>100){
        prof1b->Fill(y_e, 100.*(y_e-y_e_nm_s)/y_e);
        prof4b->Fill(Theta_e_true,100.*(y_e-y_e_nm_s)/y_e);
        prof7b->Fill(Theta_h_calc,100.*(y_e-y_e_nm_s)/y_e);
      }
      if( fabs(100.*(y_e-y_e_nm_s)/y_e)<300. )
        prof10->Fill(x_e,Q2_e,100.*(y_e-y_e_nm_s)/y_e);

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
        
        //3.1) Using Smeared DA Method
        if(detected_elec){
          Theta_h_nm_s = 2.* TMath::ATan( (Etot_jet_s - pztot_jet_s)/pttot_jet_s);    
          Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(theta_e_s/2.) )*
                   ( 1./(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
          y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
          x_da_s = Q2_da_s/(s_nm*y_da_s);  
        }

      }//jets.empty()

      //2.2) Using Smeared JB Method (Summing Over all particles)
      y_jb_sumh_s = (Etot_sumh_s - pztot_sumh_s)/(2*Ei_e);
      Q2_jb_sumh_s = (pttot_sumh_s*pttot_sumh_s)/(1. - y_jb_sumh_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s);

      prof2->Fill(y_e, 100.*(y_e-y_jb_sumh_s)/y_e);
      prof5->Fill(Theta_e_true,100.*(y_e-y_jb_sumh_s)/y_e);
      prof8->Fill(Theta_h_calc,100.*(y_e-y_jb_sumh_s)/y_e);

      if(Q2_e>10){
        prof2a->Fill(y_e, 100.*(y_e-y_jb_sumh_s)/y_e);
        prof5a->Fill(Theta_e_true,100.*(y_e-y_jb_sumh_s)/y_e);
        prof8a->Fill(Theta_h_calc,100.*(y_e-y_jb_sumh_s)/y_e);
      }
      if(Q2_e>100){
        prof2b->Fill(y_e, 100.*(y_e-y_jb_sumh_s)/y_e);
        prof5b->Fill(Theta_e_true,100.*(y_e-y_jb_sumh_s)/y_e);
        prof8b->Fill(Theta_h_calc,100.*(y_e-y_jb_sumh_s)/y_e);
      }

      if( fabs(100.*(y_e-y_jb_sumh_s)/y_e)<300. )
        prof11->Fill(x_e,Q2_e,100.*(y_e-y_jb_sumh_s)/y_e);

      //2.3) Using 4-Vector of outgoing X
      eh_s_tot.SetPxPyPzE(pxtot_sumh_s,pytot_sumh_s,pztot_sumh_s,Etot_sumh_s);

      q_h_s= (eh_s_tot - pni);
      Q2_jb_sumh_s = -1.0*(q_h_s*q_h_s);
      y_jb_sumh_s = (pni*q_h_s)/(pni*ei);
      //x_jb_sumh_s = Q2_jb_sumh_s/(2*pni*q_h_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s); //Using Q2 = xsy, instead of 4-vector product Q2_jb_sumh_s/(2*pni*q_h_s); 

      //3.2) Using Smeared DA Method (Summing Over Hadrons)
      if(detected_elec){

        Theta_h_nm_s = 2.* TMath::ATan( (Etot_sumh_s - pztot_sumh_s)/pttot_sumh_s);
        Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(theta_e_s/2.) )*
          ( 1./(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
        y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
        x_da_s = Q2_da_s/(s_nm*y_da_s);

        prof3->Fill(y_e, 100.*(y_e-y_da_s)/y_e);
        prof6->Fill(Theta_e_true,100.*(y_e-y_da_s)/y_e);
        prof9->Fill(Theta_h_calc,100.*(y_e-y_da_s)/y_e);

        if(Q2_e>10){
          prof3a->Fill(y_e, 100.*(y_e-y_da_s)/y_e);
          prof6a->Fill(Theta_e_true,100.*(y_e-y_da_s)/y_e);
          prof9a->Fill(Theta_h_calc,100.*(y_e-y_da_s)/y_e);
        }
        if(Q2_e>100){
          prof3b->Fill(y_e, 100.*(y_e-y_da_s)/y_e);
          prof6b->Fill(Theta_e_true,100.*(y_e-y_da_s)/y_e);
          prof9b->Fill(Theta_h_calc,100.*(y_e-y_da_s)/y_e);
        }

        if( fabs(100.*(y_e-y_da_s)/y_e)<300. )
          prof12->Fill(x_e,Q2_e,100.*(y_e-y_da_s)/y_e);

      }

    }//holdParticles.empty()

    //reset variables
    detected_elec = false;

  }//Finished Event Loop

  //Make Plots
  //Extract Resolutions from TProfile
  TGraph *gr1 = new TGraph();
  gr1->SetMarkerColor(kRed);gr1->SetLineColor(kRed);gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.5);gr1->SetLineWidth(2);
  TGraph *gr2 = new TGraph();
  gr2->SetMarkerColor(kBlue);gr2->SetLineColor(kBlue);gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(1.5);gr2->SetLineWidth(2);
  TGraph *gr3 = new TGraph();
  gr3->SetMarkerColor(kGreen);gr3->SetLineColor(kGreen);gr3->SetMarkerStyle(22);
  gr3->SetMarkerSize(1.5);gr3->SetLineWidth(2);

  TGraph *gr1a = new TGraph();
  gr1a->SetMarkerColor(kRed);gr1a->SetLineColor(kRed);gr1a->SetMarkerStyle(20);
  gr1a->SetMarkerSize(1.5);gr1a->SetLineWidth(2);
  TGraph *gr2a = new TGraph();
  gr2a->SetMarkerColor(kBlue);gr2a->SetLineColor(kBlue);gr2a->SetMarkerStyle(21);
  gr2a->SetMarkerSize(1.5);gr2a->SetLineWidth(2);
  TGraph *gr3a = new TGraph();
  gr3a->SetMarkerColor(kGreen);gr3a->SetLineColor(kGreen);gr3a->SetMarkerStyle(22);
  gr3a->SetMarkerSize(1.5);gr3a->SetLineWidth(2);

  TGraph *gr1b = new TGraph();
  gr1b->SetMarkerColor(kRed);gr1b->SetLineColor(kRed);gr1b->SetMarkerStyle(20);
  gr1b->SetMarkerSize(1.5);gr1b->SetLineWidth(2);
  TGraph *gr2b = new TGraph();
  gr2b->SetMarkerColor(kBlue);gr2b->SetLineColor(kBlue);gr2b->SetMarkerStyle(21);
  gr2b->SetMarkerSize(1.5);gr2b->SetLineWidth(2);
  TGraph *gr3b = new TGraph();
  gr3b->SetMarkerColor(kGreen);gr3b->SetLineColor(kGreen);gr3b->SetMarkerStyle(22);
  gr3b->SetMarkerSize(1.5);gr3b->SetLineWidth(2);

  int counter1(0),counter1a(0),counter1b(0);
  int counter2(0),counter2a(0),counter2b(0);
  int counter3(0),counter3a(0),counter3b(0);

  for(int ibin=1;ibin<=prof1->GetNbinsX();ibin++){

    if(prof1->GetBinError(ibin)>0){
      gr1->SetPoint(counter1,prof1->GetBinCenter(ibin),prof1->GetBinError(ibin));
      counter1++;
    }
    
    if(prof2->GetBinError(ibin)>0){
      gr2->SetPoint(counter2,prof2->GetBinCenter(ibin),prof2->GetBinError(ibin));
      counter2++;
    }
    
    if(prof3->GetBinError(ibin)>0){
      gr3->SetPoint(counter3,prof3->GetBinCenter(ibin),prof3->GetBinError(ibin));
      counter3++;
    }

    if(prof1a->GetBinError(ibin)>0){
      gr1a->SetPoint(counter1a,prof1a->GetBinCenter(ibin),prof1a->GetBinError(ibin));
      counter1a++;
    }
    
    if(prof2a->GetBinError(ibin)>0){
      gr2a->SetPoint(counter2a,prof2a->GetBinCenter(ibin),prof2a->GetBinError(ibin));
      counter2a++;
    }
    
    if(prof3a->GetBinError(ibin)>0){
      gr3a->SetPoint(counter3a,prof3a->GetBinCenter(ibin),prof3a->GetBinError(ibin));
      counter3a++;
    }

    if(prof1b->GetBinError(ibin)>0){
      gr1b->SetPoint(counter1b,prof1b->GetBinCenter(ibin),prof1b->GetBinError(ibin));
      counter1b++;
    }
    
    if(prof2b->GetBinError(ibin)>0){
      gr2b->SetPoint(counter2b,prof2b->GetBinCenter(ibin),prof2b->GetBinError(ibin));
      counter2b++;
    }
    
    if(prof3b->GetBinError(ibin)>0){
      gr3b->SetPoint(counter3b,prof3b->GetBinCenter(ibin),prof3b->GetBinError(ibin));
      counter3b++;
    }
 
  }

  TGraph *gr4 = new TGraph();
  gr4->SetMarkerColor(kRed);gr4->SetLineColor(kRed);gr4->SetMarkerStyle(20);
  gr4->SetMarkerSize(1.5);gr4->SetLineWidth(2);
  TGraph *gr5 = new TGraph();
  gr5->SetMarkerColor(kBlue);gr5->SetLineColor(kBlue);gr5->SetMarkerStyle(21);
  gr5->SetMarkerSize(1.5);gr5->SetLineWidth(2);
  TGraph *gr6 = new TGraph();
  gr6->SetMarkerColor(kGreen);gr6->SetLineColor(kGreen);gr6->SetMarkerStyle(22);
  gr6->SetMarkerSize(1.5);gr6->SetLineWidth(2);

  TGraph *gr4a = new TGraph();
  gr4a->SetMarkerColor(kRed);gr4a->SetLineColor(kRed);gr4a->SetMarkerStyle(20);
  gr4a->SetMarkerSize(1.5);gr4a->SetLineWidth(2);
  TGraph *gr5a = new TGraph();
  gr5a->SetMarkerColor(kBlue);gr5a->SetLineColor(kBlue);gr5a->SetMarkerStyle(21);
  gr5a->SetMarkerSize(1.5);gr5a->SetLineWidth(2);
  TGraph *gr6a = new TGraph();
  gr6a->SetMarkerColor(kGreen);gr6a->SetLineColor(kGreen);gr6a->SetMarkerStyle(22);
  gr6a->SetMarkerSize(1.5);gr6a->SetLineWidth(2);

  TGraph *gr4b = new TGraph();
  gr4b->SetMarkerColor(kRed);gr4b->SetLineColor(kRed);gr4b->SetMarkerStyle(20);
  gr4b->SetMarkerSize(1.5);gr4b->SetLineWidth(2);
  TGraph *gr5b = new TGraph();
  gr5b->SetMarkerColor(kBlue);gr5b->SetLineColor(kBlue);gr5b->SetMarkerStyle(21);
  gr5b->SetMarkerSize(1.5);gr5b->SetLineWidth(2);
  TGraph *gr6b = new TGraph();
  gr6b->SetMarkerColor(kGreen);gr6b->SetLineColor(kGreen);gr6b->SetMarkerStyle(22);
  gr6b->SetMarkerSize(1.5);gr6b->SetLineWidth(2);

  int counter4(0),counter4a(0),counter4b(0);
  int counter5(0),counter5a(0),counter5b(0);
  int counter6(0),counter6a(0),counter6b(0);

  for(int ibin=1;ibin<=prof4->GetNbinsX();ibin++){
    if(prof4->GetBinError(ibin)>0){
      gr4->SetPoint(counter4,prof4->GetBinCenter(ibin),prof4->GetBinError(ibin));
      counter4++;
    }
    
    if(prof5->GetBinError(ibin)>0){
      gr5->SetPoint(counter5,prof5->GetBinCenter(ibin),prof5->GetBinError(ibin));
      counter5++;
    }
    
    if(prof6->GetBinError(ibin)>0){
      gr6->SetPoint(counter6,prof6->GetBinCenter(ibin),prof6->GetBinError(ibin));
      counter6++;
    }

    if(prof4a->GetBinError(ibin)>0){
      gr4a->SetPoint(counter4a,prof4a->GetBinCenter(ibin),prof4a->GetBinError(ibin));
      counter4a++;
    }
    
    if(prof5a->GetBinError(ibin)>0){
      gr5a->SetPoint(counter5a,prof5a->GetBinCenter(ibin),prof5a->GetBinError(ibin));
      counter5a++;
    }
    
    if(prof6a->GetBinError(ibin)>0){
      gr6a->SetPoint(counter6a,prof6a->GetBinCenter(ibin),prof6a->GetBinError(ibin));
      counter6a++;
    }

    if(prof4b->GetBinError(ibin)>0){
      gr4b->SetPoint(counter4b,prof4b->GetBinCenter(ibin),prof4b->GetBinError(ibin));
      counter4b++;
    }
    
    if(prof5b->GetBinError(ibin)>0){
      gr5b->SetPoint(counter5b,prof5b->GetBinCenter(ibin),prof5b->GetBinError(ibin));
      counter5b++;
    }
    
    if(prof6b->GetBinError(ibin)>0){
      gr6b->SetPoint(counter6b,prof6b->GetBinCenter(ibin),prof6b->GetBinError(ibin));
      counter6b++;
    }
  }

  TGraph *gr7 = new TGraph();
  gr7->SetMarkerColor(kRed);gr7->SetLineColor(kRed);gr7->SetMarkerStyle(20);
  gr7->SetMarkerSize(1.5);gr7->SetLineWidth(2);
  TGraph *gr8 = new TGraph();
  gr8->SetMarkerColor(kBlue);gr8->SetLineColor(kBlue);gr8->SetMarkerStyle(21);
  gr8->SetMarkerSize(1.5);gr8->SetLineWidth(2);
  TGraph *gr9 = new TGraph();
  gr9->SetMarkerColor(kGreen);gr9->SetLineColor(kGreen);gr9->SetMarkerStyle(22);
  gr9->SetMarkerSize(1.5);gr9->SetLineWidth(2);

  TGraph *gr7a = new TGraph();
  gr7a->SetMarkerColor(kRed);gr7a->SetLineColor(kRed);gr7a->SetMarkerStyle(20);
  gr7a->SetMarkerSize(1.5);gr7a->SetLineWidth(2);
  TGraph *gr8a = new TGraph();
  gr8a->SetMarkerColor(kBlue);gr8a->SetLineColor(kBlue);gr8a->SetMarkerStyle(21);
  gr8a->SetMarkerSize(1.5);gr8a->SetLineWidth(2);
  TGraph *gr9a = new TGraph();
  gr9a->SetMarkerColor(kGreen);gr9a->SetLineColor(kGreen);gr9a->SetMarkerStyle(22);
  gr9a->SetMarkerSize(1.5);gr9a->SetLineWidth(2);

  TGraph *gr7b = new TGraph();
  gr7b->SetMarkerColor(kRed);gr7b->SetLineColor(kRed);gr7b->SetMarkerStyle(20);
  gr7b->SetMarkerSize(1.5);gr7b->SetLineWidth(2);
  TGraph *gr8b = new TGraph();
  gr8b->SetMarkerColor(kBlue);gr8b->SetLineColor(kBlue);gr8b->SetMarkerStyle(21);
  gr8b->SetMarkerSize(1.5);gr8b->SetLineWidth(2);
  TGraph *gr9b = new TGraph();
  gr9b->SetMarkerColor(kGreen);gr9b->SetLineColor(kGreen);gr9b->SetMarkerStyle(22);
  gr9b->SetMarkerSize(1.5);gr9b->SetLineWidth(2);

  int counter7(0),counter7a(0),counter7b(0);
  int counter8(0),counter8a(0),counter8b(0);
  int counter9(0),counter9a(0),counter9b(0);

  for(int ibin=1;ibin<=prof7->GetNbinsX();ibin++){
    if(prof7->GetBinError(ibin)>0){
      gr7->SetPoint(counter7,prof7->GetBinCenter(ibin),prof7->GetBinError(ibin));
      counter7++;
    }
    
    if(prof8->GetBinError(ibin)>0){
      gr8->SetPoint(counter8,prof8->GetBinCenter(ibin),prof8->GetBinError(ibin));
      counter8++;
    }
    
    if(prof9->GetBinError(ibin)>0){
      gr9->SetPoint(counter9,prof9->GetBinCenter(ibin),prof9->GetBinError(ibin));
      counter9++;
    }

    if(prof7a->GetBinError(ibin)>0){
      gr7a->SetPoint(counter7a,prof7a->GetBinCenter(ibin),prof7a->GetBinError(ibin));
      counter7a++;
    }
    
    if(prof8a->GetBinError(ibin)>0){
      gr8a->SetPoint(counter8a,prof8a->GetBinCenter(ibin),prof8a->GetBinError(ibin));
      counter8a++;
    }
    
    if(prof9a->GetBinError(ibin)>0){
      gr9a->SetPoint(counter9a,prof9a->GetBinCenter(ibin),prof9a->GetBinError(ibin));
      counter9a++;
    }

    if(prof7b->GetBinError(ibin)>0){
      gr7b->SetPoint(counter7b,prof7b->GetBinCenter(ibin),prof7b->GetBinError(ibin));
      counter7b++;
    }
    
    if(prof8b->GetBinError(ibin)>0){
      gr8b->SetPoint(counter8b,prof8b->GetBinCenter(ibin),prof8b->GetBinError(ibin));
      counter8b++;
    }
    
    if(prof9b->GetBinError(ibin)>0){
      gr9b->SetPoint(counter9b,prof9b->GetBinCenter(ibin),prof9b->GetBinError(ibin));
      counter9b++;
    }
  }

  //Plot 1
  PadMxN *pad2x2_1 = new PadMxN("c1",500,500,150,150,100,125,2,2);
  pad2x2_1->Draw();
  TPad *mypad = {0};
	TH1 *hframe_1 = new TH2F("hframe_1","",100,0.01,0.99,100,0.01,55);

  for(int iCan=0;iCan<4;iCan++){
    mypad = pad2x2_1->GetPad(iCan+1);
    //gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    hframe_1->Draw();
    hframe_1->GetXaxis()->SetTitle("y");hframe_1->GetYaxis()->SetTitle("#sigma_{y}/y [%]");
		hframe_1->GetXaxis()->SetLabelFont(63);hframe_1->GetYaxis()->SetLabelFont(63);
    hframe_1->GetXaxis()->SetLabelSize(25);hframe_1->GetYaxis()->SetLabelSize(25);
    hframe_1->GetXaxis()->SetLabelOffset(0.01);hframe_1->GetYaxis()->SetLabelOffset(0.01);
    hframe_1->GetXaxis()->CenterTitle(1);hframe_1->GetYaxis()->CenterTitle(1);
    hframe_1->GetXaxis()->SetTitleSize(40);hframe_1->GetXaxis()->SetTitleOffset(2.5); 
    hframe_1->GetYaxis()->SetTitleSize(40);hframe_1->GetYaxis()->SetTitleOffset(3.0);

    if(iCan==0){ 
      gr1->Draw("PL same");gr2->Draw("PL same");gr3->Draw("PL same");
    }
    else if(iCan==1){
      gr1a->Draw("PL same");gr2a->Draw("PL same");gr3a->Draw("PL same");
    }
    else if(iCan==2){
      gr1b->Draw("PL same");gr2b->Draw("PL same");gr3b->Draw("PL same");
    }
    
  }

  // get the text pad
  pad2x2_1->GetPad(5);

  TLegend *leg1 = new TLegend(0.525,0.2,0.825,0.45);
  leg1->SetBorderSize(0);leg1->SetTextSize(0.03);
  leg1->AddEntry((TObject*)0, "18 GeV e^{-} on 275 GeV p", "");
  leg1->AddEntry(gr1,"Electron Method","p");
  leg1->AddEntry(gr2,"JB Method","p");
  leg1->AddEntry(gr3,"DA Method","p");
  leg1->Draw();

  TPaveText* pave1 = new TPaveText(0.1,0.8,0.35,0.9,"NDCNB");
  pave1->AddText("Q^{2} > 0.5 GeV^{2}");
	pave1->SetFillStyle(4000);
  pave1->SetTextFont(63);pave1->SetTextSize(30);
  pave1->SetTextColor(kBlack);
  pave1->Draw();

  TPaveText* pave2 = new TPaveText(0.5,0.8,0.8,0.9,"NDCNB");
  pave2->AddText("Q^{2} > 10 GeV^{2}");
	pave2->SetFillStyle(4000);
  pave2->SetTextFont(63);pave2->SetTextSize(30);
  pave2->SetTextColor(kBlack);
  pave2->Draw();

  TPaveText* pave3 = new TPaveText(0.1,0.35,0.35,0.55,"NDCNB");
  pave3->AddText("Q^{2} > 100 GeV^{2}");
	pave3->SetFillStyle(4000);
  pave3->SetTextFont(63);pave3->SetTextSize(30);
  pave3->SetTextColor(kBlack);
  pave3->Draw();

  TPaveText* pave4 = new TPaveText(0.90,0.09,0.94,0.12,"NDCNB");
  pave4->AddText("1");
	pave4->SetFillStyle(4000);
  pave4->SetTextFont(63);pave4->SetTextSize(25);
  pave4->Draw();

  TPaveText* pave5 = new TPaveText(0.1,0.09,0.14,0.12,"NDCNB");
  pave5->AddText("0");
	pave5->SetFillStyle(4000);
  pave5->SetTextFont(63);pave5->SetTextSize(25);
  pave5->Draw();

  TPaveText* pave6 = new TPaveText(0.5,0.09,0.54,0.12,"NDCNB");
  pave6->AddText("0");
	pave6->SetFillStyle(4000);
  pave6->SetTextFont(63);pave6->SetTextSize(25);
  pave6->Draw();

  TPaveText* pave7 = new TPaveText(0.1,0.1,0.12,0.14,"NDCNB");
  pave7->AddText("0");
	pave7->SetFillStyle(4000);
  pave7->SetTextFont(63);pave7->SetTextSize(25);
  pave7->Draw();

  TPaveText* pave8 = new TPaveText(0.1,0.5,0.12,0.53,"NDCNB");
  pave8->AddText("0");
	pave8->SetFillStyle(4000);
  pave8->SetTextFont(63);pave8->SetTextSize(25);
  pave8->Draw();

  //Plot 2
  PadMxN *pad2x2_2 = new PadMxN("c2",500,500,150,150,100,125,2,2);
  pad2x2_2->Draw();
	TH1 *hframe_2 = new TH2F("hframe_2","",100,0.01,179.99,100,0.01,55);

  for(int iCan=0;iCan<4;iCan++){
    mypad = pad2x2_2->GetPad(iCan+1);
    //gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    hframe_2->Draw();
    hframe_2->GetXaxis()->SetTitle("#theta_{e} [deg]");hframe_2->GetYaxis()->SetTitle("#sigma_{y}/y [%]");
		hframe_2->GetXaxis()->SetLabelFont(63);hframe_2->GetYaxis()->SetLabelFont(63);
    hframe_2->GetXaxis()->SetLabelSize(25);hframe_2->GetYaxis()->SetLabelSize(25);
    hframe_2->GetXaxis()->SetLabelOffset(0.01);hframe_2->GetYaxis()->SetLabelOffset(0.01);
    hframe_2->GetXaxis()->CenterTitle(1);hframe_2->GetYaxis()->CenterTitle(1);
    hframe_2->GetXaxis()->SetTitleSize(40);hframe_2->GetXaxis()->SetTitleOffset(2.5); 
    hframe_2->GetYaxis()->SetTitleSize(40);hframe_2->GetYaxis()->SetTitleOffset(3.0);

    if(iCan==0){ 
      gr4->Draw("PL same");gr5->Draw("PL same");gr6->Draw("PL same");
    }
    else if(iCan==1){
      gr4a->Draw("PL same");gr5a->Draw("PL same");gr6a->Draw("PL same");
    }
    else if(iCan==2){
      gr4b->Draw("PL same");gr5b->Draw("PL same");gr6b->Draw("PL same");
    }
    
  }

  // get the text pad
  pad2x2_2->GetPad(5);
  leg1->Draw();
  pave1->Draw();
  pave2->Draw();
  pave3->Draw();
  pave5->Draw();
  pave6->Draw();
  pave7->Draw();
  pave8->Draw();

  TPaveText* pave4a = new TPaveText(0.90,0.09,0.94,0.12,"NDCNB");
  pave4a->AddText("180");
	pave4a->SetFillStyle(4000);
  pave4a->SetTextFont(63);pave4a->SetTextSize(25);
  pave4a->Draw();

  //Plot 2
  PadMxN *pad2x2_3 = new PadMxN("c3",500,500,150,150,100,125,2,2);
  pad2x2_3->Draw();
  TH1 *hframe_3 = new TH2F("hframe_3","",100,0.01,179.99,100,0.01,55);

  for(int iCan=0;iCan<4;iCan++){
    mypad = pad2x2_3->GetPad(iCan+1);
    //gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    hframe_3->Draw();
    hframe_3->GetXaxis()->SetTitle("#theta_{h} [deg]");hframe_3->GetYaxis()->SetTitle("#sigma_{y}/y [%]");
		hframe_3->GetXaxis()->SetLabelFont(63);hframe_3->GetYaxis()->SetLabelFont(63);
    hframe_3->GetXaxis()->SetLabelSize(25);hframe_3->GetYaxis()->SetLabelSize(25);
    hframe_3->GetXaxis()->SetLabelOffset(0.01);hframe_3->GetYaxis()->SetLabelOffset(0.01);
    hframe_3->GetXaxis()->CenterTitle(1);hframe_3->GetYaxis()->CenterTitle(1);
    hframe_3->GetXaxis()->SetTitleSize(40);hframe_3->GetXaxis()->SetTitleOffset(2.5); 
    hframe_3->GetYaxis()->SetTitleSize(40);hframe_3->GetYaxis()->SetTitleOffset(3.0);

    if(iCan==0){ 
      gr7->Draw("PL same");gr8->Draw("PL same");gr9->Draw("PL same");
    }
    else if(iCan==1){
      gr7a->Draw("PL same");gr8a->Draw("PL same");gr9a->Draw("PL same");
    }
    else if(iCan==2){
      gr7b->Draw("PL same");gr8b->Draw("PL same");gr9b->Draw("PL same");
    }
    
  }

  // get the text pad
  pad2x2_3->GetPad(5);
  leg1->Draw();
  pave1->Draw();
  pave2->Draw();
  pave3->Draw();
  pave4a->Draw();
  pave5->Draw();
  pave6->Draw();
  pave7->Draw();
  pave8->Draw();

  //Make 2d plots
  for(int ibinx=1;ibinx<=prof10->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=prof10->GetNbinsY();ibiny++){
      
      hh1->SetBinContent(ibinx,ibiny,prof10->GetBinError(ibinx,ibiny) );
      hh2->SetBinContent(ibinx,ibiny,prof11->GetBinError(ibinx,ibiny) );
      hh3->SetBinContent(ibinx,ibiny,prof12->GetBinError(ibinx,ibiny) );

      cout<<hh1->GetXaxis()->GetBinCenter(ibinx)<<"   "<<hh1->GetYaxis()->GetBinCenter(ibiny)<<"   "
          <<hh1->GetBinContent(ibinx,ibiny)<<"   "<<hh2->GetBinContent(ibinx,ibiny)<<"   "<<hh3->GetBinContent(ibinx,ibiny)
          <<endl;

    }
  }

  PadMxN *pad3x1 = new PadMxN("c4",500,500,175,50,175,50,3,1);
  pad3x1->Draw();
  //TH1 *hframe_4 = new TH2F("hframe_4","",100,1e-5,1,100,1e-1,1e4);

  const Double_t min = 0.;
  const Double_t max = 55.;
  
  const Int_t nLevels = 999;
  Double_t levels[nLevels];
  gStyle->SetNumberContours(nLevels);

  for(int i = 1; i < nLevels; i++) {
    levels[i] = min + (max - min) / (nLevels - 1) * (i);
  }
  //levels[0] = 0.01;
  levels[0] = -1; //This also works as we want

  for(int iCan=0;iCan<3;iCan++){
    mypad = pad3x1->GetPad(iCan+1);
    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    /*
    hframe_4->Draw();
    hframe_4->GetXaxis()->SetTitle("x");hframe_4->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
	  hframe_4->GetXaxis()->SetLabelFont(63);hframe_4->GetYaxis()->SetLabelFont(63);
    hframe_4->GetXaxis()->SetLabelSize(25);hframe_4->GetYaxis()->SetLabelSize(25);
    hframe_4->GetXaxis()->SetLabelOffset(0.01);hframe_4->GetYaxis()->SetLabelOffset(0.01);
    hframe_4->GetXaxis()->CenterTitle(1);hframe_4->GetYaxis()->CenterTitle(1);
    hframe_4->GetXaxis()->SetTitleSize(40);hframe_4->GetXaxis()->SetTitleOffset(2.5); 
    hframe_4->GetYaxis()->SetTitleSize(40);hframe_4->GetYaxis()->SetTitleOffset(3.0);
    */

    if(iCan==0){
      hh1->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh1->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
      hh1->Draw("col");
    }
    else if(iCan==1){
      hh2->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh2->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
      hh2->Draw("col");
    }
    else if(iCan==2){
      hh3->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      hh3->GetZaxis()->SetRangeUser(min, max); // ... set the range ..
      hh3->Draw("colz");
    }

  }

  // get the text pad
  pad3x1->GetPad(4);

  TPaveText* pave9 = new TPaveText(0.45,0.0,0.55,0.125,"NDCNB");
  pave9->AddText("x");
	pave9->SetFillStyle(4000);
  pave9->SetTextFont(63);pave9->SetTextSize(45);
  pave9->Draw();

  TPaveText* pave10 = new TPaveText(0.0,0.45,0.1,0.55,"NDCNB");
  auto text10 = pave10->AddText("Q^{2} [GeV^{2}]");
	pave10->SetFillStyle(4000);
  pave10->SetTextFont(63);pave10->SetTextSize(45);
  text10->SetTextAngle(90.);
  pave10->Draw();

  TPaveText* pave11 = new TPaveText(0.92,0.70,1.0,0.85,"NDCNB");
  auto text11 = pave11->AddText("#sigma_{y}/y [%]");
	pave11->SetFillStyle(4000);
  pave11->SetTextFont(63);pave11->SetTextSize(35);
  text11->SetTextAngle(90.);
  pave11->Draw();

  TPaveText* pave12 = new TPaveText(0.075,0.65,0.25,0.85,"NDCNB");
  pave12->AddText("Electron Method");
	pave12->SetFillStyle(4000);
  pave12->SetTextFont(42);pave12->SetTextSize(0.055);
  pave12->Draw();

  TPaveText* pave13 = new TPaveText(0.335,0.65,0.485,0.85,"NDCNB");
  pave13->AddText("JB Method");
	pave13->SetFillStyle(4000);
  pave13->SetTextFont(42);pave13->SetTextSize(0.055);
  pave13->Draw();
  
  TPaveText* pave14 = new TPaveText(0.615,0.65,0.765,0.85,"NDCNB");
  pave14->AddText("DA Method");
	pave14->SetFillStyle(4000);
  pave14->SetTextFont(42);pave14->SetTextSize(0.055);
  pave14->Draw();

  TPaveText* pave15 = new TPaveText(0.11,0.80,0.26,0.9,"NDCNB");
  pave15->AddText("18 GeV e^{-} on 275 GeV p");
	pave15->SetFillStyle(4000);
  pave15->SetTextFont(42);pave15->SetTextSize(0.055);
  pave15->Draw();

  TPaveText* pave16 = new TPaveText(0.3605,0.1,0.3705,0.155,"NDCNB");
  pave16->AddText("10^{-5}");
  pave16->SetFillStyle(1001);pave16->SetFillColor(0);
	//pave16->SetFillStyle(4000);
  pave16->SetTextFont(63);pave16->SetTextSize(25);
  pave16->Draw();

  TPaveText* pave17 = new TPaveText(0.631,0.1,0.641,0.155,"NDCNB");
  pave17->AddText("10^{-5}");
  pave17->SetFillStyle(1001);pave17->SetFillColor(kWhite);
	//pave17->SetFillStyle(4000);
  pave17->SetTextFont(63);pave17->SetTextSize(25);
  pave17->Draw();

  //Print to File
  if(energy_set == 1){
    gROOT->ProcessLine("c1->Print(\"./plots/resolution_comparison_5_41.pdf[\");");
    gROOT->ProcessLine("c1->Print(\"./plots/resolution_comparison_5_41.pdf\");");
    gROOT->ProcessLine("c2->Print(\"./plots/resolution_comparison_5_41.pdf\");");
    gROOT->ProcessLine("c3->Print(\"./plots/resolution_comparison_5_41.pdf\");");
    gROOT->ProcessLine("c4->Print(\"./plots/resolution_comparison_5_41.pdf\");");
    gROOT->ProcessLine("c4->Print(\"./plots/resolution_comparison_5_41.pdf]\");");
  }
  if(energy_set == 2){
    gROOT->ProcessLine("c1->Print(\"./plots/resolution_comparison_18_275.pdf[\");");
    gROOT->ProcessLine("c1->Print(\"./plots/resolution_comparison_18_275.pdf\");");
    gROOT->ProcessLine("c2->Print(\"./plots/resolution_comparison_18_275.pdf\");");
    gROOT->ProcessLine("c3->Print(\"./plots/resolution_comparison_18_275.pdf\");");
    gROOT->ProcessLine("c4->Print(\"./plots/resolution_comparison_18_275.pdf\");");
    gROOT->ProcessLine("c4->Print(\"./plots/resolution_comparison_18_275.pdf]\");");
  }
}