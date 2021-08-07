R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

void resolution_ybins(){

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

  //True y bins
  const int nybins = 5;
  double y_low[nybins] = {0.5,0.2,0.1,0.05,0.01};
  double y_hi[nybins]  = {0.8,0.5,0.2,0.1,0.05};

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

  //Electron Method using ecal energy (NOT USING HERE)
  //JB Method -- using jet (NOT USING HERE)
  //JB 4-Vector Method (NOT USING HERE)
  //DA Method -- using jet (NOT USING HERE)

  //Electron Method using momentum as energy
  TH1 *h1a[nybins]; //Q2 Resolution
  TH1 *h1b[nybins]; //y Resolution
  TH1 *h1c[nybins]; //x Resolution

  for(int ibin=0;ibin<nybins;ibin++){
    h1a[ibin] = new TH1D(Form("h1a[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-25,25);
    h1a[ibin]->SetLineColor(kBlue);h1a[ibin]->SetLineWidth(2);
    h1a[ibin]->GetXaxis()->SetTitle("Difference from true Q^{2} [%]");h1a[ibin]->GetXaxis()->CenterTitle();

    h1b[ibin] = new TH1D(Form("h1b[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h1b[ibin]->SetLineColor(kBlue);h1b[ibin]->SetLineWidth(2);
    h1b[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h1b[ibin]->GetXaxis()->CenterTitle();

    h1c[ibin] = new TH1D(Form("h1c[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h1c[ibin]->SetLineColor(kBlue);h1c[ibin]->SetLineWidth(2);
    h1c[ibin]->GetXaxis()->SetTitle("Difference from true x [%]");h1c[ibin]->GetXaxis()->CenterTitle();
  }
  
  //JB Method -- Summing over all particles
  TH1 *h2a[nybins]; //Q2 Resolution
  TH1 *h2b[nybins]; //y Resolution
  TH1 *h2c[nybins]; //x Resolution

  for(int ibin=0;ibin<nybins;ibin++){
    h2a[ibin] = new TH1D(Form("h2a[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h2a[ibin]->SetLineColor(kBlue);h2a[ibin]->SetLineWidth(2);
    h2a[ibin]->GetXaxis()->SetTitle("Difference from true Q^{2} [%]");h2a[ibin]->GetXaxis()->CenterTitle();

    h2b[ibin] = new TH1D(Form("h2b[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h2b[ibin]->SetLineColor(kBlue);h2b[ibin]->SetLineWidth(2);
    h2b[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h2b[ibin]->GetXaxis()->CenterTitle();

    h2c[ibin] = new TH1D(Form("h2c[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h2c[ibin]->SetLineColor(kBlue);h2c[ibin]->SetLineWidth(2);
    h2c[ibin]->GetXaxis()->SetTitle("Difference from true x [%]");h2c[ibin]->GetXaxis()->CenterTitle();
  }

  //DA Method -- Summing over all particles
  TH1 *h3a[nybins]; //Q2 Resolution
  TH1 *h3b[nybins]; //y Resolution
  TH1 *h3c[nybins]; //x Resolution

  for(int ibin=0;ibin<nybins;ibin++){
    h3a[ibin] = new TH1D(Form("h3a[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h3a[ibin]->SetLineColor(kBlue);h3a[ibin]->SetLineWidth(2);
    h3a[ibin]->GetXaxis()->SetTitle("Difference from true Q^{2} [%]");h3a[ibin]->GetXaxis()->CenterTitle();

    h3b[ibin] = new TH1D(Form("h3b[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h3b[ibin]->SetLineColor(kBlue);h3b[ibin]->SetLineWidth(2);
    h3b[ibin]->GetXaxis()->SetTitle("Difference from true y [%]");h3b[ibin]->GetXaxis()->CenterTitle();

    h3c[ibin] = new TH1D(Form("h3c[%d]",ibin),Form("%.2f < y_{true} < %.2f",y_low[ibin],y_hi[ibin]),400,-100,100);
    h3c[ibin]->SetLineColor(kBlue);h3c[ibin]->SetLineWidth(2);
    h3c[ibin]->GetXaxis()->SetTitle("Difference from true x [%]");h3c[ibin]->GetXaxis()->CenterTitle();
  }

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

  Int_t nevents = tree->GetEntries();
  
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
        //Use simple 'energy-flow' algorithm
	      if(j!=electronIndex && Status_s[j]==1){

          //handle neutral particles (use ecal for photons or hcal for kaons/nuetrons and assume perfect pid)
          if(id[j]==22 || id[j]==130 || id[j]==2112){

            holdE_h = particle_s->GetE();

            double neutral_mom = 0;
            int neutral_id = id[j];

            switch(neutral_id){
              case 22:
                neutral_mom = holdE_h;
                break;
              case 130:
                if(holdE_h>0.4976)
                  neutral_mom = sqrt(holdE_h*holdE_h - 0.4976*0.4976);
                else
                  neutral_mom = 0;
                break;
              case 2112:
                if(holdE_h>0.9396)
                  neutral_mom = sqrt(holdE_h*holdE_h - 0.9396*0.9396);
                else
                  neutral_mom = 0;
                break;
            }

            holdpx_h = neutral_mom * sin( particle_s->GetTheta() ) * cos( particle_s->GetPhi() );
            holdpy_h = neutral_mom * sin( particle_s->GetTheta() ) * sin( particle_s->GetPhi() );
            holdpz_h = neutral_mom * cos( particle_s->GetTheta() );

          }else{

            //Use tracker and assume perfect pid
            holdpx_h = particle_s->GetPx();
            holdpy_h = particle_s->GetPy();
            holdpz_h = particle_s->GetPz();
            holdE_h = sqrt(particle_s->GetP()*particle_s->GetP() + particle->GetM()*particle->GetM());
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
  
    if(detected_elec){
      //1.1) Smeared scattered election using Ecal energy the final energy
      Q2_e_nm_s = 4.*Ei_e*Ef_e_s*TMath::Cos(theta_e_s/2.)*TMath::Cos(theta_e_s/2.);
      y_e_nm_s = 1. - ( (Ef_e_s/(2.*Ei_e))*(1. - TMath::Cos(theta_e_s)) );
      x_e_nm_s = Q2_e_nm_s/(s_nm*y_e_nm_s);

      //1.2) Using Smeared scattered election using track momentum as the total momentum as final energy.  
      Ef_e_p_s = TMath::Sqrt( TMath::Power(pxf_e_s,2)+TMath::Power(pyf_e_s,2)+TMath::Power(pzf_e_s,2) + (0.511e-3)*(0.511e-3) ); 
      Q2_e_nm_p_s = 4.*Ei_e*Ef_e_p_s*TMath::Cos(theta_e_s/2.)*TMath::Cos(theta_e_s/2.);
      y_e_nm_p_s = 1. - ( (Ef_e_p_s/(2.*Ei_e))*(1. - TMath::Cos(theta_e_s)) );
      x_e_nm_p_s = Q2_e_nm_p_s/(s_nm*y_e_nm_p_s);

      for(int ibin=0;ibin<nybins;ibin++){
        if( y_e>y_low[ibin] && y_e<y_hi[ibin] ){
          h1a[ibin]->Fill(100.*(Q2_e-Q2_e_nm_p_s)/Q2_e );
          h1b[ibin]->Fill(100.*(y_e-y_e_nm_p_s)/y_e );
          h1c[ibin]->Fill(100.*(x_e-x_e_nm_p_s)/x_e );
        }
      }

    } //detected electron

    if(!(holdParticles.empty())){
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

      for(int ibin=0;ibin<nybins;ibin++){
        if( y_e>y_low[ibin] && y_e<y_hi[ibin] ){
          h2a[ibin]->Fill(100.*(Q2_e-Q2_jb_sumh_s)/Q2_e );
          h2b[ibin]->Fill(100.*(y_e-y_jb_sumh_s)/y_e );
          h2c[ibin]->Fill(100.*(x_e-x_jb_sumh_s)/x_e );
        }
      }

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

        for(int ibin=0;ibin<nybins;ibin++){
          if( y_e>y_low[ibin] && y_e<y_hi[ibin] ){
            h3a[ibin]->Fill(100.*(Q2_e-Q2_da_s)/Q2_e );
            h3b[ibin]->Fill(100.*(y_e-y_da_s)/y_e );
            h3c[ibin]->Fill(100.*(x_e-x_da_s)/x_e );
          }
        }

      }

    }//holdParticles.empty()

    //reset variables
    detected_elec = false;

  }//Finished Event Loop

  //Make Latex
  TPaveText* tex_energy = new TPaveText(0.1,0.7,0.9,0.9,"NDCNB");

  if(energy_set == 1){
    tex_energy->AddText("5 GeV e^{-} on 41 GeV p, #sqrt{s}=28.6 GeV");
	  tex_energy->SetFillStyle(4000);tex_energy->SetTextFont(63);tex_energy->SetTextSize(12);
  }

  if(energy_set == 2){
    tex_energy->AddText("18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
	  tex_energy->SetFillStyle(4000);tex_energy->SetTextFont(63);tex_energy->SetTextSize(12);
  }

  TPaveText *tex1_1 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex1_1->AddText("Electron Method (using track momentum)");
  tex1_1->SetFillStyle(4000);tex1_1->SetTextFont(63);tex1_1->SetTextSize(12);

  TPaveText *tex2_1 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex2_1->AddText("J.B. Method (summing all particles)");
  tex2_1->SetFillStyle(4000);tex2_1->SetTextFont(63);tex2_1->SetTextSize(12);

  TPaveText *tex3_1 = new TPaveText(0.1,0.55,0.9,0.75,"NDCNB");
  tex3_1->AddText("D.A. Method (summing all particles)");
  tex3_1->SetFillStyle(4000);tex3_1->SetTextFont(63);tex3_1->SetTextSize(12);

  //Make Plots

  //----------------------------------
  //Electron Method using momentum as energy
  //----------------------------------
  TCanvas *c1a = new TCanvas("c1a");
  c1a->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1a->cd(ibin+1);
    h1a[ibin]->Draw();
  }
  c1a->cd(6);
  tex_energy->Draw();tex1_1->Draw();

  TCanvas *c1b = new TCanvas("c1b");
  c1b->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1b->cd(ibin+1);
    h1b[ibin]->Draw();
  }
  c1b->cd(6);
  tex_energy->Draw();tex1_1->Draw();

  TCanvas *c1c = new TCanvas("c1c");
  c1c->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c1c->cd(ibin+1);
    h1c[ibin]->Draw();
  }
  c1c->cd(6);
  tex_energy->Draw();tex1_1->Draw();

  //----------------------------------
  //JB method summing over all particles
  //----------------------------------
  TCanvas *c2a = new TCanvas("c2a");
  c2a->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c2a->cd(ibin+1);
    h2a[ibin]->Draw();
  }
  c2a->cd(6);
  tex_energy->Draw();tex2_1->Draw();

  TCanvas *c2b = new TCanvas("c2b");
  c2b->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c2b->cd(ibin+1);
    h2b[ibin]->Draw();
  }
  c2b->cd(6);
  tex_energy->Draw();tex2_1->Draw();

  TCanvas *c2c = new TCanvas("c2c");
  c2c->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c2c->cd(ibin+1);
    h2c[ibin]->Draw();
  }
  c2c->cd(6);
  tex_energy->Draw();tex2_1->Draw();

  //----------------------------------
  //DA method summing over all particles
  //----------------------------------
  TCanvas *c3a = new TCanvas("c3a");
  c3a->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c3a->cd(ibin+1);
    h3a[ibin]->Draw();
  }
  c3a->cd(6);
  tex_energy->Draw();tex3_1->Draw();

  TCanvas *c3b = new TCanvas("c3b");
  c3b->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c3b->cd(ibin+1);
    h3b[ibin]->Draw();
  }
  c3b->cd(6);
  tex_energy->Draw();tex3_1->Draw();

  TCanvas *c3c = new TCanvas("c3c");
  c3c->Divide(3,2);
  for(int ibin=0;ibin<nybins;ibin++){
    c3c->cd(ibin+1);
    h3c[ibin]->Draw();
  }
  c3c->cd(6);
  tex_energy->Draw();tex3_1->Draw();

  //Print offsets and resolutions (both in %) to screen
  //---
  cout<<endl<<"For Electron method:"<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("Q2 Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h1a[ibin]->GetMean()),h1a[ibin]->GetMeanError(),fabs(h1a[ibin]->GetRMS()),h1a[ibin]->GetRMSError());
  }
  cout<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("y Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h1b[ibin]->GetMean()),h1b[ibin]->GetMeanError(),fabs(h1b[ibin]->GetRMS()),h1b[ibin]->GetRMSError());
  }
  cout<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("x Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h1c[ibin]->GetMean()),h1c[ibin]->GetMeanError(),fabs(h1c[ibin]->GetRMS()),h1c[ibin]->GetRMSError());
  }
  cout<<endl<<"----------------------"<<endl;
  //---
  cout<<endl<<"For JB method:"<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("Q2 Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h2a[ibin]->GetMean()),h2a[ibin]->GetMeanError(),fabs(h2a[ibin]->GetRMS()),h2a[ibin]->GetRMSError());
  }
  cout<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("y Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h2b[ibin]->GetMean()),h2b[ibin]->GetMeanError(),fabs(h2b[ibin]->GetRMS()),h2b[ibin]->GetRMSError());
  }
  cout<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("x Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h2c[ibin]->GetMean()),h2c[ibin]->GetMeanError(),fabs(h2c[ibin]->GetRMS()),h2c[ibin]->GetRMSError());
  }
  cout<<endl<<"----------------------"<<endl;
  //---
  cout<<endl<<"For DA method:"<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("Q2 Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h3a[ibin]->GetMean()),h3a[ibin]->GetMeanError(),fabs(h3a[ibin]->GetRMS()),h3a[ibin]->GetRMSError());
  }
  cout<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("y Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h3b[ibin]->GetMean()),h3b[ibin]->GetMeanError(),fabs(h3b[ibin]->GetRMS()),h3b[ibin]->GetRMSError());
  }
  cout<<endl;
  for(int ibin=0;ibin<nybins;ibin++){
    printf("x Offset and Resolution (%.2f < y_{true} < %.2f): Offset = %.3f +- %.4f %%, Resolution = %.3f +- %.4f %% \n",
            y_low[ibin],y_hi[ibin],fabs(h3c[ibin]->GetMean()),h3c[ibin]->GetMeanError(),fabs(h3c[ibin]->GetRMS()),h3c[ibin]->GetRMSError());
  }
  cout<<endl<<"----------------------"<<endl;

  //Print to File
  if(energy_set == 1){
    c1a->Print("./plots/resolution_ybins_5_41.pdf[");
    c1a->Print("./plots/resolution_ybins_5_41.pdf");
    c1b->Print("./plots/resolution_ybins_5_41.pdf");
    c1c->Print("./plots/resolution_ybins_5_41.pdf");
    c2a->Print("./plots/resolution_ybins_5_41.pdf");
    c2b->Print("./plots/resolution_ybins_5_41.pdf");
    c2c->Print("./plots/resolution_ybins_5_41.pdf");
    c3a->Print("./plots/resolution_ybins_5_41.pdf");
    c3b->Print("./plots/resolution_ybins_5_41.pdf");
    c3c->Print("./plots/resolution_ybins_5_41.pdf");
    c3c->Print("./plots/resolution_ybins_5_41.pdf]");
  }
  if(energy_set == 2){
    c1a->Print("./plots/resolution_ybins_18_275.pdf[");
    c1a->Print("./plots/resolution_ybins_18_275.pdf");
    c1b->Print("./plots/resolution_ybins_18_275.pdf");
    c1c->Print("./plots/resolution_ybins_18_275.pdf");
    c2a->Print("./plots/resolution_ybins_18_275.pdf");
    c2b->Print("./plots/resolution_ybins_18_275.pdf");
    c2c->Print("./plots/resolution_ybins_18_275.pdf");
    c3a->Print("./plots/resolution_ybins_18_275.pdf");
    c3b->Print("./plots/resolution_ybins_18_275.pdf");
    c3c->Print("./plots/resolution_ybins_18_275.pdf");
    c3c->Print("./plots/resolution_ybins_18_275.pdf]");
  }
}