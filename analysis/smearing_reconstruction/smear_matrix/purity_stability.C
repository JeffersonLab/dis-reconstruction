R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

void purity_stability(){

  int energy_set(0);

  cout << "Choose Which Beam Energies for ep" << endl;
  cout << "1) e = 5 GeV, p = 41 GeV" << endl;
  cin >> energy_set;

  //Constants                                                                                                        
  const double alpha(1./137.036);
  double s_cm(0);
  const double fbgev(1./(0.3894E12));
  const double Mp(0.9383);

  if(energy_set == 1) s_cm = 4.*5*41;

  //Cross Section Bins  
  //Q2 Binning                                                                        
  double Q2_min = 1E-1;
  double Q2_max = 1E4;
  const int nbins_Q2 = 25;
  double log_bw_Q2 = (log10(Q2_max) - log10(Q2_min))/(nbins_Q2); //Determine bin width                                                         
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
  double log_bw_x = (log10(x_max) - log10(x_min))/(nbins_x); //Determine bin width                                                             
  double log_x_div;
  double x_bins[nbins_x+1];
  for(int i=0;i<nbins_x+1;i++){
    log_x_div = log10(x_min) + (i*log_bw_x);
    x_bins[i] = pow(10,log_x_div);
  }

  //Bin Yields/Kinematics/Factors                                   
  double y[nbins_x][nbins_Q2];
  double Q2_center[nbins_Q2];
  double x_center[nbins_x];
  double Q2_width[nbins_Q2];
  double x_width[nbins_x];

  //We won't use these cuts when calculating the purity/stability right now.
  //But keep them in case we want to make use of them in the future.
  //----
  //Since this is a study on resolution effects, we may want to remove the edge bins
  //which are affected by generation ranges on angle and/or momentum acceptance losses.
  bool cut_ul[nbins_x][nbins_Q2];
  bool cut_lr[nbins_x][nbins_Q2];
  double y_max = 0.98;double y_min = 0.001;double W2_min = 10;
  double y_temp, W2_temp;

  //Index holding variables
  int gen1(0), gen2(0);
  for(int i=0;i<nbins_x;i++){ //x loop                                                                                             
    for(int j=0;j<nbins_Q2;j++){ //Q2 loop

      Q2_center[j] = (Q2_bins[j+1] + Q2_bins[j])/2.;
      x_center[i] = (x_bins[i+1] + x_bins[i])/2.;
      Q2_width[j] = (Q2_bins[j+1] - Q2_bins[j]);
      x_width[i] = (x_bins[i+1] - x_bins[i]);

      y[i][j] = Q2_center[j]/(x_center[i]*s_cm);

      //Cut Bins with top-left edge > ymax and bottom-right edge < ymin and > W2min                                                                  
      y_temp = (Q2_bins[j+1])/(x_bins[i]*s_cm); //Upper left                                                                                         
      if(y_temp>y_max) //Remove bin                                                                                                                  
        cut_ul[i][j] = 0;
      else
        cut_ul[i][j] = 1;

      y_temp = (Q2_bins[j])/(x_bins[i+1]*s_cm); //Lower Right                                                                                        
      W2_temp = Mp*Mp + Q2_bins[j]*( (1./x_bins[i+1]) - 1.);
   
      if(y_temp<y_min || W2_temp<W2_min) //Remove bin                                                                                                
        cut_lr[i][j] = 0;
      else
        cut_lr[i][j] = 1;
    }
  }

  //Histograms
  TH2D *h0 = new TH2D("h0","Generated",25,x_bins,25,Q2_bins);
  h0->GetXaxis()->SetTitle("x");h0->GetXaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h0->GetYaxis()->CenterTitle();

  //Electron Method using full energy
  TH2D *h1_1 = new TH2D("h1_1","Electron Method (using ECal energy) Purity",25,x_bins,25,Q2_bins);
  h1_1->GetXaxis()->SetTitle("x");h1_1->GetXaxis()->CenterTitle();
  h1_1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h1_1->GetYaxis()->CenterTitle();
  TH2D *h1_2 = new TH2D("h1_2","Electron Method (using ECal energy) Stability",25,x_bins,25,Q2_bins);
  h1_2->GetXaxis()->SetTitle("x");h1_2->GetXaxis()->CenterTitle();
  h1_2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h1_2->GetYaxis()->CenterTitle();
  TH2D *h1_3 = new TH2D("h1_3","Electron Method (using ECal energy) Reconstructed",25,x_bins,25,Q2_bins);
  h1_3->GetXaxis()->SetTitle("x");h1_3->GetXaxis()->CenterTitle();
  h1_3->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h1_3->GetYaxis()->CenterTitle();

  //Electron Method using momentum as energy
  TH2D *h2_1 = new TH2D("h2_1","Electron Method (using track momentum) Purity",25,x_bins,25,Q2_bins);
  h2_1->GetXaxis()->SetTitle("x");h2_1->GetXaxis()->CenterTitle();
  h2_1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h2_1->GetYaxis()->CenterTitle();
  TH2D *h2_2 = new TH2D("h2_2","Electron Method (using track momentum) Stability",25,x_bins,25,Q2_bins);
  h2_2->GetXaxis()->SetTitle("x");h2_2->GetXaxis()->CenterTitle();
  h2_2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h2_2->GetYaxis()->CenterTitle();
  TH2D *h2_3 = new TH2D("h2_3","Electron Method (using track momentum) Reconstructed",25,x_bins,25,Q2_bins);
  h2_3->GetXaxis()->SetTitle("x");h2_3->GetXaxis()->CenterTitle();
  h2_3->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h2_3->GetYaxis()->CenterTitle();

  //JB Method
  TH2D *h3_1 = new TH2D("h3_1","JB Method Purity",25,x_bins,25,Q2_bins);
  h3_1->GetXaxis()->SetTitle("x");h3_1->GetXaxis()->CenterTitle();
  h3_1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_1->GetYaxis()->CenterTitle();
  TH2D *h3_2 = new TH2D("h3_2","JB Method Stability",25,x_bins,25,Q2_bins);
  h3_2->GetXaxis()->SetTitle("x");h3_2->GetXaxis()->CenterTitle();
  h3_2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_2->GetYaxis()->CenterTitle();
  TH2D *h3_3 = new TH2D("h3_3","JB Method Reconstructed",25,x_bins,25,Q2_bins);
  h3_3->GetXaxis()->SetTitle("x");h3_3->GetXaxis()->CenterTitle();
  h3_3->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_3->GetYaxis()->CenterTitle();

  TH2D *h3_4 = new TH2D("h3_4","JB Method Purity (Summing over all particles)",25,x_bins,25,Q2_bins);
  h3_4->GetXaxis()->SetTitle("x");h3_4->GetXaxis()->CenterTitle();
  h3_4->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_4->GetYaxis()->CenterTitle();
  TH2D *h3_5 = new TH2D("h3_5","JB Method Stability (Summing over all particles)",25,x_bins,25,Q2_bins);
  h3_5->GetXaxis()->SetTitle("x");h3_5->GetXaxis()->CenterTitle();
  h3_5->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_5->GetYaxis()->CenterTitle();
  TH2D *h3_6 = new TH2D("h3_6","JB Method Reconstructed (Summing over all particles)",25,x_bins,25,Q2_bins);
  h3_6->GetXaxis()->SetTitle("x");h3_6->GetXaxis()->CenterTitle();
  h3_6->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_6->GetYaxis()->CenterTitle();

  TH2D *h3_7 = new TH2D("h3_7","Purity: using 4-momentum of X",25,x_bins,25,Q2_bins);
  h3_7->GetXaxis()->SetTitle("x");h3_7->GetXaxis()->CenterTitle();
  h3_7->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_7->GetYaxis()->CenterTitle();
  TH2D *h3_8 = new TH2D("h3_8","Stability: using 4-momentum of X",25,x_bins,25,Q2_bins);
  h3_8->GetXaxis()->SetTitle("x");h3_8->GetXaxis()->CenterTitle();
  h3_8->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_8->GetYaxis()->CenterTitle();
  TH2D *h3_9 = new TH2D("h3_9","Reconstructed: using 4-momentum of X",25,x_bins,25,Q2_bins);
  h3_9->GetXaxis()->SetTitle("x");h3_9->GetXaxis()->CenterTitle();
  h3_9->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h3_9->GetYaxis()->CenterTitle();

  //DA Method
  TH2D *h4_1 = new TH2D("h4_1","DA Method Purity",25,x_bins,25,Q2_bins);
  h4_1->GetXaxis()->SetTitle("x");h4_1->GetXaxis()->CenterTitle();
  h4_1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h4_1->GetYaxis()->CenterTitle();
  TH2D *h4_2 = new TH2D("h4_2","DA Method Stability",25,x_bins,25,Q2_bins);
  h4_2->GetXaxis()->SetTitle("x");h4_2->GetXaxis()->CenterTitle();
  h4_2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h4_2->GetYaxis()->CenterTitle();
  TH2D *h4_3 = new TH2D("h4_3","DA Method Reconstructed",25,x_bins,25,Q2_bins);
  h4_3->GetXaxis()->SetTitle("x");h4_3->GetXaxis()->CenterTitle();
  h4_3->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h4_3->GetYaxis()->CenterTitle();

  TH2D *h4_4 = new TH2D("h4_4","DA Method Purity (Summing over all particles)",25,x_bins,25,Q2_bins);
  h4_4->GetXaxis()->SetTitle("x");h4_4->GetXaxis()->CenterTitle();
  h4_4->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h4_4->GetYaxis()->CenterTitle();
  TH2D *h4_5 = new TH2D("h4_5","DA Method Stability (Summing over all particles)",25,x_bins,25,Q2_bins);
  h4_5->GetXaxis()->SetTitle("x");h4_5->GetXaxis()->CenterTitle();
  h4_5->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h4_5->GetYaxis()->CenterTitle();
  TH2D *h4_6 = new TH2D("h4_6","DA Method Reconstructed (Summing over all particles)",25,x_bins,25,Q2_bins);
  h4_6->GetXaxis()->SetTitle("x");h4_6->GetXaxis()->CenterTitle();
  h4_6->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h4_6->GetYaxis()->CenterTitle();

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

      //if(!particle_s and Status[j]==1) cout<<"Missing particle "<<j<<" in event "<<i<<" with id = "<< id[j] <<endl;
      
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

    gen1 = h0->Fill(x_e,Q2_e);

    //-------------Calculate *Smeared* invariants using 4-vectors------------------//
  
    if(detected_elec){
      //1.1) Smeared scattered election using Ecal energy the final energy
      Q2_e_nm_s = 4.*Ei_e*Ef_e_s*TMath::Cos(theta_e_s/2.)*TMath::Cos(theta_e_s/2.);
      y_e_nm_s = 1. - ( (Ef_e_s/(2.*Ei_e))*(1. - TMath::Cos(theta_e_s)) );
      x_e_nm_s = Q2_e_nm_s/(s_nm*y_e_nm_s);

      gen2 = h1_3->Fill(x_e_nm_s,Q2_e_nm_s);
      if(gen1 == gen2){
        h1_1->Fill(x_e_nm_s,Q2_e_nm_s);
        h1_2->Fill(x_e_nm_s,Q2_e_nm_s);
      }
      /*else{
        cout<<"Event = "<<i<<" x_gen = "<<x_e<<" Q2_gen = "<<Q2_e<<" xe_rec = "<<x_e_nm_s<<" Q2e_rec = "
            <<Q2_e_nm_s<<" gen bin = "<<gen1<<" rec bin = "<<gen2<<endl;
      }*/

      //1.2) Using Smeared scattered election using track momentum as the total momentum as final energy.  
      Ef_e_p_s = TMath::Sqrt( TMath::Power(pxf_e_s,2)+TMath::Power(pyf_e_s,2)+TMath::Power(pzf_e_s,2) + (0.511e-3)*(0.511e-3) ); 
      Q2_e_nm_p_s = 4.*Ei_e*Ef_e_p_s*TMath::Cos(theta_e_s/2.)*TMath::Cos(theta_e_s/2.);
      y_e_nm_p_s = 1. - ( (Ef_e_p_s/(2.*Ei_e))*(1. - TMath::Cos(theta_e_s)) );
      x_e_nm_p_s = Q2_e_nm_p_s/(s_nm*y_e_nm_p_s);

      gen2 = h2_3->Fill(x_e_nm_p_s,Q2_e_nm_p_s);
      if(gen1 == gen2){
        h2_1->Fill(x_e_nm_p_s,Q2_e_nm_p_s);
        h2_2->Fill(x_e_nm_p_s,Q2_e_nm_p_s);
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

        gen2 = h3_3->Fill(x_jb_jet_s,Q2_jb_jet_s);
        if(gen1 == gen2){
          h3_1->Fill(x_jb_jet_s,Q2_jb_jet_s);
          h3_2->Fill(x_jb_jet_s,Q2_jb_jet_s);
        }
        
        //3.1) Using Smeared DA Method
        if(detected_elec){
          Theta_h_nm_s = 2.* TMath::ATan( (Etot_jet_s - pztot_jet_s)/pttot_jet_s);    
          Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(theta_e_s/2.) )*
                   ( 1./(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
          y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
          x_da_s = Q2_da_s/(s_nm*y_da_s);  

          gen2 = h4_3->Fill(x_da_s,Q2_da_s);
          if(gen1 == gen2){
            h4_1->Fill(x_da_s,Q2_da_s);
            h4_2->Fill(x_da_s,Q2_da_s);
          }
        }

      }//jets.empty()

      //2.2) Using Smeared JB Method (Summing Over Hadrons)
      y_jb_sumh_s = (Etot_sumh_s - pztot_sumh_s)/(2*Ei_e);
      Q2_jb_sumh_s = (pttot_sumh_s*pttot_sumh_s)/(1. - y_jb_sumh_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s);

      gen2 = h3_6->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      if(gen1 == gen2){
        h3_4->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
        h3_5->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      }

      //2.3) Using 4-Vector of outgoing X
      eh_s_tot.SetPxPyPzE(pxtot_sumh_s,pytot_sumh_s,pztot_sumh_s,Etot_sumh_s);

      q_h_s= (eh_s_tot - pni);
      Q2_jb_sumh_s = -1.0*(q_h_s*q_h_s);
      y_jb_sumh_s = (pni*q_h_s)/(pni*ei);
      //x_jb_sumh_s = Q2_jb_sumh_s/(2*pni*q_h_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s); //Using Q2 = xsy, instead of 4-vector product Q2_jb_sumh_s/(2*pni*q_h_s); 

      gen2 = h3_9->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      if(gen1 == gen2){
        h3_7->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
        h3_8->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      }

      /*
      if(fabs(x_e-x_jb_sumh_s)/x_e>1e-3 || fabs(Q2_e-Q2_jb_sumh_s)/Q2_e>1e-3){
        cout<<"Event = "<<i<<" x_gen = "<<x_e<<" Q2_gen = "<<Q2_e<<" xe_rec = "<<x_e_nm_s<<" Q2e_rec = "
            <<Q2_e_nm_s<<" xh_rec = "<<x_jb_sumh_s<<" Q2h_rec = "<<Q2_jb_sumh_s<<endl;
      }
      */

      //3.2) Using Smeared DA Method (Summing Over Hadrons)
      if(detected_elec){
        Theta_h_nm_s = 2.* TMath::ATan( (Etot_sumh_s - pztot_sumh_s)/pttot_sumh_s);
        Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(theta_e_s/2.) )*
          ( 1./(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
        y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(theta_e_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
        x_da_s = Q2_da_s/(s_nm*y_da_s);

        gen2 = h4_6->Fill(x_da_s,Q2_da_s);
        if(gen1 == gen2){
          h4_4->Fill(x_da_s,Q2_da_s);
          h4_5->Fill(x_da_s,Q2_da_s);
        }
      }

    }//holdParticles.empty()

    //reset variables
    gen1 = 0;
    gen2 = 0;
    detected_elec = false;

  }//Finished Event Loop

  //Divide histograms
  h1_1->Divide(h1_3);
  h1_2->Divide(h0);

  h2_1->Divide(h2_3);
  h2_2->Divide(h0);

  h3_1->Divide(h3_3);
  h3_2->Divide(h0);

  h3_4->Divide(h3_6);
  h3_5->Divide(h0);

  h3_7->Divide(h3_9);
  h3_8->Divide(h0);

  h4_1->Divide(h4_3);
  h4_2->Divide(h0);

  h4_4->Divide(h4_6);
  h4_5->Divide(h0);

  //Make Latex
  TLatex *tex_energy = new TLatex();
  if(energy_set == 1){
    tex_energy->SetText(1e-4,1e3,"5 GeV e^{-} on 41 GeV p, #sqrt{s}=28.6 GeV");
    tex_energy->SetTextColor(kBlack);
    tex_energy->SetTextSize(0.035);
  }  

  //Make Plots
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1");
  c1->Divide(2,1);
  c1->cd(1);gPad->SetLogx();gPad->SetLogy();h1_1->Draw("colz");
  c1->cd(2);gPad->SetLogx();gPad->SetLogy();h1_2->Draw("colz");
  tex_energy->Draw();

  TCanvas *c2 = new TCanvas("c2");
  c2->Divide(2,1);
  c2->cd(1);gPad->SetLogx();gPad->SetLogy();h2_1->Draw("colz");
  c2->cd(2);gPad->SetLogx();gPad->SetLogy();h2_2->Draw("colz");
  tex_energy->Draw();

  TCanvas *c3 = new TCanvas("c3");
  c3->Divide(2,1);
  c3->cd(1);gPad->SetLogx();gPad->SetLogy();h3_1->Draw("colz");
  c3->cd(2);gPad->SetLogx();gPad->SetLogy();h3_2->Draw("colz");
  tex_energy->Draw();

  TCanvas *c4 = new TCanvas("c4");
  c4->Divide(2,1);
  c4->cd(1);gPad->SetLogx();gPad->SetLogy();h3_4->Draw("colz");
  c4->cd(2);gPad->SetLogx();gPad->SetLogy();h3_5->Draw("colz");
  tex_energy->Draw();

  TCanvas *c5 = new TCanvas("c5");
  c5->Divide(2,1);
  c5->cd(1);gPad->SetLogx();gPad->SetLogy();h3_7->Draw("colz");
  c5->cd(2);gPad->SetLogx();gPad->SetLogy();h3_8->Draw("colz");
  tex_energy->Draw();

  TCanvas *c6 = new TCanvas("c6");
  c6->Divide(2,1);
  c6->cd(1);gPad->SetLogx();gPad->SetLogy();h4_1->Draw("colz");
  c6->cd(2);gPad->SetLogx();gPad->SetLogy();h4_2->Draw("colz");
  tex_energy->Draw();

  TCanvas *c7 = new TCanvas("c7");
  c7->Divide(2,1);
  c7->cd(1);gPad->SetLogx();gPad->SetLogy();h4_4->Draw("colz");
  c7->cd(2);gPad->SetLogx();gPad->SetLogy();h4_5->Draw("colz");
  tex_energy->Draw();

  //Print to File
  if(energy_set == 1){
    c1->Print("./plots/purity_stability_5_41.pdf[");
    c1->Print("./plots/purity_stability_5_41.pdf");    
    c2->Print("./plots/purity_stability_5_41.pdf");
    c3->Print("./plots/purity_stability_5_41.pdf");
    c4->Print("./plots/purity_stability_5_41.pdf");
    c5->Print("./plots/purity_stability_5_41.pdf");
    c6->Print("./plots/purity_stability_5_41.pdf");
    c7->Print("./plots/purity_stability_5_41.pdf");
    c7->Print("./plots/purity_stability_5_41.pdf]");
  }
}