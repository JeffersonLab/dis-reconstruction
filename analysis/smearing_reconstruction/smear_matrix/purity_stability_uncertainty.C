R__LOAD_LIBRARY(libeicsmear);
R__LOAD_LIBRARY(libfastjet);
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "includes/energy_resolution_optimization.h"

void purity_stability_uncertainty(){

  int energy_set(0);

  cout << "Choose Which Beam Energies for ep" << endl;
  cout << "1) e = 5 GeV, ep = 41 GeV" << endl;
  cout << "2) e = 18 GeV, ep = 275 GeV" << endl;
  cin >> energy_set;

  //Constants                                                                                                                                        
  const double alpha(1./137.036);
  double s_cm(0);
  const double fbgev(1./(0.3894E12));
  const double Mp(0.9383);

  if(energy_set == 1){s_cm = 4.*5*41;}
  if(energy_set == 2){s_cm = 4.*18*275;}

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

  double Q2_center[nbins_Q2];
  double x_center[nbins_x];
  double error[10][2][nbins_x][nbins_Q2];
 
  for(int k=0;k<10;k++){
    for(int l=0;l<2;l++){
      for(int i=0;i<nbins_x;i++){    //x loop                                                                                                               
	for(int j=0;j<nbins_Q2;j++){ //Q2 loop                                                                                                              
	  Q2_center[j] = (Q2_bins[j+1] + Q2_bins[j])/2.;
	  x_center[i] = (x_bins[i+1] + x_bins[i])/2.;
	  error[k][l][i][j] = 0;
	}
      }
    }
  }

  //Index holding variables                                                                                                                           
  int gen1(0), gen2(0);

  //Histograms
  TH2D *h[10][3];

  //Set the bins for the histograms
  for(int i=0;i<10;i++){
    for(int j=0;j<3;j++){
      h[i][j] = new TH2D("","",25,x_bins,25,Q2_bins);
    }
  }

  /*
    LEGEND for h[][] and error[][][][]
    First Index is Method, Second Index is Type
    
    First Index:
    0 -> Electron Method (Using optimized calculation)
    1 -> Electron Method (Using Track Momentum)
    2 -> JB Method (Using Jets Scattering off Nucleon)
    3 -> JB Method (Using Jets Scattering off Nucleus with New Term)
    4 -> JB Method (Summing Over Hadrons Scattering off Nucleon)
    5 -> JB Method (Summing Over Hadrons Scattering off Nucleus with New Term)
    6 -> DA Method (Using Jets Scattering off Nucleon)
    7 -> DA Method (Using Jets Scattering off Nucleus with New Term)  
    8 -> DA Method (Summing Over Hadrons Scattering off Nucleon)
    9 -> DA Method (Summing Over Hadrons Scattering off Nucleus with New Term)
  
    Second Index:
    0 -> Purity
    1 -> Stability
    2 -> Reconstructed

    Third Index and Fourth Index are the x and Q^2 bins
  */

  TH2D *h0 = new TH2D("h0","Generated",25,x_bins,25,Q2_bins);
  h0->GetXaxis()->SetTitle("x");h0->GetXaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h0->GetYaxis()->CenterTitle();

  //--------------------------------//
  //   Analyse BeAGLE Simulation    //
  //--------------------------------//

  //Load ROOT Files
  erhic::EventBeagle *event(NULL);    //Event Class
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


  tree->AddFriend(tree_s, "Smeared"); //Smeared Tree
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
  TLorentzVector Ai;    //Initial Nucleus 4-Vector
  TLorentzVector q_e; //q 4-vector using scattered electron
  Double_t Q2_e(0),y_e(0),x_e(0); //Scattered electron (no assumptions)
  Double_t s_nm(0); //CM Energy Squared (massless e,p)
  Int_t electronIndex;

  //Smeared Variables
  Int_t Status_s[500];
  Int_t id_s[500];
  Double_t Ef_e_s(0), pxf_e_s(0), pyf_e_s(0), pzf_e_s(0), ptf_e_s(0) ; //Scattered Electron
  Double_t Ef_e_p_s(0); //Scattered electron, using total momentum as final energy
  Double_t pztot_jet_s(0),pttot_jet_s(0),Etot_jet_s(0); // Total added hadronic variables via. FastJet
  Double_t pxtot_sumh_s(0),pytot_sumh_s(0),pztot_sumh_s(0),pttot_sumh_s(0),Etot_sumh_s(0); // Total added hadronic variables by summing over hadrons without Remnant
  Double_t holdpx_h(0),holdpy_h(0),holdpz_h(0),holdE_h(0),holdP_h(0); //Hold variables for summing over hadrons
  TLorentzVector ef_s;
  TLorentzVector q_e_s; //q 4-vector using scattered electron
  TLorentzVector X4_s_tot; //4-vector for sum over hadrons
  TLorentzVector q_h_s; //q 4-vector for sum over hadrons
  Double_t Theta_e_nm_s(0),Q2_e_nm_s(0),y_e_nm_s(0),x_e_nm_s(0); //Scattered electron (massless e,e',p)
  Double_t Q2_e_nm_p_s(0),y_e_nm_p_s(0),x_e_nm_p_s(0); //Scattered electron using total momentum as final energy(massless e,e',p)
  Double_t Q2_jb_jet_s(0),y_jb_jet_s(0),x_jb_jet_s(0); // JB Method (Using Jets)
  Double_t Q2_jb_sumh_s(0),y_jb_sumh_s(0),x_jb_sumh_s(0); // JB Method (Summing Over Hadrons) 
  Double_t Theta_h_nm_s(0), Q2_da_s(0), y_da_s(0), x_da_s(0); //DA Method
  Double_t mass[500];
  bool detected_elec(false);

  //FastJet Variables
  vector<fastjet::PseudoJet> holdParticles;
  valarray<double> fourvec(4);
  Double_t R = 1.0;     //For now, going to use the paramters and cuts from help-script
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm,R,fastjet::E_scheme,fastjet::Best); 
  Double_t etaMax = 4.0;
  Double_t etaMin = -etaMax;

  //Not going to apply pt cuts to jet. Don't see a reason to yet
  //auto selectJetEta = fastjet::SelectorEtaRange(etaMin+R,etaMax+R);
 
  Int_t nevents = tree->GetEntries();

  cout << nevents << endl;
  
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

    //Loop to get unsmeard particle information 
    for(int j=0;j<nParticles;j++){
      
      particle = event->GetTrack(j);
      Status[j] = (Int_t) particle->GetStatus();
      id[j] = (Int_t) particle->Id();
      mass[j] = (Double_t) particle->GetM();
      orig[j] = (Int_t) particle->GetParentIndex();  

      // Get Scattered Electron
      if(Status[j]==1 && id[j]==11 && orig[j]==3){
	      Ef_e = particle->GetE();
	      pxf_e = particle->GetPx();
	      pyf_e = particle->GetPy();
	      pzf_e = particle->GetPz();
	      ef.SetPxPyPzE(pxf_e,pyf_e,pzf_e,Ef_e);
	      electronIndex = j;	  
      }	

      // Get initial state election
      if(j==0 && Status[j]==21 && id[j]==11){
	      Ei_e = particle->GetE();
	      pxi_e = particle->GetPx();
	      pyi_e = particle->GetPy();
	      pzi_e = particle->GetPz();
	      ei.SetPxPyPzE(pxi_e,pyi_e,pzi_e,Ei_e);
      }
	 
      // Get initial state proton
      if( j==1 && id[j]==2212 && Status[j]==21 ){
	        E_pn = particle->GetE(); 
	        px_pn = particle->GetPx();
	        py_pn = particle->GetPy();
	        pz_pn = particle->GetPz();
	        Mpn = particle->GetM();
	        pni.SetPxPyPzE(px_pn,py_pn,pz_pn,E_pn);
	    }  
    }

    // Loop to get rest of particle information
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
	
	      // Get Smeared-Scattered Electron
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

            Theta_e_nm_s = particle_s->GetTheta();
	          ef_s.SetPxPyPzE(pxf_e_s,pyf_e_s,pzf_e_s,Ef_e_s);
            detected_elec = true;

          } else if(sigma_best == 1){
            pxf_e_s = particle_s->GetPx();
            pyf_e_s = particle_s->GetPy();
            pzf_e_s = particle_s->GetPz();
            Ef_e_s = sqrt( pxf_e_s*pxf_e_s + pyf_e_s* pyf_e_s + pzf_e_s*pzf_e_s + mass[j]*mass[j] );

            Theta_e_nm_s = particle_s->GetTheta();
	          ef_s.SetPxPyPzE(pxf_e_s,pyf_e_s,pzf_e_s,Ef_e_s);
            detected_elec = true;

          }else{
            detected_elec = false;
          }
	      }

	// Put hadrons into PseudoJet object
	if(j != electronIndex && Status_s[j]==1){     
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

	  fourvec[0] = holdpx_h;
	  fourvec[1] = holdpy_h;
	  fourvec[2] = holdpz_h;
	  fourvec[3] = holdE_h;
	    
	  fastjet::PseudoJet finalStateParticle(fourvec);
	  finalStateParticle.set_user_index(j);
	  holdParticles.push_back(finalStateParticle);

	}
      }        
    }//Finished Loop Over Particles

    //-------------Calculate invariants using 4-vectors------------------//
    //Since we are comparing to smeared tree here, only calculate...
    //...variables using scattered electron method.
    
    //1) Using Scattered Electron
    s_nm = 4*Ei_e*E_pn;
    q_e = (ei - ef);
    Q2_e = -1.0*(q_e*q_e);
    y_e = (pni*q_e)/(pni*ei);
    x_e =  Q2_e/(s_nm*y_e); // Q2_e/(2*pni*q_e); 
        
    gen1 = h0->Fill(x_e,Q2_e);
    //-------------Calculate *Smeared* invariants using 4-vectors------------------//
    
    if(detected_elec){
      //1.1) Using Smeared scattered election with optimized calculation  
      Q2_e_nm_s = 4.*Ei_e*Ef_e_s*TMath::Cos(Theta_e_nm_s/2.)*TMath::Cos(Theta_e_nm_s/2.);
      y_e_nm_s = 1. - ( (Ef_e_s/(2.*Ei_e))*(1. - TMath::Cos(Theta_e_nm_s)) );
      x_e_nm_s = Q2_e_nm_s/(s_nm*y_e_nm_s);
 
      gen2 = h[0][2]->Fill(x_e_nm_s,Q2_e_nm_s);
      if(gen1 == gen2){
	      h[0][0]->Fill(x_e_nm_s,Q2_e_nm_s);
	      h[0][1]->Fill(x_e_nm_s,Q2_e_nm_s);
      }
      
      //Using smeared theta value obtained from track, NOT DERIVED
      //1.2) Using Smeared scattered election (w/ assumptions) using the track momentum as final energy. 
      Ef_e_p_s = TMath::Sqrt( TMath::Power(pxf_e_s,2)+TMath::Power(pyf_e_s,2)+TMath::Power(pzf_e_s,2) + (0.511e-3)*(0.511e-3)); 
      Q2_e_nm_p_s = 4.*Ei_e*Ef_e_p_s*TMath::Cos(Theta_e_nm_s/2.)*TMath::Cos(Theta_e_nm_s/2.);
      y_e_nm_p_s = 1. - ( (Ef_e_p_s/(2.*Ei_e))*(1. - TMath::Cos(Theta_e_nm_s)) );
      x_e_nm_p_s = Q2_e_nm_p_s/(s_nm*y_e_nm_p_s);
      
      gen2 = h[1][2]->Fill(x_e_nm_p_s,Q2_e_nm_p_s);
      if(gen1 == gen2){
	      h[1][0]->Fill(x_e_nm_p_s,Q2_e_nm_p_s);
	      h[1][1]->Fill(x_e_nm_p_s,Q2_e_nm_p_s);
      }
    } //detected_elec
    
    if(!(holdParticles.empty())){
      pttot_sumh_s = TMath::Sqrt( TMath::Power(pxtot_sumh_s,2)+TMath::Power(pytot_sumh_s,2) ); 
  
      fastjet::ClusterSequence cs(holdParticles, jetDef);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

      if(!(jets.empty())){
	//For now, I am going to calculate the DIS Kinematics via. JB Method by using the highest pt jet
	pztot_jet_s = jets[0].pz();
	pttot_jet_s = jets[0].pt();
	Etot_jet_s = jets[0].E();
      
	//2.1) Using Smeared JB Method (Using Jets Scattering off Nucleon)                                                                              
	y_jb_jet_s = (Etot_jet_s - pztot_jet_s)/(2*Ei_e);
	Q2_jb_jet_s = (pttot_jet_s*pttot_jet_s)/(1. - y_jb_jet_s);
	x_jb_jet_s = Q2_jb_jet_s/(s_nm*y_jb_jet_s);  
	
	gen2 = h[2][2]->Fill(x_jb_jet_s,Q2_jb_jet_s);
	if(gen1 == gen2){
	  h[2][0]->Fill(x_jb_jet_s,Q2_jb_jet_s);
	  h[2][1]->Fill(x_jb_jet_s,Q2_jb_jet_s);
	}
	
	//2.2) Using Smeared JB Method (Using Jets Scattering off Nucleus with New Term)                                                                                                   
	y_jb_jet_s = (Etot_jet_s - pztot_jet_s)/(2*Ei_e) - (pow(Mp,2.0)/(4.0*Ei_e*E_pn));
	Q2_jb_jet_s = (pttot_jet_s*pttot_jet_s)/(1. - y_jb_jet_s);
	x_jb_jet_s = Q2_jb_jet_s/(s_nm*y_jb_jet_s);  
	
	gen2 = h[3][2]->Fill(x_jb_jet_s,Q2_jb_jet_s);
	if(gen1 == gen2){
	  h[3][0]->Fill(x_jb_jet_s,Q2_jb_jet_s);
	  h[3][1]->Fill(x_jb_jet_s,Q2_jb_jet_s);
	}

	if(detected_elec){
	  //3.1) Using Smeared DA Method Scattering off Nucleon
	  //For the angle, we can use jets or summing over all hadrons, 
	  //I will first use jets then use the same varialbes for summing over hadrons
	  Theta_h_nm_s = 2.* TMath::ATan( (Etot_jet_s - pztot_jet_s)/pttot_jet_s);    	
	  Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(Theta_e_nm_s/2.) )*
	    ( 1./(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
	  y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
	  x_da_s = Q2_da_s/(s_nm*y_da_s);      
	  
	  gen2 = h[6][2]->Fill(x_da_s,Q2_da_s);
	  if(gen1 == gen2){
	    h[6][0]->Fill(x_da_s,Q2_da_s);
	    h[6][1]->Fill(x_da_s,Q2_da_s);
	  }	
	
	  //3.2) Using Smeared DA Method Scattering off Nucleus with New Term
	  //For the angle, we can use jets or summing over all hadrons, 
	  //I will first use jets then use the same varialbes for summing over hadrons
	  Theta_h_nm_s = 2.* TMath::ATan( (Etot_jet_s - pztot_jet_s)/pttot_jet_s);    	
	  Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(Theta_e_nm_s/2.) )*
	    ( 1./(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.)-((E_pn - pz_pn)/ptf_e_s) ) );
	  y_da_s = (TMath::Tan(Theta_h_nm_s/2.)-((E_pn - pz_pn)/ptf_e_s))/(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.) -((E_pn - pz_pn)/ptf_e_s) );
	  x_da_s = Q2_da_s/(s_nm*y_da_s);      
	  
	  gen2 = h[7][2]->Fill(x_da_s,Q2_da_s);
	  if(gen1 == gen2){
	    h[7][0]->Fill(x_da_s,Q2_da_s);
	    h[7][1]->Fill(x_da_s,Q2_da_s);
	  }	

	}
      } //jets.empty()

      //2.3) Using Smeared JB Method (Summing Over Hadrons Scattering Off Nucleon)         
      y_jb_sumh_s = (Etot_sumh_s - pztot_sumh_s)/(2*Ei_e);
      Q2_jb_sumh_s = (pttot_sumh_s*pttot_sumh_s)/(1. - y_jb_sumh_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s);  

      gen2 = h[4][2]->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      if(gen1 == gen2){
	h[4][0]->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
	h[4][1]->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      }
    
      //2.4) Using Smeared JB Method (Summing Over Hadrons Scattering Off Nucleus)         
      y_jb_sumh_s = (Etot_sumh_s - pztot_sumh_s)/(2*Ei_e) - (pow(Mp,2.0)/(4.0*Ei_e*E_pn)) ;
      Q2_jb_sumh_s = (pttot_sumh_s*pttot_sumh_s)/(1. - y_jb_sumh_s);
      x_jb_sumh_s = Q2_jb_sumh_s/(s_nm*y_jb_sumh_s);  

      gen2 = h[5][2]->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      if(gen1 == gen2){
	h[5][0]->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
	h[5][1]->Fill(x_jb_sumh_s,Q2_jb_sumh_s);
      }
      
      if(detected_elec){
	//3.3) Using Smeared DA Method Scattering off Nucleon
	//Summing over Hadrons angle
	Theta_h_nm_s = 2.* TMath::ATan( (Etot_sumh_s - pztot_sumh_s)/pttot_sumh_s);    
	Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(Theta_e_nm_s/2.) )*
	  ( 1./(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.)) );
	y_da_s = TMath::Tan(Theta_h_nm_s/2.)/(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.));
	x_da_s = Q2_da_s/(s_nm*y_da_s);      
	
	gen2 = h[8][2]->Fill(x_da_s,Q2_da_s);
	if(gen1 == gen2){
	  h[8][0]->Fill(x_da_s,Q2_da_s);
	  h[8][1]->Fill(x_da_s,Q2_da_s);
	}
      
	//3.4) Using Smeared DA Method Scattering off Nucleus with New Term
	//Summing over Hadrons angle
	Theta_h_nm_s = 2.* TMath::ATan( (Etot_sumh_s - pztot_sumh_s)/pttot_sumh_s);    
	Q2_da_s = 4.*Ei_e*Ei_e*( 1./TMath::Tan(Theta_e_nm_s/2.) )*
	  ( 1./(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.) -((E_pn - pz_pn)/ptf_e_s) ) );
	y_da_s = (TMath::Tan(Theta_h_nm_s/2.) - ((E_pn - pz_pn)/ptf_e_s))/(TMath::Tan(Theta_e_nm_s/2.)+TMath::Tan(Theta_h_nm_s/2.) - ((E_pn - pz_pn)/ptf_e_s)) ;
	x_da_s = Q2_da_s/(s_nm*y_da_s);      
	
	gen2 = h[9][2]->Fill(x_da_s,Q2_da_s);
	if(gen1 == gen2){
	  h[9][0]->Fill(x_da_s,Q2_da_s);
	  h[9][1]->Fill(x_da_s,Q2_da_s);
	}
      

      }
    } //holdParticles.empty()
    
    //reset variables   
    gen1 = 0;
    gen2 = 0;
    detected_elec = false;
    
  }//Finished Event Loop

  //At this point we have all the information we need to make the uncertainty plots. I kept the above because we need to know k/n for the uncertainties

  //Calculate Errors via. Binomial Distribution modification in paper
  for(int k=0;k<10;k++){
    for(int l=0;l<2;l++){
      for(int i=0;i<nbins_x;i++){    //x loop                                                                                                               
	for(int j=0;j<nbins_Q2;j++){ //Q2 loop           

	  //There are underflow and overflow bins, bin numbers range from 0-nbins+1, so a given plot has a total of 27 bins. We need to not plots the 0th and 26th bin. Thus the reason for i+1 and j+1
	  
	  if(l==0){                                                                   
	    error[k][l][i][j] = TMath::Sqrt( (h[k][l]->GetBinContent(i+1,j+1) + 1)*(h[k][l]->GetBinContent(i+1,j+1) + 2) / ( (h[k][2]->GetBinContent(i+1,j+1) + 2)*(h[k][2]->GetBinContent(i+1,j+1) + 3) ) - (h[k][l]->GetBinContent(i+1,j+1) + 1)*(h[k][l]->GetBinContent(i+1,j+1) + 1) / ( (h[k][2]->GetBinContent(i+1,j+1) + 2)*(h[k][2]->GetBinContent(i+1,j+1) + 2) ));
	  }

	  if(l==1){                                                                                                   
	    error[k][l][i][j] = TMath::Sqrt( (h[k][l]->GetBinContent(i+1,j+1) + 1)*(h[k][l]->GetBinContent(i+1,j+1) + 2) / ( (h0->GetBinContent(i+1,j+1) + 2)*(h0->GetBinContent(i+1,j+1) + 3) ) - (h[k][l]->GetBinContent(i+1,j+1) + 1)*(h[k][l]->GetBinContent(i+1,j+1) + 1) / ( (h0->GetBinContent(i+1,j+1) + 2)*(h0->GetBinContent(i+1,j+1) + 2) ));
	  }	    									

	}
      }
    }
  }

  //Divide histograms to get purity-stability plots
  for(int i=0;i<10;i++){
    for(int j=0;j<2;j++){
    
      if(j==0){
	h[i][j]->Divide(h[i][2]);
      }

      if(j==1){
	h[i][j]->Divide(h0);
      }

    }
  }

  //Create and fill TGraphs                                                                                                                             
  TGraphErrors *gr[10][2][nbins_Q2]; //Indicies read: methods,purity(0)-stability(1),the Q^2 bins being plotted

  for(int j=0;j<10;j++){
    for(int k=0;k<2;k++){
      for(int l=0;l<nbins_Q2;l++){
	gr[j][k][l] = new TGraphErrors();
	gr[j][k][l]->SetName(Form("gr%d%d%d",j,k,l));
	gr[j][k][l]->SetLineWidth(2);
	gr[j][k][l]->SetLineColor(kRed);
	gr[j][k][l]->SetMarkerStyle(20);
	gr[j][k][l]->SetMarkerColor(kRed);
      }
    }
  }
  
  int counter(0);

  for(int k=0;k<10;k++){
    for(int l=0;l<2;l++){
      for(int j=0;j<nbins_Q2;j++){
	for(int i=0;i<nbins_x;i++){

	  //There are underflow and overflow bins that we dont need to plot. bin numbers range from 0-nbins+1, so a given plot has a total of 27 bins. We need to not plots the 0th and 26th bin.

	    gr[k][l][j]->SetPoint(counter,x_center[i],h[k][l]->GetBinContent(i+1,j+1));
	    gr[k][l][j]->SetPointError(counter,0,error[k][l][i][j]);
	    counter++;
	  
	}
	counter=0;
      }
    }
  }
  
  TLatex *tex_energy = new TLatex();
  
  if(energy_set == 1){
  tex_energy->SetText(0.2*1E-2,0.8,"5 GeV e^{-} on 41 GeV p, #sqrt{s}=28.6 GeV");
  tex_energy->SetTextColor(kBlack);
  tex_energy->SetTextSize(0.04);
  }  

  if(energy_set == 2){
  tex_energy->SetText(0.2*1E-2,0.8,"18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
  tex_energy->SetTextColor(kBlack);
  tex_energy->SetTextSize(0.04);
  }  

  TLatex *tex1[nbins_Q2];
  for(int i=0;i<nbins_Q2;i++){

    if(i >= 12){
      tex1[i] = new TLatex(1.0E-3,0.9,Form(" %.2f #leq Q^{2} < %.2f GeV^{2}",Q2_bins[i],Q2_bins[i+1])); 
    } else{
      tex1[i] = new TLatex(0.75E-2,0.9,Form(" %.2f #leq Q^{2} < %.2f GeV^{2}",Q2_bins[i],Q2_bins[i+1]));    
    }

    tex1[i]->SetTextFont(42);
    tex1[i]->SetTextSize(0.06);

  }

  TLatex *tex2[10];
  tex2[0] = new TLatex(0.15,0.95,"Electron Method (Using Optimized Calc.)");tex2[0]->SetNDC();tex2[0]->SetTextSize(0.05);
  tex2[1] = new TLatex(0.15,0.95,"Electron Method (Using Track Momentum)");tex2[1]->SetNDC();tex2[1]->SetTextSize(0.05);
  tex2[2] = new TLatex(0.15,0.95,"JB Method (Using Jets Scattering off Nucleon)");tex2[2]->SetNDC();tex2[2]->SetTextSize(0.05);
  tex2[3] = new TLatex(0.05,0.95,"JB Method (Using Jets Scattering off Nucleus with New Term)");tex2[3]->SetNDC();tex2[3]->SetTextSize(0.05);
  tex2[4] = new TLatex(0.05,0.95,"JB Method (Summing Over Hadrons Scattering off Nucleon)");tex2[4]->SetNDC();tex2[4]->SetTextSize(0.05);
  tex2[5] = new TLatex(0.1,0.95,"JB Method (Summing Over Hadrons Scattering off Nucleus with New Term)");tex2[5]->SetNDC();tex2[5]->SetTextSize(0.035);
  tex2[6] = new TLatex(0.15,0.95,"DA Method (Using Jets Scattering off Nucleon)");tex2[6]->SetNDC();tex2[6]->SetTextSize(0.05);
  tex2[7] = new TLatex(0.05,0.95,"DA Method (Using Jets Scattering off Nucleus with New Term)");tex2[7]->SetNDC();tex2[7]->SetTextSize(0.05);
  tex2[8] = new TLatex(0.05,0.95,"DA Method (Summing Over Hadrons Scattering off Nucleon)");tex2[8]->SetNDC();tex2[8]->SetTextSize(0.05);
  tex2[9] = new TLatex(0.1,0.95,"DA Method (Summing Over Hadrons Scattering off Nucleus with New Term)");tex2[9]->SetNDC();tex2[9]->SetTextSize(0.035);
				
  //Draw plots                                                                                                                                        
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelOffset(0.01,"X");
  gStyle->SetLabelOffset(0.01,"Y");
  gStyle->SetTitleXSize(0.07);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYSize(0.07);
  gStyle->SetTitleYOffset(0.8);

  float small = 1E-20;
  TCanvas *c[10][2][7]; //Index indicies read: method, purity(0)-stability(1), and number of canvases needed to plot all 25 Q^2 bins with 4 plots on
                        //each canvas.

  //The first index in the loop is method, second is purity(0)-stability(1), and the last index is referring to the number of canvases one method needs
  //For instance, the electron-method using the track momentum as energy purity plot has 25 Q^2 bins, plotting all 25 Q^2 bins and placing 4 plots in 
  //one canvas means tha we need a total of 7 canvases, thus the k range goes from 0-6 inclusive.
  for(int i=0;i<10;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<7;k++){

	c[i][j][k] = new TCanvas("");
	c[i][j][k]->Divide(2,2,small,small);

	for(int l=0;l<4;l++){
	  c[i][j][k]->cd(l+1);
	  gPad->SetLogx();
	  gPad->SetTickx();
	  gPad->SetTicky();

	  if(l==0){
	    gPad->SetBottomMargin(small);
	    gPad->SetRightMargin(small);
	    gPad->SetLeftMargin(0.25);
	    gPad->SetTopMargin(0.15);
	  } else if(l==1){
	    gPad->SetBottomMargin(small);
	    gPad->SetRightMargin(0.25);
	    gPad->SetLeftMargin(small);
	    gPad->SetTopMargin(0.15);
	  } else if(l==2){
	    gPad->SetTopMargin(small);
	    gPad->SetRightMargin(small);
	    gPad->SetLeftMargin(0.25);
	    gPad->SetBottomMargin(0.15);
	  } else{
	    gPad->SetTopMargin(small);
	    gPad->SetRightMargin(0.25);
	    gPad->SetLeftMargin(small);
	    gPad->SetBottomMargin(0.15);
	  }
	      
	  //There is a break statement because I want the 25th plot will reside alone in a canvas. The k*4+l is to ensure that 
	  //the plots being placed on the same canvas are 0-3, 4-7, 8-11, etc.
	      if(j==0){
		if(k==6){ 
		  TH1 *hframe1 = gPad->DrawFrame(1E-4,0,1,1);
		  hframe1->GetXaxis()->SetTitle("x");hframe1->GetXaxis()->CenterTitle();
		  hframe1->GetYaxis()->SetTitle("Purity");hframe1->GetYaxis()->CenterTitle();
		  gr[i][j][k*4+l]->Draw("P Same");
		  tex1[k*4+l]->Draw();
		  break;
		} else{
		  TH1 *hframe1 = gPad->DrawFrame(1E-4,0,1,1);
		  hframe1->GetXaxis()->SetTitle("x");hframe1->GetXaxis()->CenterTitle();
		  hframe1->GetYaxis()->SetTitle("Purity");hframe1->GetYaxis()->CenterTitle();
		  gr[i][j][k*4+l]->Draw("P Same");
		  tex1[k*4+l]->Draw();	     
		}
	      } else{
		if(k==6){ 
		  TH1 *hframe2 = gPad->DrawFrame(1E-4,0,1,1);
		  hframe2->GetXaxis()->SetTitle("x");hframe2->GetXaxis()->CenterTitle();
		  hframe2->GetYaxis()->SetTitle("Stability");hframe2->GetYaxis()->CenterTitle();
	          gr[i][j][k*4+l]->Draw("P Same");
		  tex1[k*4+l]->Draw();
		  break;
		} else{
		  TH1 *hframe2 = gPad->DrawFrame(1E-4,0,1,1);
		  hframe2->GetXaxis()->SetTitle("x");hframe2->GetXaxis()->CenterTitle();
		  hframe2->GetYaxis()->SetTitle("Stability");hframe2->GetYaxis()->CenterTitle();
		  gr[i][j][k*4+l]->Draw("P Same");
		  tex1[k*4+l]->Draw();
		}	  
	      }
	      
	      if(l==0){
		tex_energy->Draw();
	      }

	}
	
	c[i][j][k]->cd();  
	TPad *pad5 = new TPad("","",0,0,1,1);
	pad5->SetFillStyle(4000);  // transparent
	pad5->Draw();
	pad5->cd();
	tex2[i]->Draw();

      }
    }
  }


  //Print to File
  if(energy_set == 1){

    TFile *fout = new TFile("root_out/purity_stability_5_41_uncertainty_out.root","RECREATE");
    
    for(int j=0;j<10;j++){
      for(int k=0;k<2;k++){
	for(int l=0;l<nbins_Q2;l++){
	  gr[j][k][l]->Write();
	}
      }
    }

    fout->Write();fout->Close();

    
    for(int i=0;i<10;i++){
      for(int j=0;j<2;j++){
	for(int k=0;k<7;k++){

	  if(i==0 && j==0 & k==0){
	    c[i][j][k]->Print("./plots/purity_stability_5_41_uncertainty.pdf[");
	    c[i][j][k]->Print("./plots/purity_stability_5_41_uncertainty.pdf");
	  } else if(i==9 && j==1 && k==6){
	    c[i][j][k]->Print("./plots/purity_stability_5_41_uncertainty.pdf");
	    c[i][j][k]->Print("./plots/purity_stability_5_41_uncertainty.pdf]");
	  } else {
	    c[i][j][k]->Print("./plots/purity_stability_5_41_uncertainty.pdf");
	  }
	}
      }
    }

  }
 
 if(energy_set == 2){

    TFile *fout = new TFile("root_out/purity_stability_18_275_uncertainty_out.root","RECREATE");
    
    for(int j=0;j<10;j++){
      for(int k=0;k<2;k++){
	for(int l=0;l<nbins_Q2;l++){
	  gr[j][k][l]->Write();
	}
      }
    }
    
    fout->Write();fout->Close();

    for(int i=0;i<10;i++){
      for(int j=0;j<2;j++){
	for(int k=0;k<7;k++){

	  if(i==0 && j==0 & k==0){
	    c[i][j][k]->Print("./plots/purity_stability_18_275_uncertainty.pdf[");
	    c[i][j][k]->Print("./plots/purity_stability_18_275_uncertainty.pdf");
	  } else if(i==9 && j==1 && k==6){
	    c[i][j][k]->Print("./plots/purity_stability_18_275_uncertainty.pdf");
	    c[i][j][k]->Print("./plots/purity_stability_18_275_uncertainty.pdf]");
	  } else {
	    c[i][j][k]->Print("./plots/purity_stability_18_275_uncertainty.pdf");
	  }
	}
      }
    }

  }
  
}


