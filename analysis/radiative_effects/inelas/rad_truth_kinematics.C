R__LOAD_LIBRARY(libeicsmear);

void rad_truth_kinematics(){

  //For DJANGOH with QED radiation:
  //Using version 4.6.20
  //
  //Incoming electron is first in particle list. It has status==21.
  //Incoming proton is second in particle list. It also has status==21.
  //
  //Scattered electron is fourth in particle list. It has status==1, 
  //and has origin==1 (meaning incoming electron). In the case of FSR,
  //this should be the electron post radiation (even though it is the parent of 
  //the radiated photon)
  //
  //Radiated photon is fifth in particle list. It has status==1. 
  //For IChannel==6 events, origin==1 (meaning incoming electron).
  //For IChannel==7 events, origin==4 (meaning scattered electron).
  //For IChannel==8 events, origin==1 (meaning incoming electron).
  //
  //We want to calculate the true kinmematics using the post-ISR incoming electron
  //or the pre-ISR outgoing electron and compare that to what Djangoh gives for truth
  //quantities

  //constants
  const double Mp = 0.93827;

  //Load ROOT Files
  erhic::EventDjangoh *event(NULL); //Event Class
  erhic::ParticleMC *particle(NULL); //Particle Class
  
  TChain *tree = new TChain("EICTree");
  //Djangoh 4.6.20
  tree->Add("/gpfs02/eic/baraks/djangoh/djangoh_local/outfiles/djangoh.NC.Rad_inelas.10x100_evt.root");

  int nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  //Variables
  int channel(0);
  int nParticles(0);
  int status(0);
  int orig(0);
  int id(0);

  bool found_rad,found_scat,found_pbeam,found_ebeam;
  TLorentzVector rad_photon;
  TLorentzVector scat_electron;
  TLorentzVector beam_prot;
  TLorentzVector beam_elec;
  TLorentzVector beam_elec_preISR;

  //True kinematic quantities from Djangoh
  double Q2_true; 
  double x_true; 
  double y_true;
  double W_true;
  double nu_true;

  //True kinematic quantities from 4-momentum formulae
  double Q2_true_calc; 
  double x_true_calc; 
  double y_true_calc;
  double y_had_calc;
  double W_true_calc;
  double nu_true_calc;

  //Define Histograms
  TH2 *h1 = new TH2D("h1","",100,0,100,100,0,100);
  h1->GetXaxis()->SetTitle("Djangoh True Q^{2} [GeV^{2}]");h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle("Calc. True Q^{2} [GeV^{2}]");h1->GetYaxis()->CenterTitle();

  //Q2 Binning                                                         
  double Q2_min = 1E-10;
  double Q2_max = 1E4;
  const int nbins_Q2 = 200;
  double log_bw_Q2 = (log10(Q2_max) - log10(Q2_min))/(nbins_Q2); //Determine bin width
  double Q2_bins[nbins_Q2+1];
  for(int i=0;i<nbins_Q2+1;i++){
    double log_Q2_div = log10(Q2_min) + (i*log_bw_Q2);
    Q2_bins[i] = pow(10,log_Q2_div);
  }

  TH2 *h1a = new TH2D("h1a","",nbins_Q2,Q2_bins,nbins_Q2,Q2_bins);
  h1a->GetXaxis()->SetTitle("Djangoh True Q^{2} [GeV^{2}]");h1a->GetXaxis()->CenterTitle();
  h1a->GetYaxis()->SetTitle("Calc. True Q^{2} [GeV^{2}]");h1a->GetYaxis()->CenterTitle();

  TH2 *h2 = new TH2D("h2","",100,0,70,100,0,70);
  h2->GetXaxis()->SetTitle("Djangoh True W [GeV]");h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->SetTitle("Calc. True W [GeV]");h2->GetYaxis()->CenterTitle();

  TH2 *h2a = new TH2D("h2a","",100,0,10,100,0,10);
  h2a->GetXaxis()->SetTitle("Djangoh True W [GeV]");h2a->GetXaxis()->CenterTitle();
  h2a->GetYaxis()->SetTitle("Calc. True W [GeV]");h2a->GetYaxis()->CenterTitle();

  //x Binning                                                         
  double x_min = 1e-4;
  double x_max = 10;
  const int nbins_x = 200;
  double log_bw_x = (log10(x_max) - log10(x_min))/(nbins_x); //Determine bin width
  double x_bins[nbins_x+1];
  for(int i=0;i<nbins_x+1;i++){
    double log_x_div = log10(x_min) + (i*log_bw_x);
    x_bins[i] = pow(10,log_x_div);
  }

  TH2 *h3 = new TH2D("h3","",nbins_x,x_bins,nbins_x,x_bins);
  h3->GetXaxis()->SetTitle("Djangoh True x");h3->GetXaxis()->CenterTitle();
  h3->GetYaxis()->SetTitle("Calc. True x");h3->GetYaxis()->CenterTitle();

  TH2 *h3a = new TH2D("h3a","Channel 1",nbins_x,x_bins,nbins_x,x_bins);
  h3a->GetXaxis()->SetTitle("Djangoh True x");h3a->GetXaxis()->CenterTitle();
  h3a->GetYaxis()->SetTitle("Calc. True x");h3a->GetYaxis()->CenterTitle();

  TH2 *h3b = new TH2D("h3b","Channel 6",nbins_x,x_bins,nbins_x,x_bins);
  h3b->GetXaxis()->SetTitle("Djangoh True x");h3b->GetXaxis()->CenterTitle();
  h3b->GetYaxis()->SetTitle("Calc. True x");h3b->GetYaxis()->CenterTitle();

  TH2 *h3c = new TH2D("h3c","Channel 7",nbins_x,x_bins,nbins_x,x_bins);
  h3c->GetXaxis()->SetTitle("Djangoh True x");h3c->GetXaxis()->CenterTitle();
  h3c->GetYaxis()->SetTitle("Calc. True x");h3c->GetYaxis()->CenterTitle();

  TH2 *h3d = new TH2D("h3d","Channel 8",nbins_x,x_bins,nbins_x,x_bins);
  h3d->GetXaxis()->SetTitle("Djangoh True x");h3d->GetXaxis()->CenterTitle();
  h3d->GetYaxis()->SetTitle("Calc. True x");h3d->GetYaxis()->CenterTitle();

  TH2 *h4 = new TH2D("h4","",100,0,100,100,0,100);
  h4->GetXaxis()->SetTitle("Djangoh True #nu [GeV]");h4->GetXaxis()->CenterTitle();
  h4->GetYaxis()->SetTitle("Calc. True #nu [GeV]");h4->GetYaxis()->CenterTitle();

  //y Binning                                                         
  double y_min = 1e-6;
  double y_max = 2;
  const int nbins_y = 200;
  double log_bw_y = (log10(y_max) - log10(y_min))/(nbins_y); //Determine bin width
  double y_bins[nbins_y+1];
  for(int i=0;i<nbins_y+1;i++){
    double log_y_div = log10(y_min) + (i*log_bw_y);
    y_bins[i] = pow(10,log_y_div);
  }

  TH2 *h5 = new TH2D("h5","",nbins_y,y_bins,nbins_y,y_bins);
  h5->GetXaxis()->SetTitle("Djangoh True y");h5->GetXaxis()->CenterTitle();
  h5->GetYaxis()->SetTitle("Calc. True y");h5->GetYaxis()->CenterTitle();

  TH2 *h5a = new TH2D("h5a","",nbins_y,y_bins,nbins_y,y_bins);
  h5a->GetXaxis()->SetTitle("Djangoh True y");h5a->GetXaxis()->CenterTitle();
  h5a->GetYaxis()->SetTitle("Calc. Hadron y");h5a->GetYaxis()->CenterTitle();

  //Loop Over Events
  for(int iEvent=0;iEvent<nEntries;iEvent++){
    
    if(iEvent%10000==0) cout<<"Events Analysed = "<<iEvent<<"!"<<endl;
    tree->GetEntry(iEvent);

    //Reset variables
    found_rad = false; found_scat = false; found_pbeam = false; found_ebeam = false;

    channel = event->IChannel; //Event QED channel
    
    //Djangoh true quantities
    Q2_true = event->dtrueQ2;
    x_true = event->dtrueX;
    y_true = event->dtrueY;
    W_true = sqrt(event->dtrueW2);
    nu_true = event->dtrueNu;

    nParticles = event->GetNTracks();

    // Loop over all particles in event 
    for(int iParticle=0;iParticle<nParticles;iParticle++){

      particle = event->GetTrack(iParticle);
      
      id = (int) particle->Id();  //Particle Id Number
      status = (int) particle->GetStatus();  //Particle status
      orig = (int) particle->GetParentIndex(); //Particle origin (i.e. parent)

      //Get proton beam (same momentum for each event)
      if(id==2212 && status==21 && orig==0){ //N.B. could also just use second particle in list
          beam_prot.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
          found_pbeam = true;
      }

      //Get electron beam (same momentum for each event)
      if(id==11 && status==21 && orig==0){ //N.B. could also just use first particle in list
          beam_elec.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
          beam_elec_preISR = beam_elec;
          found_ebeam = true;
      }

      //Get scattered electron
      if(id==11 && status==1 && orig==1){ //N.B. could also just use fourth particle in list
	        scat_electron.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
	        found_scat = true;
      }	

      //Get photons radiating off incoming electron
      if(channel==6 && id==22 && status==1 && orig==1){ //N.B. could also just use fifth particle in list
	        rad_photon.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
	        found_rad = true;

          //Adjust 4-momentum of incoming electron
          beam_elec = beam_elec - rad_photon;
      }

      //Get photons radiating off scattered electron
      if(channel==7 && id==22 && status==1 && orig==4){ //N.B. could also just use fifth particle in list

          rad_photon.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
 	        found_rad = true;

          //Adjust 4-momentum of outgoin electron
          scat_electron = scat_electron + rad_photon;
      }

      //Get 'Compton' photons radiating off incoming electron
      if(channel==8 && id==22 && status==1 && orig==1){ //N.B. could also just use fifth particle in list

	        rad_photon.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
	        found_rad = true;

	        //Adjust 4-momentum of incoming electron
          beam_elec = beam_elec - rad_photon;
      }

    }//Finished Loop Over Particles

    //Calculate truth kinematic quantities
    auto q = (beam_elec - scat_electron);
    Q2_true_calc = -1.0*(q*q);
    y_true_calc = (q*beam_prot) / (beam_elec*beam_prot);
    y_had_calc = (q*beam_prot) / (beam_elec_preISR*beam_prot);
    x_true_calc = Q2_true_calc / (2*beam_prot*q);
    W_true_calc = sqrt(Mp*Mp + 2*q*beam_prot - Q2_true_calc);
    nu_true_calc = (q*beam_prot) / Mp;

    //Fill histograms
    h1->Fill(Q2_true,Q2_true_calc);h1a->Fill(Q2_true,Q2_true_calc);
    h2->Fill(W_true,W_true_calc);h2a->Fill(W_true,W_true_calc);
    h3->Fill(x_true,x_true_calc);
    h4->Fill(nu_true,nu_true_calc);
    h5->Fill(y_true,y_true_calc);h5a->Fill(y_true,y_had_calc);

    //Individual channels
    int channel = event->IChannel; //Event QED channel
    if(channel==1) h3a->Fill(x_true,x_true_calc);
    if(channel==6) h3b->Fill(x_true,x_true_calc);
    if(channel==7) h3c->Fill(x_true,x_true_calc);
    if(channel==8) h3d->Fill(x_true,x_true_calc);

  } //Finished loop over entries

  //Make plots
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);

  TCanvas *c1 = new TCanvas("c1");
  h1->Draw("colz");

  TCanvas *c1a = new TCanvas("c1a");
  c1a->SetLogx();c1a->SetLogy();
  h1a->Draw("colz");

  TCanvas *c2 = new TCanvas("c2");
  h2->Draw("colz");

  TCanvas *c2a = new TCanvas("c2a");
  h2a->Draw("colz");

  TCanvas *c3 = new TCanvas("c3");
  c3->SetLogx();c3->SetLogy();
  h3->Draw("colz");

  TCanvas *c3a = new TCanvas("c3a");
  c3a->Divide(2,2);
  c3a->cd(1);
  gPad->SetLogx();gPad->SetLogy();
  h3a->Draw("colz");
  c3a->cd(2);
  gPad->SetLogx();gPad->SetLogy();
  h3b->Draw("colz");
  c3a->cd(3);
  gPad->SetLogx();gPad->SetLogy();
  h3c->Draw("colz");
  c3a->cd(4);
  gPad->SetLogx();gPad->SetLogy();
  h3d->Draw("colz");

  TCanvas *c4 = new TCanvas("c4");
  h4->Draw("colz");

  TCanvas *c5 = new TCanvas("c5");
  c5->SetLogx();c5->SetLogy();
  h5->Draw("colz");

  TCanvas *c5a = new TCanvas("c5a");
  c5a->SetLogx();c5a->SetLogy();
  h5a->Draw("colz");

  //Print plots to pdf file
  c1->Print("plots/rad_truth_kinematics.pdf[");
  c1->Print("plots/rad_truth_kinematics.pdf");
  c1a->Print("plots/rad_truth_kinematics.pdf");
  c2->Print("plots/rad_truth_kinematics.pdf");
  c2a->Print("plots/rad_truth_kinematics.pdf");
  c3->Print("plots/rad_truth_kinematics.pdf");
  c3a->Print("plots/rad_truth_kinematics.pdf");
  c4->Print("plots/rad_truth_kinematics.pdf");
  c5->Print("plots/rad_truth_kinematics.pdf");
  c5a->Print("plots/rad_truth_kinematics.pdf");
  c5a->Print("plots/rad_truth_kinematics.pdf]");

  cout<<"-------------------------------"<<endl;

}
