R__LOAD_LIBRARY(libeicsmear);

void rad_photons(){

  //For DJANGOH with QED radiation:
  //
  //Incoming electron is first in particle list. It usually has status==21; 
  //but for some low W_true events, it has status==201 (this may be a bug).
  //
  //Scattered electron is fourth in particle list. It has status==1, 
  //and has origin==1 (meaning incoming electron)
  //
  //Radiated photon is fifth in particle list. It has status==1. 
  //For IChannel==6 events, origin==1 (meaning incoming electron); for
  //IChannel==7 events, origin==4 (meaning scattered electron)
  //
  //There are also IChannel==8 events, which I need to understand better
  //

  //Load ROOT Files
  erhic::EventDjangoh *event(NULL); //Event Class
  erhic::ParticleMC *particle(NULL); //Particle Class
  
  TChain *tree = new TChain("EICTree");
  for(int i=0;i<15;i++){
  	tree->Add(Form("/eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/10_100/Rad/djangoh.NC.10x100_%d.root",i));
  }

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

  bool found_rad,found_scat;
  TLorentzVector rad_photon;
  TLorentzVector scat_electron;
  double theta_open(0);

  //Define Histograms
  TH1 *h1 = new TH1D("h1","IChannel=6 Events",200,-15,15);
  h1->GetXaxis()->SetTitle("#eta_{rad. photon}");h1->GetXaxis()->CenterTitle();
  h1->SetLineWidth(3);h1->SetLineColor(kBlue);

  TH1 *h2 = new TH1D("h2","IChannel=7 Events",200,-15,15);
  h2->GetXaxis()->SetTitle("#eta_{rad. photon}");h2->GetXaxis()->CenterTitle();
  h2->SetLineWidth(3);h2->SetLineColor(kBlue);

  TH1 *h3 = new TH1D("h3","Angle between radiated photon and scattered electron",200,-1,15);
  h3->GetXaxis()->SetTitle("#theta_{#gammae} [deg]");h3->GetXaxis()->CenterTitle();
  h3->SetLineWidth(3);h3->SetLineColor(kBlue);

  //Loop Over Events
  for(int iEvent=0;iEvent<nEntries;iEvent++){
    
    if(iEvent%10000==0) cout<<"Events Analysed = "<<iEvent<<"!"<<endl;
    tree->GetEntry(iEvent);

    //Reset variables
    found_rad = false; found_scat = false;

    channel = event->IChannel; //Event QED channel

    nParticles = event->GetNTracks();

    // Loop over all particles in event 
    for(int iParticle=0;iParticle<nParticles;iParticle++){

      particle = event->GetTrack(iParticle);
      
      id = (int) particle->Id();  //Particle Id Number
      status = (int) particle->GetStatus();  //Particle status
      orig = (int) particle->GetParentIndex(); //Particle origin (i.e. parent)

      //Get scattered electron
      if(id==11 && status==1 && orig==1){ //N.B. could also just use fourth particle in list
	scat_electron.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
	found_scat = true;
      }	

      //Get photons radiating off incoming electron
      if(channel==6 && id==22 && status==1 && orig==1){ //N.B. could also just use fifth particle in list
	rad_photon.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
	found_rad = true;
	h1->Fill(rad_photon.Eta());
      }

      //Get photons radiating off scattered electron
      if(channel==7 && id==22 && status==1 && orig==4){ //N.B. could also just use fifth particle in list
        rad_photon.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
 	found_rad = true;
        h2->Fill(rad_photon.Eta());
      }

    }//Finished Loop Over Particles

    //Fill additional histograms
    if(channel==7 && found_rad && found_scat){
    	theta_open = TMath::RadToDeg()*scat_electron.Angle(rad_photon.Vect());
   	h3->Fill(theta_open);
    }

  } //Finished loop over entries

  //Make plots
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);

  TCanvas *c1 = new TCanvas("c1");
  h1->Draw();

  TCanvas *c2 = new TCanvas("c2");
  h2->Draw();

  TCanvas *c3 = new TCanvas("c3");
  h3->Draw();
  //Add text
  TPaveText* tex1 = new TPaveText(0.6,0.6,0.9,0.9,"NDCNB");
  tex1->AddText("IChannel=7 Events");tex1->SetFillStyle(4000);tex1->SetTextFont(63);tex1->SetTextSize(20);
  tex1->Draw();

  //Print plots to pdf file
  c1->Print("plots/rad_photons.pdf[");
  c1->Print("plots/rad_photons.pdf");
  c2->Print("plots/rad_photons.pdf");
  c3->Print("plots/rad_photons.pdf");
  c3->Print("plots/rad_photons.pdf]");

  cout<<"-------------------------------"<<endl;

}
