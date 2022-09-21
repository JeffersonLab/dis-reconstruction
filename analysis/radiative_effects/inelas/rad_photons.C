R__LOAD_LIBRARY(libeicsmear);

void rad_photons(){

  //For DJANGOH with QED radiation:
  //Updated for version 4.6.20
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

  bool found_rad,found_scat,found_prot,found_ebeam;
  TLorentzVector rad_photon;
  TLorentzVector scat_electron;
  TLorentzVector beam_prot;
  TLorentzVector beam_elec;
  TVector3 boost_vec;
  double theta_open(0);

  //Define Histograms
  TH1 *h1 = new TH1D("h1","IChannel 6 Events",200,-15,15);
  h1->GetXaxis()->SetTitle("#eta_{rad. photon}");h1->GetXaxis()->CenterTitle();
  h1->SetLineWidth(3);h1->SetLineColor(kBlue);

  TH1 *h2 = new TH1D("h2","IChannel 7 Events",200,-15,15);
  h2->GetXaxis()->SetTitle("#eta_{rad. photon}");h2->GetXaxis()->CenterTitle();
  h2->SetLineWidth(3);h2->SetLineColor(kBlue);

  TH1 *h3 = new TH1D("h3","IChannel 8 Events",200,-15,15);
  h3->GetXaxis()->SetTitle("#eta_{rad. photon}");h3->GetXaxis()->CenterTitle();
  h3->SetLineWidth(3);h3->SetLineColor(kBlue);

  TH1 *h4 = new TH1D("h4","Angle between radiated photon and scattered electron",200,-1,15);
  h4->GetXaxis()->SetTitle("#theta_{#gammae} [deg]");h4->GetXaxis()->CenterTitle();
  h4->SetLineWidth(3);h4->SetLineColor(kBlue);

  TH2 *etaE6 = new TH2D("etaE6", "Photon energy vs. direction (IChannel 6)", 100, -15, 15, 100, 0, 12);
  etaE6->GetXaxis()->SetTitle("#eta_{#gamma}");etaE6->GetXaxis()->CenterTitle();
  etaE6->GetYaxis()->SetTitle("E_{#gamma} [GeV]");etaE6->GetYaxis()->CenterTitle();

  TH2 *etaE7 = new TH2D("etaE7", "Photon energy vs. direction (IChannel 7)", 100, -15, 15, 100, 0, 11);
  etaE7->GetXaxis()->SetTitle("#eta_{#gamma}");etaE7->GetXaxis()->CenterTitle();
  etaE7->GetYaxis()->SetTitle("E_{#gamma} [GeV]");etaE7->GetYaxis()->CenterTitle();

  TH2 *etaE8 = new TH2D("etaE8", "Photon energy vs. direction (IChannel 8)", 100, -15, 15, 100, 0, 11);
  etaE8->GetXaxis()->SetTitle("#eta_{#gamma}");etaE8->GetXaxis()->CenterTitle();
  etaE8->GetYaxis()->SetTitle("E_{#gamma} [GeV]");etaE8->GetYaxis()->CenterTitle();

  TH2 *etaEall = new TH2D("etaEall", "Photon energy vs. direction (all Rad channels)", 100, -15, 15, 100, 0, 11);
  etaEall->GetXaxis()->SetTitle("#eta_{#gamma}");etaEall->GetXaxis()->CenterTitle();
  etaEall->GetYaxis()->SetTitle("E_{#gamma} [GeV]");etaEall->GetYaxis()->CenterTitle();

  TH1 *Eall = new TH1D("Eall","Photon energy (all Rad channels)", 100, -0.2, 10.5);
  Eall->GetXaxis()->SetTitle("E_{#gamma} [GeV]");Eall->GetXaxis()->CenterTitle();
  Eall->SetLineWidth(3);Eall->SetLineColor(kBlue);

  TH1 *Eall_rest = new TH1D("Eall_rest", "Photon energy (all Rad channels) in proton rest frame", 105, -5, 100);
  Eall_rest->GetXaxis()->SetTitle("E_{#gamma} [GeV]");Eall_rest->GetXaxis()->CenterTitle();
  Eall_rest->SetLineWidth(3);Eall_rest->SetLineColor(kBlue);

  TH1 *E_rest_scale = new TH1D("E_rest_scale", "(Photon energy / Electron beam energy) in proton rest frame",100,-0.1,1.1);
  E_rest_scale->GetXaxis()->SetTitle("E_{#gamma}/E_{e}");E_rest_scale->GetXaxis()->CenterTitle();
  E_rest_scale->SetLineWidth(3);E_rest_scale->SetLineColor(kBlue);

  TH1 *Eall_Q2 = new TH2D("Eall_Q2","Photon energy propagator dependence", 100, -0.2, 10.5,100,1e-1,1e4);
  Eall_Q2->GetXaxis()->SetTitle("E_{#gamma} [GeV]");Eall_Q2->GetXaxis()->CenterTitle();
  Eall_Q2->GetYaxis()->SetTitle("Q^{4}_{elec.} / Q^{4}_{true}");Eall_Q2->GetYaxis()->CenterTitle();
  Eall_Q2->SetTitleOffset(1.15,"Y");

  TH1 *Eall_weight = new TH1D("Eall_weight","Photon energy w/ propagator weighting", 100, -0.2, 10.5);
  Eall_weight->GetXaxis()->SetTitle("E_{#gamma} [GeV]");Eall_weight->GetXaxis()->CenterTitle();
  Eall_weight->SetLineWidth(3);Eall_weight->SetLineColor(kBlue);

  //Loop Over Events
  for(int iEvent=0;iEvent<nEntries;iEvent++){
    
    if(iEvent%10000==0) cout<<"Events Analysed = "<<iEvent<<"!"<<endl;
    tree->GetEntry(iEvent);

    //Reset variables
    found_rad = false; found_scat = false; found_prot = false; found_ebeam = false;

    channel = event->IChannel; //Event QED channel
    double prop_weight = (event->GetQ2()*event->GetQ2())/(event->dtrueQ2*event->dtrueQ2);

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
          boost_vec = -1.*beam_prot.BoostVector(); //To Boost to rest frame of the proton
          found_prot = true;

          //Write out for first event
          if(iEvent==0){
            cout<<"Proton Beam:"<<endl;
            beam_prot.Print();
          }
      }

      //Get electron beam (same momentum for each event)
      if(id==11 && status==21 && orig==0){ //N.B. could also just use first particle in list

          beam_elec.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
          found_ebeam = true;

          //Write out for first event
          if(iEvent==0){
            cout<<"Elec Beam:"<<endl;
            beam_elec.Print();
          }
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

	        h1->Fill(rad_photon.Eta());
          etaE6->Fill(rad_photon.Eta(),rad_photon.E());
          etaEall->Fill(rad_photon.Eta(),rad_photon.E()); //Fill for all rad channels (6,7,8)
          Eall->Fill(rad_photon.E()); //Fill for all rad channels (6,7,8)
          Eall_Q2->Fill( rad_photon.E(),prop_weight); //Fill for all rad channels (6,7,8)
          Eall_weight->Fill(rad_photon.E(),1./prop_weight); //Fill for all rad channels (6,7,8)
      }

      //Get photons radiating off scattered electron
      if(channel==7 && id==22 && status==1 && orig==4){ //N.B. could also just use fifth particle in list

          rad_photon.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
 	        found_rad = true;

          h2->Fill(rad_photon.Eta());
          etaE7->Fill(rad_photon.Eta(),rad_photon.E());
          etaEall->Fill(rad_photon.Eta(),rad_photon.E()); //Fill for all rad channels (6,7,8)
          Eall->Fill(rad_photon.E()); //Fill for all rad channels (6,7,8)
          Eall_Q2->Fill( rad_photon.E(),prop_weight); //Fill for all rad channels (6,7,8)
          Eall_weight->Fill(rad_photon.E(),1./prop_weight); //Fill for all rad channels (6,7,8)
      }

      //Get 'Compton' photons radiating off incoming electron
      if(channel==8 && id==22 && status==1 && orig==1){ //N.B. could also just use fifth particle in list

	        rad_photon.SetPxPyPzE(particle->GetPx(),particle->GetPy(),particle->GetPz(),particle->GetE());
	        found_rad = true;

	        h3->Fill(rad_photon.Eta());
          etaE8->Fill(rad_photon.Eta(),rad_photon.E());
          etaEall->Fill(rad_photon.Eta(),rad_photon.E()); //Fill for all rad channels (6,7,8)
          Eall->Fill(rad_photon.E()); //Fill for all rad channels (6,7,8)
          Eall_Q2->Fill( rad_photon.E(),prop_weight); //Fill for all rad channels (6,7,8)
          Eall_weight->Fill(rad_photon.E(),1./prop_weight); //Fill for all rad channels (6,7,8)
      }

    }//Finished Loop Over Particles

    //Fill additional histograms
    if(channel==7 && found_rad && found_scat){
    	  theta_open = TMath::RadToDeg()*scat_electron.Angle(rad_photon.Vect());
      	h4->Fill(theta_open);
    }

    if(found_prot && found_rad && found_ebeam){
        TLorentzVector rad_photon_rest = rad_photon;
        rad_photon_rest.Boost(boost_vec);
        Eall_rest->Fill(rad_photon_rest.E());

        TLorentzVector beam_elec_rest = beam_elec;
        beam_elec_rest.Boost(boost_vec);

        //Write out for first event
          if(iEvent==0){
            cout<<"Elec Beam in proton rest frame:"<<endl;
            beam_elec_rest.Print();
          }

        E_rest_scale->Fill(rad_photon_rest.E()/beam_elec_rest.E());
    }
    //Finished additional histograms

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

  TCanvas *c4 = new TCanvas("c4");
  h4->Draw();
  //Add text
  TPaveText* tex1 = new TPaveText(0.6,0.6,0.9,0.9,"NDCNB");
  tex1->AddText("IChannel=7 Events");tex1->SetFillStyle(4000);tex1->SetTextFont(63);tex1->SetTextSize(20);
  tex1->Draw();

  TCanvas *c5a = new TCanvas("c5a");
  etaE6->Draw("colz");

  TCanvas *c5b = new TCanvas("c5b");
  etaE7->Draw("colz");

  TCanvas *c5c = new TCanvas("c5c");
  etaE8->Draw("colz");

  TCanvas *c5d = new TCanvas("c5d");
  gPad->SetLogz();
  etaEall->Draw("colz");

  TCanvas *c6a = new TCanvas("c6a");
  Eall->Draw();

  TCanvas *c6b = new TCanvas("c6b");
  Eall_rest->Draw();

  TCanvas *c6c = new TCanvas("c6c");
  E_rest_scale->Draw();

  TCanvas *c7 = new TCanvas("c7");
  gPad->SetLogz();
  Eall_Q2->Draw("colz");
  //Add text
  TPaveText* tex2 = new TPaveText(0.1,0.6,0.4,0.9,"NDCNB");
  tex2->AddText("Q^{2}_{elec.} > 0.5 GeV^{2}");tex2->SetFillStyle(4000);tex2->SetTextFont(63);tex2->SetTextSize(20);
  tex2->Draw();

  TCanvas *c8 = new TCanvas("c8");
  Eall_weight->Draw("hist");
  //Add text
  TPaveText* tex3 = new TPaveText(0.6,0.6,0.6,0.9,"NDCNB");
  tex3->AddText("Event weight = #frac{Q^{4}_{elec}}{Q^{4}_{true}}");
  tex3->SetFillStyle(4000);tex3->SetTextFont(63);tex3->SetTextSize(20);
  tex3->Draw();

  //Print plots to pdf file
  c1->Print("plots/rad_photons.pdf[");
  c1->Print("plots/rad_photons.pdf");
  c2->Print("plots/rad_photons.pdf");
  c3->Print("plots/rad_photons.pdf");
  c4->Print("plots/rad_photons.pdf");
  c5a->Print("plots/rad_photons.pdf");
  c5b->Print("plots/rad_photons.pdf");
  c5c->Print("plots/rad_photons.pdf");
  c5d->Print("plots/rad_photons.pdf");
  c6a->Print("plots/rad_photons.pdf");
  c6b->Print("plots/rad_photons.pdf");
  c6c->Print("plots/rad_photons.pdf");
  c7->Print("plots/rad_photons.pdf");
  c8->Print("plots/rad_photons.pdf");
  c8->Print("plots/rad_photons.pdf]");

  cout<<"-------------------------------"<<endl;

}
