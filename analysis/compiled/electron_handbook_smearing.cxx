#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/smear/EventSmear.h"

#include "Riostream.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TRint.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLine.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

int main(int argc, char **argv){

  #ifdef WITHRINT
  TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
  #else
  TApplication *myapp = new TApplication("myapp",0,0);
  #endif
  
  //Event Class
  erhic::EventPythia *event(NULL);
  
  //Particle Class
  erhic::ParticleMC *particle(NULL);

  //Smeared Event Class
  Smear::Event *event_s(NULL);
  
  //Smeared Particle Class
  Smear::ParticleMCS *particle_s(NULL);
  
  TFile *f = new TFile("/eic/data/baraks/pythiaeRHIC/outfiles/new_tune/ep_10_100_newtune_0.root");

  TTree *tree = (TTree*)f->Get("EICTree");
  tree->AddFriend("Smeared","/eic/data/baraks/pythiaeRHIC/outfiles/new_tune/ep_10_100_newtune_0_handbook_smeared.root");

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("eventS",&event_s);
  Int_t nevents = tree->GetEntries();
  
  //MC Variables
  Int_t nParticles(0);
  Int_t Status[500];
  Int_t id[500];
  Int_t orig[500];

  //Smeared Variables
  //Int_t nParticles_s(0);
  Int_t Status_s[500];
  //Int_t id_s[500];

  //Histograms
  TH2 *h1 = new TH2D("h1","Smeared Electron Energy vs. True Energy",100,0,20,100,0,20);
  h1->GetXaxis()->SetTitle("True Electron Energy [GeV]");h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->SetTitle("Smeared Electron Energy [GeV]");h1->GetYaxis()->CenterTitle();
  
  TH2 *h2 = new TH2D("h2","Smeared Electron Momentum vs. True Momentum",100,0,20,100,0,20);
  h2->GetXaxis()->SetTitle("True Electron Momentum [GeV/c]");h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->SetTitle("Smeared Electron Momentum [GeV/c]");h2->GetYaxis()->CenterTitle();

  TH2 *h3 = new TH2D("h3","Smeared Electron Polar Angle vs. True Polar Angle",100,0,180,100,0,180);
  h3->GetXaxis()->SetTitle("True Electron Angle [deg]");h3->GetXaxis()->CenterTitle();
  h3->GetYaxis()->SetTitle("Smeared Electron Angle [deg]");h3->GetYaxis()->CenterTitle(); 

  //Begin Analysis
  cout<<"Starting Analysis..."<<endl;
  cout<<"Total Number of Events = "<<nevents<<endl;
  
  // Loop over all events
  for(int i=0;i<nevents;i++){
    if(i%10000==0) cout<<"Events Analysed = "<<i<<"!"<<endl;
    
    tree->GetEntry(i);    

    nParticles = event->GetNTracks(); //Same as event_s->GetNTracks();
   
    // Loop over all particles in event 
    for(int j=0;j<nParticles;j++){
      
      particle = event->GetTrack(j);
      Status[j] = (Int_t) particle->GetStatus();
      id[j] = (Int_t) particle->Id();
      orig[j] = (Int_t) particle->orig;

      particle_s = event_s->GetTrack(j);

      if(particle_s){ //make sure not null pointer
        Status_s[j] = (Int_t) particle_s->GetStatus();
        //id_s[j] = (Int_t) particle_s->Id();

        if(Status_s[j]==1 /*&& id_s[j]==11*/ && 
	   id[j]==11 && Status[j]==1 && orig[j]==3){
      		h1->Fill(particle->E,particle_s->E);
        	h2->Fill(particle->p,particle_s->p);
                h3->Fill(particle->theta*TMath::RadToDeg(),particle_s->theta*TMath::RadToDeg());
      	}//Scattered Electron
      }//Particle NULL check 
    }//Finished Loop Over Particles
  }//Finished Event Loop


  //Make Plots
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1");
  h1->Draw("colz");

  TCanvas *c2 = new TCanvas("c2");
  h2->Draw("colz");

  TCanvas *c3 = new TCanvas("c3");
  h3->Draw("colz");

  //Print Plots to File
  c1->Print("plots/electron_handbook_smearing.pdf[");
  c1->Print("plots/electron_handbook_smearing.pdf");
  c2->Print("plots/electron_handbook_smearing.pdf");
  c3->Print("plots/electron_handbook_smearing.pdf");
  c3->Print("plots/electron_handbook_smearing.pdf]"); 
 
  myapp->Run();
  return 0;

}

