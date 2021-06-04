#include "eicsmear/erhic/EventPythia.h"

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
#include "TGraph.h"
#include "TF1.h"
#include "TChain.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

//Globals
const double Mp(0.9383);
const double Me(0.511E-3);
const double s_cm = 4.*5.*41; //beam energy combination

//---------------------
//W_min = 2 GeV limit (for Q2 vs. y plot)
double W_limit_y(Double_t *x, Double_t *par){

    double y = x[0];
    double W_limit = 2;

    double Q2_limit = (Mp*Mp - W_limit*W_limit + s_cm*y);
    return Q2_limit;
}

//---------------------
//Q2_min = me**2 * y**2 / (1-y) limit  (for Q2 vs. y plot)
double photo_limit_y(Double_t *x, Double_t *par){

    double y = x[0];

    double Q2_limit = Me*Me*y*y / (1.-y);
    return Q2_limit;
}

//---------------------
//W_min = 2 GeV limit (for Q2 vs. x plot)
double W_limit_x(Double_t *x, Double_t *par){

    double xb = x[0];
    double W_limit = 2;

    double Q2_limit = (W_limit*W_limit - Mp*Mp )/( 1./xb - 1.);
    return Q2_limit;
}

//---------------------
//Q2_max = x**2 * s**2 / (me**2 + xs) limit (for Q2 vs. x plot)
double photo_limit_x(Double_t *x, Double_t *par){

    double xb = x[0];

    double Q2_limit = xb*xb*s_cm*s_cm / (Me*Me + xb*s_cm);
    return Q2_limit;
}

//---------------------
//Main function: 1) plot Q2 vs. y and Q2 vs. x
//               2) plot electron, positron, pi-minus yields as function...
//                  ...of momentum for eta bins
//               3) plot selection variables for scattered electron and pi-minus
int main(int argc, char **argv){

    #ifdef WITHRINT
    TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
    #else
    TApplication *myapp = new TApplication("myapp",0,0);
    #endif

    //Define TGraphs and Histograms
    TGraph *gr1 = new TGraph(); //Q2 vs. y
    TGraph *gr2 = new TGraph(); //Q2 vs. x

    int nEta = 6; //Number of eta bins
    double Eta_low[] = {-3.5, -2.0, -1.0, 0.0, 1.0, 2.0};
    double Eta_hi[] =  {-2.0, -1.0,  0.0, 1.0, 2.0, 3.5};

    int nEta_fine = 6; //Number of finer eta bins for eta<0;
    double Eta_fine_low[] = {-4.0, -3.0, -2.0, -1.5, -1.0, -0.5};
    double Eta_fine_hi[] =  {-3.0, -2.0, -1.5, -1.0, -0.5,  0.0};

    int nEta_finer = 3; //Number of finer eta bins for -1<eta<0;
    double Eta_finer_low[] = {-1.0, -0.6, -0.3};
    double Eta_finer_hi[] =  {-0.6, -0.3,  0.0};

    TH1* h_se[nEta];
    TH1* h_ae[nEta];
    TH1* h_pos[nEta];
    TH1* h_npi[nEta];

    //Four 'topological' distributions for scattered electron and negative pions
    //1)  Number of particles with R<0.7
    //2)  Number of additional particles that have |eta|<3.5
    //3)  Phi difference with sum of all particles that have |eta|<3.5
    //3a) Phi difference with sum of all particles that have |eta|<4.0
    //4)  |Pt| difference with sum of all particles that have |eta|<3.5
    //4a) |Pt| difference with sum of all particles that have |eta|<4.0
    //5)  Total E-P_z - 2E_e for all particles that have |eta|<3.5
    //6)  Momentum distributions w/ cuts on 4) and 5)
    TH1* h_se_1[nEta];
    TH1* h_npi_1[nEta];
    TH1* h_se_2[nEta];
    TH1* h_npi_2[nEta];
    TH1* h_se_3[nEta];
    TH1* h_npi_3[nEta];
    TH1* h_se_3a[nEta];
    TH1* h_npi_3a[nEta];
    TH1* h_se_4[nEta];
    TH1* h_npi_4[nEta];
    TH1* h_se_4a[nEta];
    TH1* h_npi_4a[nEta];
    TH1* h_se_5[nEta];
    TH1* h_npi_5[nEta];
    TH1* h_se_6[nEta];
    TH1* h_npi_6[nEta];

    for(int ihist=0;ihist<nEta;ihist++){
        h_se[ihist] = new TH1D(Form("h_se[%d]",ihist),"",1000,0,50);
        h_se[ihist]->SetLineColor(kGreen);h_se[ihist]->SetLineWidth(2);

        h_ae[ihist] = new TH1D(Form("h_ae[%d]",ihist),"",1000,0,50);
        h_ae[ihist]->SetLineColor(kBlack);h_ae[ihist]->SetLineWidth(2);

        h_pos[ihist] = new TH1D(Form("h_pos[%d]",ihist),"",1000,0,50);
        h_pos[ihist]->SetLineColor(kRed);h_pos[ihist]->SetLineWidth(2);

        h_npi[ihist] = new TH1D(Form("h_npi[%d]",ihist),"",1000,0,50);
        h_npi[ihist]->SetLineColor(kBlue);h_npi[ihist]->SetLineWidth(2);

        h_se_1[ihist] = new TH1D(Form("h_se_1[%d]",ihist),"",10,0,10);
        h_se_1[ihist]->SetLineColor(kGreen);h_se_1[ihist]->SetLineWidth(2);

        h_npi_1[ihist] = new TH1D(Form("h_npi_1[%d]",ihist),"",10,0,10);
        h_npi_1[ihist]->SetLineColor(kBlue);h_npi_1[ihist]->SetLineWidth(2);

        h_se_2[ihist] = new TH1D(Form("h_se_2[%d]",ihist),"",20,0,20);
        h_se_2[ihist]->SetLineColor(kGreen);h_se_2[ihist]->SetLineWidth(2);

        h_npi_2[ihist] = new TH1D(Form("h_npi_2[%d]",ihist),"",20,0,20);
        h_npi_2[ihist]->SetLineColor(kBlue);h_npi_2[ihist]->SetLineWidth(2);

        h_se_3[ihist] = new TH1D(Form("h_se_3[%d]",ihist),"",100,-0.1,3.25);
        h_se_3[ihist]->SetLineColor(kGreen);h_se_3[ihist]->SetLineWidth(2);

        h_npi_3[ihist] = new TH1D(Form("h_npi_3[%d]",ihist),"",100,-0.1,3.25);
        h_npi_3[ihist]->SetLineColor(kBlue);h_npi_3[ihist]->SetLineWidth(2);

        h_se_3a[ihist] = new TH1D(Form("h_se_3a[%d]",ihist),"",100,-0.1,3.25);
        h_se_3a[ihist]->SetLineColor(kGreen);h_se_3a[ihist]->SetLineWidth(2);

        h_npi_3a[ihist] = new TH1D(Form("h_npi_3a[%d]",ihist),"",100,-0.1,3.25);
        h_npi_3a[ihist]->SetLineColor(kBlue);h_npi_3a[ihist]->SetLineWidth(2);

        h_se_4[ihist] = new TH1D(Form("h_se_4[%d]",ihist),"",100,-3,3);
        h_se_4[ihist]->SetLineColor(kGreen);h_se_4[ihist]->SetLineWidth(2);

        h_npi_4[ihist] = new TH1D(Form("h_npi_4[%d]",ihist),"",100,-3,3);
        h_npi_4[ihist]->SetLineColor(kBlue);h_npi_4[ihist]->SetLineWidth(2);

        h_se_4a[ihist] = new TH1D(Form("h_se_4a[%d]",ihist),"",100,-3,3);
        h_se_4a[ihist]->SetLineColor(kGreen);h_se_4a[ihist]->SetLineWidth(2);

        h_npi_4a[ihist] = new TH1D(Form("h_npi_4a[%d]",ihist),"",100,-3,3);
        h_npi_4a[ihist]->SetLineColor(kBlue);h_npi_4a[ihist]->SetLineWidth(2);

        h_se_5[ihist] = new TH1D(Form("h_se_5[%d]",ihist),"",100,-11,1);
        h_se_5[ihist]->SetLineColor(kGreen);h_se_5[ihist]->SetLineWidth(2);

        h_npi_5[ihist] = new TH1D(Form("h_npi_5[%d]",ihist),"",100,-11,1);
        h_npi_5[ihist]->SetLineColor(kBlue);h_npi_5[ihist]->SetLineWidth(2);

        h_se_6[ihist] = new TH1D(Form("h_se_6[%d]",ihist),"",1000,0,50);
        h_se_6[ihist]->SetLineColor(kGreen);h_se_6[ihist]->SetLineWidth(2);

        h_npi_6[ihist] = new TH1D(Form("h_npi_6[%d]",ihist),"",1000,0,50);
        h_npi_6[ihist]->SetLineColor(kBlue);h_npi_6[ihist]->SetLineWidth(2);
    }

    TH1* h_se_fine[nEta_fine];
    TH1* h_npi_fine[nEta_fine];

    for(int ihist=0;ihist<nEta_fine;ihist++){
        h_se_fine[ihist] = new TH1D(Form("h_se_fine[%d]",ihist),"",1000,0,50);
        h_se_fine[ihist]->SetLineColor(kGreen);h_se_fine[ihist]->SetLineWidth(2);

        h_npi_fine[ihist] = new TH1D(Form("h_npi_fine[%d]",ihist),"",1000,0,50);
        h_npi_fine[ihist]->SetLineColor(kBlue);h_npi_fine[ihist]->SetLineWidth(2);
    }

    TH1* h_se_finer[nEta_finer];
    TH1* h_npi_finer[nEta_finer];

    for(int ihist=0;ihist<nEta_finer;ihist++){
        h_se_finer[ihist] = new TH1D(Form("h_se_finer[%d]",ihist),"",1000,0,50);
        h_se_finer[ihist]->SetLineColor(kGreen);h_se_finer[ihist]->SetLineWidth(2);

        h_npi_finer[ihist] = new TH1D(Form("h_npi_finer[%d]",ihist),"",1000,0,50);
        h_npi_finer[ihist]->SetLineColor(kBlue);h_npi_finer[ihist]->SetLineWidth(2);
    }

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class
    erhic::ParticleMC *part_se(NULL);

    TChain *t = new TChain("EICTree");

   //Using files created with my installation
   for(int i=0;i<500;i++){
        t->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/other_studies/5_41/ep_minbias_%d.root",i));
    }

    t->SetBranchAddress("event",&event);

    //We need to use the 'true' kinematic quantities given by pythia (e.g. trueQ2)
    //since the electron ones (e.g. QSquared) seem to be calculated assuming 
    //zero electron mass -- a bad assumption at very low Q2. At higher Q2, both
    //sets of quantities are the same when the QED effects are turned OFF.
    //When using the true quantities, it seems like they have any accuracy limit
    //of 10^-11; so events generated at lower x or Q2 will be put at the limiting value.
    //This is just a problem with the writing out of those variables, rather 
    //than an issue with the event generation for the very low x and Q2 events.
    //I fixed this issue in my personal version of pythiaeRHIC, and now I see the
    //events generated all the way down to the minimum allowed values.
    double Q2_true;
    double x_true;
    double y_true;

    //Particle variables
    int nParticles(0);
    int id(0);
    int status(0);
    int orig(0);

    //'Topological' variables
    double R(0);
    double phi_diff(0),eta_diff(0);
    int index_se(0),bin_se(0),count_se(0);
    int num_npi(0);
    int index_npi[50],bin_npi[50],count_npi[50];
    double phi_npi[50],eta_npi[50],p_npi[50],pt_npi[50];
    int ntotal(0);
    double pxsum_se(0),pysum_se(0);
    double pxsum_npi[50],pysum_npi[50];
    double pxsum_se_a(0),pysum_se_a(0);
    double pxsum_npi_a[50],pysum_npi_a[50];
    double empztot(0);

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    double cross_tot = 80E9; //Total Cross Section in fb
    double lum = ( (double) nevents)/cross_tot; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"PYTHIA6 Simulation:"<<endl;
    cout<<"Total Number of Events = "<<nevents<<endl;
    cout<<"Integrated Luminosity = "<<lum<<" fb^-1"<<endl<<endl;

    //Loop over events
    for(int iEvent=0;iEvent<nevents;iEvent++){
        //for(int iEvent=0;iEvent<1E5;iEvent++){ //for debugging
        
        if(iEvent%10000==0) cout<<"Events Analysed = "<<iEvent<<"!"<<endl;
        t->GetEntry(iEvent);

        Q2_true = event->GetTrueQ2();
        x_true = event->GetTrueX();
        y_true = event->GetTrueY();

        //Fill kinematic graphs
        //If we draw all events, the output pdf file size gets really big
        //...so just draw the first 200k events here
        //...since we just need to see the cuts
        if(iEvent<2E5){
            gr1->SetPoint(iEvent,y_true,Q2_true);
            gr2->SetPoint(iEvent,x_true,Q2_true);        
        }

        //Get Number of tracks
        nParticles = event->GetNTracks();

        //Reset variables
        num_npi = 0;count_se = 0; part_se = NULL;
        ntotal = 0;
        pxsum_se = 0; pysum_se = 0;
        pxsum_se_a = 0; pysum_se_a = 0;
        empztot = 0;

        // Loop over all particles in event
        for(int iParticle=0;iParticle<nParticles;iParticle++){
      
            particle = event->GetTrack(iParticle);

            id = (int) particle->Id();
            status = (int) particle->GetStatus();
            orig = (int) particle->GetParentIndex();

            if(id==11 && status==1 && orig==3){ //Scattered Electron
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                        h_se[iEta]->Fill(particle->GetP());
                        index_se = iParticle;
                        bin_se = iEta;
                        part_se = event->GetTrack(iParticle);
                    }
                }
                for(int iEta=0;iEta<nEta_fine;iEta++){
                    if(particle->GetEta()>Eta_fine_low[iEta] && particle->GetEta()<Eta_fine_hi[iEta]){
                        h_se_fine[iEta]->Fill(particle->GetP());
                    }
                }
                for(int iEta=0;iEta<nEta_finer;iEta++){
                    if(particle->GetEta()>Eta_finer_low[iEta] && particle->GetEta()<Eta_finer_hi[iEta]){
                        h_se_finer[iEta]->Fill(particle->GetP());
                    }
                }
            }

            if(id==11 && status==1){ //All Electrons
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta])
                        h_ae[iEta]->Fill(particle->GetP());
                }
            }

            if(id==-11 && status==1){ //Positrons
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta])
                        h_pos[iEta]->Fill(particle->GetP());
                }
            }

            if(id==-211 && status==1){ //Negative Pions
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                        h_npi[iEta]->Fill(particle->GetP());

                        index_npi[num_npi] = iParticle;
                        bin_npi[num_npi] = iEta;
                        phi_npi[num_npi] = particle->GetPhi();
                        eta_npi[num_npi] = particle->GetEta();
                        p_npi[num_npi] = particle->GetP();
                        pt_npi[num_npi] = particle->GetPt();
                        num_npi++;
                    }
                }
                for(int iEta=0;iEta<nEta_fine;iEta++){
                    if(particle->GetEta()>Eta_fine_low[iEta] && particle->GetEta()<Eta_fine_hi[iEta]){
                        h_npi_fine[iEta]->Fill(particle->GetP());
                    }
                }
                for(int iEta=0;iEta<nEta_finer;iEta++){
                    if(particle->GetEta()>Eta_finer_low[iEta] && particle->GetEta()<Eta_finer_hi[iEta]){
                        h_npi_finer[iEta]->Fill(particle->GetP());
                    }
                }
            }

        }

        //Reset pion variables
        for(int iPion=0;iPion<num_npi;iPion++){
            count_npi[iPion] = 0;
            pxsum_npi[iPion] = 0;pysum_npi[iPion] = 0;
            pxsum_npi_a[iPion] = 0;pysum_npi_a[iPion] = 0;
        }

        //Second Loop over all particles
        for(int iParticle=0;iParticle<nParticles;iParticle++){
      
            particle = event->GetTrack(iParticle);
            status = (int) particle->GetStatus();

            //Scattered electron
            if(status==1 && part_se!=NULL){ 

                phi_diff = fabs(particle->GetPhi() - part_se->GetPhi());
                if(phi_diff>TMath::Pi())
                    phi_diff = TMath::TwoPi() - phi_diff;
                
                eta_diff = particle->GetEta() - part_se->GetEta();
                R = sqrt( pow(phi_diff,2) + pow(eta_diff,2) );
                
                if( R<0.7 && iParticle!=index_se ){
                    count_se++;
                }

                if( iParticle!=index_se && fabs(particle->GetEta())<3.5 ){
                    pxsum_se+=particle->GetPx();
                    pysum_se+=particle->GetPy();

                }
                if( iParticle!=index_se && fabs(particle->GetEta())<4.0 ){
                    pxsum_se_a+=particle->GetPx();
                    pysum_se_a+=particle->GetPy();
                }       
            }

            //Negative pions
            for(int iPion=0;iPion<num_npi;iPion++){
                if(status==1){

                    phi_diff = fabs(particle->GetPhi() - phi_npi[iPion]);
                    if(phi_diff>TMath::Pi())
                        phi_diff = TMath::TwoPi() - phi_diff;
                    
                    eta_diff = particle->GetEta() - eta_npi[iPion];
                    R = sqrt( pow(phi_diff,3) + pow(eta_diff,2) );
                    
                    if(R<0.7 && iParticle!=index_npi[iPion]){
                        count_npi[iPion]++;
                    }

                    if( iParticle!=index_npi[iPion] && fabs(particle->GetEta())<3.5 ){
                        pxsum_npi[iPion]+=particle->GetPx();
                        pysum_npi[iPion]+=particle->GetPy();
                    }
                    if( iParticle!=index_npi[iPion] && fabs(particle->GetEta())<4.0 ){
                        pxsum_npi_a[iPion]+=particle->GetPx();
                        pysum_npi_a[iPion]+=particle->GetPy();
                    }
                }
            }

            //Total mulitiplicity for |eta|<3.5
            if(status==1 && fabs(particle->GetEta())<3.5 ){
                ntotal++;
            }

            //Total E-pz for event (particles with |eta|<3.5)
            if(status==1 && fabs(particle->GetEta())<3.5 ){
                empztot+=( particle->GetE()-particle->GetPz() );
            }
        }

        if(part_se!=NULL && part_se->GetP()>0.1 ){ //Only consider electrons that have P>0.1 GeV/c
            h_se_1[bin_se]->Fill(count_se);
            h_se_2[bin_se]->Fill(ntotal-1); //Total number of particles other than scattered electron

            phi_diff = fabs( part_se->GetPhi() - atan2(pysum_se,pxsum_se) );
            if(phi_diff>TMath::Pi())
                phi_diff = TMath::TwoPi() - phi_diff;
            h_se_3[bin_se]->Fill(phi_diff);

            //Filled here to use |eta|<3.5 phi_diff
            if( (empztot - (2.*5.))>-5 && phi_diff>2 ){
                h_se_6[bin_se]->Fill(part_se->GetP());
            }

            phi_diff = fabs( part_se->GetPhi() - atan2(pysum_se_a,pxsum_se_a) );
            if(phi_diff>TMath::Pi())
                phi_diff = TMath::TwoPi() - phi_diff;
            h_se_3a[bin_se]->Fill(phi_diff);

            h_se_4[bin_se]->Fill( part_se->GetPt() - sqrt(pxsum_se*pxsum_se - pysum_se*pysum_se) );
            h_se_4a[bin_se]->Fill( part_se->GetPt() - sqrt(pxsum_se_a*pxsum_se_a - pysum_se_a*pysum_se_a) );

            h_se_5[bin_se]->Fill( empztot - (2.*5.) );       
        }

        for(int iPion=0;iPion<num_npi;iPion++){
            if( p_npi[iPion]>0.1 ){ //Only consider pions that have P>0.1 GeV/c
                h_npi_1[bin_npi[iPion]]->Fill(count_npi[iPion]);
                h_npi_2[bin_npi[iPion]]->Fill(ntotal-1); //Total number of particles other than this specific pi-minus

                phi_diff = fabs( phi_npi[iPion] - atan2(pysum_npi[iPion],pxsum_npi[iPion]) );
                if(phi_diff>TMath::Pi())
                    phi_diff = TMath::TwoPi() - phi_diff;
                h_npi_3[bin_npi[iPion]]->Fill(phi_diff);

                //Filled here to use |eta|<3.5 phi_diff
                if( (empztot - (2.*5.))>-5 && phi_diff>2 ){
                    h_npi_6[bin_npi[iPion]]->Fill(p_npi[iPion]);
                }

                phi_diff = fabs( phi_npi[iPion] - atan2(pysum_npi_a[iPion],pxsum_npi_a[iPion]) );
                if(phi_diff>TMath::Pi())
                    phi_diff = TMath::TwoPi() - phi_diff;
                h_npi_3a[bin_npi[iPion]]->Fill(phi_diff);

                h_npi_4[bin_npi[iPion]]->Fill( pt_npi[iPion] - sqrt(pxsum_npi[iPion]*pxsum_npi[iPion] - pysum_npi[iPion]*pysum_npi[iPion]) );
                h_npi_4a[bin_npi[iPion]]->Fill( pt_npi[iPion] - sqrt(pxsum_npi_a[iPion]*pxsum_npi_a[iPion] - pysum_npi_a[iPion]*pysum_npi_a[iPion]) );

                h_npi_5[bin_npi[iPion]]->Fill( empztot - (2.*5.) );
            }
        }

    }

    //Scale down pion momentum histograms by 10^4 suppression
    TH1* h_npi_sup[nEta];
    TH1* h_npi_6_sup[nEta];

    for(int ihist=0;ihist<nEta;ihist++){
        h_npi_sup[ihist] = (TH1*) h_npi[ihist]->Clone(Form("h_npi_sup[%d]",ihist));
        h_npi_sup[ihist]->SetLineColor(kMagenta);

        h_npi_sup[ihist]->Scale(1E-4);

        h_npi_6_sup[ihist] = (TH1*) h_npi_6[ihist]->Clone(Form("h_npi_6_sup[%d]",ihist));
        h_npi_6_sup[ihist]->SetLineColor(kMagenta);

        h_npi_6_sup[ihist]->Scale(1E-4);
    }

    TH1* h_npi_fine_sup[nEta_fine];

    for(int ihist=0;ihist<nEta_fine;ihist++){
        h_npi_fine_sup[ihist] = (TH1*) h_npi_fine[ihist]->Clone(Form("h_npi_fine_sup[%d]",ihist));
        h_npi_fine_sup[ihist]->SetLineColor(kMagenta);

        h_npi_fine_sup[ihist]->Scale(1E-4);
    }

    TH1* h_npi_finer_sup[nEta_finer];
    TH1* h_npi_finer_sup_a[nEta_finer];

    for(int ihist=0;ihist<nEta_finer;ihist++){
        h_npi_finer_sup[ihist] = (TH1*) h_npi_finer[ihist]->Clone(Form("h_npi_finer_sup[%d]",ihist));
        h_npi_finer_sup[ihist]->SetLineColor(kMagenta);
        h_npi_finer_sup[ihist]->Scale(1E-4);

        h_npi_finer_sup_a[ihist] = (TH1*) h_npi_finer[ihist]->Clone(Form("h_npi_finer_sup_a[%d]",ihist));
        h_npi_finer_sup_a[ihist]->SetLineColor(kRed);
        h_npi_finer_sup_a[ihist]->Scale(1E-5);
    }

    //Create vertical lines at minimum momentum values --
    //defined as electron minimum electron momentum in that eta range... 
    //satifying both Q2>1 and y<0.95
    double mom_low[] = {2.780,0.419,0.325,0.500,2.097,13.90};
    TLine *linea[nEta];

    for(int iline=0;iline<nEta;iline++){
        linea[iline] = new TLine(mom_low[iline],1E-1,mom_low[iline],2E6);
        linea[iline]->SetLineColor(kOrange);
        linea[iline]->SetLineWidth(2);
        linea[iline]->SetLineStyle(2);
    }

    double mom_low_fine[] = {0.,2.789,1.054,0.419,0.325,0.342}; //First one is unphysical
    TLine *linea_fine[nEta_fine];

    for(int iline=0;iline<nEta_fine;iline++){
        linea_fine[iline] = new TLine(mom_low_fine[iline],1E-1,mom_low_fine[iline],2E6);
        linea_fine[iline]->SetLineColor(kOrange);
        linea_fine[iline]->SetLineWidth(2);
        linea_fine[iline]->SetLineStyle(2);
    }

    double mom_low_finer[] = {0.325,0.325,0.342};
    TLine *linea_finer[nEta_finer];

    for(int iline=0;iline<nEta_finer;iline++){
        linea_finer[iline] = new TLine(mom_low_finer[iline],1E-1,mom_low_finer[iline],2E6);
        linea_finer[iline]->SetLineColor(kOrange);
        linea_finer[iline]->SetLineWidth(2);
        linea_finer[iline]->SetLineStyle(2);
    }

    //Draw plots
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.055);
    gStyle->SetTitleXOffset(0.85);
    gStyle->SetTitleYSize(0.055);
    gStyle->SetTitleYOffset(0.85);

    TCanvas *c1 = new TCanvas("c1");
    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();
    TH1 *hframe1 = gPad->DrawFrame(5E-4,1E-13,1.1,1E4);
    hframe1->SetTitle("5 GeV e^{-} on 41 GeV p");
    hframe1->GetXaxis()->SetTitle("y");hframe1->GetXaxis()->CenterTitle();
    hframe1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");hframe1->GetYaxis()->CenterTitle();
    gr1->Draw("P Same");

    c1->Print("plots/raw_yields_norad.pdf[");
    c1->Print("plots/raw_yields_norad.pdf");

    //Redraw with limit functions
    TF1 *f1a = new TF1("f1a",W_limit_y,1E-4,1,0);
    f1a->SetLineColor(kRed);f1a->SetLineWidth(3);
    f1a->Draw("Same");

    TF1 *f1b = new TF1("f1b",photo_limit_y,1E-4,1,0);
    f1b->SetLineColor(kGreen);f1b->SetLineWidth(3);
    f1b->Draw("same");

    TLatex *tex1a = new TLatex(1E-3,1E2,"W > 2 GeV Limit");
    tex1a->SetTextColor(kRed);tex1a->SetTextFont(42);
    tex1a->Draw();

    TLatex *tex1b = new TLatex(5E-2,1E-11,"Q^{2} > #frac{m_{e}^{2} y^{2}}{1-y} Limit");
    tex1b->SetTextColor(kGreen);tex1b->SetTextFont(42);
    tex1b->Draw();
    c1->Print("plots/raw_yields_norad.pdf");

    TCanvas *c2 = new TCanvas("c2");
    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();
    TH1 *hframe2 = gPad->DrawFrame(1E-13,1E-13,1.1,1E4);
    hframe2->SetTitle("5 GeV e^{-} on 41 GeV p");
    hframe2->GetXaxis()->SetTitle("x");hframe2->GetXaxis()->CenterTitle();
    hframe2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");hframe2->GetYaxis()->CenterTitle();
    gr2->Draw("P Same");

    c2->Print("plots/raw_yields_norad.pdf");

    //Redraw with limit functions
    TF1 *f2a = new TF1("f2a",W_limit_x,1E-14,1,0);
    f2a->SetLineColor(kRed);f2a->SetLineWidth(3);
    f2a->Draw("Same");

    TF1 *f2b = new TF1("f1b",photo_limit_x,1E-14,1,0);
    f2b->SetLineColor(kGreen);f2b->SetLineWidth(3);
    f2b->Draw("same");

    TLatex *tex2a = new TLatex(1E-5,1E-6,"W > 2 GeV Limit");
    tex2a->SetTextColor(kRed);tex2a->SetTextFont(42);
    tex2a->Draw();

    TLatex *tex2b = new TLatex(1E-10,1E1,"Q^{2} < #frac{x^{2} s^{2}}{m_{e}^{2} + xs} (Q^{2} > #frac{m_{e}^{2} y^{2}}{1-y}) Limit");
    tex2b->SetTextColor(kGreen);tex2b->SetTextFont(42);
    tex2b->Draw();

    c2->Print("plots/raw_yields_norad.pdf");
    
    TCanvas *c3 = new TCanvas("c3");
    c3->Divide(3,2);

    TH1 *hframe3[nEta];

    TLatex *tex3a = new TLatex(0.15,4E6,"5 GeV e^{-} on 41 GeV p");
    tex3a->SetTextColor(kBlack);tex3a->SetTextFont(42);

    TLatex *tex3b = new TLatex(2,5E6,"Scattered Electron");
    tex3b->SetTextColor(kGreen);tex3b->SetTextFont(42);

    TLatex *tex3c = new TLatex(2,1.5E6,"Negative Pions");
    tex3c->SetTextColor(kBlue);tex3c->SetTextFont(42);
    
    for(int iCan=0;iCan<nEta;iCan++){
        c3->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3[iCan] = gPad->DrawFrame(0.1,0.1,50,2E7);
        hframe3[iCan]->GetXaxis()->SetTitle("Momentum [GeV/c]");hframe3[iCan]->GetXaxis()->CenterTitle();
        hframe3[iCan]->GetYaxis()->SetTitle("Yield");hframe3[iCan]->GetYaxis()->CenterTitle();
        hframe3[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_low[iCan],Eta_hi[iCan]));
        
        h_se[iCan]->Draw("same");
        h_npi[iCan]->Draw("same");
        
        if(iCan==0){
            tex3a->Draw();
        }
        if(iCan==2){
            tex3b->Draw();tex3c->Draw();
        }
    }
    c3->Print("plots/raw_yields_norad.pdf");

    TCanvas *c3a = new TCanvas("c3a");
    c3a->Divide(3,2);

    TLatex *tex3d = new TLatex(2,5E5,"w/10^{4} suppression");
    tex3d->SetTextColor(kMagenta);tex3d->SetTextFont(42);

    TLatex *tex3e = new TLatex(2,1E6,"P_{min.} Value");
    tex3e->SetTextColor(kOrange);tex3e->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c3a->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3[iCan]->Draw();
        h_se[iCan]->Draw("same");
        h_npi[iCan]->Draw("same");
        h_npi_sup[iCan]->Draw("hist same");

        linea[iCan]->Draw();

        if(iCan==0){
            tex3a->Draw();
        }
        if(iCan==2){
            tex3b->Draw();tex3c->Draw();tex3d->Draw();
        }
        if(iCan==3){
            tex3e->Draw();
        }
    }
    c3a->Print("plots/raw_yields_norad.pdf");

    TCanvas *c3b = new TCanvas("c3b");
    c3b->Divide(3,2);

    TH1 *hframe3b[nEta_fine];

    for(int iCan=0;iCan<nEta_fine;iCan++){
        c3b->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3b[iCan] = gPad->DrawFrame(0.1,0.1,50,2E7);
        hframe3b[iCan]->GetXaxis()->SetTitle("Momentum [GeV/c]");hframe3b[iCan]->GetXaxis()->CenterTitle();
        hframe3b[iCan]->GetYaxis()->SetTitle("Yield");hframe3b[iCan]->GetYaxis()->CenterTitle();
        hframe3b[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_fine_low[iCan],Eta_fine_hi[iCan]));

        h_se_fine[iCan]->Draw("same");
        h_npi_fine[iCan]->Draw("same");

        if(iCan==0){
            tex3a->Draw();
        }
        if(iCan==2){
            tex3b->Draw();tex3c->Draw();
        }
    }
    c3b->Print("plots/raw_yields_norad.pdf");

    TCanvas *c3c = new TCanvas("c3c");
    c3c->Divide(3,2);
    
    for(int iCan=0;iCan<nEta_fine;iCan++){
        c3c->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3b[iCan]->Draw();
        h_se_fine[iCan]->Draw("same");
        h_npi_fine[iCan]->Draw("same");
        h_npi_fine_sup[iCan]->Draw("hist same");

        if(iCan!=0)linea_fine[iCan]->Draw();

        if(iCan==0){
            tex3a->Draw();
        }
        if(iCan==2){
            tex3b->Draw();tex3c->Draw();tex3d->Draw();
        }
        if(iCan==3){
            tex3e->Draw();
        }
    }
    c3c->Print("plots/raw_yields_norad.pdf");

    TCanvas *c3d = new TCanvas("c3d");
    c3d->Divide(3,2);

    TH1 *hframe3d[nEta_finer];

    TLatex *tex3f = new TLatex(2,2E5,"+DIRC suppression");
    tex3f->SetTextColor(kRed);tex3f->SetTextFont(42);

    for(int iCan=0;iCan<nEta_finer;iCan++){
        c3d->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3d[iCan] = gPad->DrawFrame(0.1,0.1,50,2E7);
        hframe3d[iCan]->GetXaxis()->SetTitle("Momentum [GeV/c]");hframe3d[iCan]->GetXaxis()->CenterTitle();
        hframe3d[iCan]->GetYaxis()->SetTitle("Yield");hframe3d[iCan]->GetYaxis()->CenterTitle();
        hframe3d[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_finer_low[iCan],Eta_finer_hi[iCan]));

        h_se_finer[iCan]->Draw("same");
        h_npi_finer[iCan]->Draw("same");

        if(iCan==0){
            tex3a->Draw();
        }
        if(iCan==2){
            tex3b->Draw();tex3c->Draw();
        }
    }
    for(int iCan=0;iCan<nEta_finer;iCan++){
        c3d->cd(iCan+4);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3d[iCan]->Draw();
        h_se_finer[iCan]->Draw("same");
        h_npi_finer[iCan]->Draw("same");
        h_npi_finer_sup[iCan]->Draw("hist same");
        h_npi_finer_sup_a[iCan]->Draw("hist same");

        linea_finer[iCan]->Draw();

        if(iCan==0){
            tex3a->Draw();tex3e->Draw();
        }
        if(iCan==2){
            tex3b->Draw();tex3c->Draw();tex3d->Draw();tex3f->Draw();
        }
    }
    c3d->Print("plots/raw_yields_norad.pdf");

    TCanvas *c4 = new TCanvas("c4");
    c4->Divide(3,2);

    TLatex *tex4a = new TLatex(2,5E6,"All Electrons");
    tex4a->SetTextColor(kBlack);tex4a->SetTextFont(42);

    TLatex *tex4b = new TLatex(2,1.5E6,"All Positrons");
    tex4b->SetTextColor(kRed);tex4b->SetTextFont(42);
    
    TLatex *tex4c = new TLatex(2,6E5,"Negative Pions");
    tex4c->SetTextColor(kBlue);tex4c->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c4->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();
        
        hframe3[iCan]->Draw();
        h_ae[iCan]->Draw("same");
        h_pos[iCan]->Draw("same");
        h_npi[iCan]->Draw("same");

        if(iCan==0){
            tex3a->Draw();
        }
        if(iCan==2){
            tex4a->Draw();tex4b->Draw();tex4c->Draw();
        }    
    }
    c4->Print("plots/raw_yields_norad.pdf");

    TCanvas *c5 = new TCanvas("c5");
    c5->Divide(3,2);
    TH1 *hframe5[nEta];

    TLatex *tex5a = new TLatex(0.5,3E7,"5 GeV e^{-} on 41 GeV p");
    tex5a->SetTextColor(kBlack);tex5a->SetTextFont(42);

    TLatex *tex5b = new TLatex(1.5,3E7,"Scattered Electron, p>0.1 GeV/c");
    tex5b->SetTextColor(kGreen);tex5b->SetTextFont(42);

    TLatex *tex5c = new TLatex(1.5,8E6,"Negative Pions, p>0.1 GeV/c");
    tex5c->SetTextColor(kBlue);tex5c->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c5->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe5[iCan] = gPad->DrawFrame(0,0.1,10,2E8);
        hframe5[iCan]->GetXaxis()->SetTitle("Number of Particles with R<0.7");hframe5[iCan]->GetXaxis()->CenterTitle();
        hframe5[iCan]->GetXaxis()->CenterLabels();
        hframe5[iCan]->GetYaxis()->SetTitle("Yield");hframe5[iCan]->GetYaxis()->CenterTitle();
        hframe5[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_low[iCan],Eta_hi[iCan]));

        h_se_1[iCan]->Draw("same");
        h_npi_1[iCan]->Draw("same");

        if(iCan==0){
            tex5a->Draw();
        }
        if(iCan==2){
            tex5b->Draw();tex5c->Draw();
        }
    }
    c5->Print("plots/raw_yields_norad.pdf");

    TCanvas *c6 = new TCanvas("c6");
    c6->Divide(3,2);
    TH1 *hframe6[nEta];

    TLatex *tex6a = new TLatex(1,3E7,"5 GeV e^{-} on 41 GeV p");
    tex6a->SetTextColor(kBlack);tex6a->SetTextFont(42);

    TLatex *tex6b = new TLatex(3,3E7,"Scattered Electron, p>0.1 GeV/c");
    tex6b->SetTextColor(kGreen);tex6b->SetTextFont(42);

    TLatex *tex6c = new TLatex(3,8E6,"Negative Pions, p>0.1 GeV/c");
    tex6c->SetTextColor(kBlue);tex6c->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c6->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe6[iCan] = gPad->DrawFrame(0,0.1,20,2E8);
        hframe6[iCan]->GetXaxis()->SetTitle("Number of additional particles w/|#eta|<3.5");hframe6[iCan]->GetXaxis()->CenterTitle();
        hframe6[iCan]->GetYaxis()->SetTitle("Yield");hframe6[iCan]->GetYaxis()->CenterTitle();
        hframe6[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_low[iCan],Eta_hi[iCan]));

        h_se_2[iCan]->Draw("same");
        h_npi_2[iCan]->Draw("same");

        if(iCan==0){
            tex6a->Draw();
        }
        if(iCan==2){
            tex6b->Draw();tex6c->Draw();
        }
    }
    c6->Print("plots/raw_yields_norad.pdf");

    TCanvas *c7 = new TCanvas("c7");
    c7->Divide(3,2);
    TH1 *hframe7[nEta];

    TLatex *tex7a = new TLatex(0.1,3E7,"5 GeV e^{-} on 41 GeV p");
    tex7a->SetTextColor(kBlack);tex7a->SetTextFont(42);

    TLatex *tex7b = new TLatex(0.4,3E7,"Scattered Electron, p>0.1 GeV/c");
    tex7b->SetTextColor(kGreen);tex7b->SetTextFont(42);

    TLatex *tex7c = new TLatex(0.4,9E6,"Negative Pions, p>0.1 GeV/c");
    tex7c->SetTextColor(kBlue);tex7c->SetTextFont(42);

    TLatex *tex7d = new TLatex(0.4,2.5E6,"Compare particles w/|#eta|<3.5");
    tex7d->SetTextColor(kBlack);tex7d->SetTextFont(42);

    TLatex *tex7e = new TLatex(0.4,2.5E6,"Compare particles w/|#eta|<4.0");
    tex7e->SetTextColor(kBlack);tex7e->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c7->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe7[iCan] = gPad->DrawFrame(-0.1,0.1,3.25,2E8);
        hframe7[iCan]->GetXaxis()->SetTitle("#phi Difference");hframe7[iCan]->GetXaxis()->CenterTitle();
        hframe7[iCan]->GetYaxis()->SetTitle("Yield");hframe7[iCan]->GetYaxis()->CenterTitle();
        hframe7[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_low[iCan],Eta_hi[iCan]));

        h_se_3[iCan]->Draw("same");
        h_npi_3[iCan]->Draw("same");

        if(iCan==0){
            tex7a->Draw();
        }
        if(iCan==2){
            tex7b->Draw();tex7c->Draw();tex7d->Draw();
        }
    }
    c7->Print("plots/raw_yields_norad.pdf");

    TCanvas *c7a = new TCanvas("c7a");
    c7a->Divide(3,2);

    for(int iCan=0;iCan<nEta;iCan++){
        c7a->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe7[iCan]->Draw();

        h_se_3a[iCan]->Draw("same");
        h_npi_3a[iCan]->Draw("same");

        if(iCan==0){
            tex7a->Draw();
        }
        if(iCan==2){
            tex7b->Draw();tex7c->Draw();tex7e->Draw();
        }
    }
    c7a->Print("plots/raw_yields_norad.pdf");

    TCanvas *c8 = new TCanvas("c8");
    c8->Divide(3,2);
    TH1 *hframe8[nEta];

    TLatex *tex8a = new TLatex(-2.8,3E7,"5 GeV e^{-} on 41 GeV p");
    tex8a->SetTextColor(kBlack);tex8a->SetTextFont(42);

    TLatex *tex8b = new TLatex(-2.2,3E7,"Scattered Electron, p>0.1 GeV/c");
    tex8b->SetTextColor(kGreen);tex8b->SetTextFont(42);

    TLatex *tex8c = new TLatex(-2.2,9E6,"Negative Pions, p>0.1 GeV/c");
    tex8c->SetTextColor(kBlue);tex8c->SetTextFont(42);

    TLatex *tex8d = new TLatex(-2.2,2.5E6,"Compare particles w/|#eta|<3.5");
    tex8d->SetTextColor(kBlack);tex8d->SetTextFont(42);

    TLatex *tex8e = new TLatex(-2.2,2.5E6,"Compare particles w/|#eta|<4.0");
    tex8e->SetTextColor(kBlack);tex8e->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c8->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe8[iCan] = gPad->DrawFrame(-3,0.1,3,2E8);
        hframe8[iCan]->GetXaxis()->SetTitle("|P_{T}| Difference");hframe8[iCan]->GetXaxis()->CenterTitle();
        hframe8[iCan]->GetYaxis()->SetTitle("Yield");hframe8[iCan]->GetYaxis()->CenterTitle();
        hframe8[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_low[iCan],Eta_hi[iCan]));

        h_se_4[iCan]->Draw("same");
        h_npi_4[iCan]->Draw("same");

        if(iCan==0){
            tex8a->Draw();
        }
        if(iCan==2){
            tex8b->Draw();tex8c->Draw();tex8d->Draw();
        }
    }
    c8->Print("plots/raw_yields_norad.pdf");

    TCanvas *c8a = new TCanvas("c8a");
    c8a->Divide(3,2);

    for(int iCan=0;iCan<nEta;iCan++){
        c8a->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe8[iCan]->Draw();

        h_se_4a[iCan]->Draw("same");
        h_npi_4a[iCan]->Draw("same");

        if(iCan==0){
            tex8a->Draw();
        }
        if(iCan==2){
            tex8b->Draw();tex8c->Draw();tex8e->Draw();
        }
    }
    c8a->Print("plots/raw_yields_norad.pdf");

    TCanvas *c9 = new TCanvas("c9");
    c9->Divide(3,2);
    TH1 *hframe9[nEta];

    TLatex *tex9a = new TLatex(-20,3E7,"5 GeV e^{-} on 41 GeV p");
    tex9a->SetTextColor(kBlack);tex9a->SetTextFont(42);

    TLatex *tex9b = new TLatex(-20,3E7,"Scattered Electron, p>0.1 GeV/c");
    tex9b->SetTextColor(kGreen);tex9b->SetTextFont(42);

    TLatex *tex9c = new TLatex(-20,9E6,"Negative Pions, p>0.1 GeV/c");
    tex9c->SetTextColor(kBlue);tex9c->SetTextFont(42);

    TLatex *tex9d = new TLatex(-20,2.5E6,"Sum particles w/|#eta|<3.5");
    tex9d->SetTextColor(kBlack);tex9d->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c9->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe9[iCan] = gPad->DrawFrame(-11,0.1,1,2E8);
        hframe9[iCan]->GetXaxis()->SetTitle("(E - P_{z})^{Tot.} - 2E_{e}");hframe9[iCan]->GetXaxis()->CenterTitle();
        hframe9[iCan]->GetXaxis()->SetTitleSize(0.045);
        hframe9[iCan]->GetYaxis()->SetTitle("Yield");hframe9[iCan]->GetYaxis()->CenterTitle();
        hframe9[iCan]->SetTitle(Form("%0.1f < #eta < %0.1f",Eta_low[iCan],Eta_hi[iCan]));

        h_se_5[iCan]->Draw("same");
        h_npi_5[iCan]->Draw("same");

        if(iCan==0){
            tex9a->Draw();
        }
        if(iCan==2){
            tex9b->Draw();tex9c->Draw();tex9d->Draw();
        }
    }
    c9->Print("plots/raw_yields_norad.pdf");

    TCanvas *c10 = new TCanvas("c10");
    c10->Divide(3,2);

    TLatex *tex10a = new TLatex(2,1.5E5,"all w/kin. cuts");
    tex10a->SetTextColor(kBlack);tex10a->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c10->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3[iCan]->Draw();
        h_se_6[iCan]->Draw("same");
        h_npi_6[iCan]->Draw("same");
        h_npi_6_sup[iCan]->Draw("hist same");

        linea[iCan]->Draw();
        
        if(iCan==0){
            tex3a->Draw();
        }
        if(iCan==2){
            tex3b->Draw();tex3c->Draw();tex3d->Draw();tex10a->Draw();
        }
        if(iCan==3){
            tex3e->Draw();
        }
    }
    c10->Print("plots/raw_yields_norad.pdf");
    c10->Print("plots/raw_yields_norad.pdf]");

    myapp->Run();
    return 0;
}
