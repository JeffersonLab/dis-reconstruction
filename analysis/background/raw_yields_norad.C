R__LOAD_LIBRARY(libeicsmear);

//Globals
const double Mp(0.9383);
const double Me(0.511E-3);
const double s_cm = 4.*10.*100; //beam energy combination

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
void raw_yields_norad(){

    //Define TGraphs and Histograms
    TGraph *gr1 = new TGraph(); //Q2 vs. y
    TGraph *gr2 = new TGraph(); //Q2 vs. x

    int nEta = 6; //Number of eta bins
    double Eta_low[] = {-3.5, -2.0, -1.0, 0.0, 1.0, 2.0};
    double Eta_hi[] =  {-2.0, -1.0,  0.0, 1.0, 2.0, 3.5};

    TH1* h_se[nEta];
    TH1* h_ae[nEta];
    TH1* h_pos[nEta];
    TH1* h_npi[nEta];

    //Four 'topological' distributions for scattered electron and negative pions
    //1) Number of particles with R<0.7
    //2) Number of additional particles that have |eta|<3.5
    //3) Phi difference with sum of all particles that have |eta|<3.5
    //3a) Phi difference with sum of all particles that have |eta|<4.0
    //4) |Pt| difference with sum of all particles that have |eta|<3.5
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
    }

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class
    erhic::ParticleMC *part_se(NULL);

    TChain *t = new TChain("EICTree");
    for(int i=0;i<15;i++){
        t->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/other_studies/10_100/ep_minbias_%d.root",i));
    }

    t->SetBranchAddress("event",&event);

    //We need to use the 'true' kinematic quantities given by pythia (e.g. trueQ2)
    //since the electron ones (e.g. QSquared) seem to be calculated assuming 
    //zero electron mass -- a bad assumption at very low Q2. At higher Q2, both
    //sets of quantities are the same when the QED effects are turned OFF.
    //When using the true quantities, it seems like they have any accuracy limit
    //of 10^-11; so events generated at lower x or Q2 will be put at the limiting value.
    //This is probably be just a problem with the writing out of those variables, rather 
    //than an issue with the event generation for the very low x and Q2 events.
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

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    double cross_tot = 125.8E9; //Total Cross Section in fb
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
                    R = sqrt( pow(phi_diff,2) + pow(eta_diff,2) );
                    
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

            //Total mulitiplicity for eta<|3.5|
            if(status==1 && fabs(particle->GetEta())<3.5 )
                ntotal++;
        }

        if(part_se!=NULL && part_se->GetP()>0.1 ){ //Only consider electrons that have P>0.1 GeV/c
            h_se_1[bin_se]->Fill(count_se);
            h_se_2[bin_se]->Fill(ntotal-1); //Total number of particles other than scattered electron

            phi_diff = fabs( part_se->GetPhi() - atan2(pysum_se,pxsum_se) );
            if(phi_diff>TMath::Pi())
                phi_diff = TMath::TwoPi() - phi_diff;
            h_se_3[bin_se]->Fill(phi_diff);

            phi_diff = fabs( part_se->GetPhi() - atan2(pysum_se_a,pxsum_se_a) );
            if(phi_diff>TMath::Pi())
                phi_diff = TMath::TwoPi() - phi_diff;
            h_se_3a[bin_se]->Fill(phi_diff);

            h_se_4[bin_se]->Fill( part_se->GetPt() - sqrt(pxsum_se*pxsum_se - pysum_se*pysum_se) );
            h_se_4a[bin_se]->Fill( part_se->GetPt() - sqrt(pxsum_se_a*pxsum_se_a - pysum_se_a*pysum_se_a) );
        }

        for(int iPion=0;iPion<num_npi;iPion++){
            if( p_npi[iPion]>0.1 ){ //Only consider pions that have P>0.1 GeV/c
                h_npi_1[bin_npi[iPion]]->Fill(count_npi[iPion]);
                h_npi_2[bin_npi[iPion]]->Fill(ntotal-1); //Total number of particles other than this specific pi-minus

                phi_diff = fabs( phi_npi[iPion] - atan2(pysum_npi[iPion],pxsum_npi[iPion]) );
                if(phi_diff>TMath::Pi())
                    phi_diff = TMath::TwoPi() - phi_diff;
                h_npi_3[bin_npi[iPion]]->Fill(phi_diff);

                phi_diff = fabs( phi_npi[iPion] - atan2(pysum_npi_a[iPion],pxsum_npi_a[iPion]) );
                if(phi_diff>TMath::Pi())
                    phi_diff = TMath::TwoPi() - phi_diff;
                h_npi_3a[bin_npi[iPion]]->Fill(phi_diff);

                h_npi_4[bin_npi[iPion]]->Fill( pt_npi[iPion] - sqrt(pxsum_npi[iPion]*pxsum_npi[iPion] - pysum_npi[iPion]*pysum_npi[iPion]) );
                h_npi_4a[bin_npi[iPion]]->Fill( pt_npi[iPion] - sqrt(pxsum_npi_a[iPion]*pxsum_npi_a[iPion] - pysum_npi_a[iPion]*pysum_npi_a[iPion]) );
            }
        }

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
    TH1 *hframe1 = gPad->DrawFrame(1E-4,1E-14,1.1,1E4);
    hframe1->SetTitle("10 GeV e^{-} on 100 GeV p");
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
    TH1 *hframe2 = gPad->DrawFrame(1E-14,1E-14,1.1,1E4);
    hframe2->SetTitle("10 GeV e^{-} on 100 GeV p");
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

    TLatex *tex2b = new TLatex(1E-10,1E2,"Q^{2} < #frac{x^{2} s^{2}}{M_{e}^{2} + xs} (Q^{2} > #frac{m_{e}^{2} y^{2}}{1-y}) Limit");
    tex2b->SetTextColor(kGreen);tex2b->SetTextFont(42);
    tex2b->Draw();

    c2->Print("plots/raw_yields_norad.pdf");
    
    TCanvas *c3 = new TCanvas("c3");
    c3->Divide(3,2);

    TH1 *hframe3[nEta];

    TLatex *tex3a = new TLatex(0.15,2E5,"10 GeV e^{-} on 100 GeV p");
    tex3a->SetTextColor(kBlack);tex3a->SetTextFont(42);

    TLatex *tex3b = new TLatex(2,3E5,"Scattered Electron");
    tex3b->SetTextColor(kGreen);tex3b->SetTextFont(42);

    TLatex *tex3c = new TLatex(2,1E5,"Negative Pions");
    tex3c->SetTextColor(kBlue);tex3c->SetTextFont(42);
    
    for(int iCan=0;iCan<nEta;iCan++){
        c3->cd(iCan+1);
        gPad->SetLogx();gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe3[iCan] = gPad->DrawFrame(0.1,0.1,50,1E6);
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

    TCanvas *c4 = new TCanvas("c4");
    c4->Divide(3,2);

    TLatex *tex4a = new TLatex(2,3E5,"All Electrons");
    tex4a->SetTextColor(kBlack);tex4a->SetTextFont(42);

    TLatex *tex4b = new TLatex(2,1E5,"All Positrons");
    tex4b->SetTextColor(kRed);tex4b->SetTextFont(42);
    
    TLatex *tex4c = new TLatex(2,3E4,"Negative Pions");
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

    TLatex *tex5a = new TLatex(0.5,2E6,"10 GeV e^{-} on 100 GeV p");
    tex5a->SetTextColor(kBlack);tex5a->SetTextFont(42);

    TLatex *tex5b = new TLatex(1.5,2E6,"Scattered Electron, p>0.1 GeV/c");
    tex5b->SetTextColor(kGreen);tex5b->SetTextFont(42);

    TLatex *tex5c = new TLatex(1.5,7E5,"Negative Pions, p>0.1 GeV/c");
    tex5c->SetTextColor(kBlue);tex5c->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c5->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe5[iCan] = gPad->DrawFrame(0,0.1,10,1E7);
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

    TLatex *tex6a = new TLatex(1,2E6,"10 GeV e^{-} on 100 GeV p");
    tex6a->SetTextColor(kBlack);tex6a->SetTextFont(42);

    TLatex *tex6b = new TLatex(3,2E6,"Scattered Electron, p>0.1 GeV/c");
    tex6b->SetTextColor(kGreen);tex6b->SetTextFont(42);

    TLatex *tex6c = new TLatex(3,7E5,"Negative Pions, p>0.1 GeV/c");
    tex6c->SetTextColor(kBlue);tex6c->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c6->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe6[iCan] = gPad->DrawFrame(0,0.1,20,1E7);
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

    TLatex *tex7a = new TLatex(0.1,2E6,"10 GeV e^{-} on 100 GeV p");
    tex7a->SetTextColor(kBlack);tex7a->SetTextFont(42);

    TLatex *tex7b = new TLatex(0.4,2E6,"Scattered Electron, p>0.1 GeV/c");
    tex7b->SetTextColor(kGreen);tex7b->SetTextFont(42);

    TLatex *tex7c = new TLatex(0.4,7E5,"Negative Pions, p>0.1 GeV/c");
    tex7c->SetTextColor(kBlue);tex7c->SetTextFont(42);

    TLatex *tex7d = new TLatex(0.4,2E5,"Compare particles w/|#eta|<3.5");
    tex7d->SetTextColor(kBlack);tex7d->SetTextFont(42);

    TLatex *tex7e = new TLatex(0.4,2E5,"Compare particles w/|#eta|<4.0");
    tex7e->SetTextColor(kBlack);tex7e->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c7->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe7[iCan] = gPad->DrawFrame(-0.1,0.1,3.25,1E7);
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

    TLatex *tex8a = new TLatex(-2.8,2E6,"10 GeV e^{-} on 100 GeV p");
    tex8a->SetTextColor(kBlack);tex8a->SetTextFont(42);

    TLatex *tex8b = new TLatex(-2.2,2E6,"Scattered Electron, p>0.1 GeV/c");
    tex8b->SetTextColor(kGreen);tex8b->SetTextFont(42);

    TLatex *tex8c = new TLatex(-2.2,7E5,"Negative Pions, p>0.1 GeV/c");
    tex8c->SetTextColor(kBlue);tex8c->SetTextFont(42);

    TLatex *tex8d = new TLatex(-2.2,2E5,"Compare particles w/|#eta|<3.5");
    tex8d->SetTextColor(kBlack);tex8d->SetTextFont(42);

    TLatex *tex8e = new TLatex(-2.2,2E5,"Compare particles w/|#eta|<4.0");
    tex8e->SetTextColor(kBlack);tex8e->SetTextFont(42);

    for(int iCan=0;iCan<nEta;iCan++){
        c8->cd(iCan+1);
        gPad->SetLogy();
        gPad->SetTickx();gPad->SetTicky();

        hframe8[iCan] = gPad->DrawFrame(-3,0.1,3,1E7);
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
    c8a->Print("plots/raw_yields_norad.pdf]");

}