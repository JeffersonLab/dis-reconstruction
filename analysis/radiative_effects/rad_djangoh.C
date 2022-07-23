R__LOAD_LIBRARY(libeicsmear);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void rad_djangoh(){

    //Constants
    const double alpha(1./137.036);
    const double s_cm = 4.*10*100;
    const double fbgev(1./(0.3894E12));
    const double Mp(0.9383);

    //Cross Section Bins
    double Q2_min = 1E-1;
    double Q2_max = 1E4;
	const int nbins_Q2 = 25;
	double log_bw_Q2 = (log10(Q2_max) - log10(Q2_min))/nbins_Q2;
    double log_low_Q2,log_hi_Q2;
    double Q2_low[nbins_Q2],Q2_hi[nbins_Q2];
    for(int i=0;i<nbins_Q2;i++){
		log_low_Q2 = log10(Q2_min) + (i*log_bw_Q2);
		log_hi_Q2 = log10(Q2_min) + (i+1)*log_bw_Q2;
		Q2_low[i] = pow(10,log_low_Q2);
		Q2_hi[i] = pow(10,log_hi_Q2);
	}
   
    double x_min = 1E-5;
    double x_max = 1;
	const int nbins = 25;
	double log_bw = (log10(x_max) - log10(x_min))/nbins;
    double log_low,log_hi;
    double x_low[nbins],x_hi[nbins];
    for(int i=0;i<nbins;i++){
		log_low = log10(x_min) + (i*log_bw);
		log_hi = log10(x_min) + (i+1)*log_bw;
		x_low[i] = pow(10,log_low);
		x_hi[i] = pow(10,log_hi);
	}

    //Define Radiative kinematics hisotgrams
    TH2 *h1 = new TH2D("h1","",140,0,70,140,0,70);
    h1->GetXaxis()->SetTitle("W_{True} [GeV]");h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("W_{Elec.} [GeV]");h1->GetYaxis()->CenterTitle();

    TH2 *h1a = new TH2D("h1a","Q^{2}_{True} < 50 GeV^{2}",140,0,70,140,0,70);
    h1a->GetXaxis()->SetTitle("W_{True} [GeV]");h1a->GetXaxis()->CenterTitle();
    h1a->GetYaxis()->SetTitle("W_{Elec.} [GeV]");h1a->GetYaxis()->CenterTitle();

    TH2 *h1b = new TH2D("h1b","Q^{2}_{True} > 50 GeV^{2}",140,0,70,140,0,70);
    h1b->GetXaxis()->SetTitle("W_{True} [GeV]");h1b->GetXaxis()->CenterTitle();
    h1b->GetYaxis()->SetTitle("W_{Elec.} [GeV]");h1b->GetYaxis()->CenterTitle();

    TH2 *h1c = new TH2D("h1c","W^{2}_{Elec.} > 10 GeV^{2}, y_{Elec.}>0.01",200,0,200,140,0,70);
    h1c->GetXaxis()->SetTitle("Q^{2}_{Elec.} [GeV^{2}]");h1c->GetXaxis()->CenterTitle();
    h1c->GetYaxis()->SetTitle("W_{True} [GeV]");h1c->GetYaxis()->CenterTitle();

    TH2 *h1d = new TH2D("h1d","",200,0,200,140,0,70);
    h1d->GetXaxis()->SetTitle("Q^{2}_{True.} [GeV^{2}]");h1d->GetXaxis()->CenterTitle();
    h1d->GetYaxis()->SetTitle("W_{True} [GeV]");h1d->GetYaxis()->CenterTitle();

    TH1 *h2 = new TH1D("h2","",100,0,10);
    h2->GetXaxis()->SetTitle("W_{True} [GeV]");h2->GetXaxis()->CenterTitle();
    h2->SetLineColor(kBlue);h2->SetLineWidth(2);

    TH1 *h2a = new TH1D("h2a","",200,-2,80);
    h2a->GetXaxis()->SetTitle("W_{True} [GeV]");h2a->GetXaxis()->CenterTitle();
    h2a->SetLineColor(kBlue);h2a->SetLineWidth(2);

    TH1 *h2b = new TH1D("h2b","",200,-2,80);
    h2b->GetXaxis()->SetTitle("W_{Elec.} [GeV]");h2b->GetXaxis()->CenterTitle();
    h2b->SetLineColor(kBlue);h2b->SetLineWidth(2);

    TH2 *h3 = new TH2D("h3","",100,0,100,100,0,100);
    h3->GetXaxis()->SetTitle("Q^{2}_{True} [GeV^{2}]");h3->GetXaxis()->CenterTitle();
    h3->GetYaxis()->SetTitle("Q^{2}_{Elec.} [GeV^{2}]");h3->GetYaxis()->CenterTitle();

    //Use (same) variable binning for x and y
    const int nbins_hist = 100;
    double hist_min = 1E-4; double hist_max = 1;
    double hist_log_bw = (log10(hist_max) - log10(hist_min))/nbins_hist;
    double hist_bin[nbins_hist+1];

    for(int i=0;i<=nbins_hist;i++){
		double log_hist = log10(hist_min) + (i*hist_log_bw);
		hist_bin[i] = pow(10,log_hist);
    }

    TH2 *h4 = new TH2D("h4","",nbins_hist,hist_bin,nbins_hist,hist_bin);
    h4->GetXaxis()->SetTitle("y_{True}");h4->GetXaxis()->CenterTitle();
    h4->GetYaxis()->SetTitle("y_{Elec.}");h4->GetYaxis()->CenterTitle();

    TH2 *h5 = new TH2D("h5","",nbins_hist,hist_bin,nbins_hist,hist_bin);
    h5->GetXaxis()->SetTitle("x_{True}");h5->GetXaxis()->CenterTitle();
    h5->GetYaxis()->SetTitle("x_{Elec.}");h5->GetYaxis()->CenterTitle();

    TH2 *h5a = new TH2D("h5a","Q^{2}_{True} < 50 GeV^{2}",nbins_hist,hist_bin,nbins_hist,hist_bin);
    h5a->GetXaxis()->SetTitle("x_{True}");h5a->GetXaxis()->CenterTitle();
    h5a->GetYaxis()->SetTitle("x_{Elec.}");h5a->GetYaxis()->CenterTitle();

    TH2 *h5b = new TH2D("h5b","Q^{2}_{True} > 50 GeV^{2}",nbins_hist,hist_bin,nbins_hist,hist_bin);
    h5b->GetXaxis()->SetTitle("x_{True}");h5b->GetXaxis()->CenterTitle();
    h5b->GetYaxis()->SetTitle("x_{Elec.}");h5b->GetYaxis()->CenterTitle();

    TH1 *h6 = new TH1D("h6","",10,0,10);
    h6->GetXaxis()->SetTitle("IChannel");h6->GetXaxis()->CenterTitle();
    h6->GetXaxis()->CenterLabels();
    h6->SetLineColor(kBlue);h6->SetLineWidth(2);

    TH1 *h7a = new TH1D("h7a","Sum over final-state particles w/ -4<#eta<4",100,-10,22); //Channel 1
    h7a->GetXaxis()->SetTitle("(E-p_{z})_{Tot.}");h7a->GetXaxis()->CenterTitle();
    h7a->SetLineColor(kBlack);h7a->SetLineWidth(2);

    TH1 *h7b = new TH1D("h7b","",100,-10,22); //Channel 6
    h7b->GetXaxis()->SetTitle("(E-p_{z})_{Tot.}");h7b->GetXaxis()->CenterTitle();
    h7b->SetLineColor(kBlue);h7b->SetLineWidth(2);

    TH1 *h7c = new TH1D("h7c","",100,-10,22); //Channel 7
    h7c->GetXaxis()->SetTitle("(E-p_{z})_{Tot.}");h7c->GetXaxis()->CenterTitle();
    h7c->SetLineColor(kGreen);h7c->SetLineWidth(2);

    TH1 *h7d = new TH1D("h7d","",100,-10,22); //Channel 8
    h7d->GetXaxis()->SetTitle("(E-p_{z})_{Tot.}");h7d->GetXaxis()->CenterTitle();
    h7d->SetLineColor(kRed);h7d->SetLineWidth(2);

    //Bin Yields/Kinematics/Factors
    //int nbins_tot = nbins_Q2*nbins;
    double yield[nbins][nbins_Q2];
    double yield_true[nbins][nbins_Q2];
    double error[nbins][nbins_Q2];
    double rcs[nbins][nbins_Q2];
    double y[nbins][nbins_Q2];
    double Q2_center[nbins_Q2];
    double x_center[nbins];
    double Q2_width[nbins_Q2];
    double x_width[nbins];
    double rcs_factor[nbins][nbins_Q2];
    bool cut_ul[nbins][nbins_Q2];
    bool cut_lr[nbins][nbins_Q2];
    double ymax = 1;double ymin = 0.001; double W2min = 10;
    double y_temp, W2_temp;

    for(int i=0;i<nbins;i++){ //x loop
        for(int j=0;j<nbins_Q2;j++){ //Q2 loop
            yield[i][j]=0;
            yield_true[i][j]=0;
            error[i][j]=0;
            rcs[i][j]=0;  

            Q2_center[j] = (Q2_hi[j]+Q2_low[j])/2.;
            x_center[i] = (x_hi[i]+x_low[i])/2.;
            Q2_width[j] = (Q2_hi[j]-Q2_low[j]);
            x_width[i] = (x_hi[i]-x_low[i]);

            y[i][j] = Q2_center[j]/(x_center[i]*s_cm);
            rcs_factor[i][j] = fbgev * ( pow(Q2_center[j],2)*x_center[i] ) / 
                            ( TMath::TwoPi() * pow(alpha,2) * ( 1.+pow((1.-y[i][j]),2) ) );

            //Cut Bins with top-left edge > ymax and bottom-right edge < ymin and > W2min
            y_temp = (Q2_hi[j])/(x_low[i]*s_cm); //Upper left
            if(y_temp>ymax) //Remove bin
                cut_ul[i][j] = 0;
            else
                cut_ul[i][j] = 1;

            y_temp = (Q2_low[j])/(x_hi[i]*s_cm); //Lower Right
            W2_temp = Mp*Mp + Q2_low[j]*( (1./x_hi[i]) - 1.);
            if(y_temp<ymin || W2_temp<W2min) //Remove bin
                cut_lr[i][j] = 0;
            else
                cut_lr[i][j] = 1;

        }
    }

    //Q2 bins (center) to plot; 8 bins
    const int num_plot = 8;
    int Q2_bin_plot[num_plot] = {6,8,9,10,11,13,15,18};

    //-----------------------------------
    //Analyse DJANGOH with QED radiation
    //-----------------------------------

    //Load ROOT Files
    erhic::EventDjangoh *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *t = new TChain("EICTree");
    
    //Djangoh 4.6.20
    t->Add("/gpfs02/eic/baraks/djangoh/djangoh_local/outfiles/djangoh.NC.Rad.10x100_evt.root");

    //Djangoh 4.6.10
    //for(int i=0;i<15;i++){
        //t->Add(Form("/eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/10_100/Rad/djangoh.NC.10x100_%d.root",i));
    //}

    t->SetBranchAddress("event",&event);
    double Q2_event; //Get this from scattered electron
    double x_event; //''
    double y_event; //''
    double W2_event; //''
    double Q2_true; //Get this from exchanged virtual boson
    double x_true; //''
    double y_true; //''
    double W2_true; //''

    // There is some problem with the calculation of nu_true. 
    // The scattered electron one (nu_event) is calculated fine.
    //-----------------
    //double nu_event;
    //double nu_true;

    //Other variables
    double Eout1,pzout1;
    double Eout6,pzout6;
    double Eout7,pzout7;
    double Eout8,pzout8;

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    
    //Djangoh 4.6.20
    double cross_tot = 1.51E9; //Total Cross Section in fb
    //Djangoh 4.6.10
    //double cross_tot = 1.25E9; //Total Cross Section in fb
    
    double lum = ( (double) nevents)/cross_tot; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"DJANGOH (QED Rad ON) Simulation:"<<endl;
    cout<<"Total Number of Events = "<<nevents<<endl;
    cout<<"Integrated Luminosity = "<<lum<<" fb^-1"<<endl<<endl;

    //Loop over events
    for(int k=0;k<nevents;k++){

        //Reset Variables
        Eout1=0;pzout1=0;
        Eout6=0;pzout6=0;
        Eout7=0;pzout7=0;
        Eout8=0;pzout8=0;

        if(k%10000==0) cout<<"Events Analysed = "<<k<<"!"<<endl;
        t->GetEntry(k);

        Q2_event = event->GetQ2();
        x_event = event->GetX();
        y_event = event->GetY();
        W2_event = event->GetW2();

        //Doesn't seem to be methods for these in DJANGOH
        Q2_true = event->dtrueQ2;
        x_true = event->dtrueX;
        y_true = event->dtrueY;
        W2_true = event->dtrueW2;

        //Using scattered electron quantities
        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_event>x_low[i] && x_event<x_hi[i] && Q2_event>Q2_low[j] && Q2_event<Q2_hi[j] &&
                    y_event>ymin && W2_event>W2min && y_event<ymax)
                    yield[i][j]+=1.;
            }
        }

        //Using true quantities
        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_true>x_low[i] && x_true<x_hi[i] && Q2_true>Q2_low[j] && Q2_true<Q2_hi[j] &&
                    y_true>ymin && W2_true>W2min && y_true<ymax)
                    yield_true[i][j]+=1.;
            }
        }

        //Fill Histograms
        h1->Fill( sqrt(W2_true), sqrt(W2_event) );
        if(Q2_true<50) h1a->Fill( sqrt(W2_true), sqrt(W2_event) );
        if(Q2_true>50) h1b->Fill( sqrt(W2_true), sqrt(W2_event) );
        if(W2_event>10 && y_event>0.01) h1c->Fill(Q2_event,sqrt(W2_true));
        h1d->Fill(Q2_true,sqrt(W2_true));

        h2->Fill( sqrt(W2_true) );
        h2a->Fill( sqrt(W2_true) );
        h2b->Fill( sqrt(W2_event) );
        h3->Fill(Q2_true,Q2_event);
        h4->Fill(y_true,y_event);
        h5->Fill(x_true,x_event);
        if(Q2_true<50)h5a->Fill(x_true,x_event);
        if(Q2_true>50)h5b->Fill(x_true,x_event);
        h6->Fill(event->IChannel);


        // Loop over all particles in event 
        for(int iParticle=0;iParticle<event->GetNTracks();iParticle++){

            particle = event->GetTrack(iParticle);
        
            //Get total outgoing energy and z momentum component
            if( particle->GetStatus()==1 && fabs(particle->GetEta())<4 ){

                if(event->IChannel==1){
                    Eout1+=particle->GetE();
                    pzout1+=particle->GetPz();
                }
                else if(event->IChannel==6){
                    Eout6+=particle->GetE();
                    pzout6+=particle->GetPz();
                }
                else if(event->IChannel==7){
                    Eout7+=particle->GetE();
                    pzout7+=particle->GetPz();
                }
                else if(event->IChannel==8){
                    Eout8+=particle->GetE();
                    pzout8+=particle->GetPz();
                }
            }
        }

        //Fill additional histograms
        if( event->IChannel==1 ) h7a->Fill(Eout1-pzout1);
        if( event->IChannel==6 ) h7b->Fill(Eout6-pzout6);
        if( event->IChannel==7 ) h7c->Fill(Eout7-pzout7);
        if( event->IChannel==8 ) h7d->Fill(Eout8-pzout8);

    } //End event loop

    //Convert yields to reduced cross sections
    //First do for scattered electron kinematics
    for(int i=0;i<nbins;i++){
        for(int j=0;j<nbins_Q2;j++){
            rcs[i][j] = rcs_factor[i][j] * (yield[i][j]) / (lum*x_width[i]*Q2_width[j]);
            if(yield[i][j]>0) error[i][j] = (rcs[i][j]) / (sqrt(yield[i][j]));
        }
    }

    //Create and fill TGraphs
    TGraphErrors *gr1[nbins_Q2];
    for(int j=0;j<nbins_Q2;j++){
        gr1[j] = new TGraphErrors();
        gr1[j]->SetLineWidth(2);
        gr1[j]->SetLineColor(kBlue);
        gr1[j]->SetMarkerStyle(20);
        gr1[j]->SetMarkerColor(kBlue);
    }
    
    int counter(0);
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){
            if(yield[i][j]>50 && cut_lr[i][j] && cut_ul[i][j]){
                gr1[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1[j]->SetPointError(counter,0,error[i][j]);
                counter++;
            }
        }
        counter=0;
    }

    //Next do for true kinematics
    for(int i=0;i<nbins;i++){
        for(int j=0;j<nbins_Q2;j++){
            rcs[i][j] = rcs_factor[i][j] * (yield_true[i][j]) / (lum*x_width[i]*Q2_width[j]);
            if(yield_true[i][j]>0) error[i][j] = (rcs[i][j]) / (sqrt(yield_true[i][j]));
        }
    }

    //Create and fill TGraphs
    TGraphErrors *gr1t[nbins_Q2];
    for(int j=0;j<nbins_Q2;j++){
        gr1t[j] = new TGraphErrors();
        gr1t[j]->SetLineWidth(2);
        gr1t[j]->SetLineColor(kGreen);
        gr1t[j]->SetMarkerStyle(21);
        gr1t[j]->SetMarkerColor(kGreen);
    }

    counter = 0;
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){
            if(yield_true[i][j]>50 && cut_lr[i][j] && cut_ul[i][j]){
                gr1t[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1t[j]->SetPointError(counter,0,error[i][j]);
                counter++;
            }
        }
        counter=0;
    }

    //------------------------------------------------
    //Analyse DJANGOH with no QED radiation simulation
    //------------------------------------------------

    //Load ROOT Files
    erhic::EventDjangoh *event_d(NULL); //Event Class
    erhic::ParticleMC *particle_d(NULL); //Particle Class

    TChain *t_d = new TChain("EICTree");
    
    //Djangoh 4.6.20
    t_d->Add("/gpfs02/eic/baraks/djangoh/djangoh_local/outfiles/djangoh.NC.10x100_evt.root");

    //Djangoh 4.6.10
    //for(int i=0;i<10;i++){
        //t_d->Add(Form("/eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/10_100/djangoh.NC.10x100_%d.root",i));
    //}

    t_d->SetBranchAddress("event",&event_d);

    //Calculate Generated Luminosity
    int nevents_d = t_d->GetEntries();
    
    //Djangoh 4.6.20
    double cross_tot_d = 1.18E9; //Total Cross Section in fb

    //Djangoh 4.6.10
    //double cross_tot_d = 1.06E9; //Total Cross Section in fb
    double lum_d = ( (double) nevents_d)/cross_tot_d; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"DJANGOH (QED Rad OFF) Simulation:"<<endl;
    cout<<"Total Number of Events = "<<nevents_d<<endl;
    cout<<"Integrated Luminosity = "<<lum_d<<" fb^-1"<<endl<<endl;

    //Reset Variables
    for(int i=0;i<nbins;i++){ //x loop
        for(int j=0;j<nbins_Q2;j++){ //Q2 loop
            yield[i][j]=0;  
            error[i][j]=0;
            rcs[i][j]=0;
        }
    }

    //Loop over events
    for(int k=0;k<nevents_d;k++){
        if(k%10000==0) cout<<"Events Analysed = "<<k<<"!"<<endl;
        t_d->GetEntry(k);

        Q2_event = event_d->GetQ2();
        x_event = event_d->GetX();
        y_event = event_d->GetY();
        W2_event = event_d->GetW2();

        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_event>x_low[i] && x_event<x_hi[i] && Q2_event>Q2_low[j] && Q2_event<Q2_hi[j] &&
                    y_event>ymin && W2_event>W2min && y_event<ymax)
                    yield[i][j]+=1.;
            }
        }
    } //End event loop

    //Convert Yields to reduced cross sections
    for(int i=0;i<nbins;i++){
        for(int j=0;j<nbins_Q2;j++){
            rcs[i][j] = rcs_factor[i][j] * (yield[i][j]) / (lum_d*x_width[i]*Q2_width[j]);
            if(yield[i][j]>0) error[i][j] = (rcs[i][j]) / (sqrt(yield[i][j]));
        }
    }

    //Create and fill TGraphs
    TGraphErrors *gr1_d[nbins_Q2];
    for(int j=0;j<nbins_Q2;j++){
        gr1_d[j] = new TGraphErrors();
        gr1_d[j]->SetLineWidth(2);
        gr1_d[j]->SetLineColor(kMagenta);
        gr1_d[j]->SetMarkerStyle(23);
        gr1_d[j]->SetMarkerColor(kMagenta);
    }
    
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){
            if(yield[i][j]>50 && cut_lr[i][j] && cut_ul[i][j]){
                gr1_d[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1_d[j]->SetPointError(counter,0,error[i][j]);
                counter++;
            }
        }
        counter=0;
    }

    //Latex Labels
    TLatex *tex1[num_plot];
    for(int i=0;i<num_plot;i++){
        tex1[i] = new TLatex(2.5E-2,1.8,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]]));
        tex1[i]->SetTextFont(42);
        tex1[i]->SetTextSize(0.06);
    }

    TLatex *tex2 = new TLatex(4E-4,1.6,"10 GeV e^{-} on 100 GeV p, #sqrt{s}=63.2 GeV");
    tex2->SetTextColor(kBlack);tex2->SetTextFont(42);

    TLatex *tex3 = new TLatex(4E-4,1.45,"Djangoh, y>0.001 and W^{2}>10GeV^{2}");
    tex3->SetTextColor(kBlack);tex3->SetTextFont(42);

    TLatex *tex4 = new TLatex(4E-4,1.45,"#sigma_{Rad.} using electron quantities");
    tex4->SetTextColor(kBlue);tex4->SetTextFont(62);

    TLatex *tex5 = new TLatex(4E-4,1.30,"#sigma_{Rad.} using true quantities");
    tex5->SetTextColor(kGreen);tex5->SetTextFont(62);

    TLatex *tex6 = new TLatex(4E-4,1.15,"#sigma_{Born.}");
    tex6->SetTextColor(kMagenta);tex6->SetTextFont(62);

    //Draw plots
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    //gStyle->SetNdivisions(405,"X");
    //gStyle->SetNdivisions(405,"Y");
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetLabelOffset(0.01,"X");
    gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.07);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYSize(0.07);
    gStyle->SetTitleYOffset(0.8);

    float small = 1E-20;

    TCanvas *c1 = new TCanvas("c1");
    c1->Divide(2,2,small,small);

    c1->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky(); 
    TH1 *hframe = gPad->DrawFrame(1E-4,0,1,2);
    hframe->GetXaxis()->SetTitle("x");hframe->GetXaxis()->CenterTitle();
    hframe->GetYaxis()->SetTitle("#sigma_{r,NC}^{-}");hframe->GetYaxis()->CenterTitle();
    counter = 0;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");  
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();

    c1->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex4->Draw();tex5->Draw();tex6->Draw();

    c1->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    
    c1->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    TCanvas *c2 = new TCanvas("c2");
    c2->Divide(2,2,small,small);

    c2->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();

    c2->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex4->Draw();tex5->Draw();tex6->Draw();
    
    c2->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    
    c2->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    hframe->Draw();
    counter++;
    gr1[Q2_bin_plot[counter]]->Draw("P Same");
    gr1t[Q2_bin_plot[counter]]->Draw("P Same");
    gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    //Plot Histograms
    TCanvas *cc1 = new TCanvas("cc1");
    cc1->SetLogz();
    h1->Draw("colz");

    TCanvas *cc1a = new TCanvas("cc1a");
    cc1a->SetLogz();
    h1a->Draw("colz");

    TCanvas *cc1b = new TCanvas("cc1b");
    cc1b->SetLogz();
    h1b->Draw("colz");

    TCanvas *cc1c = new TCanvas("cc1c");
    cc1c->SetLogz();
    h1c->Draw("colz");

    TCanvas *cc1d = new TCanvas("cc1d");
    cc1d->SetLogz();
    h1d->Draw("colz");

    TCanvas *cc2 = new TCanvas("cc2");
    h2->Draw();

    TCanvas *cc2a = new TCanvas("cc2a");
    h2a->Draw();

    TCanvas *cc2b = new TCanvas("cc2b");
    h2b->Draw();

    TCanvas *cc3 = new TCanvas("cc3");
    cc3->SetLogz();
    h3->Draw("colz");

    TCanvas *cc4 = new TCanvas("cc4");
    cc4->SetLogx();cc4->SetLogy();cc4->SetLogz();
    h4->Draw("colz");

    TCanvas *cc5 = new TCanvas("cc5");
    cc5->SetLogx();cc5->SetLogy();cc5->SetLogz();
    h5->Draw("colz");

    TCanvas *cc5a = new TCanvas("cc5a");
    cc5a->SetLogx();cc5a->SetLogy();cc5a->SetLogz();
    h5a->Draw("colz");

    TCanvas *cc5b = new TCanvas("cc5b");
    cc5b->SetLogx();cc5b->SetLogy();cc5b->SetLogz();
    h5b->Draw("colz");

    TCanvas *cc6 = new TCanvas("cc6");
    h6->Draw();

    TCanvas *cc7 = new TCanvas("cc7");
    gPad->SetLogy();
    h7a->Draw();
    h7b->Draw("same");
    h7c->Draw("same");
    h7d->Draw("same");

    TLegend *leg1 = new TLegend(0.125,0.6,0.325,0.8);
    leg1->SetBorderSize(0);leg1->SetTextSize(0.035);
    leg1->AddEntry(h7a,"Channel 1","l");
    leg1->AddEntry(h7b,"Channel 6","l");
    leg1->AddEntry(h7c,"Channel 7","l");
    leg1->AddEntry(h7d,"Channel 8","l");
    leg1->Draw();

    //Print to File
    c1->Print("plots/rad_djangoh.pdf[");
    c1->Print("plots/rad_djangoh.pdf");
    c2->Print("plots/rad_djangoh.pdf");
    cc1->Print("plots/rad_djangoh.pdf");
    cc1a->Print("plots/rad_djangoh.pdf");
    cc1b->Print("plots/rad_djangoh.pdf");
    cc1c->Print("plots/rad_djangoh.pdf");
    cc1d->Print("plots/rad_djangoh.pdf");
    cc2->Print("plots/rad_djangoh.pdf");
    cc2a->Print("plots/rad_djangoh.pdf");
    cc2b->Print("plots/rad_djangoh.pdf");
    cc3->Print("plots/rad_djangoh.pdf");
    cc4->Print("plots/rad_djangoh.pdf");
    cc5->Print("plots/rad_djangoh.pdf");
    cc5a->Print("plots/rad_djangoh.pdf");
    cc5b->Print("plots/rad_djangoh.pdf");
    cc6->Print("plots/rad_djangoh.pdf");
    cc7->Print("plots/rad_djangoh.pdf");
    cc7->Print("plots/rad_djangoh.pdf]");

}