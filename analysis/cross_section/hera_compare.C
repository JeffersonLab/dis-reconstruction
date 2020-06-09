R__LOAD_LIBRARY(libeicsmear);

void hera_compare(){

    //Constants
    const double alpha(1./137.036);
    const double s_cm = 4.*27.5*920.;
    const double fbgev(1./(0.3894E12));

    //Cross Section Bins
    const int nbins_Q2=8;
    double Q2_low[] = {1.8, 3.25, 6.0,  9.0, 18, 55, 110, 175};
    double Q2_hi[] =  {2.2, 3.75, 7.0, 11.0, 26, 65, 130, 225};
            
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

    //Remove edge x points
    double xcut_low[] = {2.5E-5, 4E-5,   6E-5, 1E-4, 2.5E-4, 6E-4, 1.5E-3, 2.5E-3};
    double xcut_hi[] =  {1.5E-1, 2.5E-1, 7E-1, 7E-1, 7E-1,   1.,   1.,     1.};

    //Bin Yields/Kinematics/Factors
    //int nbins_tot = nbins_Q2*nbins;
    double yield[nbins][nbins_Q2];
    double error[nbins][nbins_Q2];
    double rcs[nbins][nbins_Q2];
    double y[nbins][nbins_Q2];
    double Q2_center[nbins_Q2];
    double x_center[nbins];
    double Q2_width[nbins_Q2];
    double x_width[nbins];
    double rcs_factor[nbins][nbins_Q2];

    for(int i=0;i<nbins;i++){ //x loop
        for(int j=0;j<nbins_Q2;j++){ //Q2 loop
            yield[i][j]=0;  
            error[i][j]=0;
            rcs[i][j]=0;  

            Q2_center[j] = (Q2_hi[j]+Q2_low[j])/2.;
            x_center[i] = (x_hi[i]+x_low[i])/2.;
            Q2_width[j] = (Q2_hi[j]-Q2_low[j]);
            x_width[i] = (x_hi[i]-x_low[i]);

            y[i][j] = Q2_center[j]/(x_center[i]*s_cm);
            rcs_factor[i][j] = fbgev * ( pow(Q2_center[j],2)*x_center[i] ) / 
                            ( TMath::TwoPi() * pow(alpha,2) * ( 1.+pow((1.-y[i][j]),2) ) );

        }
    }

    //---------------------------------
    //Analyse PYTHIA6 simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *t = new TChain("EICTree");
    for(int i=0;i<=9;i++){
        t->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/hera/ep_hera_%d.root",i));
    }

    t->SetBranchAddress("event",&event);
    double Q2_event; //Get this from scattered electron (may differ from true if rad. effect included)
    double x_event; //''

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    double cross_tot = 0.8876E9; //Total Cross Section in fb
    double lum = ( (double) nevents)/cross_tot; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"PYTHIA6 Simulation:"<<endl;
    cout<<"Total Number of Events = "<<nevents<<endl;
    cout<<"Integrated Luminosity = "<<lum<<" fb^-1"<<endl<<endl;

    //Loop over events
    for(int k=0;k<nevents;k++){
        if(k%10000==0) cout<<"Events Analysed = "<<k<<"!"<<endl;
        t->GetEntry(k);

        Q2_event = event->GetQ2();
        x_event = event->GetX();

        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_event>x_low[i] && x_event<x_hi[i] && Q2_event>Q2_low[j] && Q2_event<Q2_hi[j])
                    yield[i][j]+=1.;
            }
        }
    }

    //Convert Yields to reduced cross sections
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
        gr1[j]->SetLineColor(kRed);
        gr1[j]->SetMarkerStyle(20);
        gr1[j]->SetMarkerColor(kRed);
    }
    
    int counter(0);
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){
            if(yield[i][j]>50 && x_center[i]>xcut_low[j] && x_center[i]<xcut_hi[j]){
                gr1[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1[j]->SetPointError(counter,0,error[i][j]);
                counter++;
            }
        }
        counter=0;
    }


    //---------------------------------
    //Analyse DJANGOH simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventDjangoh *event_d(NULL); //Event Class
    erhic::ParticleMC *particle_d(NULL); //Particle Class

    TChain *t_d = new TChain("EICTree");
    for(int i=0;i<=9;i++){
        t_d->Add(Form("/eic/data/baraks/dis-reconstruction/djangoh/outfiles/hera/djangoh_hera_%d.root",i));
    }

    t_d->SetBranchAddress("event",&event_d);

    //Calculate Generated Luminosity
    int nevents_d = t_d->GetEntries();
    double cross_tot_d = 0.724E9; //Total Cross Section in fb
    double lum_d = ( (double) nevents_d)/cross_tot_d; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"DJANGOH Simulation:"<<endl;
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

        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_event>x_low[i] && x_event<x_hi[i] && Q2_event>Q2_low[j] && Q2_event<Q2_hi[j])
                    yield[i][j]+=1.;
            }
        }
    }

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
            if(yield[i][j]>50 && x_center[i]>xcut_low[j] && x_center[i]<xcut_hi[j]){
                gr1_d[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1_d[j]->SetPointError(counter,0,error[i][j]);
                counter++;
            }
        }
        counter=0;
    }

    //---------------------------------
    //Read-in and fill HERA data
    //---------------------------------
    ifstream in;
    in.open("hera_data.txt",ios::in);
  
    double read_q2,read_x,read_sr,read_pe,read_s;
    int line1(0);

    int hera_counter[nbins_Q2];
    TGraphErrors *gr2[nbins_Q2];
    for(int j=0;j<nbins_Q2;j++){

        hera_counter[j]=0;

        gr2[j] = new TGraphErrors();
        gr2[j]->SetLineWidth(2);
        gr2[j]->SetLineColor(kBlue);
        gr2[j]->SetMarkerStyle(22);
        gr2[j]->SetMarkerColor(kBlue);
    }

    while(in.good()){
        in >> read_q2 >> read_x >> read_sr >> read_pe >> read_s;

        if(in.good()){
            for(int j=0;j<nbins_Q2;j++){
                if( fabs(read_q2-Q2_center[j])<0.1 && fabs(read_s-319)<0.1 ){
                    gr2[j]->SetPoint(hera_counter[j],read_x,read_sr);
                    gr2[j]->SetPointError(hera_counter[j],0, (read_pe/100.)*read_sr );
                    hera_counter[j]++;
                }
            }
            line1++;
        }
    }

    cout << "We read " << line1 << " HERA data lines."<<endl;
    in.close();

    //Latex Labels
    TLatex *tex1[nbins_Q2];
    for(int j=0;j<nbins_Q2;j++){
        tex1[j] = new TLatex(1E-2,1.8,Form("Q^{2} = %.0f GeV^{2}",Q2_center[j]));
        tex1[j]->SetTextFont(42);
        tex1[j]->SetTextSize(0.06);
    }

    TLatex *tex2 = new TLatex(2E-3,1.6,"HERA e^{+}p, #sqrt{s}=319 GeV");
    tex2->SetTextColor(kBlue);tex2->SetTextFont(42);

    TLatex *tex3 = new TLatex(2E-3,1.45,"Pythia6 e^{+}p, #sqrt{s}=319 GeV");
    tex3->SetTextColor(kRed);tex3->SetTextFont(42);

    TLatex *tex4 = new TLatex(2E-3,1.3,"Djangoh e^{+}p, #sqrt{s}=319 GeV");
    tex4->SetTextColor(kMagenta);tex4->SetTextFont(42);

    TLatex *tex2a = new TLatex(3E-5,0.6,"HERA e^{+}p, #sqrt{s}=319 GeV");
    tex2a->SetTextColor(kBlue);tex2a->SetTextFont(42);

    TLatex *tex3a = new TLatex(3E-5,0.45,"Pythia6 e^{+}p, #sqrt{s}=319 GeV");
    tex3a->SetTextColor(kRed);tex3a->SetTextFont(42);

    TLatex *tex4a = new TLatex(3E-5,0.3,"Djangoh e^{+}p, #sqrt{s}=319 GeV");
    tex4a->SetTextColor(kMagenta);tex4a->SetTextFont(42);

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
    TH1 *hframe = gPad->DrawFrame(1E-5,0,1,2);
    hframe->GetXaxis()->SetTitle("x");hframe->GetXaxis()->CenterTitle();
    hframe->GetYaxis()->SetTitle("#sigma_{r,NC}^{+}");hframe->GetYaxis()->CenterTitle();
    gr2[0]->Draw("PL Same");gr1[0]->Draw("P Same");gr1_d[0]->Draw("P Same");
    tex1[0]->Draw();

    c1->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    gr2[1]->Draw("PL Same");gr1[1]->Draw("P Same");gr1_d[1]->Draw("P Same");
    tex1[1]->Draw();
    
    c1->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    gr2[2]->Draw("PL Same");gr1[2]->Draw("P Same");gr1_d[2]->Draw("P Same");
    tex1[2]->Draw();
    
    c1->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    gr2[3]->Draw("PL Same");gr1[3]->Draw("P Same");gr1_d[3]->Draw("P Same");
    tex1[3]->Draw();tex2->Draw();tex3->Draw();tex4->Draw();

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
    gr2[4]->Draw("PL Same");gr1[4]->Draw("P Same");gr1_d[4]->Draw("P Same");
    tex1[4]->Draw();

    c2->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    gr2[5]->Draw("PL Same");gr1[5]->Draw("P Same");gr1_d[5]->Draw("P Same");
    tex1[5]->Draw();
    
    c2->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    gr2[6]->Draw("PL Same");gr1[6]->Draw("P Same");gr1_d[6]->Draw("P Same");
    tex1[6]->Draw();
    
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
    gr2[7]->Draw("PL Same");gr1[7]->Draw("P Same");gr1_d[7]->Draw("P Same");
    tex1[7]->Draw();tex2a->Draw();tex3a->Draw();tex4a->Draw();

    //Print to File
    c1->Print("plots/hera_compare.pdf[");
    c1->Print("plots/hera_compare.pdf");
    c2->Print("plots/hera_compare.pdf");
    c2->Print("plots/hera_compare.pdf]");

}