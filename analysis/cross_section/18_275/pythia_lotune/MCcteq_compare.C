//Only want to show cteq LO (cteq6l1, 10042) pdf that was used in the MC
//We show 3 sets of plots:
//1) Average MC diff. cs in bin vs. cteq diff. cs
//2) Bin-centered MC vs. cteq diff. cs
//3) Average MC diff. cs in bin vs. average cteq diff. cs in bin

R__LOAD_LIBRARY(libeicsmear);

//Globals
std::vector<double> vec_Q2bc_a, vec_xbc_a, vec_sigcenter_a, vec_sigaverage_a, vec_bc_a; //cteq6l1

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
TGraph* Theory_Graph(const char* infile, double Q2_set=1, 
                        int color=1, int style=20, int lw=2, bool abs=0){

  //Create graph
  TGraph *gr = new TGraph();
  gr->SetMarkerStyle(style);gr->SetMarkerColor(color);
  gr->SetLineWidth(lw);gr->SetLineColor(color);

  //Read input file
  ifstream in;
  in.open(infile,ios::in);

  //Write to Graph
  double Q2(0),x(0),y(0);
  int counter(0);

  while(in.good()){
    in >> Q2 >> x >> y;

    if( in.good() && fabs(Q2-Q2_set)<0.5 ){

      if(abs)
        gr->SetPoint(counter,x,fabs(y));
      else
        gr->SetPoint(counter,x,y);

      counter++;
    }
  }

  cout<<"For Q2 = "<<Q2_set<<" GeV^2... ";
  cout<<"We read "<<counter<<" matching lines from data file "<<infile<<" !"<<endl;
  in.close();
  return gr;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void read_bc(const char* infile,int model){

    //Read input file
    ifstream in;
    in.open(infile,ios::in);

    //Write to Graph
    double Q2_temp,x_temp,cs_center_temp,cs_average_temp,bc_temp;
    int counter(0);

    while(in.good()){
        in >> Q2_temp >> x_temp >> cs_center_temp >> cs_average_temp >>bc_temp;

        if( in.good() ){

            if(model==1){
                vec_Q2bc_a.push_back(Q2_temp);
                vec_xbc_a.push_back(x_temp);
                vec_sigcenter_a.push_back(cs_center_temp);
                vec_sigaverage_a.push_back(cs_average_temp);
                vec_bc_a.push_back(bc_temp);
            }
            counter++;
        }
    }

  cout<<"We read "<<counter<<" lines from bin-centering data file "<<infile<<" !"<<endl;
  in.close();

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
double get_bc_fac(double x_center,double Q2_center,int model){
    
    double bc_fac = 0;
    int size = 0;

    if(model==1){
        size = (int) vec_Q2bc_a.size();

        for(int i=0;i<size;i++){
            if( fabs( (vec_xbc_a[i]-x_center)/x_center ) < 0.02 && fabs( (vec_Q2bc_a[i]-Q2_center)/Q2_center ) < 0.02 ){
                bc_fac = vec_bc_a[i];
                //cout <<x_center<<" "<<Q2_center<<" "<<bc_fac <<endl;
                break;
            }
        }
    }

    return bc_fac;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
double get_model_average(double x_center,double Q2_center,int model){
    double average_cs = 0;
    int size = 0;

    if(model==1){
        size = (int) vec_Q2bc_a.size();

        for(int i=0;i<size;i++){
            if( fabs( (vec_xbc_a[i]-x_center)/x_center ) < 0.02 && fabs( (vec_Q2bc_a[i]-Q2_center)/Q2_center ) < 0.02 ){
                average_cs = vec_sigaverage_a[i];
                
                break;
            }
        }
    }

    return average_cs;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

void MCcteq_compare(){

    //Constants
    const double alpha(1./137.036);
    const double s_cm = 4.*18*275;
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

    //Bin Yields/Kinematics/Factors
    //int nbins_tot = nbins_Q2*nbins;
    double yield[nbins][nbins_Q2];
    double error[nbins][nbins_Q2];
    double error_10[nbins][nbins_Q2];
    double rcs[nbins][nbins_Q2];
    double y[nbins][nbins_Q2];
    double Q2_center[nbins_Q2];
    double x_center[nbins];
    double Q2_width[nbins_Q2];
    double x_width[nbins];
    double rcs_factor[nbins][nbins_Q2];
    bool cut_ul[nbins][nbins_Q2];
    bool cut_lr[nbins][nbins_Q2];
    double bc_fac_a[nbins][nbins_Q2];
    double model_average[nbins][nbins_Q2];
    double ymax = 0.98;double ymin = 0.001; double W2min = 10;
    double y_temp, W2_temp;

    for(int i=0;i<nbins;i++){ //x loop
        for(int j=0;j<nbins_Q2;j++){ //Q2 loop
            yield[i][j]=0;  
            error[i][j]=0;
            error_10[i][j]=0;
            rcs[i][j]=0;  
            bc_fac_a[i][j]=0;
            model_average[i][j]=0;

            Q2_center[j] = (Q2_hi[j]+Q2_low[j])/2.;
            x_center[i] = (x_hi[i]+x_low[i])/2.;
            Q2_width[j] = (Q2_hi[j]-Q2_low[j]);
            x_width[i] = (x_hi[i]-x_low[i]);

            y[i][j] = Q2_center[j]/(x_center[i]*s_cm);
            rcs_factor[i][j] = fbgev * ( pow(Q2_center[j],2)*x_center[i] ) / 
                            ( TMath::TwoPi() * pow(alpha,2) * ( 1.+pow((1.-y[i][j]),2) ) );

            //Print out bins and centers
            //cout<<"i, j = "<< i <<" , "<< j << endl;
            //cout<<"x, Q2 = "<< x_center[i] <<" , "<< Q2_center[j] << endl;

            //Cut Bins with top-left edge > ymax and bottom-right edge < ymin and > W2min
            y_temp = (Q2_hi[j])/(x_low[i]*s_cm); //Upper left
            if(y_temp>ymax) //Remove bin
                cut_ul[i][j] = 0;
            else
                cut_ul[i][j] = 1;

            //cout<<"y_temp = "<<y_temp<<endl;

            y_temp = (Q2_low[j])/(x_hi[i]*s_cm); //Lower Right
            W2_temp = Mp*Mp + Q2_low[j]*( (1./x_hi[i]) - 1.);
            if(y_temp<ymin || W2_temp<W2min) //Remove bin
                cut_lr[i][j] = 0;
            else
                cut_lr[i][j] = 1;

            //cout<<"y_temp = "<<y_temp<<endl;
            //cout<<"W2_temp = "<<W2_temp<<endl;

        }
    }

    //Q2 bins (center) to plot; 8 bins
    const int num_plot = 8;
    int Q2_bin_plot[num_plot] = {6,8,9,10,11,13,15,18};

    //Read in Theory for RCS 
    TGraph *cteq6l1[num_plot];

    for(int i=0;i<num_plot;i++){
        cteq6l1[i] = Theory_Graph("input_files/cteq6l1_SF_rcs.txt",Q2_center[Q2_bin_plot[i]],kBlue,20,2);
    }

    //Read in Bin-Centering Correction Files
    read_bc("input_files/cteq6l1_SF_bin_centering.txt",1);

    //---------------------------------
    //Analyse PYTHIA6 simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *t = new TChain("EICTree");
    for(int i=0;i<15;i++){
        t->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/yellow_report/18_275/ep_18_275_newtune_%d.root",i));
    }

    t->SetBranchAddress("event",&event);
    double Q2_event; //Get this from scattered electron (may differ from true if rad. effect included)
    double x_event; //''
    double y_event; //''
    double W2_event; //''

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    double cross_tot = 1.67E9; //Total Cross Section in fb
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
        y_event = event->GetY();
        W2_event = event->GetW2();

        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_event>x_low[i] && x_event<x_hi[i] && Q2_event>Q2_low[j] && Q2_event<Q2_hi[j] &&
                    y_event>ymin && W2_event>W2min && y_event<ymax) //pythia has a hard cut W>2GeV applied
                    yield[i][j]+=1.;
            }
        }
    }

    //Convert Yields to reduced cross sections
    for(int i=0;i<nbins;i++){
        for(int j=0;j<nbins_Q2;j++){
            rcs[i][j] = rcs_factor[i][j] * (yield[i][j]) / (lum*x_width[i]*Q2_width[j]);

            bc_fac_a[i][j] = get_bc_fac(x_center[i],Q2_center[j],1);
            model_average[i][j] = rcs_factor[i][j] * get_model_average(x_center[i],Q2_center[j],1);
            
            if(yield[i][j]>0) { 
                error[i][j] = (rcs[i][j]) / (sqrt(yield[i][j]));
                error_10[i][j] = 100.*( 1./ sqrt(yield[i][j]*(10./lum)) ); //% Error for 10 fb^-1
            }
        }
    }

    //Create and fill TGraphs
    TGraphErrors *gr1[nbins_Q2]; //Without Bin-Centering Correction
    TGraphErrors *gr2[nbins_Q2]; //With Bin-Centering Correction
    TGraph *gbc1[nbins_Q2]; //Bin-Centering Correction Factor
    TGraph *ge1[nbins_Q2]; //Expected Statistical Errors (%) for L=10 fb^-1
    TGraph *cteq6l1_a[nbins_Q2]; //Average model cross section (divided by RCS factor at bin center) 

    for(int j=0;j<nbins_Q2;j++){
        gr1[j] = new TGraphErrors();
        gr1[j]->SetLineWidth(2);
        gr1[j]->SetLineColor(kRed);
        gr1[j]->SetMarkerStyle(20);
        gr1[j]->SetMarkerColor(kRed);

        gr2[j] = new TGraphErrors();
        gr2[j]->SetLineWidth(2);
        gr2[j]->SetLineColor(kRed);
        gr2[j]->SetMarkerStyle(20);
        gr2[j]->SetMarkerColor(kRed);

        gbc1[j] = new TGraph(); //Correction comes from cteq6l1
        gbc1[j]->SetLineWidth(2);
        gbc1[j]->SetLineColor(kGreen);
        gbc1[j]->SetMarkerStyle(21);
        gbc1[j]->SetMarkerColor(kGreen);

        ge1[j] = new TGraph();
        ge1[j]->SetLineWidth(2);
        ge1[j]->SetLineColor(kRed);
        ge1[j]->SetMarkerStyle(20);
        ge1[j]->SetMarkerColor(kRed);

        cteq6l1_a[j] = new TGraph();
        cteq6l1_a[j]->SetLineWidth(2);
        cteq6l1_a[j]->SetLineColor(kBlue);
        cteq6l1_a[j]->SetMarkerStyle(21);
        cteq6l1_a[j]->SetMarkerColor(kBlue);
    }
    
    int counter(0);
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){
            if(yield[i][j]>50 && cut_lr[i][j] && cut_ul[i][j]){
                gr1[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1[j]->SetPointError(counter,0,error[i][j]);

                gbc1[j]->SetPoint(counter,x_center[i],bc_fac_a[i][j]);

                ge1[j]->SetPoint(counter,x_center[i],error_10[i][j]);

                gr2[j]->SetPoint(counter,x_center[i],rcs[i][j]*bc_fac_a[i][j]);
                gr2[j]->SetPointError(counter,0,error[i][j]*bc_fac_a[i][j]);

                cteq6l1_a[j]->SetPoint(counter,x_center[i],model_average[i][j]);
                
                counter++;
            }
        }
        counter=0;
    }

    //---------------------------------
    //Analyse Pythia6 tuned simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventPythia *event_d(NULL); //Event Class
    erhic::ParticleMC *particle_d(NULL); //Particle Class

    TChain *t_d = new TChain("EICTree");
    for(int i=0;i<15;i++){
        t_d->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/other_studies/18_275/ep_minbias_%d.root",i));
    }

    t_d->SetBranchAddress("event",&event_d);

    //Calculate Generated Luminosity
    int nevents_d = t_d->GetEntries();
    double cross_tot_d = 1.67E9; //Total Cross Section in fb
    double lum_d = ( (double) nevents_d)/cross_tot_d; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"PYTHIA6 Tuned Simulation:"<<endl;
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
                    y_event>ymin && W2_event>W2min && y_event<ymax && 
                    event_d->GetProcess()==99) //pythia has a hard cut W>2GeV applied
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
    TGraphErrors *gr1_d[nbins_Q2]; //Without Bin-Centering Correction
    TGraphErrors *gr2_d[nbins_Q2]; //With Bin-Centering Correction

    for(int j=0;j<nbins_Q2;j++){
        gr1_d[j] = new TGraphErrors();
        gr1_d[j]->SetLineWidth(2);
        gr1_d[j]->SetLineColor(kMagenta);
        gr1_d[j]->SetMarkerStyle(23);
        gr1_d[j]->SetMarkerColor(kMagenta);

        gr2_d[j] = new TGraphErrors();
        gr2_d[j]->SetLineWidth(2);
        gr2_d[j]->SetLineColor(kMagenta);
        gr2_d[j]->SetMarkerStyle(23);
        gr2_d[j]->SetMarkerColor(kMagenta);
    }
    
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){
            if(yield[i][j]>50 && cut_lr[i][j] && cut_ul[i][j]){
                gr1_d[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1_d[j]->SetPointError(counter,0,error[i][j]);
                
                gr2_d[j]->SetPoint(counter,x_center[i],rcs[i][j]*bc_fac_a[i][j]);
                gr2_d[j]->SetPointError(counter,0,error[i][j]*bc_fac_a[i][j]);

                counter++;
            }
        }
        counter=0;
    }

    //Latex Labels
    TLatex *tex1[num_plot], *tex1a[num_plot];
    TLatex *tex1b[num_plot], *tex1c[num_plot]; //Error plots
    TLatex *tex1d[num_plot], *tex1e[num_plot]; //Ratio plots
    TLatex *tex1f[num_plot]; //Model BC Ratiow
    for(int i=0;i<num_plot;i++){
        tex1[i] = new TLatex(2.5E-2,1.8,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]]));
        tex1[i]->SetTextFont(42);
        tex1[i]->SetTextSize(0.06);

        tex1a[i] = new TLatex(2.5E-2,1.05,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]]));
        tex1a[i]->SetTextFont(42);
        tex1a[i]->SetTextSize(0.06);

        tex1b[i] = new TLatex(1E-2,0.09,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]])); //Only 0-3 drawn
        tex1b[i]->SetTextFont(42);
        tex1b[i]->SetTextSize(0.06);

        tex1c[i] = new TLatex(1E-2,0.3,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]])); //Only 4-7 drawn
        tex1c[i]->SetTextFont(42);
        tex1c[i]->SetTextSize(0.06);

        tex1d[i] = new TLatex(3E-2,1.15,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]])); //Only 0-3 drawn
        tex1d[i]->SetTextFont(42);
        tex1d[i]->SetTextSize(0.06);

        tex1e[i] = new TLatex(2E-4,1.15,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]])); //Only 4-7 drawn
        tex1e[i]->SetTextFont(42);
        tex1e[i]->SetTextSize(0.06);

        tex1f[i] = new TLatex(2.5E-2,1.008,Form("Q^{2} = %.1f GeV^{2}",Q2_center[Q2_bin_plot[i]]));
        tex1f[i]->SetTextFont(42);
        tex1f[i]->SetTextSize(0.06);
    }

    TLatex *tex2 = new TLatex(3E-4,1.6,"18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
    tex2->SetTextColor(kBlack);tex2->SetTextFont(42);

    TLatex *tex2a = new TLatex(5E-4,1.02,"18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
    tex2a->SetTextColor(kBlack);tex2a->SetTextFont(42);

    TLatex *tex2b = new TLatex(5E-4,0.0825,"18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
    tex2b->SetTextColor(kBlack);tex2b->SetTextFont(42);

    TLatex *tex2c = new TLatex(5E-4,0.275,"18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
    tex2c->SetTextColor(kBlack);tex2c->SetTextFont(42);

    TLatex *tex3 = new TLatex(3E-4,1.45,Form("Pythia6 e^{-}p, %.0e<y<%.2f and W^{2}>%.0fGeV^{2}",ymin,ymax,W2min));
    tex3->SetTextColor(kRed);tex3->SetTextFont(62);

    TLatex *tex4 = new TLatex(3E-4,1.3,"Pythia6 e^{-}p, w/parameter tuning");
    tex4->SetTextColor(kMagenta);tex4->SetTextFont(62);

    TLatex *tex6 = new TLatex(1E-3,1.6,"cteq6l1 (10042)");
    tex6->SetTextColor(kBlue);tex6->SetTextFont(62);

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
    gStyle->SetTitleYOffset(0.95);

    float small = 1E-20;

    //Average MC diff. cs in bin vs. cteq diff. cs
    //------------------------------------
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
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();

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
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex6->Draw();

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
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();

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
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex6->Draw();
    
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
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    //Bin-centered MC vs. cteq diff. cs
    //------------------------------------
    TCanvas *c1a = new TCanvas("c1a");
    c1a->Divide(2,2,small,small);

    c1a->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky(); 
    hframe->Draw();
    counter = 0;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();

    c1a->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex6->Draw();

    c1a->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    
    c1a->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    TCanvas *c2a = new TCanvas("c2a");
    c2a->Divide(2,2,small,small);

    c2a->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();

    c2a->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex6->Draw();
    
    c2a->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    
    c2a->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    //Average MC diff. cs in bin vs. average cteq diff. cs in bin
    //------------------------------------
    TCanvas *c1b = new TCanvas("c1b");
    c1b->Divide(2,2,small,small);

    c1b->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky(); 
    hframe->Draw();
    counter = 0;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();

    c1b->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex6->Draw();

    c1b->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    
    c1b->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    TCanvas *c2b = new TCanvas("c2b");
    c2b->Divide(2,2,small,small);

    c2b->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();

    c2b->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex6->Draw();
    
    c2b->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    
    c2b->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    counter++;
    cteq6l1_a[Q2_bin_plot[counter]]->Draw("L Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    //Print to File
    c1->Print("plots/MCcteq_compare.pdf[");
    c1->Print("plots/MCcteq_compare.pdf");
    c2->Print("plots/MCcteq_compare.pdf");
    c1a->Print("plots/MCcteq_compare.pdf");
    c2a->Print("plots/MCcteq_compare.pdf");
    c1b->Print("plots/MCcteq_compare.pdf");
    c2b->Print("plots/MCcteq_compare.pdf");
    c2b->Print("plots/MCcteq_compare.pdf]");

}