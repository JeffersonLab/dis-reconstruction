R__LOAD_LIBRARY(libeicsmear);

//Globals
std::vector<double> vec_Q2bc_a, vec_xbc_a, vec_sigcenter_a, vec_sigaverage_a, vec_bc_a; //JAM
std::vector<double> vec_Q2bc_b, vec_xbc_b, vec_sigcenter_b, vec_sigaverage_b, vec_bc_b; //CTEQ
std::vector<double> vec_Q2bc_c, vec_xbc_c, vec_sigcenter_c, vec_sigaverage_c, vec_bc_c; //NNPDF

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
            else if(model==2){
                vec_Q2bc_b.push_back(Q2_temp);
                vec_xbc_b.push_back(x_temp);
                vec_sigcenter_b.push_back(cs_center_temp);
                vec_sigaverage_b.push_back(cs_average_temp);
                vec_bc_b.push_back(bc_temp);
            }
            else if(model==3){
                vec_Q2bc_c.push_back(Q2_temp);
                vec_xbc_c.push_back(x_temp);
                vec_sigcenter_c.push_back(cs_center_temp);
                vec_sigaverage_c.push_back(cs_average_temp);
                vec_bc_c.push_back(bc_temp);
            }

            counter++;
        }
    }

  cout<<"We read "<<counter<<" lines from bin-centering data file "<<infile<<" !"<<endl;
  in.close();

}

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
    else if(model==2){
        size = (int) vec_Q2bc_b.size();

        for(int i=0;i<size;i++){
            if( fabs( (vec_xbc_b[i]-x_center)/x_center ) < 0.02 && fabs( (vec_Q2bc_b[i]-Q2_center)/Q2_center ) < 0.02 ){
                bc_fac = vec_bc_b[i];
                //cout <<x_center<<" "<<Q2_center<<" "<<bc_fac <<endl;
                break;
            }
        }
    }
    else if(model==3){
        size = (int) vec_Q2bc_c.size();

        for(int i=0;i<size;i++){
            if( fabs( (vec_xbc_c[i]-x_center)/x_center ) < 0.02 && fabs( (vec_Q2bc_c[i]-Q2_center)/Q2_center ) < 0.02 ){
                bc_fac = vec_bc_c[i];
                //cout <<x_center<<" "<<Q2_center<<" "<<bc_fac <<endl;
                break;
            }
        }
    }

    return bc_fac;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
TGraph* Ratio_Theory(TGraph *gr_bottom, TGraph *gr_top, int color=1, int style=20, int lw=2){

  //Create graph
  TGraph *gr = new TGraph();
  gr->SetMarkerStyle(style);gr->SetMarkerColor(color);
  gr->SetLineWidth(lw);gr->SetLineColor(color);

  double xup(0),yup(0),xdown(0),ydown(0);

  for(int i=0;i<gr_top->GetN();i++){
    gr_top->GetPoint(i,xup,yup);
    gr_bottom->GetPoint(i,xdown,ydown); //xup should always equal xdown as we generate the theory grids at the same x

    gr->SetPoint(i,xup,(yup/ydown));
  }

  return gr;

}

TGraphErrors* Ratio_Data(TGraph *gr_theory, TGraphErrors *gr_data, int color=1, int style=20, int lw=2){

  //Create graph
  TGraphErrors *gr = new TGraphErrors();
  gr->SetMarkerStyle(style);gr->SetMarkerColor(color);
  gr->SetLineWidth(lw);gr->SetLineColor(color);

  double xdata(0),ydata(0),yerr(0);
  double ytheory(0);

  for(int i=0;i<gr_data->GetN();i++){
    gr_data->GetPoint(i,xdata,ydata);
    yerr = gr_data->GetErrorY(i);

    ytheory = gr_theory->Eval(xdata,0,"S"); //Theory grid points aren't same as center of data bins

    gr->SetPoint(i,xdata,(ydata/ytheory));
    gr->SetPointError(i,0,(yerr/ytheory));
  }

  return gr;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

void theory_compare(){

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
    double bc_fac_b[nbins][nbins_Q2];
    double bc_fac_c[nbins][nbins_Q2];
    double ymax = 0.98;double ymin = 0.001; double W2min = 10;
    double y_temp, W2_temp;

    for(int i=0;i<nbins;i++){ //x loop
        for(int j=0;j<nbins_Q2;j++){ //Q2 loop
            yield[i][j]=0;  
            error[i][j]=0;
            error_10[i][j]=0;
            rcs[i][j]=0;  
            bc_fac_a[i][j]=0;
            bc_fac_b[i][j]=0;
            bc_fac_c[i][j]=0;

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
    TGraph *CT18ptxg[num_plot];
    TGraph *NNPDF31_nnlo_pch_as_0118_SF[num_plot];
    TGraph *JAM4EIC[num_plot];

    TGraph *ratio_CT18[num_plot];  //Ratio to JAM
    TGraph *ratio_NNPDF[num_plot]; //''
    TGraph *ratio_JAM[num_plot];   //''

    for(int i=0;i<num_plot;i++){
        CT18ptxg[i] = Theory_Graph("input_files/CT18ptxg_rcs.txt",Q2_center[Q2_bin_plot[i]],kBlack,20,2);
        NNPDF31_nnlo_pch_as_0118_SF[i] = Theory_Graph("input_files/NNPDF31_nnlo_pch_as_0118_SF_rcs.txt",Q2_center[Q2_bin_plot[i]],kBlue,20,2);
        JAM4EIC[i] = Theory_Graph("input_files/JAM4EIC_p_rcs.txt",Q2_center[Q2_bin_plot[i]],kGreen,20,2);

        ratio_CT18[i] = Ratio_Theory(JAM4EIC[i],CT18ptxg[i],kBlack,20,2);
        ratio_NNPDF[i] = Ratio_Theory(JAM4EIC[i],NNPDF31_nnlo_pch_as_0118_SF[i],kBlue,20,2);
        ratio_JAM[i] = Ratio_Theory(JAM4EIC[i],JAM4EIC[i],kGreen,20,2);
    }

    //Read in Bin-Centering Correction Files
    read_bc("input_files/JAM4EIC_p_bin_centering.txt",1);
    read_bc("input_files/CT18ptxg_bin_centering.txt",2);
    read_bc("input_files/NNPDF31_nnlo_pch_as_0118_SF_bin_centering.txt",3);

    //---------------------------------
    //Analyse PYTHIA6 simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *t = new TChain("EICTree");
    for(int i=0;i<15;i++){
        t->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/yellow_report/18_275_nlopdf/ep_18_275_newtune_%d.root",i));
    }

    t->SetBranchAddress("event",&event);
    double Q2_event; //Get this from scattered electron (may differ from true if rad. effect included)
    double x_event; //''
    double y_event; //''
    double W2_event; //''

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    double cross_tot = 1.52E9; //Total Cross Section in fb
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
            bc_fac_b[i][j] = get_bc_fac(x_center[i],Q2_center[j],2);
            bc_fac_c[i][j] = get_bc_fac(x_center[i],Q2_center[j],3);
            if(yield[i][j]>0) { 
                error[i][j] = (rcs[i][j]) / (sqrt(yield[i][j]));
                error_10[i][j] = 100.*( 1./ sqrt(yield[i][j]*(10./lum)) ); //% Error for 10 fb^-1
            }
        }
    }

    //Create and fill TGraphs
    TGraphErrors *gr1[nbins_Q2]; //Without Bin-Centering Correction
    TGraphErrors *gr2[nbins_Q2]; //With Bin-Centering Correction
    TGraph *gbc1[nbins_Q2],*gbc2[nbins_Q2],*gbc3[nbins_Q2]; //Bin-Centering Correction Factor
    TGraph *ge1[nbins_Q2]; //Expected Statistical Errors (%) for L=10 fb^-1

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

        gbc1[j] = new TGraph(); //Correction comes from JAMEIC
        gbc1[j]->SetLineWidth(2);
        gbc1[j]->SetLineColor(kGreen);
        gbc1[j]->SetMarkerStyle(21);
        gbc1[j]->SetMarkerColor(kGreen);

        gbc2[j] = new TGraph(); //Correction comes from CTEQ
        gbc2[j]->SetLineWidth(2);
        gbc2[j]->SetLineColor(kBlack);
        gbc2[j]->SetMarkerStyle(20);
        gbc2[j]->SetMarkerColor(kBlack);

        gbc3[j] = new TGraph(); //Correction comes from NNPDF
        gbc3[j]->SetLineWidth(2);
        gbc3[j]->SetLineColor(kBlue);
        gbc3[j]->SetMarkerStyle(22);
        gbc3[j]->SetMarkerColor(kBlue);

        ge1[j] = new TGraph();
        ge1[j]->SetLineWidth(2);
        ge1[j]->SetLineColor(kRed);
        ge1[j]->SetMarkerStyle(20);
        ge1[j]->SetMarkerColor(kRed);
    }
    
    int counter(0);
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){
            if(yield[i][j]>50 && cut_lr[i][j] && cut_ul[i][j]){
                gr1[j]->SetPoint(counter,x_center[i],rcs[i][j]);
                gr1[j]->SetPointError(counter,0,error[i][j]);

                gbc1[j]->SetPoint(counter,x_center[i],bc_fac_a[i][j]);
                gbc2[j]->SetPoint(counter,x_center[i],bc_fac_b[i][j]);
                gbc3[j]->SetPoint(counter,x_center[i],bc_fac_c[i][j]);

                ge1[j]->SetPoint(counter,x_center[i],error_10[i][j]);

                gr2[j]->SetPoint(counter,x_center[i],rcs[i][j]*bc_fac_a[i][j]);
                gr2[j]->SetPointError(counter,0,error[i][j]*bc_fac_a[i][j]);
                
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
    for(int i=0;i<10;i++){
        t_d->Add(Form("/eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/18_275/djangoh.NC.18x275_%d.root",i));
    }

    t_d->SetBranchAddress("event",&event_d);

    //Calculate Generated Luminosity
    int nevents_d = t_d->GetEntries();
    double cross_tot_d = 1.51E9; //Total Cross Section in fb
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
        y_event = event_d->GetY();
        W2_event = event_d->GetW2();

        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_event>x_low[i] && x_event<x_hi[i] && Q2_event>Q2_low[j] && Q2_event<Q2_hi[j] &&
                    y_event>ymin && W2_event>W2min && y_event<ymax) //djangoh has a hard cut W>3GeV applied for Rad=0 unpol. NC)
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

    TLatex *tex4 = new TLatex(3E-4,1.3,Form("Djangoh e^{-}p, %.0e<y<%.2f and W^{2}>%.0fGeV^{2}",ymin,ymax,W2min));
    tex4->SetTextColor(kMagenta);tex4->SetTextFont(62);

    TLatex *tex5 = new TLatex(1E-3,1.45,"CT18ptxg");
    tex5->SetTextColor(kBlack);tex5->SetTextFont(62);

    TLatex *tex6 = new TLatex(1E-3,1.6,"NNPDF31_nnlo_pch_as_0118_SF");
    tex6->SetTextColor(kBlue);tex6->SetTextFont(62);

    TLatex *tex7 = new TLatex(1E-3,1.3,"JAM4EIC");
    tex7->SetTextColor(kGreen);tex7->SetTextFont(62);

    TLatex *tex7a = new TLatex(2.5E-2,0.99,"JAM4EIC");
    tex7a->SetTextColor(kGreen);tex7a->SetTextFont(62);

    TLatex *tex8 = new TLatex(1E-2,0.075,"Pythia6");
    tex8->SetTextColor(kRed);tex8->SetTextFont(62);

    TLatex *tex8a = new TLatex(1E-2,0.25,"Pythia6");
    tex8a->SetTextColor(kRed);tex8a->SetTextFont(62);

    TLatex *tex9a = new TLatex(2E-4,0.94,"Pythia6");
    tex9a->SetTextColor(kRed);tex9a->SetTextFont(62);

    TLatex *tex9b = new TLatex(2E-4,0.92,"Djangoh");
    tex9b->SetTextColor(kMagenta);tex9a->SetTextFont(62);

    TLatex *tex9c = new TLatex(2E-4,0.90,"JAM4EIC");
    tex9c->SetTextColor(kGreen);tex9c->SetTextFont(62);

    TLatex *tex9d = new TLatex(2E-4,0.88,"CT18ptxg");
    tex9d->SetTextColor(kBlack);tex9d->SetTextFont(62);

    TLatex *tex9e = new TLatex(2E-4,0.86,"NNPDF31_nnlo_pch_as_0118_SF");
    tex9e->SetTextColor(kBlue);tex9e->SetTextFont(62);

    TLatex *tex10a = new TLatex(2E-4,0.995,"JAM4EIC");
    tex10a->SetTextColor(kGreen);tex10a->SetTextFont(62);

    TLatex *tex10b = new TLatex(2E-4,0.994,"CT18ptxg");
    tex10b->SetTextColor(kBlack);tex10b->SetTextFont(62);

    TLatex *tex10c = new TLatex(2E-4,0.993,"NNPDF31_nnlo_pch_as_0118_SF");
    tex10c->SetTextColor(kBlue);tex10c->SetTextFont(62);

    //Data to Theory (JAM) Ratio
    TGraph *ratio_pythia[num_plot];  //Ratio (with bin-centering corr.) to JAM
    TGraph *ratio_djangoh[num_plot]; //''

    for(int i=0;i<num_plot;i++){
        ratio_pythia[i] = Ratio_Data(JAM4EIC[i],gr2[Q2_bin_plot[i]],kRed,20,2);
        ratio_djangoh[i] = Ratio_Data(JAM4EIC[i],gr2_d[Q2_bin_plot[i]],kMagenta,23,2);
    }

    //Bin-Centering Ratios
    TGraph *ratio_BC_CT18[num_plot];  //Ratio to JAM
    TGraph *ratio_BC_NNPDF[num_plot]; //''
    TGraph *ratio_BC_JAM[num_plot];   //''

    for(int i=0;i<num_plot;i++){
        ratio_BC_CT18[i] = Ratio_Theory(gbc1[Q2_bin_plot[i]],gbc2[Q2_bin_plot[i]],kBlack,22,2);
        ratio_BC_NNPDF[i] = Ratio_Theory(gbc1[Q2_bin_plot[i]],gbc3[Q2_bin_plot[i]],kBlue,23,2);
        ratio_BC_JAM[i] = Ratio_Theory(gbc1[Q2_bin_plot[i]],gbc1[Q2_bin_plot[i]],kGreen,20,2);
    }

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

    //Cross Section (no bin-centering)
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex5->Draw();tex6->Draw();tex7->Draw();

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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex5->Draw();tex6->Draw();tex7->Draw();
    
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    //Bin-centering correction factor (from JAM4EIC)
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
    TH1 *hframe_1a = gPad->DrawFrame(1E-4,0.7,1,1.1);
    hframe_1a->GetXaxis()->SetTitle("x");hframe_1a->GetXaxis()->CenterTitle();
    hframe_1a->GetYaxis()->SetTitle("Bin-Centering Correction");hframe_1a->GetYaxis()->CenterTitle();
    counter = 0;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();
    tex2a->Draw();tex7a->Draw();

    c1a->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1a->Draw();
    counter++;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();

    c1a->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1a->Draw();
    counter++;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();
    
    c1a->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1a->Draw();
    counter++;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();

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
    hframe_1a->Draw();
    counter++;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();
    tex2a->Draw();tex7a->Draw();

    c2a->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1a->Draw();
    counter++;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();
    
    c2a->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1a->Draw();
    counter++;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();
    
    c2a->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1a->Draw();
    counter++;
    gbc1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1a[counter]->Draw();

    //Bin-Centering Model Ratios
    //------------------------------------
    TCanvas *c1ar = new TCanvas("c1ar");
    c1ar->Divide(2,2,small,small);

    c1ar->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky(); 
    TH1 *hframe_1ar = gPad->DrawFrame(1E-4,0.99,1,1.01);
    hframe_1ar->GetXaxis()->SetTitle("x");hframe_1ar->GetXaxis()->CenterTitle();
    hframe_1ar->GetYaxis()->SetTitle("B.C. Correction Ratio");hframe_1ar->GetYaxis()->CenterTitle();
    hframe_1ar->GetYaxis()->SetTitleOffset(1.1);
    counter = 0;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();

    c1ar->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1ar->Draw();
    counter++;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();

    c1ar->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1ar->Draw();
    counter++;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();
    
    c1ar->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1ar->Draw();
    counter++;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();
    tex10a->Draw();tex10b->Draw();tex10c->Draw();

    TCanvas *c2ar = new TCanvas("c2ar");
    c2ar->Divide(2,2,small,small);

    c2ar->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1ar->Draw();
    counter++;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();

    c2ar->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1ar->Draw();
    counter++;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();
    
    c2ar->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1ar->Draw();
    counter++;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();

    c2ar->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1ar->Draw();
    counter++;
    ratio_BC_JAM[counter]->Draw("L Same");
    ratio_BC_CT18[counter]->Draw("P Same");ratio_BC_NNPDF[counter]->Draw("P Same");
    tex1f[counter]->Draw();
    tex10a->Draw();tex10b->Draw();tex10c->Draw();

    //Cross section including bin-centering correction
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex5->Draw();tex6->Draw();tex7->Draw();

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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();
    tex5->Draw();tex6->Draw();tex7->Draw();
    
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
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
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr2[Q2_bin_plot[counter]]->Draw("P Same");gr2_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    //Ratio Plots
    //------------------------------------
    TCanvas *c1r = new TCanvas("c1r");
    c1r->Divide(2,2,small,small);

    c1r->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky(); 
    TH1 *hframe_1r = gPad->DrawFrame(1E-4,0.8,1,1.2);
    hframe_1r->GetXaxis()->SetTitle("x");hframe_1r->GetXaxis()->CenterTitle();
    hframe_1r->GetYaxis()->SetTitle("Cross Section Ratio");hframe_1r->GetYaxis()->CenterTitle();
    counter = 0;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1d[counter]->Draw();

    c1r->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1r->Draw();
    counter++;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1d[counter]->Draw();

    c1r->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1r->Draw();
    counter++;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1d[counter]->Draw();
    
    c1r->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1r->Draw();
    counter++;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1d[counter]->Draw();
    tex9a->Draw();tex9b->Draw();tex9c->Draw();tex9d->Draw();tex9e->Draw();

    TCanvas *c2r = new TCanvas("c2r");
    c2r->Divide(2,2,small,small);

    c2r->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1r->Draw();
    counter++;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1e[counter]->Draw();

    c2r->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1r->Draw();
    counter++;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1e[counter]->Draw();
    
    c2r->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1r->Draw();
    counter++;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1e[counter]->Draw();
    
    c2r->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1r->Draw();
    counter++;
    ratio_JAM[counter]->Draw("L Same");ratio_CT18[counter]->Draw("L Same");ratio_NNPDF[counter]->Draw("L Same");
    ratio_pythia[counter]->Draw("P Same");ratio_djangoh[counter]->Draw("P Same");
    tex1e[counter]->Draw();
    tex9a->Draw();tex9b->Draw();tex9c->Draw();tex9d->Draw();tex9e->Draw();
    
    //Expected Statistical uncertainties (%) for 10 fb^-1 (pythia only)
    //------------------------------------
    TCanvas *c1c = new TCanvas("c1c");
    c1c->Divide(2,2,small,small);

    c1c->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky(); 
    TH1 *hframe_1c = gPad->DrawFrame(1E-4,0,1,0.1);
    hframe_1c->GetXaxis()->SetTitle("x");hframe_1c->GetXaxis()->CenterTitle();
    hframe_1c->GetYaxis()->SetTitle("#intL=10 fb^{-1} Stat. Uncer. [%]");hframe_1c->GetYaxis()->CenterTitle();
    counter = 0;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1b[counter]->Draw();
    tex2b->Draw();tex8->Draw();
    
    c1c->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1c->Draw();
    counter++;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1b[counter]->Draw();

    c1c->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1c->Draw();
    counter++;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1b[counter]->Draw();
    
    c1c->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1c->Draw();
    counter++;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1b[counter]->Draw();

    TCanvas *c2c = new TCanvas("c2c");
    c2c->Divide(2,2,small,small);

    c2c->cd(1);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    TH1 *hframe_1d = gPad->DrawFrame(1E-4,0,1,0.35);
    hframe_1d->GetXaxis()->SetTitle("x");hframe_1d->GetXaxis()->CenterTitle();
    hframe_1d->GetYaxis()->SetTitle("#intL=10 fb^{-1} Stat. Uncer. [%]");hframe_1d->GetYaxis()->CenterTitle();
    counter++;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1c[counter]->Draw();
    tex2c->Draw();tex8a->Draw();

    c2c->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1d->Draw();
    counter++;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1c[counter]->Draw();
    
    c2c->cd(3);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(0.25);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1d->Draw();
    counter++;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1c[counter]->Draw();
    
    c2c->cd(4);
    gPad->SetLogx();
    gPad->SetTopMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe_1d->Draw();
    counter++;
    ge1[Q2_bin_plot[counter]]->Draw("P Same");
    tex1c[counter]->Draw();

    //Print to File
    c1->Print("plots/theory_compare.pdf[");
    c1->Print("plots/theory_compare.pdf");
    c2->Print("plots/theory_compare.pdf");
    c1a->Print("plots/theory_compare.pdf");
    c2a->Print("plots/theory_compare.pdf");
    c1ar->Print("plots/theory_compare.pdf");
    c2ar->Print("plots/theory_compare.pdf");
    c1b->Print("plots/theory_compare.pdf");
    c2b->Print("plots/theory_compare.pdf");
    c1r->Print("plots/theory_compare.pdf");
    c2r->Print("plots/theory_compare.pdf");
    c1c->Print("plots/theory_compare.pdf");
    c2c->Print("plots/theory_compare.pdf");
    c2c->Print("plots/theory_compare.pdf]");

}