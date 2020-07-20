R__LOAD_LIBRARY(libeicsmear);

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

void theory_compare(){

    //Constants
    const double alpha(1./137.036);
    const double s_cm = 4.*5*41;
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
            error[i][j]=0;
            rcs[i][j]=0;  

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
    TGraph *CT18ptxg[nbins_Q2];
    TGraph *NNPDF31_nnlo_pch_as_0118_SF[nbins_Q2];
    TGraph *JAM4EIC[nbins_Q2];

    for(int i=0;i<num_plot;i++){
        CT18ptxg[i] = Theory_Graph("input_files/CT18ptxg_rcs.txt",Q2_center[Q2_bin_plot[i]],kBlack,20,2);
        NNPDF31_nnlo_pch_as_0118_SF[i] = Theory_Graph("input_files/NNPDF31_nnlo_pch_as_0118_SF_rcs.txt",Q2_center[Q2_bin_plot[i]],kBlue,20,2);
        JAM4EIC[i] = Theory_Graph("input_files/JAM4EIC_p_rcs.txt",Q2_center[Q2_bin_plot[i]],kGreen,20,2);
    }

    //---------------------------------
    //Analyse PYTHIA6 simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *t = new TChain("EICTree");
    for(int i=0;i<15;i++){
        t->Add(Form("/eic/data/baraks/pythiaeRHIC/outfiles/yellow_report/5_41/ep_5_41_newtune_%d.root",i));
    }

    t->SetBranchAddress("event",&event);
    double Q2_event; //Get this from scattered electron (may differ from true if rad. effect included)
    double x_event; //''
    double y_event; //''
    double W2_event; //''

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    double cross_tot = 0.716E9; //Total Cross Section in fb
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
            if(yield[i][j]>50 && cut_lr[i][j] && cut_ul[i][j]){
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
    for(int i=0;i<10;i++){
        t_d->Add(Form("/eic/data/baraks/dis-reconstruction/djangoh/outfiles/yellow/5_41/djangoh.NC.5x41_%d.root",i));
    }

    t_d->SetBranchAddress("event",&event_d);

    //Calculate Generated Luminosity
    int nevents_d = t_d->GetEntries();
    double cross_tot_d = 0.741E9; //Total Cross Section in fb
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

    TLatex *tex2 = new TLatex(7E-4,1.6,"5 GeV e^{-} on 41 GeV p, #sqrt{s}=28.6 GeV");
    tex2->SetTextColor(kBlack);tex2->SetTextFont(42);

    TLatex *tex3 = new TLatex(7E-4,1.45,"Pythia6 e^{-}p, y>0.001 and W^{2}>10GeV^{2}");
    tex3->SetTextColor(kRed);tex3->SetTextFont(62);

    TLatex *tex4 = new TLatex(7E-4,1.3,"Djangoh e^{-}p, y>0.001 and W^{2}>10GeV^{2}");
    tex4->SetTextColor(kMagenta);tex4->SetTextFont(62);

    TLatex *tex5 = new TLatex(1E-3,1.45,"CT18ptxg");
    tex5->SetTextColor(kBlack);tex5->SetTextFont(62);

    TLatex *tex6 = new TLatex(1E-3,1.6,"NNPDF31_nnlo_pch_as_0118_SF");
    tex6->SetTextColor(kBlue);tex6->SetTextFont(62);

    TLatex *tex7 = new TLatex(1E-3,1.3,"JAM4EIC");
    tex7->SetTextColor(kGreen);tex7->SetTextFont(62);

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
    hframe->Draw();
    counter++;
    CT18ptxg[counter]->Draw("C Same");JAM4EIC[counter]->Draw("C Same");
    NNPDF31_nnlo_pch_as_0118_SF[counter]->Draw("C Same");
    gr1[Q2_bin_plot[counter]]->Draw("P Same");gr1_d[Q2_bin_plot[counter]]->Draw("P Same");
    tex1[counter]->Draw();

    //Print to File
    c1->Print("plots/theory_compare.pdf[");
    c1->Print("plots/theory_compare.pdf");
    c2->Print("plots/theory_compare.pdf");
    c2->Print("plots/theory_compare.pdf]");

}