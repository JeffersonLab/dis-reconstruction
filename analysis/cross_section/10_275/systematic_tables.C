R__LOAD_LIBRARY(libeicsmear);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

void systematic_tables(){

    //Constants
    const double alpha(1./137.036);
    const double s_cm = 4.*10*275;
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
    double yield[nbins][nbins_Q2];
    double error[nbins][nbins_Q2];
    double error_10[nbins][nbins_Q2];
    double rcs[nbins][nbins_Q2];
    double y_center[nbins][nbins_Q2];
    double Q2_center[nbins_Q2];
    double x_center[nbins];
    double Q2_width[nbins_Q2];
    double x_width[nbins];
    double W2_center[nbins][nbins_Q2];
    double rcs_factor[nbins][nbins_Q2];
    bool cut_ul[nbins][nbins_Q2];
    bool cut_lr[nbins][nbins_Q2];
    bool cut_ul_bad[nbins][nbins_Q2];
    bool cut_lr_bad[nbins][nbins_Q2];
    double ptp_err[nbins][nbins_Q2];
    double norm_err[nbins][nbins_Q2];
    double ptp_err_bad[nbins][nbins_Q2];
    double norm_err_bad[nbins][nbins_Q2];

    double ymax = 0.98; double ymin = 0.001; double W2min = 0.; //10;
    double ymax_bad = 0.95; double ymin_bad = 0.01; double W2min_bad = 0.; double Q2_min_bad = 1.;
    double y_temp, W2_temp;

    for(int i=0;i<nbins;i++){ //x loop
        for(int j=0;j<nbins_Q2;j++){ //Q2 loop
            yield[i][j]=0;  
            error[i][j]=0;
            error_10[i][j]=0;
            rcs[i][j]=0;

            Q2_center[j] = (Q2_hi[j]+Q2_low[j])/2.;
            x_center[i] = (x_hi[i]+x_low[i])/2.;
            Q2_width[j] = (Q2_hi[j]-Q2_low[j]);
            x_width[i] = (x_hi[i]-x_low[i]);

            y_center[i][j] = Q2_center[j]/(x_center[i]*s_cm);
            W2_center[i][j] = Mp*Mp + Q2_center[j]*( (1./x_center[i]) - 1.);

            rcs_factor[i][j] = fbgev * ( pow(Q2_center[j],2)*x_center[i] ) / 
                            ( TMath::TwoPi() * pow(alpha,2) * ( 1.+pow((1.-y_center[i][j]),2) ) );

            //Cut Bins with top-left edge > ymax and bottom-right edge < ymin and > W2min
            y_temp = (Q2_hi[j])/(x_low[i]*s_cm); //Upper left
            if(y_temp>ymax) //Remove bin
                cut_ul[i][j] = 0;
            else
                cut_ul[i][j] = 1;

            if(y_temp>ymax_bad) //Remove bin
                cut_ul_bad[i][j] = 0;
            else
                cut_ul_bad[i][j] = 1;

            y_temp = (Q2_low[j])/(x_hi[i]*s_cm); //Lower Right
            W2_temp = Mp*Mp + Q2_low[j]*( (1./x_hi[i]) - 1.);

            if(y_temp<ymin || W2_temp<W2min) //Remove bin
                cut_lr[i][j] = 0;
            else
                cut_lr[i][j] = 1;

            if(y_temp<ymin_bad || W2_temp<W2min_bad) //Remove bin
                cut_lr_bad[i][j] = 0;
            else
                cut_lr_bad[i][j] = 1;

            //Set systematic uncertainties
            if(y_center[i][j]<0.01){
                ptp_err[i][j] = 2.5;
                ptp_err_bad[i][j] = 3.0;
            }
            else{
                ptp_err[i][j] = 1.5;
                ptp_err_bad[i][j] = 2.3;
            }
            norm_err[i][j] = 2.5;
            norm_err_bad[i][j] = 4.3;
        }

    }

    //Output files
    ofstream myfile_good;
    myfile_good.open("data_optimistic.out");

    ofstream myfile_bad;
    myfile_bad.open("data_pessimistic.out");

    //---------------------------------
    //Analyse PYTHIA6 simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *t = new TChain("EICTree");
    for(int i=0;i<50;i++){
        t->Add(Form("/gpfs02/eic/baraks/pythia/outfiles/other_studies/10_275/ep_10_275_newtune_%d.root",i));
    }

    t->SetBranchAddress("event",&event);
    double Q2_event; //Get this from scattered electron (may differ from true if rad. effect included)
    double x_event; //''
    double y_event; //''
    double W2_event; //''

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    double cross_tot = 0.841E9; //Total Cross Section in fb
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
    int shift(0);

    for(int i=0;i<nbins;i++){
        for(int j=0;j<nbins_Q2;j++){
            rcs[i][j] = rcs_factor[i][j] * (yield[i][j]) / (lum*x_width[i]*Q2_width[j]);
            shift=0;

            if(yield[i][j]>0) { 
                error[i][j] = (rcs[i][j]) / (sqrt(yield[i][j]));
                //Percent error for 100 fb^-1. Note we do 100 fb^-1 here since it's 
                //the maximum luminosity and what we use for ATHENA. For the other energies,
                //we extracted for 10fb^-1
                error_10[i][j] = 100.*( 1./ sqrt(yield[i][j]*(100./lum)) );

                if(x_center[i]>1e-3)
                    shift++;

                if(x_center[i]>1e-2)
                    shift++;
                
                if(x_center[i]>1e-1)
                    shift++;

                //Print to files
                if(yield[i][j]>1 && cut_lr[i][j] && cut_ul[i][j]){
                    myfile_good << setprecision(3) << x_center[i] << setw(15+shift) << setprecision(3) << Q2_center[j] 
                    << setw(15) << setprecision(3) << y_center[i][j] << setw(15) << setprecision(4)
                    << W2_center[i][j] << setw(15) << setprecision(3) << error_10[i][j] 
                    << setw(10) << ptp_err[i][j] << setw(10) << norm_err[i][j] << endl;
                }
                if(yield[i][j]>1 && cut_lr_bad[i][j] && cut_ul_bad[i][j] && Q2_center[j]>Q2_min_bad){
                    myfile_bad << setprecision(3) << x_center[i] << setw(15+shift) << setprecision(3) << Q2_center[j] 
                    << setw(15) << setprecision(3) << y_center[i][j] << setw(15) << setprecision(4)
                    << W2_center[i][j] << setw(15) << setprecision(3) << error_10[i][j] 
                    << setw(10) << ptp_err_bad[i][j] << setw(10) << norm_err_bad[i][j] << endl;
                }
            }
        }
    }
    //Close output files
    myfile_good.close();
    myfile_bad.close();
}