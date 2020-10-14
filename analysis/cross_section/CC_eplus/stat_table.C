R__LOAD_LIBRARY(libeicsmear);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

void stat_table(){

    //Constants
    const double alpha(1./137.036);
    const double s_cm = 4.*18*275;
    const double fbgev(1./(0.3894E12));
    const double Mp(0.9383);

    //Cross Section Bins (Xiaoxuan's binning)
	const int nbins_Q2 = 10;
	double Q2_low[] =   {100, 120, 190, 390, 690, 990,  1400, 2900, 4800, 8000};
	double Q2_hi[] =    {120, 140, 210, 410, 710, 1010, 1600, 3100, 5200, 11000};
    double Q2_print[] = {110, 130, 200, 400, 700, 1000, 1400, 3000, 5000, 9000};
   
    double x_min = 1E-3;
    double x_max = 1;
	const int nbins = 20;
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
    double error_10[nbins][nbins_Q2];
    double Q2_center[nbins_Q2];
    double x_center[nbins];
    double Q2_width[nbins_Q2];
    double x_width[nbins];

    for(int i=0;i<nbins;i++){ //x loop
        for(int j=0;j<nbins_Q2;j++){ //Q2 loop
            yield[i][j]=0;  
            error_10[i][j]=0;

            Q2_center[j] = (Q2_hi[j]+Q2_low[j])/2.;
            x_center[i] = (x_hi[i]+x_low[i])/2.;
            Q2_width[j] = (Q2_hi[j]-Q2_low[j]);
            x_width[i] = (x_hi[i]-x_low[i]);
        }
    }

    //Output files
    ofstream myfile_good;
    myfile_good.open("data_CC_eplus.out");

    //---------------------------------
    //Analyse Djangoh simulation
    //---------------------------------

    //Load ROOT Files
    erhic::EventDjangoh *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *t = new TChain("EICTree");
    //t->Add("/gpfs02/eic/baraks/djangoh/outfiles/djangoh.CC.Rad.18x275_evt.root"); //Rad=1 File
    t->Add("/gpfs02/eic/baraks/djangoh/outfiles/djangoh.CC.18x275_evt.root"); //Rad=0 File
    

    t->SetBranchAddress("event",&event);
    double Q2_event; //Xiaoxuan used the true value (not scattered electron); so we'll do that
    double x_event; //''

    //Calculate Generated Luminosity
    int nevents = t->GetEntries();
    //double cross_tot = 8.376E3; //Total Cross Section in fb (Rad=1 File)
    double cross_tot = 9.057E3; //Total Cross Section in fb (Rad=0 File)
    double lum = ( (double) nevents)/cross_tot; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"DJANGOH Simulation:"<<endl;
    cout<<"Total Number of Events = "<<nevents<<endl;
    cout<<"Integrated Luminosity = "<<lum<<" fb^-1"<<endl<<endl;

    //Loop over events
    for(int k=0;k<nevents;k++){
        if(k%10000==0) cout<<"Events Analysed = "<<k<<"!"<<endl;
        t->GetEntry(k);

        Q2_event = event->dtrueQ2;
        x_event = event->dtrueX;

        for(int i=0;i<nbins;i++){
            for(int j=0;j<nbins_Q2;j++){
                if(x_event>x_low[i] && x_event<x_hi[i] && Q2_event>Q2_low[j] && Q2_event<Q2_hi[j]) 
                    yield[i][j]+=1.;
            }
        }
    }

    //Print table
    for(int j=0;j<nbins_Q2;j++){
        for(int i=0;i<nbins;i++){

            if(yield[i][j]>0) { 
                error_10[i][j] = 100.*( 1./ sqrt( yield[i][j]*(10./lum) ) ); //% Error for 10 fb^-1

                myfile_good << setprecision(3) << Q2_print[j] << setw(15) << setprecision(3) << x_center[i]
                << setw(15) << setprecision(3) << error_10[i][j] << endl;

            }
        }
    }
    //Close output files
    myfile_good.close();
}
