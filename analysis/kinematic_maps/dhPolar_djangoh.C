R__LOAD_LIBRARY(libeicsmear);

void dhPolar_djangoh(){

	//Label indicating the particle ID. Normalize your counts to an integrated luminosity of 10 fb-1
  	//TString title("Electron Momentum vs Theta (#intL=10 fb^{-1})"); // Title of your plot
	//TString title("Proton Momentum vs Theta (#intL=10 fb^{-1})"); // Title of your plot
	//TString title("Photon Momentum vs Theta (#intL=10 fb^{-1})"); // Title of your plot
	//TString title("Neutron Momentum vs Theta (#intL=10 fb^{-1})"); // Title of your plot
	TString title("Hadron Momentum vs Theta (#intL=10 fb^{-1})"); // Title of your plot

	//Load ROOT Files
	TChain *tree = new TChain("EICTree");

	//Read in Beam Energies
  	Int_t energy_set(0);

  	cout<<"Enter Beam Energy Combination:"<<endl;
  	cout<<"1) e:5 GeV; p:41 GeV"<<endl;
  	cout<<"2) e:5 GeV; p:100 GeV"<<endl;
  	cout<<"3) e:10 GeV; p:100 GeV"<<endl;
  	cout<<"4) e:18 GeV; p:275 GeV"<<endl;
  	cin>>energy_set;

	//5x41
	if(energy_set==1){
		for(Int_t i=0;i<10;i++){
			tree->Add(Form("../../djangoh/outfiles/yellow/5_41/djangoh.NC.5x41_%d.root",i));
    	}
	}

	//5x100
	if(energy_set==2){
		for(Int_t i=0;i<10;i++){
			tree->Add(Form("../../djangoh/outfiles/yellow/5_100/djangoh.NC.5x100_%d.root",i));
		}
	}

	//10x100
	if(energy_set==3){
		for(Int_t i=0;i<10;i++){
    		tree->Add(Form("../../djangoh/outfiles/yellow/10_100/djangoh.NC.10x100_%d.root",i));
  		}
	}

	//18x275
	if(energy_set==4){
		for(Int_t i=0;i<10;i++){
			tree->Add(Form("../../djangoh/outfiles/yellow/18_275/djangoh.NC.18x275_%d.root",i));
    	}
	}

	//Print out total number of events
	int nevents = tree->GetEntries();
	cout << "Total Number of events = "<< nevents <<endl;

	//const float rmax = 50; // momentum range (for electrons;5x41 protons,photons,neutrons,hadrons)
	//int nbins = (int) rmax*2;

	//const float rmax = 110; // momentum range (5x100 and 10x100 for protons,neutrons,hadrons)
	//int nbins = (int) rmax;

	//const float rmax = 100; // momentum range (5x100 and 10x100 for photons)
	//int nbins = (int) rmax;

	const float rmax = 300; // momentum range (18x275 for protons,photons,neutrons,hadrons)
	int nbins = (int) rmax/3.;
	
	int rmax_int = (int) rmax;

	//Please keep the following binning so that we can compare across different WG
	TH2D *h=new TH2D("h",title.Data(),36,0.,TMath::Pi(),nbins,0.,rmax_int);
 
	//Fill your 2D histogram with Theta (in deg) and Momenta (in GeV) for a given particle ID
	TCanvas *ctemp = new TCanvas("ctemp");
	//tree->Draw("particles.p:particles.theta>>h","particles.id==11 && particles.KS==1"); //All Final-State Electrons
	//tree->Draw("particles.p:particles.theta>>h","particles.id==11 && particles.orig==1 && particles.KS==1"); //Scattered Electron Only
	//tree->Draw("particles.p:particles.theta>>h","particles.id==2212 && particles.KS==1"); //All Final-State Protons
	//tree->Draw("particles.p:particles.theta>>h","particles.id==22 && particles.KS==1"); //All Final-State Photons
	//tree->Draw("particles.p:particles.theta>>h","particles.id==2112 && particles.KS==1"); //All Final-State Neutrons
	tree->Draw("particles.p:particles.theta>>h","fabs(particles.id)>110 && particles.KS==1"); //All Final-State Hadrons

	//10 GeV^2 < Q2 < 100 GeV^2
	//tree->Draw("particles.p:particles.theta>>h","particles.id==11 && particles.KS==1 && QSquared>10 && QSquared<100"); //All Final-State Electrons
	//tree->Draw("particles.p:particles.theta>>h","particles.id==11 && particles.orig==1 && particles.KS==1 && QSquared>10 && QSquared<100"); //Scattered Electron Only
	//tree->Draw("particles.p:particles.theta>>h","particles.id==2212 && particles.KS==1 && QSquared>10 && QSquared<100"); //All Final-State Protons
	//tree->Draw("particles.p:particles.theta>>h","particles.id==22 && particles.KS==1 && QSquared>10 && QSquared<100"); //All Final-State Photons
	//tree->Draw("particles.p:particles.theta>>h","particles.id==2112 && particles.KS==1 && QSquared>10 && QSquared<100"); //All Final-State Neutrons
	//tree->Draw("particles.p:particles.theta>>h","fabs(particles.id)>110 && particles.KS==1 && QSquared>10 && QSquared<100"); //All Final-State Hadrons

	double cross_section(0); //In fb
	double sf_10(0); //For 10 fb^-1 integrated luminosity

	//5x41
	if(energy_set==1){
		cross_section = 0.741E9;
	}
	
	//5x100
	if(energy_set==2){
		cross_section = 0.908E9;
	}
	
	//10x100
	if(energy_set==3){
		cross_section = 1.06E9;
	}
 
	//18x275
	if(energy_set==4){
		cross_section = 1.51E9;
	}

	sf_10 = 10. / ( (double) nevents / cross_section );
	cout<<"Scale Factor is "<<sf_10<<endl;
	h->Scale(sf_10);

	// --------- Plotting ----------------
	//const float cmax = 400; // pad height
	//TCanvas *c1 = new TCanvas("c1","", 2.4*cmax, 1.25*cmax); 
	TCanvas *c1 = new TCanvas("c1");
	c1->SetFrameLineColor(0);
	gPad->SetBottomMargin(0.15); // leave more space for labels
	gPad->SetTopMargin(0.1);  // uniform up, left/right
	gPad->SetRightMargin(0.1); 
	gPad->SetLeftMargin(0.1); 
	gPad->SetBorderSize(0);  // no border

	// prepare a transparent frame
	TH1F *frame = gPad->DrawFrame(-1.*rmax,-0., rmax, rmax); 
	frame->GetXaxis()->SetAxisColor(0); 
	frame->GetYaxis()->SetAxisColor(0); 
	frame->GetXaxis()->SetLabelOffset(999); 
	frame->GetXaxis()->SetLabelColor(0); 
	frame->GetXaxis()->SetTickLength(0); 
	frame->GetYaxis()->SetLabelOffset(999); 
	frame->GetYaxis()->SetLabelColor(0); 
	frame->GetYaxis()->SetTickLength(0); 
	gPad->SetLogz();

	//polar frame 
	TEllipse *el3 = new TEllipse(0.,0.,rmax,rmax,0,180);
	el3->SetFillColor(0); // empty 
	el3->SetLineColor(1);
	el3->SetLineWidth(2);
	el3->Draw("");

	//draw theta label 
	TLatex *latex = new TLatex(); 
	latex->SetTextSize(0.06);
	latex->SetTextFont(42);
	latex->DrawLatex(1.01*0.5*sqrt(3)*rmax,1.05*0.5*sqrt(1)*rmax,"30#circ");
	latex->DrawLatex(0.98*0.5*sqrt(1)*rmax,1.05*0.5*sqrt(3)*rmax,"60#circ");
	latex->DrawLatex(-2.5,1.04*rmax,"90#circ");
	latex->DrawLatex(-1.2*0.5*sqrt(1)*rmax,1.05*0.5*sqrt(3)*rmax,"120#circ");
	latex->DrawLatex(-1.18*0.5*sqrt(3)*rmax,1.05*0.5*sqrt(1)*rmax,"150#circ");
	latex->SetTextSize(0.06);
	latex->DrawLatex(-1.1*rmax,0.02*rmax,"#theta");

	TLatex *latex_a = new TLatex();
	latex_a->SetTextSize(0.035);
	latex_a->SetTextFont(42);
	latex_a->DrawLatex(-1.25*rmax,1.025*rmax,title.Data());
	if(energy_set==1)
		latex_a->DrawLatex(-1.25*rmax,0.95*rmax,"5 GeV e^{-} on 41 GeV p");
	if(energy_set==2)
		latex_a->DrawLatex(-1.25*rmax,0.95*rmax,"5 GeV e^{-} on 100 GeV p");
	if(energy_set==3)
		latex_a->DrawLatex(-1.25*rmax,0.95*rmax,"10 GeV e^{-} on 100 GeV p");
	if(energy_set==4)
		latex_a->DrawLatex(-1.25*rmax,0.95*rmax,"18 GeV e^{-} on 275 GeV p");

	latex_a->DrawLatex(-1.25*rmax,0.875*rmax,"Q^{2} > 0.5 GeV^{2}");
	//latex_a->DrawLatex(-1.25*rmax,0.875*rmax,"10 GeV^{2} < Q^{2} < 100 GeV^{2}");

	//draw theta ticks 
	TLine *line = new TLine(); 
	line->DrawLine(0.5*sqrt(3)*rmax, 0.5*sqrt(1)*rmax, 1.03*0.5*sqrt(3)*rmax, 1.03*0.5*sqrt(1)*rmax); 
	line->DrawLine(0.5*sqrt(1)*rmax, 0.5*sqrt(3)*rmax, 1.03*0.5*sqrt(1)*rmax, 1.03*0.5*sqrt(3)*rmax); 
	line->DrawLine(0,rmax,0,1.02*rmax); 
	line->DrawLine(-0.5*sqrt(1)*rmax, 0.5*sqrt(3)*rmax, -1.03*0.5*sqrt(1)*rmax, 1.03*0.5*sqrt(3)*rmax); 
	line->DrawLine(-0.5*sqrt(3)*rmax, 0.5*sqrt(1)*rmax, -1.03*0.5*sqrt(3)*rmax, 1.03*0.5*sqrt(1)*rmax); 

	// positive momentum axis
	TGaxis *gaxis = new TGaxis(0,0,rmax,0,0,rmax,505,"-+");
	gaxis->SetLabelSize(0.06);
	gaxis->SetTitle("Momentum (GeV)"); 
	gaxis->SetTitleSize(0.06);
	gaxis->SetTitleFont(42);
	gaxis->SetLabelFont(42);
	gaxis->Draw();

	//negative momentum axis
	TF1 *f1=new TF1("f1","-x",0,rmax);
	TGaxis *A1 = new TGaxis(-1.*rmax,0,0,0,"f1",505,"-+");
 	A1->SetLabelSize(0.06);
    	A1->SetTitleSize(0.06);
    	A1->SetTitleFont(42);
    	A1->SetLabelFont(42);
   	A1->Draw("same");
	
	// proton/ion going direction
	//TArrow *arrow = new TArrow(-1.*rmax, -0.12*rmax, -0.5*rmax, -0.12*rmax, 0.02,"|>");
	//arrow->SetLineWidth(2);
	//arrow->SetAngle(30);
	//arrow->Draw();
	//latex->SetTextSize(0.06);
	//latex->DrawLatex(-0.8*rmax, -0.09*rmax, "p/A");
 
  	h->Draw("POL COLZ SAME");
	if(energy_set==1)
		c1->Print("yellow_ps/dhPolar_5_41.pdf");
	if(energy_set==2)
		c1->Print("yellow_ps/dhPolar_5_100.pdf");
	if(energy_set==3)
		c1->Print("yellow_ps/dhPolar_10_100.pdf");
	if(energy_set==4)
		c1->Print("yellow_ps/dhPolar_18_275.pdf");
	
	//Write out text file
	ofstream myfile;
	if(energy_set==1)
		myfile.open("yellow_ps/dhPolar_5_41.txt");
	if(energy_set==2)
		myfile.open("yellow_ps/dhPolar_5_100.txt");
	if(energy_set==3)
		myfile.open("yellow_ps/dhPolar_10_100.txt");
	if(energy_set==4)
		myfile.open("yellow_ps/dhPolar_18_275.txt");


	myfile <<"Angle [Rad]  |  Bin Width [Rad]  |  Momentum [GeV/c]  |  Bin Width[GeV/c]  |  Events/10fb^-1"<<endl;
	myfile << "============================================"<<endl<<endl;
	for(int i=1;i<=h->GetNbinsX();i++){
		for(int j=1;j<=h->GetNbinsY();j++){
	 		if((h->GetBinContent(i,j))!=0){
				myfile << setw(5) << setprecision(4) << h->GetXaxis()->GetBinCenter(i) << 
				setw(10) << setprecision(4) << h->GetXaxis()->GetBinWidth(i) << 
				setw(10) << setprecision(4) << h->GetYaxis()->GetBinCenter(j) << 
				setw(10) << setprecision(4) << h->GetYaxis()->GetBinWidth(i) << 
				setw(15) << setprecision(4) << h->GetBinContent(i,j) << endl;
	 		}
		}
	}
	myfile << "============================================"<<endl;

	myfile.close();

}
