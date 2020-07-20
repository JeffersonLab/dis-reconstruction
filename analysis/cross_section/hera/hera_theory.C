//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

//Read in data file and put into graph (TGraphErrors)
//Assumes file is in format x,y,yerr
TGraphErrors* Data_Graph(const char* infile, double rs_set=0,double Q2_set=1,
                        int color=1, int style=20, int lw=2, bool abs=0){

  //Create graph
  TGraphErrors *gr = new TGraphErrors();
  gr->SetMarkerStyle(style);gr->SetMarkerColor(color);
  gr->SetLineWidth(lw);gr->SetLineColor(color);

  //Read input file
  ifstream in;
  in.open(infile,ios::in);

  //Write to Graph
  double x(0),y(0),ey(0);
  double rs(0),Q2(0);
  int counter(0);

  while(in.good()){
    in >> Q2 >> x >> y >> ey >> rs;

    if(in.good() && fabs(rs-rs_set)<1E-2 && fabs(Q2-Q2_set)<1E-2){

      if(abs)
        gr->SetPoint(counter,x,fabs(y));
      else
        gr->SetPoint(counter,x,y);
      
      gr->SetPointError(counter,0,y*(ey/100.));

      counter++;
    }
  }

  cout<<"For Q2 = "<<Q2_set<<" GeV^2... ";
  cout<<"We read "<<counter<<" matching lines from data file "<<infile<<" !"<<endl;
  in.close();
  return gr;

}

//Read in data file and put into graph (TGraph)
//Assumes file is in format Q2,x,y
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

    if( in.good() && fabs(Q2-Q2_set)<1E-2 ){

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

void hera_theory(){

    //Cross Section Bins
    const int nbins_Q2=8;
    double Q2_low[] = {1.8, 3.25, 6.0,  9.0, 18, 55, 110, 175};
    double Q2_hi[] =  {2.2, 3.75, 7.0, 11.0, 26, 65, 130, 225};
    double Q2_center[nbins_Q2];
    for(int j=0;j<nbins_Q2;j++){
        Q2_center[j] = (Q2_hi[j]+Q2_low[j])/2.;
    }

    //Read in HERA data and Theory for RCS 
    TGraphErrors *gr_hera[nbins_Q2];
    TGraph *CT18ptxg[nbins_Q2];
    TGraph *NNPDF31_lo_as_0118_SF[nbins_Q2];
    TGraph *NNPDF31_nnlo_pch_as_0118_SF[nbins_Q2];
    TGraph *JAM4EIC[nbins_Q2];

    for(int i=0;i<nbins_Q2;i++){
        gr_hera[i] = Data_Graph("input_files/hera_data.txt",319,Q2_center[i],kBlue,22,2);
        CT18ptxg[i] = Theory_Graph("input_files/CT18ptxg_rcs.txt",Q2_center[i],kBlack,20,2);
        NNPDF31_lo_as_0118_SF[i] = Theory_Graph("input_files/NNPDF31_lo_as_0118_SF_rcs.txt",Q2_center[i],kOrange,20,2);
        NNPDF31_nnlo_pch_as_0118_SF[i] = Theory_Graph("input_files/NNPDF31_nnlo_pch_as_0118_SF_rcs.txt",Q2_center[i],kRed,20,2);
        JAM4EIC[i] = Theory_Graph("input_files/JAM4EIC_p_rcs.txt",Q2_center[i],kGreen,20,2);
    }

    //Latex Labels
    TLatex *tex1[nbins_Q2];
    for(int j=0;j<nbins_Q2;j++){
        tex1[j] = new TLatex(1E-2,1.8,Form("Q^{2} = %.1f GeV^{2}",Q2_center[j]));
        tex1[j]->SetTextFont(42);
        tex1[j]->SetTextSize(0.06);
    }

    TLatex *tex2 = new TLatex(4E-3,1.2,"HERA e^{+}p, #sqrt{s}=319 GeV");
    tex2->SetTextColor(kBlue);tex2->SetTextFont(62);

    TLatex *tex3 = new TLatex(2E-3,1.35,"CT18ptxg");
    tex3->SetTextColor(kBlack);tex3->SetTextFont(62);

    TLatex *tex4 = new TLatex(1E-3,1.5,"NNPDF31_lo_as_0118_SF");
    tex4->SetTextColor(kOrange);tex4->SetTextFont(62);

    TLatex *tex5 = new TLatex(4E-4,1.65,"NNPDF31_nnlo_pch_as_0118_SF");
    tex5->SetTextColor(kRed);tex5->SetTextFont(62);

    TLatex *tex6 = new TLatex(2E-2,1.35,"JAM4EIC");
    tex6->SetTextColor(kGreen);tex6->SetTextFont(62);

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
    CT18ptxg[0]->Draw("C Same");JAM4EIC[0]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[0]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[0]->Draw("C Same");
    gr_hera[0]->Draw("P Same");
    tex1[0]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();tex5->Draw();tex6->Draw();

    c1->cd(2); 
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    CT18ptxg[1]->Draw("C Same");JAM4EIC[1]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[1]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[1]->Draw("C Same");
    gr_hera[1]->Draw("P Same");
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
    CT18ptxg[2]->Draw("C Same");JAM4EIC[2]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[2]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[2]->Draw("C Same");
    gr_hera[2]->Draw("P Same");
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
    CT18ptxg[3]->Draw("C Same");JAM4EIC[3]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[3]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[3]->Draw("C Same");
    gr_hera[3]->Draw("P Same");
    tex1[3]->Draw();

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
    CT18ptxg[4]->Draw("C Same");JAM4EIC[4]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[4]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[4]->Draw("C Same");
    gr_hera[4]->Draw("P Same");
    tex1[4]->Draw();
    tex2->Draw();tex3->Draw();tex4->Draw();tex5->Draw();tex6->Draw();

    c2->cd(2);
    gPad->SetLogx();
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(small);
    gPad->SetTopMargin(0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    hframe->Draw();
    CT18ptxg[5]->Draw("C Same");JAM4EIC[5]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[5]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[5]->Draw("C Same");
    gr_hera[5]->Draw("P Same");
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
    CT18ptxg[6]->Draw("C Same");JAM4EIC[6]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[6]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[6]->Draw("C Same");
    gr_hera[6]->Draw("P Same");
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
    CT18ptxg[7]->Draw("C Same");JAM4EIC[7]->Draw("C Same");
    NNPDF31_lo_as_0118_SF[7]->Draw("C Same");NNPDF31_nnlo_pch_as_0118_SF[7]->Draw("C Same");
    gr_hera[7]->Draw("P Same");
    tex1[7]->Draw();

    //Print to File
    c1->Print("plots/hera_theory.pdf[");
    c1->Print("plots/hera_theory.pdf");
    c2->Print("plots/hera_theory.pdf");
    c2->Print("plots/hera_theory.pdf]");

}
