//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

//Read in data file and put into graph (TGraph)
//Assumes file is in format grid number,Q2,x,y
TGraph* Theory_Graph(const char* infile, int grid_set=0, double Q2_set=1, 
                        int color=1, int style=20, int lw=2, bool abs=0){

  //Create graph
  TGraph *gr = new TGraph();
  gr->SetMarkerStyle(style);gr->SetMarkerColor(color);
  gr->SetLineWidth(lw);gr->SetLineColor(color);

  //Read input file
  ifstream in;
  in.open(infile,ios::in);

  //Write to Graph
  int grid(0);
  double Q2(0),x(0),y(0);
  int counter(0);

  while(in.good()){
    in >> grid >> Q2 >> x >> y;

    if( in.good() && grid==grid_set && fabs(Q2-Q2_set)<1E-2 ){

      if(abs)
        gr->SetPoint(counter,x,fabs(y));
      else
        gr->SetPoint(counter,x,y);

      counter++;
    }
  }

  cout<<"For Grid #"<<grid_set<<" and Q2 = "<<Q2_set<<" GeV^2... ";
  cout<<"We read "<<counter<<" matching lines from data file "<<infile<<" !"<<endl;
  in.close();
  return gr;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

void sf_theory(){
    

    //Read Theory Files
    TGraph *CT18ptxg_F2[4],*CT18ptxg_FL[4],*CT18ptxg_F3[4];
    TGraph *NNPDF31_lo_as_0118_SF_F2[4],*NNPDF31_lo_as_0118_SF_FL[4],*NNPDF31_lo_as_0118_SF_F3[4];
    TGraph *NNPDF31_nnlo_pch_as_0118_SF_F2[4],*NNPDF31_nnlo_pch_as_0118_SF_FL[4],*NNPDF31_nnlo_pch_as_0118_SF_F3[4];
    TGraph *JAM4EIC_F2[4],*JAM4EIC_FL[4],*JAM4EIC_F3[4];

    double Q2[] = {5,10,100,1000};

    cout<<"-------------------------------"<<endl;
    for(int i=0;i<4;i++){
        CT18ptxg_F2[i] = Theory_Graph("input_files/F2_Theory.txt",0,Q2[i],kBlack,20,3);
        CT18ptxg_FL[i] = Theory_Graph("input_files/FL_Theory.txt",0,Q2[i],kBlack,20,3);
        CT18ptxg_F3[i] = Theory_Graph("input_files/F3_Theory.txt",0,Q2[i],kBlack,20,3);
        NNPDF31_lo_as_0118_SF_F2[i] = Theory_Graph("input_files/F2_Theory.txt",1,Q2[i],kBlue,20,3);
        NNPDF31_lo_as_0118_SF_FL[i] = Theory_Graph("input_files/FL_Theory.txt",1,Q2[i],kBlue,20,3);
        NNPDF31_lo_as_0118_SF_F3[i] = Theory_Graph("input_files/F3_Theory.txt",1,Q2[i],kBlue,20,3);
        NNPDF31_nnlo_pch_as_0118_SF_F2[i] = Theory_Graph("input_files/F2_Theory.txt",2,Q2[i],kRed,20,3);
        NNPDF31_nnlo_pch_as_0118_SF_FL[i] = Theory_Graph("input_files/FL_Theory.txt",2,Q2[i],kRed,20,3);
        NNPDF31_nnlo_pch_as_0118_SF_F3[i] = Theory_Graph("input_files/F3_Theory.txt",2,Q2[i],kRed,20,3);
        JAM4EIC_F2[i] = Theory_Graph("input_files/F2_Theory.txt",3,Q2[i],kGreen,20,3);
        JAM4EIC_FL[i] = Theory_Graph("input_files/FL_Theory.txt",3,Q2[i],kGreen,20,3);
        JAM4EIC_F3[i] = Theory_Graph("input_files/F3_Theory.txt",3,Q2[i],kGreen,20,3);
    }
    cout<<"-------------------------------"<<endl;

    //Draw plots
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetNdivisions(405,"X");
    gStyle->SetNdivisions(405,"Y");
    gStyle->SetLabelSize(0.06,"X");
    gStyle->SetLabelSize(0.06,"Y");
    gStyle->SetLabelOffset(-0.02,"X");
    gStyle->SetLabelOffset(0.02,"Y");
    gStyle->SetTitleXSize(0.07);
    gStyle->SetTitleXOffset(0.7);
    gStyle->SetTitleYSize(0.07);
    gStyle->SetTitleYOffset(0.7);

    float small = 1E-20;

    TCanvas *c1[4];
    TH1 *hframe1[4],*hframe2[4],*hframe3[4];
    TLatex *tex1[4];
    
    for(int i=0;i<4;i++){

        c1[i] = new TCanvas(Form("c[%d]",i));
        c1[i]->Divide(3,1,small,small);

        c1[i]->cd(1);
        gPad->SetLogx();
        gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
        //gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        gPad->SetTickx();
        gPad->SetTicky(); 
        hframe1[i] = gPad->DrawFrame(1E-4,0,1,3.0);
        hframe1[i]->GetXaxis()->SetTitle("x");hframe1[i]->GetXaxis()->CenterTitle();
        hframe1[i]->GetYaxis()->SetTitle("F_{2,NC}");hframe1[i]->GetYaxis()->CenterTitle();
        CT18ptxg_F2[i]->Draw("C same");
        NNPDF31_lo_as_0118_SF_F2[i]->Draw("C same");
        NNPDF31_nnlo_pch_as_0118_SF_F2[i]->Draw("C same");
        JAM4EIC_F2[i]->Draw("C same");

        tex1[i] = new TLatex(2E-3,2.75,Form("Q^{2} = %.0f GeV^{2}",Q2[i]));
        tex1[i]->SetTextColor(kBlack);tex1[i]->SetTextSize(0.06);//tex1[i]->SetTextFont(42);
        tex1[i]->Draw();

        gPad->Modified();gPad->Update();

        c1[i]->cd(2); 
        gPad->SetLogx();
        gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
        //gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        gPad->SetTickx();
        gPad->SetTicky(); 
        hframe2[i] = gPad->DrawFrame(1E-4,0,1,0.5);
        hframe2[i]->GetXaxis()->SetTitle("x");hframe2[i]->GetXaxis()->CenterTitle();
        hframe2[i]->GetYaxis()->SetTitle("F_{L,NC}");hframe2[i]->GetYaxis()->CenterTitle();
        CT18ptxg_FL[i]->Draw("C same");
        NNPDF31_lo_as_0118_SF_FL[i]->Draw("C same");
        NNPDF31_nnlo_pch_as_0118_SF_FL[i]->Draw("C same");
        JAM4EIC_FL[i]->Draw("C same");
        gPad->Modified();gPad->Update();
    
        c1[i]->cd(3);
        gPad->SetLogx();
        gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
        //gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        gPad->SetTickx();
        gPad->SetTicky(); 
        hframe3[i] = gPad->DrawFrame(1E-4,0,1,1.0);
        hframe3[i]->GetXaxis()->SetTitle("x");hframe3[i]->GetXaxis()->CenterTitle();
        hframe3[i]->GetYaxis()->SetTitle("F_{3,NC}");hframe3[i]->GetYaxis()->CenterTitle();
        CT18ptxg_F3[i]->Draw("C same");
        NNPDF31_lo_as_0118_SF_F3[i]->Draw("C same");
        NNPDF31_nnlo_pch_as_0118_SF_F3[i]->Draw("C same");
        JAM4EIC_F3[i]->Draw("C same");

        if(i==0){
            TLegend *leg1 = new TLegend(0.15,0.65,0.8,0.875);
            leg1->SetBorderSize(0);leg1->SetTextSize(0.035);
            leg1->AddEntry(CT18ptxg_F3[i],"CT18ptxg","l");
            leg1->AddEntry(NNPDF31_lo_as_0118_SF_F3[i],"NNPDF31_lo_as_0118_SF","l");
            leg1->AddEntry(NNPDF31_nnlo_pch_as_0118_SF_F3[i],"NNPDF31_nnlo_pch_as_0118_SF","l");
            leg1->AddEntry(JAM4EIC_F3[i],"JAM4EIC","l");
            leg1->Draw();
        }

        gPad->Modified();gPad->Update();

        if(i==0)c1[i]->Print("sf_theory.pdf[");
        c1[i]->Print("sf_theory.pdf");
        if(i==3)c1[i]->Print("sf_theory.pdf]");
    }
}