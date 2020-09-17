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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

void replica_compare(){

    //Q2 values to plot; 8 Q2 points
    const int num_plot = 8;
    double Q2_plot[num_plot] = {2.05,5.15,8.15,12.9,20.5,51.5,129.2,514.5};
    
    //Number of Replicas (50 for JAM4EIC)
    const int num_rep = 50;

    //Read-in Central Value
    TGraph *JAM4EIC_cen[num_plot];

    for(int i=0;i<num_plot;i++){
        JAM4EIC_cen[i] = Theory_Graph("input_files/JAM4EIC_p_rcs.txt",Q2_plot[i],kGreen,20,2);
    }

    //Read-in Replicas
    TGraph *JAM4EIC_replica[num_rep][num_plot];

    for(int j=0;j<num_rep;j++){
        for(int i=0;i<num_plot;i++){
            JAM4EIC_replica[j][i] = Theory_Graph(Form("input_files/replicas/JAM4EIC_p_%d_rcs.txt",j+1),Q2_plot[i],kGreen,20,2);
        }
    }

    //Make Ratio Plots
    TGraph *ratio_JAM4EIC[num_rep][num_plot];

    for(int j=0;j<num_rep;j++){
        for(int i=0;i<num_plot;i++){
            ratio_JAM4EIC[j][i] = Ratio_Theory(JAM4EIC_cen[i],JAM4EIC_replica[j][i],kGreen,20,1);
        }
    }


    //Latex Labels
    TLatex *tex1[num_plot];
    for(int i=0;i<num_plot;i++){
        tex1[i] = new TLatex(2E-4,1.04,Form("Q^{2} = %.1f GeV^{2}",Q2_plot[i]));
        tex1[i]->SetTextFont(42);
        tex1[i]->SetTextSize(0.06);
    }

    TLatex *tex2 = new TLatex(2E-4,0.98,"18 GeV e^{-} on 275 GeV p, #sqrt{s}=141 GeV");
    tex2->SetTextColor(kBlack);tex2->SetTextFont(42);

    TLatex *tex3 = new TLatex(2E-4,0.97,"JAM4EIC replica/central");
    tex3->SetTextColor(kGreen);tex3->SetTextFont(62);


    //Make Plots
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
    int counter = 0;

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
    TH1 *hframe = gPad->DrawFrame(1E-4,0.95,1,1.05);
    hframe->GetXaxis()->SetTitle("x");hframe->GetXaxis()->CenterTitle();
    hframe->GetYaxis()->SetTitle("Cross Section Ratio");hframe->GetYaxis()->CenterTitle();
    counter = 0;

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
    tex1[counter]->Draw();

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

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
    tex1[counter]->Draw();
    
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

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
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

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
    tex1[counter]->Draw();tex2->Draw();tex3->Draw();
    
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

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
    tex1[counter]->Draw();

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

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
    tex1[counter]->Draw();
      
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

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
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

    for(int j=0;j<num_rep;j++){
        ratio_JAM4EIC[j][counter]->Draw("L Same");
    }
    tex1[counter]->Draw();tex2->Draw();tex3->Draw();
    

    //Print to file
    c1->Print("plots/replica_compare.pdf[");
    c1->Print("plots/replica_compare.pdf");
    c2->Print("plots/replica_compare.pdf");
    c2->Print("plots/replica_compare.pdf]");
    

}