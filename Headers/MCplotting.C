#include <stdlib.h>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TProfile.h"
#include "THStack.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"


void makeProbGraph( )
{
        TCanvas *c = new TCanvas("c","A Simple Graph Example",800,600);
        float x[7] = { 0.005, 0.01, 0.02 , 0.03, 0.05, 0.07, 0.1 };
        float yOne[7] = { 0.105174, 0.352006, 0.79802, 0.965229, 0.999815, 1., 1. };
        float yNo[7]  = {0.894826, 0.647994,0.20198,0.0347715,0.000184963,0};
        float yAll[7] = {8.04699e-5,0.00187541,0.024285,0.141222,0.687043,0.956035,0.999253};

        TGraph *grOne, *grNo, *grAll;
        grOne = new TGraph (7,x,yOne);
        grNo  = new TGraph (7,x,yNo );
        grAll = new TGraph (7,x,yAll);

        grOne->SetMarkerStyle(17);
        grNo->SetMarkerStyle(17);
        grAll->SetMarkerStyle(17);

        grOne->SetLineWidth(2);
        grNo->SetLineWidth(2);
        grAll->SetLineWidth(2);

        grOne->SetLineColor(kGreen);
        grNo->SetLineColor(kBlue);
        grAll->SetLineColor(kOrange);

        grOne->SetMaximum(3.);
        grOne->SetMinimum(0.01);

        c->SetLogy();
        c->SetGridy();
        c->SetLogx();
        c->SetGridx();
        grOne->Draw("AC*");
        grNo->Draw("C*");
        grAll->Draw("C*");

        c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/14Oct2019/probGraph.png");
}



void plotHist( TString path1, TString h1Name, TString tag, bool twoD = false,  bool logy = false, bool logx = false ){

        gStyle->SetOptStat("eou");

        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        file->Close();
        if ( h1_Copy->GetEntries() == 0 ) return;

        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
        c->cd();
        c->SetGrid();
        h1_Copy->GetXaxis()->SetTitle("nPE");
        h1_Copy->GetYaxis()->SetTitle("a.u.");
        if ( logy ) c->SetLogy();
        if ( logx ) c->SetLogx();
        if ( !twoD ) h1_Copy->Draw("HIST");
        if ( twoD ) h1_Copy->Draw("colz");
        c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/211019/mCPactivity/MC_"+h1Name+tag+".png");
        c->Close();

}


void compareHists(  TString path1, TString path2, TString path3, TString path4, TString path5, TString hName, bool logy = false, bool logx = false ){


        TFile* file1 = TFile::Open(path1);
        TH1D *h1_f = (TH1D*)file1->Get(hName);
        TH1D *h1 = (TH1D*)h1_f->Clone(hName);
        h1->SetDirectory(0);
	file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2_f = (TH1D*)file2->Get(hName);
        TH1D *h2 = (TH1D*)h2_f->Clone(hName);
        h2->SetDirectory(0);
	file2->Close();

        TFile* file3 = TFile::Open(path3);
        TH1D *h3_f = (TH1D*)file3->Get(hName);
        TH1D *h3 = (TH1D*)h3_f->Clone(hName);
        h3->SetDirectory(0);
	file3->Close();

        TFile* file4 = TFile::Open(path4);
        TH1D *h4_f = (TH1D*)file4->Get(hName);
        TH1D *h4 = (TH1D*)h4_f->Clone(hName);
        h4->SetDirectory(0);
	file4->Close();

        TFile* file5 = TFile::Open(path5);
        TH1D *h5_f = (TH1D*)file5->Get(hName);
        TH1D *h5 = (TH1D*)h5_f->Clone(hName);
        h5->SetDirectory(0);
	file5->Close();


        h1->SetLineColor(kRed);
        h2->SetLineColor(kGreen);
        h3->SetLineColor(kOrange);
        h4->SetLineColor(kBlue);
        h5->SetLineColor(kBlack);

	h1->SetLineWidth(3);
	h2->SetLineWidth(3);
	h3->SetLineWidth(3);
	h4->SetLineWidth(3);
	h5->SetLineWidth(3);

//	h1->Scale(1./h1->GetEntries());
//	h2->Scale(1./h2->GetEntries());
//	h3->Scale(1./h3->GetEntries());
//	h4->Scale(1./h4->GetEntries());
//	h5->Scale(1./h5->GetEntries());

	if ( h5->GetMaximum() > h1->GetMaximum() ) h1->SetMaximum( 1.1*h5->GetMaximum() );
	if ( h4->GetMaximum() > h1->GetMaximum() ) h1->SetMaximum( 1.1*h4->GetMaximum() );
	if ( h3->GetMaximum() > h1->GetMaximum() ) h1->SetMaximum( 1.1*h3->GetMaximum() );
	if ( h2->GetMaximum() > h1->GetMaximum() ) h1->SetMaximum( 1.1*h2->GetMaximum() );

        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
	if ( logy ) c->SetLogy();
	if ( logx ) c->SetLogx();
	c->SetGridx();
	c->SetGridy();
        h1->GetYaxis()->SetLabelOffset(0.005);
        h1->GetYaxis()->SetTitle("a.u.    ");  // Define Y ..
        h1->GetXaxis()->SetLabelOffset(0.005);
	h1->SetStats(0);
	h1->GetXaxis()->SetRange(0, h1->GetNbinsX() + 1);
	h2->GetXaxis()->SetRange(0, h2->GetNbinsX() + 1);
	h3->GetXaxis()->SetRange(0, h3->GetNbinsX() + 1);
	h4->GetXaxis()->SetRange(0, h4->GetNbinsX() + 1);
	h5->GetXaxis()->SetRange(0, h5->GetNbinsX() + 1);

        h1->Draw("HIST");
        h2->Draw("HIST,SAME");
        h3->Draw("HIST,SAME");
        h4->Draw("HIST,SAME");
        h5->Draw("HIST,SAME");
	TLegend *legend = new TLegend(.50, .55, .85, .90);
        legend->AddEntry(h1,"m=1.0GeV, Q=0.01" );
        legend->AddEntry(h2,"m=1.0GeV, Q=0.03" );
        legend->AddEntry(h3,"m=1.0GeV, Q=0.05" );
        legend->AddEntry(h4,"m=1.0GeV, Q=0.07" );
        legend->AddEntry(h5,"m=1.0GeV, Q=0.1" );
        legend->Draw();
        c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/211019/mCPactivity/MC_"+hName+".png");
        c->Close();

}



void compareHists_v2(  TString path1, TString path2, TString path3, TString path4, TString path5, TString hName, bool logy = false ){


        TFile* file1 = TFile::Open(path1);
        TH1D *h1_f = (TH1D*)file1->Get(hName);
        TH1D *h1 = (TH1D*)h1_f->Clone(hName);
        h1->SetDirectory(0);
	file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2_f = (TH1D*)file2->Get(hName);
        TH1D *h2 = (TH1D*)h2_f->Clone(hName);
        h2->SetDirectory(0);
	file2->Close();

        TFile* file3 = TFile::Open(path3);
        TH1D *h3_f = (TH1D*)file3->Get(hName);
        TH1D *h3 = (TH1D*)h3_f->Clone(hName);
        h3->SetDirectory(0);
	file3->Close();

        TFile* file4 = TFile::Open(path4);
        TH1D *h4_f = (TH1D*)file4->Get(hName);
        TH1D *h4 = (TH1D*)h4_f->Clone(hName);
        h4->SetDirectory(0);
	file4->Close();

        TFile* file5 = TFile::Open(path5);
        TH1D *h5_f = (TH1D*)file5->Get(hName);
        TH1D *h5 = (TH1D*)h5_f->Clone(hName);
        h5->SetDirectory(0);
	file5->Close();


        h1->SetLineColor(kRed);
        h2->SetLineColor(kGreen);
        h3->SetLineColor(kOrange);
        h4->SetLineColor(kBlue);
        h5->SetLineColor(kBlack);

	h1->SetLineWidth(3);
	h2->SetLineWidth(3);
	h3->SetLineWidth(3);
	h4->SetLineWidth(3);
	h5->SetLineWidth(3);

	h1->Scale(1./h1->GetEntries());
	h2->Scale(h1->GetMaximum()/h2->GetMaximum());
	h3->Scale(h1->GetMaximum()/h3->GetMaximum());
	h4->Scale(h1->GetMaximum()/h4->GetMaximum());
	h5->Scale(h1->GetMaximum()/h5->GetMaximum());

        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
	if ( logy ) c->SetLogy();
	c->SetGridx();
	c->SetGridy();
        h1->GetYaxis()->SetLabelOffset(0.005);
        h1->GetYaxis()->SetTitle("a.u.    ");  // Define Y ..
        h1->GetXaxis()->SetLabelOffset(0.005);
	h1->SetStats(0);
	h1->GetXaxis()->SetRange(0, h1->GetNbinsX() + 1);
	h2->GetXaxis()->SetRange(0, h2->GetNbinsX() + 1);
	h3->GetXaxis()->SetRange(0, h3->GetNbinsX() + 1);
	h4->GetXaxis()->SetRange(0, h4->GetNbinsX() + 1);
	h5->GetXaxis()->SetRange(0, h5->GetNbinsX() + 1);

        h1->Draw("HIST");
        h2->Draw("HIST,SAME");
        h3->Draw("HIST,SAME");
        h4->Draw("HIST,SAME");
        h5->Draw("HIST,SAME");
	TLegend *legend = new TLegend(.75, .60, .99, .99);
        legend->AddEntry(h1,"m=1.0GeV, Q=0.01" );
        legend->AddEntry(h2,"m=1.0GeV, Q=0.03" );
        legend->AddEntry(h3,"m=1.0GeV, Q=0.05" );
        legend->AddEntry(h4,"m=1.0GeV, Q=0.07" );
        legend->AddEntry(h5,"m=1.0GeV, Q=0.1" );
        legend->Draw();
        c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/211019/mCPactivity/MC_"+hName+".png");
        c->Close();

}


void ratioHists(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString path2="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1" , TString outPath = "/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/SimMCratio_"){


        TFile* file1 = TFile::Open(path1);
        TH1D *h1_Copy = (TH1D*)file1->Get(h1Name);
        TH1D *h1 = (TH1D*)h1_Copy->Clone(h1Name);
        h1->SetDirectory(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH1D *h2_Copy = (TH1D*)file2->Get(h1Name);
        TH1D *h2 = (TH1D*)h2_Copy->Clone(h1Name);
        h2->SetDirectory(0);
        file2->Close();

        if ( h1->GetEntries() != 0 && h2->GetEntries() != 0 ){

		TCanvas *c = new TCanvas("c", "canvas", 800, 800);
		gStyle->SetOptStat("ou");

		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
		pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
		pad1->SetLogy();         // Vertical grid
		pad1->SetGridx();         // Vertical grid
		pad1->SetGridy();         // Vertical grid
		pad1->Draw();             // Draw the upper pad: pad1
		pad1->cd();               // pad1 becomes the current pad
		h1->SetStats(0);          // No statistics on upper plot
		h1->SetMinimum(0.11);
		if ( h1->GetMaximum() < h2->GetMaximum() ) h1->SetMaximum(2.1*h2->GetMaximum());
		h1->Draw("HIST");               // Draw h1
		h2->Draw("HIST,SAME");         // Draw h2 on top of h1
		h1->GetXaxis()->SetTitle("nPE");
	        TLegend *legend = new TLegend(.75, .75, .95, .95);
                legend->AddEntry(h1, "Data");
                legend->AddEntry(h2, "MC");
                legend->Draw();

		c->cd();          // Go back to the main canvas before defining pad2
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
		pad2->SetTopMargin(0.1);
		pad2->SetBottomMargin(0.2);
		pad2->SetGridy(); // vertical grid
		pad2->Draw();
		pad2->cd();       // pad2 becomes the current pad

		TH1D *h3 = (TH1D*)h1->Clone("h3");
		h3->SetLineColor(kBlack);
		h3->SetMinimum(0.);  // Define Y ..
		h3->SetMaximum(2.); // .. range
		h3->Sumw2();
		h3->SetStats(0);      // No statistics on lower plot
		h3->Divide(h2);
		h3->Draw("ep");       // Draw the ratio plot

		h1->SetLineColor(kBlue+1);
		h1->SetLineWidth(3);


		h2->SetLineColor(kRed);
		h2->SetLineWidth(3);

		h3->SetTitle(""); // Remove the ratio title


		h3->GetYaxis()->SetTitle("Data/MC   ");
		h3->GetYaxis()->SetNdivisions(505);
		h3->GetYaxis()->SetTitleSize(20);
		h3->GetYaxis()->SetTitleFont(43);
		h3->GetYaxis()->SetTitleOffset(1.55);
		h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		h3->GetYaxis()->SetLabelSize(15);

		h3->GetXaxis()->SetTitle("nPE");
		h3->GetXaxis()->SetTitleSize(20);
		h3->GetXaxis()->SetTitleFont(43);
		h3->GetXaxis()->SetTitleOffset(4.);
		h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		h3->GetXaxis()->SetLabelSize(15);

		c->SaveAs(outPath+h1Name+".png");

		c->Close();
        }

}
