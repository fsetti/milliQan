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


void comparePlots(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", TString h2Name="name_of_h2", TString h3Name="name_of_h3", TString h4Name="name_of_h4"){

        TFile* file = TFile::Open(path1);

        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        h1_Copy->SetStats(0);

        TH1D *h2 = (TH1D*)file->Get(h2Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h2Name);
        h2_Copy->SetDirectory(0);
        h2_Copy->SetStats(0);

        TH1D *h3 = (TH1D*)file->Get(h3Name);
        TH1D *h3_Copy = (TH1D*)h3->Clone(h3Name);
        h3_Copy->SetDirectory(0);
        h3_Copy->SetStats(0);

        TH1D *h4 = (TH1D*)file->Get(h4Name);
        TH1D *h4_Copy = (TH1D*)h4->Clone(h4Name);
        h4_Copy->SetDirectory(0);
        h4_Copy->SetStats(0);

        file->Close();

	if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 &&  h3_Copy->GetEntries() != 0 && h4_Copy->GetEntries() != 0 ){
//		gStyle->SetOptStat(11);
	        h1_Copy->SetLineColor(kRed);
	        h1_Copy->SetFillColor(kRed);
	        h1_Copy->SetFillStyle(3003);
	        h2_Copy->SetLineColor(kBlue);
	        h2_Copy->SetLineWidth(3);
//	        h2_Copy->SetFillColor(kBlue);
//	        h2_Copy->SetFillStyle(3002);
	        h3_Copy->SetLineColor(kGreen);
	        h3_Copy->SetLineWidth(3);
//	        h3_Copy->SetFillColor(kGreen);
//	        h3_Copy->SetFillStyle(3002);
	        h4_Copy->SetLineColor(kBlack);
	        h4_Copy->SetFillColor(kBlack);
	        h4_Copy->SetFillStyle(3002);
	
	        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
//		c->SetLogy();
		h1_Copy->GetYaxis()->SetLabelOffset(0.005);
		h1_Copy->GetYaxis()->SetTitle("[Hz]    ");  // Define Y ..
		h1_Copy->GetXaxis()->SetLabelOffset(0.005);
//		h1_Copy->GetXaxis()->SetTitle("time [ns]");  // Define Y ..
	        h1_Copy->Draw("HIST,ERR");
	        h2_Copy->Draw("HIST,ERR,SAME");
	        h3_Copy->Draw("HIST,ERR,SAME");
	        h4_Copy->Draw("HIST,ERR,SAME");
	        TLegend *legend = new TLegend(.75, .80, .95, .95);
	        legend->AddEntry(h1_Copy, "Ch18");
	        legend->AddEntry(h2_Copy, "Ch20");
	        legend->AddEntry(h4_Copy, "Ch28");
	        legend->AddEntry(h3_Copy, "Ch21");
		legend->Draw();
	        c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Cosmics/Comparison_"+h1Name+".png");
	        c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Cosmics/Comparison_"+h1Name+".pdf");
		c->Close();
	}

}


void plotHist( TString path1, TString h1Name, bool logy = false, bool twoD = false ){
	
	gStyle->SetOptStat("eou");

        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
//        h1_Copy->SetStats(0);
	file->Close();
	if ( h1_Copy->GetEntries() == 0 ) return;
	
	TCanvas *c = new TCanvas("c", "canvas", 900, 700);
	c->cd();
	if ( logy ) c->SetLogy();
	if ( !twoD ) h1_Copy->Draw("HIST");
	if ( twoD ) h1_Copy->Draw("colz");
	c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Cosmics/"+h1Name+".png");
	c->Close();

}

void compareHists(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", TString h2Name="name_of_h2" ){

	gStyle->SetOptStat("ou");


        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);

        TH1D *h2 = (TH1D*)file->Get(h2Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h2Name);
        h2_Copy->SetDirectory(0);
//        h2_Copy->SetStats(0);
        file->Close();

        if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 ){
                h1_Copy->SetLineColor(kRed);
                h1_Copy->SetFillColor(kRed);
                h1_Copy->SetFillStyle(3003);
                h2_Copy->SetLineColor(kBlue);
                h2_Copy->SetFillColor(kBlue);
                h2_Copy->SetFillStyle(3002);
                if ( h1_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h2_Copy->SetMaximum(h1_Copy->GetMaximum()*1.1);

                TCanvas *c = new TCanvas("c", "canvas", 900, 700);
//		c->SetLogy();
                h2_Copy->GetYaxis()->SetLabelOffset(0.005);
 //               h2_Copy->GetYaxis()->SetTitle("[Hz]    ");  // Define Y ..
                h2_Copy->GetXaxis()->SetLabelOffset(0.005);
                h2_Copy->Draw("HIST,ERR");
                h1_Copy->Draw("HIST,ERR,SAME");
                TLegend *legend = new TLegend(.1, .6, .38, .9);
                legend->AddEntry(h1_Copy, "slab");
                legend->AddEntry(h2_Copy, "same layer, #mu pulse");
                legend->Draw();
                c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Cosmics/Comparison_"+h1Name+"_vs_"+h2Name+".png");
                c->Close();
        }

}

void compareHists_v2(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString path2="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", bool logy = false ){

	gStyle->SetOptStat("uo");


        TFile* file1 = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file1->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        file1->Close();
	h1_Copy->Scale(1./h1_Copy->GetEntries());

        TFile* file2 = TFile::Open(path2);
        TH1D *h2 = (TH1D*)file2->Get(h1Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h1Name);
        h2_Copy->SetDirectory(0);
        file2->Close();
	h2_Copy->Scale(1./h2_Copy->GetEntries());

        if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 ){
                h1_Copy->SetLineColor(kRed);
                h1_Copy->SetFillColor(kRed);
                h1_Copy->SetFillStyle(3003);
                h2_Copy->SetLineColor(kBlue);
                h2_Copy->SetFillColor(kBlue);
                h2_Copy->SetFillStyle(3002);
                if ( h1_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h2_Copy->SetMaximum(h1_Copy->GetMaximum()*1.1);

                TCanvas *c = new TCanvas("c", "canvas", 900, 700);
		if ( logy ) c->SetLogy();
                h2_Copy->GetYaxis()->SetLabelOffset(0.005);
 //               h2_Copy->GetYaxis()->SetTitle("[Hz]    ");  // Define Y ..
                h2_Copy->GetXaxis()->SetLabelOffset(0.005);
                h2_Copy->Draw("HIST");
                h1_Copy->Draw("HIST,SAME");
//	        TLegend *legend = new TLegend(.75, .80, .95, .95);
                TLegend *legend = new TLegend(.1, .8, .38, .95);
                legend->AddEntry(h1_Copy, "BarHit");
                legend->AddEntry(h2_Copy, "noBarHit");
                legend->Draw();
                c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Cosmics/Comparison_"+h1Name+".png");
                c->Close();
        }

}


void compareProfiles(TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString path2="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString hName="name_of_h" ){

        TFile* file1 = TFile::Open(path1);
        TProfile *h1 = (TProfile*)file1->Get(hName);
        TProfile *h1_Copy = (TProfile*)h1->Clone(hName);
        h1_Copy->SetDirectory(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TProfile *h2 = (TProfile*)file2->Get(hName);
        TProfile *h2_Copy = (TProfile*)h2->Clone(hName);
        h2_Copy->SetDirectory(0);
        file2->Close();

        if ( h1_Copy->GetEntries() != 0 && h2_Copy->GetEntries() != 0 ){
                h1_Copy->SetLineColor(kRed);
                h2_Copy->SetLineColor(kBlue);
                if ( h1_Copy->GetMaximum() > h2_Copy->GetMaximum() ) h2_Copy->SetMaximum(h1_Copy->GetMaximum()*1.1);

                TCanvas *c = new TCanvas("c", "canvas", 900, 700);
		c->SetGrid();
                h2_Copy->Draw();
                h1_Copy->Draw("SAME");
                TLegend *legend = new TLegend(.1, .8, .38, .95);
                legend->AddEntry(h1_Copy, "BarHit");
                legend->AddEntry(h2_Copy, "noBarHit");
                legend->Draw();
                c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Cosmics/Comparison_"+hName+".png");
                c->Close();
        }

}




void stackHists (TString path1="/home/users/fsetti/milliQan/EfficienciesTrig_samples/Run2152TC_triggerEff.root", TString h1Name="name_of_h1", TString h2Name="name_of_h2" , TString h3Name="name_of_h3" ){


        THStack *hs = new THStack("hist","");

        TFile* file = TFile::Open(path1);
        TH1D *h1 = (TH1D*)file->Get(h1Name);
        TH1D *h1_Copy = (TH1D*)h1->Clone(h1Name);
        h1_Copy->SetDirectory(0);
        TH1D *h2 = (TH1D*)file->Get(h2Name);
        TH1D *h2_Copy = (TH1D*)h2->Clone(h2Name);
        h2_Copy->SetDirectory(0);
        h2_Copy->SetStats(0);
        TH1D *h3 = (TH1D*)file->Get(h3Name);
        TH1D *h3_Copy = (TH1D*)h3->Clone(h3Name);
        h3_Copy->SetDirectory(0);
        h3_Copy->SetStats(0);
        file->Close();

        h1_Copy->SetLineColor(kRed);
        h1_Copy->SetFillColor(kRed);
        h1_Copy->SetFillStyle(3003);
        h2_Copy->SetLineColor(kBlue);
        h2_Copy->SetFillColor(kBlue);
        h2_Copy->SetFillStyle(3002);
        h3_Copy->SetLineColor(kGreen);
        h3_Copy->SetFillColor(kGreen);
        h3_Copy->SetFillStyle(3001);

	hs->Add(h1_Copy);
	hs->Add(h2_Copy);
	hs->Add(h3_Copy);

        TCanvas *c = new TCanvas("c", "canvas", 900, 700);
//	c->SetLogy();
//        hs->GetYaxis()->SetLabelOffset(0.005);
//        hs->GetYaxis()->SetTitle("[Hz]    ");  // Define Y ..
//        hs->GetXaxis()->SetLabelOffset(0.005);
        hs->Draw("HIST,ERR");
        TLegend *legend = new TLegend(.75, .80, .95, .95);
        legend->AddEntry(h1_Copy, "Layer1");
        legend->AddEntry(h2_Copy, "Layer2");
        legend->AddEntry(h3_Copy, "Layer3");
        legend->Draw();
        c->SaveAs("/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/Cosmics/Comparison_Stack_"+h1Name+".png");
        c->Close();


}
