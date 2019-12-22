#include "pmtSettings.h"
#include <iostream>

using namespace std;


void timeDelay(TString path1, TString pmtType, float threshold, TString outputFileName ){


	double nPEarray[nBins];
	double edges[nBins+1];
	float x=minBin;
	for ( unsigned int i=0; i<nBins+1; i++){ 
	      edges[i]= pow(10,x);
		if ( i == nBins ) break;
	      nPEarray[i]= pow(10,x);
	      x+=increment;
	}  


        TFile* file = TFile::Open(path1);
        TH1D *h1_file = (TH1D*)file->Get(pmtTemplate);
        TH1D *h1 = (TH1D*)h1_file->Clone(pmtTemplate);
        h1->SetDirectory(0);
        h1->SetStats(0);
        file->Close();

	TString output_filename = outputFileName + pmtType + ".root";
	TFile *outFile = new TFile(output_filename.Data(),"RECREATE"); 
	TH1D *h = new TH1D(pmtType,"time delay as a function of nPE", nBins, edges ); 

	float spe_value;
	if ( pmtType == "ET" ) spe_value = spe_fit_mean_ET;
	if ( pmtType == "R878" ) spe_value = spe_fit_mean_R878;
	if ( pmtType == "R7725" ) spe_value = spe_fit_mean_R7725;
	h1->Scale(spe_value);
	double nPE, nPEprev;
	nPEprev = 1;

	for ( unsigned int i=0; i<nBins; i++){
		nPE = nPEarray[i];
		h1->Scale(nPE/nPEprev);
		unsigned int pulseBins = h1->GetNbinsX();
		float time = -1;
		for ( unsigned int j=1;j<pulseBins+1;j++){
			if ( h1->GetBinContent(j) > threshold ){
				time = h1->GetXaxis()->GetBinCenter( j );
				break;
			}
		}
		h->Fill( nPE, time );
		nPEprev = nPE;
	}
	
	outFile->cd(); 
	h->Write();
	outFile->Write();
	outFile->Close();


}


void plotAllPmts( TString inputFile1, TString h1Name, TString inputFile2,  TString h2Name, TString inputFile3, TString h3Name ){


        TFile* file1 = TFile::Open(inputFile1);
        TH1D *h1_file = (TH1D*)file1->Get(h1Name);
        TH1D *h1 = (TH1D*)h1_file->Clone(h1Name);
        h1->SetDirectory(0);
        h1->SetStats(0);
        file1->Close();

        TFile* file2 = TFile::Open(inputFile2);
        TH1D *h2_file = (TH1D*)file2->Get(h2Name);
        TH1D *h2 = (TH1D*)h2_file->Clone(h2Name);
        h2->SetDirectory(0);
        h2->SetStats(0);
	file2->Close();

        TFile* file3 = TFile::Open(inputFile3);
        TH1D *h3_file = (TH1D*)file3->Get(h3Name);
        TH1D *h3 = (TH1D*)h3_file->Clone(h3Name);
        h3->SetDirectory(0);
        h3->SetStats(0);
	file3->Close();


	h1->SetLineColor(kBlack);
	h1->SetLineWidth(3);
	h2->SetLineColor(kRed);
	h2->SetLineWidth(3);
	h3->SetLineColor(kBlue);
	h3->SetLineWidth(3);


	TCanvas *c = new TCanvas("canvas","",800,600);
	c->cd();
	c->SetLogx();
	h1->Draw("HIST");
	h2->Draw("HIST,SAME");
	h3->Draw("HIST,SAME");
        TLegend *legend = new TLegend(.75, .80, .95, .95);
        legend->AddEntry(h1, "R878");
        legend->AddEntry(h2, "R7725");
        legend->AddEntry(h3, "ET");
        legend->Draw();


	c->SaveAs("/home/users/fsetti/public_html/milliQan/PmtWaveforms/pmts.pdf");

	c->Close();
}	
