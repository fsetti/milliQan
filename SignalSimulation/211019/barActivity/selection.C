// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.

#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/MCntupleHandler.h"


void selection( TChain *chain, TString output_filename ) 
{

	float NmCP=0;

	TString ChStr[32];
	ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
	ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
	ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
	ChStr[30] = "30"; ChStr[31] = "31";
	
	
	TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
	outFile->cd();
	
	TH1D *h_nPE,*h_nPEn,*h_nPEnN, *h_nPEp ;

	float increment = (maxBin-minBin)/nBins;
	Double_t edges[nBins+1];
	float x=minBin;
	for ( unsigned int i=0; i<nBins+1; i++){
	      edges[i]= pow(10,x);
	      x+=increment;
	}
	
	double nPEbins[23] = { 0., 10., 20., 30., 50., 70., 90., 120., 150., 180., 210., 240., 270., 300., 350., 400., 450., 500., 550., 600., 650., 700.,750. };
	
	h_nPE = new TH1D("nPE_slabs","nPE",         nBins, edges);
	h_nPEn = new TH1D("nPE_adjBars","nPE",      nBins, edges);
	h_nPEnN = new TH1D("nPE_non_adjBars","nPE", nBins, edges);
	h_nPEp = new TH1D("nPE_Panels","nPE",  25, 0, 250 );


	InitializeChain(chain);
	//Number of events to loop over
	Int_t nentries = (Int_t)chain->GetEntries();
	
	
///////////////////////////////////////////////////////////////////////////
// Main Event Loop 
	for(int ia = 0; ia<nentries; ia++){
	
		chain->GetEntry(ia);
		

		if ( !mcTruth_threeBarLine || !mcTruth_fourSlab ) continue;

		vector<int> trigNonNch, trigNch;
		if ( chan_muDist[0] > 0 && chan_muDist[6] > 0 && chan_muDist[2] > 0  )    { trigNch = trigN062;  trigNonNch = trignonNtop;	}	
		if ( chan_muDist[1] > 0 && chan_muDist[7] > 0 && chan_muDist[3] > 0  )    { trigNch = trigN173;  trigNonNch = trignonNtop;	}
		if ( chan_muDist[24] > 0 && chan_muDist[16] > 0 && chan_muDist[22] > 0  ) { trigNch = trigN241622;	}	
		if ( chan_muDist[25] > 0 && chan_muDist[17] > 0 && chan_muDist[23] > 0  ) { trigNch = trigN251723;	}	
		if ( chan_muDist[8] > 0 && chan_muDist[12] > 0 && chan_muDist[4] > 0  )   { trigNch = trigN8124; trigNonNch = trignonNbot;	}
		if ( chan_muDist[9] > 0 && chan_muDist[13] > 0 && chan_muDist[5] > 0  )   { trigNch = trigN9135; trigNonNch = trignonNbot;	}

		NmCP++;	
		for (unsigned int i=0; i<chan->size(); i++){
//		don't fill hists if mCP went through bar
			if ( chan_muDist[chan->at(i)] > 0 ) continue;
			for (unsigned int j=0; j<trigNch.size(); j++){
				if ( chan->at(i) == trigNch.at(j) ){
					h_nPEn->Fill( nPE->at(i) );
//	if ( nPE->at(i) < 3 ) cout << "file: " << fileID << ", even: " << orig_evt << ", scale " << scale1fb << ", true NEP: " << chan_trueNPE[chan->at(i)] << ", reco nPE: " << nPE->at(i) << ", muDist:  " << chan_muDist[chan->at(i)] <<  " channel: " << chan->at(i) <<  endl;
					break;
				}
			}
			for (unsigned int j=0; j<trigNonNch.size(); j++){
				if ( chan->at(i) == trigNonNch.at(j) ){
					h_nPEnN->Fill( nPE->at(i) );
					break;
				}
			}
			for (unsigned int j=0; j<trigPanels.size(); j++){
				if ( chan->at(i) == trigPanels.at(j) ){
					h_nPEp->Fill( nPE->at(i) );
					break;
				}
			}
		}
//	end of main loop
	}
				


 h_nPEn->Scale(1./NmCP);
 h_nPEnN->Scale(1./NmCP);
 h_nPEp->Scale(1./NmCP);

 outFile->cd();

 h_nPEn->Write();
 h_nPEnN->Write();
 h_nPEp->Write();

 outFile->Write();
 outFile->Close();

}
