// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename, bool Data = true ) 
{

	TString ChStr[32];
	ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
	ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
	ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
	ChStr[30] = "30"; ChStr[31] = "31";
	
	
 	TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
	TH1D *h_nPE,*h_nPEt,*h_nPEb,*h_nPEbT;
	TH1D *h_nPEsec, *h_dT0[32], *h_dT1[32], *h_dT2[32], *h_dTno[32], *h_dT, *h_dT18, *h_dT21;
	TH2D *h_nPEdT18, *h_nPEdT21;

	float increment = (maxBin-minBin)/nBins;
	Double_t edges[nBins+1];
	float x=minBin;
	for ( unsigned int i=0; i<nBins+1; i++){
	      edges[i]= pow(10,x);
	      x+=increment;
	}
	
	h_nPE = new TH1D("MC_nPE_slabs","nPE",  30, 0,30 );
	h_nPEt = new TH1D("trueMC_nPE_slabs","nPE",30, 0,30 );
	h_nPEb = new TH1D("MC_nPE_bars","nPE",  125, 0,2500 );
	h_nPEbT = new TH1D("trueMC_nPE_bars","nPE",125, 0,2500 );

	for (unsigned int i=0; i<nChannels;i++){	
		h_dT0[i] = new TH1D("dT_slabs0_Ch"+ChStr[i],"#Delta_{t}",  40, 350, 460 );
		h_dT1[i] = new TH1D("dT_slabs1_Ch"+ChStr[i],"#Delta_{t}",  40, 350, 460 );
		h_dT2[i] = new TH1D("dT_slabs2_Ch"+ChStr[i],"#Delta_{t}",  40, 350, 460 );
		h_dTno[i]= new TH1D("dT_slabsNoBar_Ch"+ChStr[i],"#Delta_{t}",  40, 350, 460 );
	}

	h_dT18 = new TH1D("dT_slab18","#Delta_{t}",  40, 360, 480 );
	h_dT21 = new TH1D("dT_slab21","#Delta_{t}",  40, 360, 480 );
	h_nPEsec = new TH1D("nPE_secPulses","nPE",  125, 0,1500 );

	h_nPEdT18 = new TH2D("nPEdT_slab18","#Delta_{t}",30,0,30,  40, 360, 480 );
	h_nPEdT21 = new TH2D("nPEdT_slab21","#Delta_{t}",30,0,30,  40, 360, 480 );

	InitializeChain(chain);
	//Number of events to loop over
	Int_t nentries = (Int_t)chain->GetEntries();
	
	
///////////////////////////////////////////////////////////////////////////
// Main Event Loop 
	for(int ia = 0; ia<nentries; ia++){
	
		chain->GetEntry(ia);

	// 	synchronisation requirements
		bool synchBoards = false;                   // only consider events with synched boards
		if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true;
		if (!synchBoards && Data ) continue;
		if ( !beam  && Data ) continue;

	
	// slab selection
		vector<vector<float>> slabPulses;
		slabPulses = SlabSelection ( chan,height,ptime,nPE,time_module_calibrated );
		if ( slabPulses.size() != 4 ) continue;
		std::sort( slabPulses.begin(), slabPulses.end(), sort_wrt_channel );
		if ( fabs(slabPulses.at(2).at(4) - slabPulses.at(0).at(4)) > 12 ) continue;
					
	//	bar selection
		vector<vector<vector<float>>> barPulses0, barPulses1, barPulses2;
		barPulses0 = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173} );
		barPulses1 = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig241622,trig251723} );
		barPulses2 = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig8124,trig9135} );
		bool muonHit = MuonInBars ( chan,ptime,nPE, trigBars );


		for ( unsigned int i=0; i<slabPulses.size(); i++){
			if ( barPulses0.size() != 0 && barPulses1.size() == 0 && barPulses2.size() == 0 ) h_dT0[(int)slabPulses.at(i).at(0)]->Fill ( slabPulses.at(i).at(4) );
			if ( barPulses0.size() == 0 && barPulses1.size() != 0 && barPulses2.size() == 0 ) h_dT1[(int)slabPulses.at(i).at(0)]->Fill ( slabPulses.at(i).at(4) );
			if ( barPulses0.size() == 0 && barPulses1.size() == 0 && barPulses2.size() != 0 ) h_dT2[(int)slabPulses.at(i).at(0)]->Fill ( slabPulses.at(i).at(4) );
			if ( !muonHit ) 	     							 h_dTno[(int)slabPulses.at(i).at(0)]->Fill ( slabPulses.at(i).at(4) );
		}


//	end of main loop
	}
				

 outFile->cd();
 for (unsigned int i=0; i<nChannels;i++){
	if ( h_dT0[i]->GetEntries() != 0 ) h_dT0[i]->Write();
	if ( h_dT1[i]->GetEntries() != 0 ) h_dT1[i]->Write();
	if ( h_dT2[i]->GetEntries() != 0 ) h_dT2[i]->Write();
	if ( h_dTno[i]->GetEntries() != 0 ) h_dTno[i]->Write();
 }

// h_dT0->Write();
// h_dT1->Write();
// h_dT2->Write();

 outFile->Write();
 outFile->Close();

}
