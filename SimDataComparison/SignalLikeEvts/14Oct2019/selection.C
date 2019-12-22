// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename, TString Data = true ) 
{

	TString ChStr[32];
	ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
	ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
	ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
	ChStr[30] = "30"; ChStr[31] = "31";
	
	
	TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
	outFile->cd();
	
	TH1D *h_nPE,*h_nPEt,*h_nPEb,*h_nPEbT;
	TH1D *h_nPEsec, *h_dT0, *h_dT1, *h_dT2, *h_dTno, *h_dT, *h_dT18, *h_dT21;
	TH2D *h_nPEdT18, *h_nPEdT21;

	float increment = (maxBin-minBin)/nBins;
	Double_t edges[nBins+1];
	float x=minBin;
	for ( unsigned int i=0; i<nBins+1; i++){
	      edges[i]= pow(10,x);
	      x+=increment;
	}
	
	double nPEbins[23] = { 0., 10., 20., 30., 50., 70., 90., 120., 150., 180., 210., 240., 270., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750. };
	
	h_nPE = new TH1D("MC_nPE_slabs","nPE",  30, 0,30 );
	h_nPEt = new TH1D("trueMC_nPE_slabs","nPE",30, 0,30 );
	h_nPEb = new TH1D("MC_nPE_bars","nPE",  125, 0,2500 );
	h_nPEbT = new TH1D("trueMC_nPE_bars","nPE",125, 0,2500 );
	
	h_dT0 = new TH1D("dT_slabs0","#Delta_{t}",  40, -30, 50 );
	h_dT1 = new TH1D("dT_slabs1","#Delta_{t}",  40, -30, 50 );
	h_dT2 = new TH1D("dT_slabs2","#Delta_{t}",  40, -30, 50 );
	h_dTno = new TH1D("dT_slabsNoBar","#Delta_{t}",  40, -30, 50 );
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
	slabPulses = SlabSelection( chan,height,time,nPE,time_module_calibrated );
	if ( slabPulses.size() != 4 ) continue;
	std::sort( slabPulses.begin(), slabPulses.end(), sort_wrt_channel );
	if ( fabs(slabPulses.at(2).at(4) - slabPulses.at(0).at(4)) > 12 ) continue;
				
//	bar selection
	vector<vector<vector<float>>> barPulses0, barPulses1, barPulses2, noBarPulses;
	barPulses0 = pulseBarSelection_v2 ( chan,height,time,nPE,time_module_calibrated, {{1},{7},{3}} );
	barPulses1 = pulseBarSelection_v2 ( chan,height,time,nPE,time_module_calibrated, {{0,24,25},{6,16,17},{2,22,23}} );
	barPulses2 = pulseBarSelection_v2 ( chan,height,time,nPE,time_module_calibrated, {{8,9},{12,13},{4,5}} );
	noBarPulses= pulseBarSelection_v2 ( chan,height,time,nPE,time_module_calibrated, {trigBars} );
	
	if ( barPulses0.at(0).size() != 0 ) h_dT0->Fill ( slabPulse.at(0).at(4) ); 
	if ( barPulses1.at(0).size() != 0 ) h_dT1->Fill ( slabPulse.at(0).at(4) ); 
	if ( barPulses2.at(0).size() != 0 ) h_dT2->Fill ( slabPulse.at(0).at(4) ); 
	if ( noBarPulses.at(0).size()!= 0 ) h_dTno->Fill ( slabPulse.at(0).at(4) ); 
		

//	end of main loop
	}
				

 outFile->cd();

 h_dT0->Write();
 h_dT2->Write();
 h_dT3->Write();
 h_dTno->Write();

 outFile->Write();
 outFile->Close();

}
