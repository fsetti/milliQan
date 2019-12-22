// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/CmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename ) 
{

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

  TH1D *h_nPE[nChannels];
  for (unsigned int j=0; j<nChannels; j++){
		h_nPE[j] = new TH1D("nPE_Ch"+ChStr[j],"Panel Calibration nPE Ch " + ChStr[j], 100, 0, 1000 );
   }


  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);

//	skip non-synchronised events 
    bool synchBoards = false;
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true; 
    if ( !synchBoards ) continue;

/////////////////////////////////////////////////////////////////////////////
//	Muon selection I
//
    vector<vector<float>> cosmicsL1, cosmicsL2, cosmicsL3;
    cosmicsL1 = MuonPulses (  chan,height,ptime,nPE,time_module_calibrated, trigPanL1 );
    cosmicsL2 = MuonPulses (  chan,height,ptime,nPE,time_module_calibrated, trigPanL2 );
    cosmicsL3 = MuonPulses (  chan,height,ptime,nPE,time_module_calibrated, trigPanL3 );

    if ( MuonInBars( chan,ptime,nPE, {24,25} ) ){
	for (unsigned int i=0; i<cosmicsL1.size(); i++){
		h_nPE[(int)cosmicsL1.at(i).at(0)]->Fill( cosmicsL1.at(i).at(3) );
	}
    }
    if ( MuonInBars( chan,ptime,nPE, {16,17} ) ){
	for (unsigned int i=0; i<cosmicsL2.size(); i++){
		h_nPE[(int)cosmicsL2.at(i).at(0)]->Fill( cosmicsL2.at(i).at(3) );
	}
    }
    if ( MuonInBars( chan,ptime,nPE, {22,23} ) ){
	for (unsigned int i=0; i<cosmicsL3.size(); i++){
		h_nPE[(int)cosmicsL3.at(i).at(0)]->Fill( cosmicsL3.at(i).at(3) );
	}
    }
	    
//	end of main loop
  }

  for (unsigned int j=0; j<nChannels; j++){
	if ( h_nPE[j]->GetEntries() == 0 ) continue;
	h_nPE[j]->Write();
   }
  outFile->Write();
  outFile->Close();

}
