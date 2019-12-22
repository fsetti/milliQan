// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
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

  TH1D *h_nPE, *h_dT;
  TH2D *h_nPEdT, *h_dPEdT, *h_rPEdT;
  float minBin = 0;
  float maxBin = 5;
  const int nBins = 100;
  float increment = (maxBin-minBin)/nBins;
  Double_t edges[nBins+1];
  float x=minBin;
  for ( unsigned int i=0; i<nBins+1; i++){
	edges[i]= pow(10,x);
	x+=increment;
  }

  h_nPEdT = new TH2D("SimData_nPEdT","nPE vs dT", 80, -20, 60, nBins, edges );
  h_dPEdT = new TH2D("SimData_dPEdT","Delta nPE vs dT", 80, -20, 60, nBins, edges );
  h_rPEdT = new TH2D("SimData_rPEdT","Ratio nPE vs dT", 80, -20, 60, 1e2,0.00001,1 );
  h_nPE = new TH1D("SimData_nPE","largest nPE", nBins, edges );
  h_dT = new TH1D("SimData_dT","largest dT", 80, -20, 60 );

  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);

//	skip non-synchronised events && require beam On
    bool synchBoards = false;			// only consider events with synched boards
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true; 
    if (!synchBoards ) continue;
    if ( !beam ) continue;



    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > 12 ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

    
    vector<vector<float>> L1Pulses, L2Pulses, L3Pulses;
    L1Pulses = selectPulses_v2( chan,height,ptime,nPE,time_module_calibrated, trigL1 );
    L2Pulses = selectPulses_v2( chan,height,ptime,nPE,time_module_calibrated, trigL2 );
    L3Pulses = selectPulses_v2( chan,height,ptime,nPE,time_module_calibrated, trigL3 );


//	end of main loop
  }



  outFile->cd();

  h_nPEdT->Write();
  h_dPEdT->Write();
  h_rPEdT->Write();
  h_dT->Write();
  h_nPE->Write();
  outFile->Write();
  outFile->Close();

}
