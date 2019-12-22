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
  float increment = (maxBin-minBin)/nBins;
  Double_t edges[nBins+1];
  float x=minBin;
  for ( unsigned int i=0; i<nBins+1; i++){
	edges[i]= pow(10,x);
	x+=increment;
  }

  TH2D *h_nPEdTall, *h_nPEdTR878, *h_nPEdTR7725, *h_nPEdTET, *h_nPEdT[nChannels];
  h_nPEdTall = new TH2D("SimData_nPEdTall","nPE vs dT", 100, -50, 50, nBins, edges );
  h_nPEdTR878 = new TH2D("SimData_nPEdTR878","nPE vs dT", 100, -50, 50, nBins, edges );
  h_nPEdTR7725 = new TH2D("SimData_nPEdTR7725","nPE vs dT", 100, -50, 50, nBins, edges );
  h_nPEdTET = new TH2D("SimData_nPEdTET","nPE vs dT", 100, -50, 50, nBins, edges );
  for (unsigned int i=0; i<nChannels; i++){
	h_nPEdT[i] = new TH2D("SimData_nPEdT"+ChStr[i],"nPE vs dT", 100, -50, 50, nBins, edges );
  }

  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);


    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > 12 ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

   

/////////////////////////////////////////////////////////////////////////////
//	Muon Selection - hits in straight line 
//    bool muonHit = false;
//    muonHit = MuonInBars ( chan,ptime,nPE, trigBars );
//    if ( muonHit ) continue;
  
 
    vector<vector<float>> L1Pulses, L2Pulses, L3Pulses, pulses;
    L1Pulses = selectPulses_v2( chan,height,ptime,nPE,time_module_calibrated, trigL1 );
    L2Pulses = selectPulses_v2( chan,height,ptime,nPE,time_module_calibrated, trigL2 );
    L3Pulses = selectPulses_v2( chan,height,ptime,nPE,time_module_calibrated, trigL3 );


    for ( unsigned int i=0; i<L1Pulses.size(); i++){
	h_nPEdTall->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
	h_nPEdT[(int)L1Pulses.at(i).at(0)]->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
	if ( (int)L1Pulses.at(i).at(0) == 5 || (int)L1Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
	if ( (int)L1Pulses.at(i).at(0) == 24 || (int)L1Pulses.at(i).at(0) == 25 || (int)L1Pulses.at(i).at(0) == 9 || (int)L1Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
	if ( (int)L1Pulses.at(i).at(0) != 5 && (int)L1Pulses.at(i).at(0) != 22 && (int)L1Pulses.at(i).at(0) != 24 && (int)L1Pulses.at(i).at(0) != 25 &&(int)L1Pulses.at(i).at(0) != 9 && (int)L1Pulses.at(i).at(0) != 17 )  h_nPEdTR878->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
    }
    for ( unsigned int i=0; i<L2Pulses.size(); i++){
	h_nPEdTall->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
	h_nPEdT[(int)L2Pulses.at(i).at(0)]->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
	if ( (int)L2Pulses.at(i).at(0) == 5 || (int)L2Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
	if ( (int)L2Pulses.at(i).at(0) == 24 || (int)L2Pulses.at(i).at(0) == 25 || (int)L2Pulses.at(i).at(0) == 9 || (int)L2Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
	if ( (int)L2Pulses.at(i).at(0) != 5 && (int)L2Pulses.at(i).at(0) != 22 && (int)L2Pulses.at(i).at(0) != 24 && (int)L2Pulses.at(i).at(0) != 25 &&(int)L2Pulses.at(i).at(0) != 9 && (int)L2Pulses.at(i).at(0) != 17 ) h_nPEdTR878->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
    }
    for ( unsigned int i=0; i<L3Pulses.size(); i++){
	h_nPEdTall->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
	h_nPEdT[(int)L3Pulses.at(i).at(0)]->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
	if ( (int)L3Pulses.at(i).at(0) == 5 || (int)L3Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
	if ( (int)L3Pulses.at(i).at(0) == 24 || (int)L3Pulses.at(i).at(0) == 25 || (int)L3Pulses.at(i).at(0) == 9 || (int)L3Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
	if ( (int)L3Pulses.at(i).at(0) != 5 && (int)L3Pulses.at(i).at(0) != 22 && (int)L3Pulses.at(i).at(0) != 24 && (int)L3Pulses.at(i).at(0) != 25 &&(int)L3Pulses.at(i).at(0) != 9 && (int)L3Pulses.at(i).at(0) != 17 ) h_nPEdTR878->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
    }

//	end of main loop
  }



  outFile->cd();

  for (unsigned int i=0; i<nChannels;i++){
	if ( h_nPEdT[i]->GetEntries() != 0 ) h_nPEdT[i]->Write();
  }

  h_nPEdTall->Write();
  h_nPEdTR878->Write();
  h_nPEdTR7725->Write();
  h_nPEdTET->Write();
  outFile->Write();
  outFile->Close();

}
