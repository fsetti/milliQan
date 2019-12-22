// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/CmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename ) 
{

  float nCosmics=0;
  int nSecCosmics=0;

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

  TH1D *h_nPE[6][nChannels];
  for (unsigned int i=0; i<6; i++){
	for (unsigned int j=0; j<nChannels; j++){
		string cosmic = std::to_string(i);
		h_nPE[i][j] = new TH1D("nPE_cosmic"+cosmic+"_Ch"+ChStr[j],"cosmic " + cosmic + " nPE Ch " + ChStr[j], 1000, 0, 1000 );
	}
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
    if (!synchBoards ) continue;

/////////////////////////////////////////////////////////////////////////////
//	Muon selection I
//
    vector<vector<float>> pulses ;
    int evtSecCosmics=0;
    int evtCosmics=0;

    for ( unsigned int i=0; i<6; i++){
	vector<int> triggerChannels, layer1, layer2 ;
        vector<vector<float>> cosmicPulses, secCosmicPulses1, secCosmicPulses2 ;
	if ( i == 0 ){  triggerChannels =  trig100248; layer1 = trigL2; layer2 = trigL3; }
	if ( i == 1 ){  triggerChannels =  trig101259; layer1 = trigL2; layer2 = trigL3; }
	if ( i == 2 ){  triggerChannels = trig1161612; layer1 = trigL1; layer2 = trigL3; }
	if ( i == 3 ){  triggerChannels = trig1171713; layer1 = trigL1; layer2 = trigL3; }
	if ( i == 4 ){  triggerChannels =  trig142224; layer1 = trigL2; layer2 = trigL1; }
	if ( i == 5 ){  triggerChannels =  trig143235; layer1 = trigL2; layer2 = trigL1; }

	cosmicPulses = CosmicSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, triggerChannels );
	if ( cosmicPulses.size() == 0 ) continue;
	evtCosmics=1;

	secCosmicPulses1 = MuonPulses ( chan,height,ptime,nPE,time_module_calibrated, layer1 );
	secCosmicPulses2 = MuonPulses ( chan,height,ptime,nPE,time_module_calibrated, layer2 );
	if ( secCosmicPulses1.size() >= 2 || secCosmicPulses2.size() >= 2 ) evtSecCosmics=1;

	if ( pulses.size() == 0 ) pulses = SameLayerActivity_v2 ( chan,height,ptime,nPE,time_module_calibrated );

	vector<float> vec_nPE(nChannels,0);	
	for (unsigned int j=0; j<pulses.size(); j++){
		vec_nPE[(int)pulses.at(j).at(0)] += pulses.at(j).at(3);
        }

	for ( unsigned int j=0; j<nChannels; j++){
		h_nPE[i][j]->Fill ( vec_nPE[j] );
	}

     }
     nCosmics+=evtCosmics;
     nSecCosmics+=evtSecCosmics;
	
//	end of main loop
  }
cout << "ratio of secondary cosmic: " << nSecCosmics/nCosmics << endl;


  outFile->Write();
  outFile->Close();

}
