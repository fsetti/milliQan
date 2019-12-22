// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/SimDataComparison/beamMuons/muonAnalysis.cc"
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

  TH1D *h_nPEall, *h_NnPEall, *h_nNnPEall, *h_NnPE[nChannels], *h_nNnPE[nChannels], *h_nPE[nChannels];
  TH2D *h_nPEdTall, *h_nPEdTR878, *h_nPEdTET, *h_nPEdTR7725, *h_nPEdT[nChannels];
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

  h_nPEdTall = new TH2D("SimData_nPEdT","nPE vs dT", 80, -20, 60, nBins, edges );
  h_nPEdTR878 = new TH2D("SimData_nPEdT_R878","nPE vs dT", 80, -20, 60, nBins, edges );
  h_nPEdTET = new TH2D("SimData_nPEdT_ET","nPE vs dT", 80, -20, 60, nBins, edges );
  h_nPEdTR7725 = new TH2D("SimData_nPEdT_R7725","nPE vs dT", 80, -20, 60, nBins, edges );

  h_nPEall = new TH1D("SimData_nPE","nPE of small pulse",75, 0, 750); 
  h_NnPEall = new TH1D("SimData_NnPE","nPE of neigh. bars",75, 0, 750); 
  h_nNnPEall = new TH1D("SimData_nNnPE","nPE of non-neigh. bars",75, 0, 750); 
  for (unsigned int i=0; i<nChannels; i++){
	h_nPE[i] = new TH1D("SimData_nPE_Ch"+ChStr[i],"nPE of small pulses",75,0,750);
	h_NnPE[i] = new TH1D("SimData_NnPE_Ch"+ChStr[i],"nPE of neigh. bars",75,0,750);
	h_nNnPE[i] = new TH1D("SimData_nNnPE_Ch"+ChStr[i],"nPE of non-neigh. bars",75,0,750);
  	h_nPEdT[i] = new TH2D("SimData_nPEdT_Ch"+ChStr[i],"nPE vs dT", 80, -20, 60, nBins, edges );
  }

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


/////////////////////////////////////////////////////////////////////////////
//	Muon selection I
//
    vector<vector<vector<float>>> neighbouringPulses, nonNeighbouringPulses;
    nonNeighbouringPulses = SameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig08,trig09,trig18,trig19,trig612,trig613,trig712,trig713,trig24,trig25,trig34,trig35}, false );
    neighbouringPulses = SameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig01,trig024,trig025,trig124,trig125,trig2425,trig248,trig249,trig258,trig259,trig89,trig67,trig616,trig617,trig716,trig717,trig1617,trig1612,trig1613,trig1712,trig1713,trig1213,trig23,trig222,trig223,trig322,trig323,trig2223,trig224,trig225,trig234,trig235,trig45});

    for ( unsigned int i=0; i<neighbouringPulses.size(); i++){
        h_NnPE[(int)neighbouringPulses.at(i).at(1).at(0)]->Fill( neighbouringPulses.at(i).at(1).at(3) );
	h_NnPEall->Fill( neighbouringPulses.at(i).at(1).at(3) );
    }
    for ( unsigned int i=0; i<nonNeighbouringPulses.size(); i++){
        h_nNnPE[(int)nonNeighbouringPulses.at(i).at(1).at(0)]->Fill( nonNeighbouringPulses.at(i).at(1).at(3) );
	h_nNnPEall->Fill( nonNeighbouringPulses.at(i).at(1).at(3) );
    }



/////////////////////////////////////////////////////////////////////////////
//	Muon selection II

//    bool muonHits;	//	Require No Muon Hits in any Bar
//    muonHits = MuonInBars ( chan,ptime,nPE, trigBars ); 
//    if ( !muonHits ) continue;


    vector<vector<float>> smallPulses;
    smallPulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigBars ); 

    for ( unsigned int i=0; i<smallPulses.size(); i++){
        h_nPE[(int)smallPulses.at(i).at(0)]->Fill( smallPulses.at(i).at(3) );
	h_nPEall->Fill( smallPulses.at(i).at(3) );

        if ( (int)smallPulses.at(i).at(0) == 0 || (int)smallPulses.at(i).at(0) == 1 ||(int)smallPulses.at(i).at(0) == 24 ||(int)smallPulses.at(i).at(0) == 25 ||(int)smallPulses.at(i).at(0) == 8 ||(int)smallPulses.at(i).at(0) == 9 ){
		h_nPEdTall->Fill(  - slabPulses.at(1).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
		h_nPEdT[(int)smallPulses.at(i).at(0)]->Fill(  - slabPulses.at(1).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
	}
        if ( (int)smallPulses.at(i).at(0) == 6 || (int)smallPulses.at(i).at(0) == 7 ||(int)smallPulses.at(i).at(0) == 16 ||(int)smallPulses.at(i).at(0) == 17 ||(int)smallPulses.at(i).at(0) == 12 ||(int)smallPulses.at(i).at(0) == 13 ){
		h_nPEdTall->Fill( - slabPulses.at(3).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
		h_nPEdT[(int)smallPulses.at(i).at(0)]->Fill( - slabPulses.at(3).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
	}
        if ( (int)smallPulses.at(i).at(0) == 2 || (int)smallPulses.at(i).at(0) == 3 ||(int)smallPulses.at(i).at(0) == 22 ||(int)smallPulses.at(i).at(0) == 23 ||(int)smallPulses.at(i).at(0) == 4 ||(int)smallPulses.at(i).at(0) == 5 ){ 
		h_nPEdTall->Fill( - slabPulses.at(2).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
		h_nPEdT[(int)smallPulses.at(i).at(0)]->Fill( - slabPulses.at(2).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
	}
	

        if ( (int)smallPulses.at(i).at(0) != 9 && (int)smallPulses.at(i).at(0) != 17 && (int)smallPulses.at(i).at(0) != 24 && (int)smallPulses.at(i).at(0) != 25 
	   && (int)smallPulses.at(i).at(0) != 5 && (int)smallPulses.at(i).at(0) != 22 ){ 
		h_nPEdTR878->Fill( - slabPulses.at(2).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
	}
        if ( (int)smallPulses.at(i).at(0) == 9 || (int)smallPulses.at(i).at(0) == 17 || (int)smallPulses.at(i).at(0) == 24 ||(int)smallPulses.at(i).at(0) == 25 ){ 
		h_nPEdTET->Fill( - slabPulses.at(2).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
	}
        if ( (int)smallPulses.at(i).at(0) == 5 || (int)smallPulses.at(i).at(0) == 22 ){ 
		h_nPEdTR7725->Fill( - slabPulses.at(2).at(4) + smallPulses.at(i).at(4), smallPulses.at(i).at(3) );
	}


    }


    for ( unsigned int i=0; i<slabPulses.size(); i++){
        h_nPE[(int)slabPulses.at(i).at(0)]->Fill( slabPulses.at(i).at(3) );
    }


//	end of main loop
  }



  outFile->cd();

  for (unsigned int i=0; i<nChannels; i++){
	if ( h_nPE[i]->GetEntries() == 0 ) continue;
	h_nPE[i]->Write();
	h_NnPE[i]->Write();
	h_nNnPE[i]->Write();
  	h_nPEdT[i]->Write();
  }

  h_nPEdTall->Write();
  h_nPEdTR878->Write();
  h_nPEdTET->Write();
  h_nPEdTR7725->Write();
  h_nPEall->Write();
  h_NnPEall->Write();
  h_nNnPEall->Write();
  outFile->Write();
  outFile->Close();

}
