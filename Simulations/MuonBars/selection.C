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
    vector<vector<vector<float>>> muonPulses;
    muonPulses = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig8124,trig9135} );
    if ( muonPulses.size() == 0 ) continue;

    vector<vector<vector<float>>> nonNeighPulses;
    nonNeighPulses = SameLayerSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig18,trig19,trig08,trig09,trig612,trig613,trig712,trig713,trig24,trig25,trig34,trig35}, false );

    
    for ( unsigned int i=0; i<nonNeighPulses.size(); i++){
	h_nPEdTall->Fill ( nonNeighPulses.at(i).at(1).at(4) - nonNeighPulses.at(i).at(0).at(4), nonNeighPulses.at(i).at(1).at(3) );
	h_nPEdT[(int)nonNeighPulses.at(i).at(1).at(0)]->Fill ( nonNeighPulses.at(i).at(1).at(4) - nonNeighPulses.at(i).at(0).at(4), nonNeighPulses.at(i).at(1).at(3) );
	if ( nonNeighPulses.at(i).at(1).at(0) == 5 || nonNeighPulses.at(i).at(1).at(0) == 22 ){
		h_nPEdTR7725->Fill ( nonNeighPulses.at(i).at(1).at(4) - nonNeighPulses.at(i).at(0).at(4), nonNeighPulses.at(i).at(1).at(3) );
	}
	if ( nonNeighPulses.at(i).at(1).at(0) == 9 || nonNeighPulses.at(i).at(1).at(0) == 17 || nonNeighPulses.at(i).at(1).at(0) == 24 || nonNeighPulses.at(i).at(1).at(0) == 25 ){
		h_nPEdTET->Fill ( nonNeighPulses.at(i).at(1).at(4) - nonNeighPulses.at(i).at(0).at(4), nonNeighPulses.at(i).at(1).at(3) );
	}
	if ( nonNeighPulses.at(i).at(1).at(0) != 9 && nonNeighPulses.at(i).at(1).at(0) != 17 && nonNeighPulses.at(i).at(1).at(0) != 24 && nonNeighPulses.at(i).at(1).at(0) != 25 && nonNeighPulses.at(i).at(1).at(0) != 5 && nonNeighPulses.at(i).at(1).at(0) != 22 ){
		h_nPEdTR878->Fill ( nonNeighPulses.at(i).at(1).at(4) - nonNeighPulses.at(i).at(0).at(4), nonNeighPulses.at(i).at(1).at(3) );
	}
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
