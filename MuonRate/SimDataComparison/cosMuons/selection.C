// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/SimDataComparison/cosMuons/muonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename ) 
{

  int nCosmics = 0;
  float nMultipleCosmics = 0;

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

  TH1D *h_nPEall, *h_nPE[nChannels];
  TH2D *h_nPEdTall, *h_nPEdT[nChannels];
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

  h_nPEdTall = new TH2D("SimDataC_nPEdT","nPE vs dT", 100, -20, 60, nBins, edges );
  h_nPEall = new TH1D("SimDataC_nPE","nPE of small pulse",100, 0, 1000); 
  for (unsigned int i=0; i<nChannels; i++){
	h_nPE[i] = new TH1D("SimDataC_nPE_Ch"+ChStr[i],"nPE of small pulses",100,0,1000);
  	h_nPEdT[i] = new TH2D("SimDataC_nPEdT_Ch"+ChStr[i],"nPE vs dT", 100, -20, 60, nBins, edges );
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
    vector<vector<vector<float>>> cosmicPulses, secondaryPulses;
    cosmicPulses = CosmicSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig100248,trig101259,trig1161612,trig1171713,trig142224,trig143235} );
    secondaryPulses = SameLayerActivity ( chan,height,ptime,nPE,time_module_calibrated, cosmicPulses );

    if ( cosmicPulses.size() != 0 ) nCosmics++;
    if ( cosmicPulses.size()  > 1 ) nMultipleCosmics++;

    for ( unsigned int i=0; i< secondaryPulses.size(); i++){
	for (unsigned int j=0; j<secondaryPulses.at(i).size(); j++){
		h_nPE[(int)secondaryPulses.at(i).at(j).at(0)]->Fill( secondaryPulses.at(i).at(j).at(3) );
		h_nPEall->Fill( secondaryPulses.at(i).at(j).at(3) );
	}
    }

    if ( cosmicPulses.size() == secondaryPulses.size() ){
	for (unsigned int i=0; i<cosmicPulses.size(); i++){
		std::sort( cosmicPulses.at(i).begin(), cosmicPulses.at(i).end(), sort_wrt_channel );
		for (unsigned int j=0; j<secondaryPulses.at(i).size(); j++){
			h_nPEdTall->Fill( secondaryPulses.at(i).at(j).at(4) -  cosmicPulses.at(i).at(3).at(4), secondaryPulses.at(i).at(j).at(3));
			h_nPEdT[(int)secondaryPulses.at(i).at(j).at(0)]->Fill( secondaryPulses.at(i).at(j).at(4) -  cosmicPulses.at(i).at(3).at(4), secondaryPulses.at(i).at(j).at(3));
		}
	}
     }	

//	end of main loop
  }


  cout << " ratio of multiple cosmics within same event: " << nMultipleCosmics/nCosmics << endl;
  cout << "events with cosmics: " << nCosmics << endl;
  outFile->cd();

  for (unsigned int i=0; i<nChannels; i++){
	if ( h_nPE[i]->GetEntries() == 0 ) continue;
	h_nPE[i]->Write();
  	h_nPEdT[i]->Write();
  }
  h_nPEdTall->Write();
  h_nPEall->Write();
  outFile->Write();
  outFile->Close();

}
