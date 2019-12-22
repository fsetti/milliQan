// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename, bool MC = false ) 
{

  double weight = 1;
  double weight_qcd = 0.014301851 * totLumi;
  double weight_qcd_nonbc = 0.0033202620 * totLumi;
  double weight_DY = 0.00011057182 * totLumi;
  double weight_W = 0.00026728950 * totLumi;
  TFile *currentFile;
  TString currentFileName;

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

  TH1D *h_nPE, *h_nPEET, *h_nPER878, *h_nPER7725, *h_dT;
  float increment = (maxBin-minBin)/nBins;
  Double_t edges[nBins+1];
  float x=minBin;
  for ( unsigned int i=0; i<nBins+1; i++){
	edges[i]= pow(10,x);
	x+=increment;
  }

  double nPEbins[23] = { 0., 10., 20., 30., 50., 70., 90., 120., 150., 180., 210., 240., 270., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750. };

  h_nPE = new TH1D("SimData_nPE","nPE",  22, nPEbins  );
  h_nPEET = new TH1D("SimData_nPEET","nPE",  22, nPEbins  );
  h_nPER878 = new TH1D("SimData_nPER878","nPE",  22, nPEbins  );
  h_nPER7725 = new TH1D("SimData_nPER7725","nPE",  22, nPEbins  );

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

    bool synchBoards = false;                   // only consider events with synched boards
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true;
    if ( (!synchBoards  && !MC) || ( !beam && !MC ) ) continue;


    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > 12 ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

   

/////////////////////////////////////////////////////////////////////////////
//	Muon Selection - hits in straight line 
    bool muonHit = false;
    muonHit = MuonInBars ( chan,ptime,nPE, trigBars );
    if ( muonHit ) continue;
  
 
    vector<vector<float>> L1Pulses, L2Pulses, L3Pulses, pulses;
    L1Pulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigL1 );
    L2Pulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigL2 );
    L3Pulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigL3 );

    if ( !MC ){
	    for ( unsigned int i=0; i<L1Pulses.size(); i++){
		if ( L1Pulses.at(i).at(0) == 5 || L1Pulses.at(i).at(0) == 22 ) h_nPER7725->Fill( L1Pulses.at(i).at(3) );
		if ( L1Pulses.at(i).at(0) == 9 || L1Pulses.at(i).at(0) == 17 || L1Pulses.at(i).at(0) == 24 || L1Pulses.at(i).at(0) == 25 ) h_nPEET->Fill( L1Pulses.at(i).at(3) );
		if ( L1Pulses.at(i).at(0) != 5 && L1Pulses.at(i).at(0) != 22 && L1Pulses.at(i).at(0) != 9 && L1Pulses.at(i).at(0) != 17 && L1Pulses.at(i).at(0) != 24 && L1Pulses.at(i).at(0) != 25 ) h_nPER878->Fill( L1Pulses.at(i).at(3) );
		h_nPE->Fill( L1Pulses.at(i).at(3) );
		h_nPEdTall->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
		h_nPEdT[(int)L1Pulses.at(i).at(0)]->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
		if ( (int)L1Pulses.at(i).at(0) == 5 || (int)L1Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
		if ( (int)L1Pulses.at(i).at(0) == 24 || (int)L1Pulses.at(i).at(0) == 25 || (int)L1Pulses.at(i).at(0) == 9 || (int)L1Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
		if ( (int)L1Pulses.at(i).at(0) != 5 && (int)L1Pulses.at(i).at(0) != 22 && (int)L1Pulses.at(i).at(0) != 24 && (int)L1Pulses.at(i).at(0) != 25 &&(int)L1Pulses.at(i).at(0) != 9 && (int)L1Pulses.at(i).at(0) != 17 )  h_nPEdTR878->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) );
	    }
	    for ( unsigned int i=0; i<L2Pulses.size(); i++){
		if ( L2Pulses.at(i).at(0) == 5 || L2Pulses.at(i).at(0) == 22 ) h_nPER7725->Fill( L2Pulses.at(i).at(3) );
		if ( L2Pulses.at(i).at(0) == 9 || L2Pulses.at(i).at(0) == 17 || L2Pulses.at(i).at(0) == 24 || L2Pulses.at(i).at(0) == 25 ) h_nPEET->Fill( L2Pulses.at(i).at(3) );
		if ( L2Pulses.at(i).at(0) != 5 && L2Pulses.at(i).at(0) != 22 && L2Pulses.at(i).at(0) != 9 && L2Pulses.at(i).at(0) != 17 && L2Pulses.at(i).at(0) != 24 && L2Pulses.at(i).at(0) != 25 ) h_nPER878->Fill( L2Pulses.at(i).at(3) );
		h_nPE->Fill( L2Pulses.at(i).at(3) );
		h_nPEdTall->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
		h_nPEdT[(int)L2Pulses.at(i).at(0)]->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
		if ( (int)L2Pulses.at(i).at(0) == 5 || (int)L2Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
		if ( (int)L2Pulses.at(i).at(0) == 24 || (int)L2Pulses.at(i).at(0) == 25 || (int)L2Pulses.at(i).at(0) == 9 || (int)L2Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
		if ( (int)L2Pulses.at(i).at(0) != 5 && (int)L2Pulses.at(i).at(0) != 22 && (int)L2Pulses.at(i).at(0) != 24 && (int)L2Pulses.at(i).at(0) != 25 &&(int)L2Pulses.at(i).at(0) != 9 && (int)L2Pulses.at(i).at(0) != 17 ) h_nPEdTR878->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) );
	    }
	    for ( unsigned int i=0; i<L3Pulses.size(); i++){
		if ( L3Pulses.at(i).at(0) == 5 || L3Pulses.at(i).at(0) == 22 ) h_nPER7725->Fill( L3Pulses.at(i).at(3) );
		if ( L3Pulses.at(i).at(0) == 9 || L3Pulses.at(i).at(0) == 17 || L3Pulses.at(i).at(0) == 24 || L3Pulses.at(i).at(0) == 25 ) h_nPEET->Fill( L3Pulses.at(i).at(3) );
		if ( L3Pulses.at(i).at(0) != 5 && L3Pulses.at(i).at(0) != 22 && L3Pulses.at(i).at(0) != 9 && L3Pulses.at(i).at(0) != 17 && L3Pulses.at(i).at(0) != 24 && L3Pulses.at(i).at(0) != 25 ) h_nPER878->Fill( L3Pulses.at(i).at(3) );
		h_nPE->Fill( L3Pulses.at(i).at(3) );
		h_nPEdTall->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
		h_nPEdT[(int)L3Pulses.at(i).at(0)]->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
		if ( (int)L3Pulses.at(i).at(0) == 5 || (int)L3Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
		if ( (int)L3Pulses.at(i).at(0) == 24 || (int)L3Pulses.at(i).at(0) == 25 || (int)L3Pulses.at(i).at(0) == 9 || (int)L3Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
		if ( (int)L3Pulses.at(i).at(0) != 5 && (int)L3Pulses.at(i).at(0) != 22 && (int)L3Pulses.at(i).at(0) != 24 && (int)L3Pulses.at(i).at(0) != 25 &&(int)L3Pulses.at(i).at(0) != 9 && (int)L3Pulses.at(i).at(0) != 17 ) h_nPEdTR878->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) );
	    }


    }


//	repeat for MC with weights
    if ( MC ){
	currentFile = chain->GetFile();
	currentFileName = currentFile->GetName();
	if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcd_all.root" ) weight = weight_qcd;
	if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcdnonbc_all.root" ) weight = weight_qcd_nonbc;
	if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_dy_all.root" ) weight = weight_DY;
	if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_w_all.root" ) weight = weight_W;

	for ( unsigned int i=0; i<L1Pulses.size(); i++){
		if ( L1Pulses.at(i).at(0) == 5 || L1Pulses.at(i).at(0) == 22 ) h_nPER7725->Fill( L1Pulses.at(i).at(3) , weight);
		if ( L1Pulses.at(i).at(0) == 9 || L1Pulses.at(i).at(0) == 17 || L1Pulses.at(i).at(0) == 24 || L1Pulses.at(i).at(0) == 25 ) h_nPEET->Fill( L1Pulses.at(i).at(3) , weight);
		if ( L1Pulses.at(i).at(0) != 5 && L1Pulses.at(i).at(0) != 22 && L1Pulses.at(i).at(0) != 9 && L1Pulses.at(i).at(0) != 17 && L1Pulses.at(i).at(0) != 24 && L1Pulses.at(i).at(0) != 25 ) h_nPER878->Fill( L1Pulses.at(i).at(3) , weight);
	    h_nPE->Fill( L1Pulses.at(i).at(3), weight );
	    h_nPEdTall->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3), weight );
	    h_nPEdT[(int)L1Pulses.at(i).at(0)]->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) , weight );
	    if ( (int)L1Pulses.at(i).at(0) == 5 || (int)L1Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) , weight );
	    if ( (int)L1Pulses.at(i).at(0) == 24 || (int)L1Pulses.at(i).at(0) == 25 || (int)L1Pulses.at(i).at(0) == 9 || (int)L1Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3) , weight );
	    if ( (int)L1Pulses.at(i).at(0) != 5 && (int)L1Pulses.at(i).at(0) != 22 && (int)L1Pulses.at(i).at(0) != 24 && (int)L1Pulses.at(i).at(0) != 25 &&(int)L1Pulses.at(i).at(0) != 9 && (int)L1Pulses.at(i).at(0) != 17 )  h_nPEdTR878->Fill( - slabPulses.at(1).at(4) + L1Pulses.at(i).at(4), L1Pulses.at(i).at(3), weight );
	}
	for ( unsigned int i=0; i<L2Pulses.size(); i++){
		if ( L2Pulses.at(i).at(0) == 5 || L2Pulses.at(i).at(0) == 22 ) h_nPER7725->Fill( L2Pulses.at(i).at(3) , weight);
		if ( L2Pulses.at(i).at(0) == 9 || L2Pulses.at(i).at(0) == 17 || L2Pulses.at(i).at(0) == 24 || L2Pulses.at(i).at(0) == 25 ) h_nPEET->Fill( L2Pulses.at(i).at(3) , weight);
		if ( L2Pulses.at(i).at(0) != 5 && L2Pulses.at(i).at(0) != 22 && L2Pulses.at(i).at(0) != 9 && L2Pulses.at(i).at(0) != 17 && L2Pulses.at(i).at(0) != 24 && L2Pulses.at(i).at(0) != 25 ) h_nPER878->Fill( L2Pulses.at(i).at(3) , weight );
	    h_nPE->Fill( L2Pulses.at(i).at(3) , weight);
	    h_nPEdTall->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) , weight);
	    h_nPEdT[(int)L2Pulses.at(i).at(0)]->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) , weight );
	    if ( (int)L2Pulses.at(i).at(0) == 5 || (int)L2Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) , weight );
	    if ( (int)L2Pulses.at(i).at(0) == 24 || (int)L2Pulses.at(i).at(0) == 25 || (int)L2Pulses.at(i).at(0) == 9 || (int)L2Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3) , weight);
	    if ( (int)L2Pulses.at(i).at(0) != 5 && (int)L2Pulses.at(i).at(0) != 22 && (int)L2Pulses.at(i).at(0) != 24 && (int)L2Pulses.at(i).at(0) != 25 &&(int)L2Pulses.at(i).at(0) != 9 && (int)L2Pulses.at(i).at(0) != 17 ) h_nPEdTR878->Fill( - slabPulses.at(3).at(4) + L2Pulses.at(i).at(4), L2Pulses.at(i).at(3), weight );
	}
	for ( unsigned int i=0; i<L3Pulses.size(); i++){
		if ( L3Pulses.at(i).at(0) == 5 || L3Pulses.at(i).at(0) == 22 ) h_nPER7725->Fill( L3Pulses.at(i).at(3) , weight);
		if ( L3Pulses.at(i).at(0) == 9 || L3Pulses.at(i).at(0) == 17 || L3Pulses.at(i).at(0) == 24 || L3Pulses.at(i).at(0) == 25 ) h_nPEET->Fill( L3Pulses.at(i).at(3) , weight);
		if ( L3Pulses.at(i).at(0) != 5 && L3Pulses.at(i).at(0) != 22 && L3Pulses.at(i).at(0) != 9 && L3Pulses.at(i).at(0) != 17 && L3Pulses.at(i).at(0) != 24 && L3Pulses.at(i).at(0) != 25 ) h_nPER878->Fill( L3Pulses.at(i).at(3) , weight);
	    h_nPE->Fill( L3Pulses.at(i).at(3) , weight);
	    h_nPEdTall->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) , weight);
	    h_nPEdT[(int)L3Pulses.at(i).at(0)]->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) , weight);
	    if ( (int)L3Pulses.at(i).at(0) == 5 || (int)L3Pulses.at(i).at(0) == 22 ) h_nPEdTR7725->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) , weight);
	    if ( (int)L3Pulses.at(i).at(0) == 24 || (int)L3Pulses.at(i).at(0) == 25 || (int)L3Pulses.at(i).at(0) == 9 || (int)L3Pulses.at(i).at(0) == 17 ) h_nPEdTET->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) , weight);
	    if ( (int)L3Pulses.at(i).at(0) != 5 && (int)L3Pulses.at(i).at(0) != 22 && (int)L3Pulses.at(i).at(0) != 24 && (int)L3Pulses.at(i).at(0) != 25 &&(int)L3Pulses.at(i).at(0) != 9 && (int)L3Pulses.at(i).at(0) != 17 ) h_nPEdTR878->Fill( - slabPulses.at(2).at(4) + L3Pulses.at(i).at(4), L3Pulses.at(i).at(3) , weight );
	}
    }

//	end of main loop
  }



  outFile->cd();

for (unsigned int i=0; i<nChannels;i++){
      if ( h_nPEdT[i]->GetEntries() != 0 ) h_nPEdT[i]->Write();
  }
  h_nPER878->Write();
  h_nPER7725->Write();
  h_nPEET->Write();

  h_nPEdTall->Write();
  h_nPEdTR878->Write();
  h_nPEdTR7725->Write();
  h_nPEdTET->Write();
  outFile->Write();
  outFile->Close();

}
