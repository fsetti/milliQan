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

  TH1D *h_nPE,*h_lnPE, *h_dT;
  float increment = (maxBin-minBin)/nBins;
  Double_t edges[nBins+1];
  float x=minBin;
  for ( unsigned int i=0; i<nBins+1; i++){
	edges[i]= pow(10,x);
	x+=increment;
  }

  double nPEbins[23] = { 0., 10., 20., 30., 50., 70., 90., 120., 150., 180., 210., 240., 270., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750. };

  h_nPE = new TH1D("SimData_nPE","nPE",  25, 0,750 );
  h_lnPE = new TH1D("SimData_lnPE","nPE",  25, 0,750  );

  TH2D *h_nPEdTall, *h_nPEdT[nChannels], *h_MCnPE;
  h_nPEdTall = new TH2D("SimData_nPEdTall","nPE vs dT", 100, -50, 50, nBins, edges );
  h_MCnPE = new TH2D("SimData_MCnPE","reco vs MC nPW", nBins , edges, nBins, edges );
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


// check for true MC nPe

//    if ( MC ) {
//    int   MCnPE[32];
//    TBranch        *b_MCnPE;   //!
//    chain->SetBranchAddress("chan_trueNPE", &MCnPE, &b_MCnPE);
//    int selection[32];
//	for (unsigned int i=0; i<chan->size(); i++){
//		if ( nPE->at(i) < 50000  && selection[chan->at(i)] != 1  
//		 &&  chan->at(i) != 10 &&  chan->at(i) != 11&&  chan->at(i) != 14&&  chan->at(i) != 18&&  chan->at(i) != 19&&  chan->at(i) != 20&&  chan->at(i) != 21&&  chan->at(i) != 26&&  chan->at(i) != 27&&  chan->at(i) != 28&&  chan->at(i) != 29&&  chan->at(i) != 30&&  chan->at(i) != 31 && chan->at(i) != 15 ){ 
//			h_MCnPE->Fill( MCnPE[chan->at(i)], nPE->at(i) );
//			selection[chan->at(i)] = 1;
//		}
//	}
//    }






/////////////////////////////////////////////////////////////////////////////
//	Slab no Bar Selection - muon hits in panels but NOT in bars 
    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > 12 ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

    bool muonHit = false;
    muonHit = MuonInBars ( chan,ptime,nPE, trigBars );
    if ( muonHit ) continue;




/////////////////////////////////////////////////////////////////////////////
//	three hits in a line 
  
    vector<vector<vector<float>>> lpulses;
    lpulses = selectPulses_signaLike_line ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig241622,trig251723,trig8124,trig9135} );

    if ( !MC ){
	for (unsigned int i=0; i<lpulses.size(); i++){	
		std::sort( lpulses.at(i).begin(), lpulses.at(i).end(), sort_wrt_nPE );
		h_lnPE->Fill ( lpulses.at(i).at(0).at(3) );
	}
    }


    if ( MC ){
        currentFile = chain->GetFile();
        currentFileName = currentFile->GetName();
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcd.root" ) weight = weight_qcd;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcd_nonbc.root" ) weight = weight_qcd_nonbc;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_dy.root" ) weight = weight_DY;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_w.root" ) weight = weight_W;
	for (unsigned int i=0; i<lpulses.size(); i++){	
		std::sort( lpulses.at(i).begin(), lpulses.at(i).end(), sort_wrt_nPE );
		h_lnPE->Fill ( lpulses.at(i).at(0).at(3) , weight );
	}
    }


/////////////////////////////////////////////////////////////////////////////
//	exactly one hit in each layer
//
    vector<vector<float>> L1Pulses, L2Pulses, L3Pulses, pulses;
    L1Pulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigL1 );
    L2Pulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigL2 );
    L3Pulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigL3 );


    if ( L1Pulses.size() != 1 || L2Pulses.size() != 1 || L3Pulses.size() != 1 ) continue ;

    pulses = selectPulses( chan,height,ptime,nPE,time_module_calibrated, trigBars );
    std::sort( pulses.begin(), pulses.end(), sort_wrt_nPE );

    if ( !MC ) {
	h_nPE->Fill ( pulses.at(0).at(3) );
    }

    if ( MC ) {
       currentFile = chain->GetFile();
        currentFileName = currentFile->GetName();
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcd.root" ) weight = weight_qcd;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcdnonbc.root" ) weight = weight_qcd_nonbc;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_dy.root" ) weight = weight_DY;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_w.root" ) weight = weight_W;
	h_nPE->Fill ( pulses.at(0).at(3) , weight );
    }


//	end of main loop
  }



  outFile->cd();

for (unsigned int i=0; i<nChannels;i++){
      if ( h_nPEdT[i]->GetEntries() != 0 ) h_nPEdT[i]->Write();
  }

  h_nPE->Write();
  h_MCnPE->Write();
  h_lnPE->Write();
  outFile->Write();
  outFile->Close();

}
