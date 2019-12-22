// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/muonAnalysis.cc"
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/milliHists.h"
#include "TF1.h"


void muonSelectionSkim( TChain *chain, TString output_filename, TString runName, vector<float> channels, bool synch = false, bool muonBarHit = true ) 
{

  float nSelEvts = 0;
  float nSlabEvts = 0;
  float nBarEvts = 0;
  float oneSameLayerEvts = 0;
  float twoSameLayerEvts = 0;
  float threeSameLayerEvts = 0;


  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


  vector<int> trigCh;
  for (unsigned int i=0; i<channels.size(); i++){
        trigCh.push_back((int)channels[i]);
  }

  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

  TH1D *h_tDiffSlab, *h_tDiffAll, *h_MuHitsSL, *h_dTlayerL, *h_dTlayer, *h_nPE;  
  TH2D *h_tDiff;
  TProfile *profile;

  h_nPE = new TH1D("RunSkim_nPE","nPE of small pulse",41, 0.1, 1e3); 
  h_tDiffAll = new TH1D("RunSkim_tDiffAll","#Delta_{t} cal",51, -125, 75); 
  h_tDiffSlab = new TH1D("RunSkim_tDiffSlab","#Delta_{t} cal",51, -125, 75); 
  h_dTlayer = new TH1D("RunSkim_dTlayer","#Delta_{t} cal",31, -35, 75); 
  h_dTlayerL = new TH1D("RunSkim_dTlayerL","#Delta_{t} cal",31, -35, 75); 
  h_MuHitsSL = new TH1D("RunSkim_MuHitsSL","#mu hits in same layer",6, 0, 6); 
  h_tDiff = new TH2D("RunSkim_tDiff","#Delta_{t} #mu vs #Delta_{t} sec pulse" ,31,0,800, 21,-150,80 ); 
  profile = new TProfile("profile","profile of #Delta_{t} vs nPE" ,47,0.1,800, -125, 75 ); 

  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();

// Plots of the event time.
// First we need to find the min and max times and the total live time.
// We also count the number of events with synched digitizers.
  int nSynched = 0;
  int nUnSynched = 0;
  double minTime = 9.E15; double maxTime = -9.E15;
  int firstEvt = 0; int lastEvt = 0;
  int maxFile = 0;
  double SumTime = 0.; // Sum of all live time
  double SumSynchTime = 0.; // Sum of synched live time
  double prevTime = 0.;
  bool synched = false;
  string FirstTime="N/A", LastTime="N/A";
  for(int ia = 0; ia<nentries; ia++){ // Loop over all events
    chain->GetEntry(ia);
    if ( !beam ) continue;
    if (filenum>maxFile) maxFile = filenum;
    if ((event>2 || filenum>1) && event_time_fromTDC>1.5E9) { // Ignore the first two events in file1 which can be leftover from previous run; ignore bad time measurements
      if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) {
        nSynched++;
        synched = true;
      }
      else {
        nUnSynched++;
        synched = false;
      }
      if (prevTime>0 && ((filenum==1 && event>3) || (filenum>1 && event>1))) {
        SumTime += event_time_fromTDC-prevTime;
        if (synched) SumSynchTime += event_time_fromTDC-prevTime;
      }
      prevTime = event_time_fromTDC;
     if (event_time_fromTDC < minTime) { // New first event
      minTime = event_time_fromTDC;
      firstEvt = ia;
      FirstTime = *event_t_string;
     } // New first event
     if (event_time_fromTDC > maxTime) { // New last event
      maxTime = event_time_fromTDC;
      lastEvt = ia;
      LastTime = *event_t_string;
     } // New last event
   } // Ignore the first events, which may be buggy
  } // Loop over all events
  double runTime;
  if ( !synch ) runTime = SumTime;
  if (  synch ) runTime = (float)nSynched/(nSynched+nUnSynched)*SumTime;
  double weight = 1./runTime;




/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);
    vector<int> trigChStatus(trigCh.size(),0);  // set no pulses in tag channel at the beginning of every event

//	skip non-synchronised events && require beam On
    bool synchBoards = false;			// only consider events with synched boards
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true; 
    if (!synchBoards && synch) continue;
    if ( !beam ) continue;
    nSelEvts++;



/////////////////////////////////////////////////////////////////////////////
//	Slab selection
    vector<vector<vector<float>>>  totPulses;
    totPulses = pulseSlabSelection( chan,height,ptime,nPE,time_module_calibrated, {trigCh} );
    if ( totPulses.size() == 0 ) continue;
    vector<vector<float>> SlabPulses;
    SlabPulses = totPulses.at(0);
    if ( SlabPulses.size() != trigCh.size() ) continue;
    std::sort(SlabPulses.begin(), SlabPulses.end(), sort_wrt_channel);
    float dt_cal = SlabPulses.at(2)[4] - SlabPulses.at(0)[4];
// require dt_cal between -12 and +12 ns from dt_cal distribution
    if ( fabs( dt_cal ) > 12 ) continue;
    nSlabEvts++;




/////////////////////////////////////////////////////////////////////////////
//	Layer selection
//	Bar selection I, require through going particles in upper, lower channels 
    vector<vector<vector<float>>> throughPulses;
   if ( muonBarHit ) throughPulses = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig8124,trig9135,trig063,trig172,trig8125,trig9134, trig162,trig073,trig9124,trig8135} );
   if ( !muonBarHit) throughPulses = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig241622,trig251723,trig8124,trig9135, trig063,trig172,trig241623,trig251722,trig8125,trig9134, trig162,trig073,trig251622,trig241723,trig9124,trig8135} );
 
   bool throughParticle = false; 
   int nBarHits=0;
   for (unsigned int i=0; i<throughPulses.size(); i++){
	if (  nBarHits == 0 ){ 
		nBarHits++ ;
		throughParticle = true;
		break;
	}
   }
   nBarEvts+=nBarHits; 
   if ( muonBarHit && !throughParticle ) continue;
   if ( !muonBarHit && throughParticle ) continue;


/////////////////////////////////////////////////////////////////////////////
//	Layer selection
//	Bar selection II, store largest and second largest pulses in each layer 
  vector<vector<vector<float>>> L1Pulses, L2Pulses, L3Pulses;
  vector<vector<float>> L1PulsesNoHit, L2PulsesNoHit, L3PulsesNoHit;
  if ( muonBarHit ){
	  L1Pulses = sameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig08,trig18,trig09,trig19} );
	  L2Pulses = sameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig612,trig613,trig712,trig713} );
	  L3Pulses = sameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig24,trig25,trig34,trig35} );
  }

  if ( !muonBarHit ){
	  L1PulsesNoHit = sameLayerSelectionNoHit( chan,height,ptime,nPE,time_module_calibrated, {0,1,24,25,8,9} );
	  L2PulsesNoHit = sameLayerSelectionNoHit( chan,height,ptime,nPE,time_module_calibrated, {6,7,16,17,12,13} );
	  L3PulsesNoHit = sameLayerSelectionNoHit( chan,height,ptime,nPE,time_module_calibrated, {2,3,22,23,4,5} );
//	  L1Pulses = sameLayerSelectionNoHit_v2( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig063,trig172,trig162,trig073} );
//	  L2Pulses = sameLayerSelectionNoHit_v2( chan,height,ptime,nPE,time_module_calibrated, {trig241622,trig251723,trig241623,trig251722,trig251622,trig241723} );
//	  L3Pulses = sameLayerSelectionNoHit_v2( chan,height,ptime,nPE,time_module_calibrated, {trig8124,trig9135,trig8125,trig9134,trig9124,trig8135} );
//	  L1Pulses = sameLayerSelectionNoHit_v2( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173} );
//	  L2Pulses = sameLayerSelectionNoHit_v2( chan,height,ptime,nPE,time_module_calibrated, {trig241622,trig251723} );
//	  L3Pulses = sameLayerSelectionNoHit_v2( chan,height,ptime,nPE,time_module_calibrated, {trig8124,trig9135} );
   }

   int nSameLayerHitL1 = 0;
   int nSameLayerHitL2 = 0;
   int nSameLayerHitL3 = 0;
   if ( muonBarHit ){
	   for (unsigned int i=0; i<L1Pulses.size(); i++){
		nSameLayerHitL1++;
		break; 
	   }
	   for (unsigned int i=0; i<L2Pulses.size(); i++){
		nSameLayerHitL2++;
		break;
	   }
	   for (unsigned int i=0; i<L3Pulses.size(); i++){
		nSameLayerHitL3++;
		break; 
	   }
    }
    if ( !muonBarHit ){
	if ( L1PulsesNoHit.size() != 0 ) nSameLayerHitL1++;
	if ( L2PulsesNoHit.size() != 0 ) nSameLayerHitL2++;
	if ( L3PulsesNoHit.size() != 0 ) nSameLayerHitL3++;
    }

//   if ( nSameLayerHitL1 + nSameLayerHitL2 + nSameLayerHitL3 > 0 )   cout << throughPulses.size() << ", run: " << run << ", event: " << event << ", file: " << filenum << endl;
   if ( nSameLayerHitL1 + nSameLayerHitL2 + nSameLayerHitL3 > 0 )   oneSameLayerEvts++;
   if ( nSameLayerHitL1 + nSameLayerHitL2 + nSameLayerHitL3 > 1 )   twoSameLayerEvts++;
   if ( nSameLayerHitL1 + nSameLayerHitL2 + nSameLayerHitL3 > 2 ) threeSameLayerEvts++;


   if ( muonBarHit ){
	   for (unsigned int i=0; i<L1Pulses.size(); i++){
		h_nPE->Fill(  L1Pulses.at(i).at(1).at(3) );
		h_tDiffAll->Fill( L1Pulses.at(i).at(0).at(4) - L1Pulses.at(i).at(1).at(4) );
//		profile->Fill( L1Pulses.at(i).at(1).at(3), L1Pulses.at(i).at(0).at(4) - L1Pulses.at(i).at(1).at(4)  );
		h_tDiffSlab->Fill( SlabPulses.at(1).at(4) - L1Pulses.at(i).at(1).at(4) );
		profile->Fill( L1Pulses.at(i).at(1).at(3), SlabPulses.at(1).at(4) - L1Pulses.at(i).at(1).at(4) );
	   }
	   h_MuHitsSL->Fill( L1Pulses.size() );
	   for (unsigned int i=0; i<L2Pulses.size(); i++){
		h_nPE->Fill(  L2Pulses.at(i).at(1).at(3) );
		h_tDiffAll->Fill( L2Pulses.at(i).at(0).at(4) - L2Pulses.at(i).at(1).at(4) );
//		profile->Fill(  L2Pulses.at(i).at(1).at(3), L2Pulses.at(i).at(0).at(4) - L2Pulses.at(i).at(1).at(4) );
		h_tDiffSlab->Fill( SlabPulses.at(3).at(4) - L2Pulses.at(i).at(1).at(4) );
		profile->Fill( L2Pulses.at(i).at(1).at(3), SlabPulses.at(3).at(4) - L2Pulses.at(i).at(1).at(4) );
	   }
	   h_MuHitsSL->Fill( L2Pulses.size() );
	   for (unsigned int i=0; i<L3Pulses.size(); i++){
		h_nPE->Fill(  L3Pulses.at(i).at(1).at(3) );
		h_tDiffAll->Fill( L3Pulses.at(i).at(0).at(4) - L3Pulses.at(i).at(1).at(4) );
//		profile->Fill(  L3Pulses.at(i).at(1).at(3), L3Pulses.at(i).at(0).at(4) - L3Pulses.at(i).at(1).at(4) );
		h_tDiffSlab->Fill( SlabPulses.at(2).at(4) - L3Pulses.at(i).at(1).at(4) );
		profile->Fill( L3Pulses.at(i).at(1).at(3), SlabPulses.at(2).at(4) - L3Pulses.at(i).at(1).at(4) );
	    }
	    h_MuHitsSL->Fill( L3Pulses.size() );
	
//		check hits across different layers
	   if ( nSameLayerHitL1 == 1 && nSameLayerHitL2 == 1 && nSameLayerHitL3 == 0){
		   for (unsigned int i=0; i<L1Pulses.size(); i++){
			for (unsigned int j=0; j<L2Pulses.size(); j++){
				h_dTlayer->Fill( L2Pulses.at(j).at(1).at(4) - L1Pulses.at(i).at(1).at(4) );
				h_dTlayerL->Fill( L2Pulses.at(j).at(0).at(4) - L1Pulses.at(i).at(0).at(4) );
				if ( L2Pulses.at(j).at(1).at(3) > L1Pulses.at(i).at(1).at(3) ) h_tDiff->Fill(L2Pulses.at(j).at(1).at(3),L2Pulses.at(j).at(1).at(4)-L1Pulses.at(i).at(1).at(4)  );
				if ( L2Pulses.at(j).at(1).at(3) < L1Pulses.at(i).at(1).at(3) ) h_tDiff->Fill(L1Pulses.at(i).at(1).at(3),L2Pulses.at(j).at(1).at(4)-L1Pulses.at(i).at(1).at(4)  );
			}
		   }
	   }
	   if ( nSameLayerHitL1 == 1 && nSameLayerHitL3 == 1 ){
		   for (unsigned int i=0; i<L1Pulses.size(); i++){
			for (unsigned int j=0; j<L3Pulses.size(); j++){
				h_dTlayer->Fill( L3Pulses.at(j).at(1).at(4) - L1Pulses.at(i).at(1).at(4) );
				h_dTlayerL->Fill( L3Pulses.at(j).at(0).at(4) - L1Pulses.at(i).at(0).at(4) );
				if ( L3Pulses.at(j).at(1).at(3) > L1Pulses.at(i).at(1).at(3) ) h_tDiff->Fill(L3Pulses.at(j).at(1).at(3),L3Pulses.at(j).at(1).at(4)-L1Pulses.at(i).at(1).at(4)  );
				if ( L3Pulses.at(j).at(1).at(3) < L1Pulses.at(i).at(1).at(3) ) h_tDiff->Fill(L1Pulses.at(i).at(1).at(3),L3Pulses.at(j).at(1).at(4)-L1Pulses.at(i).at(1).at(4)  );
			}
		   }
	   }
	   if ( nSameLayerHitL1 == 0 && nSameLayerHitL2 == 1 && nSameLayerHitL3 == 1 ){
		   for (unsigned int i=0; i<L2Pulses.size(); i++){
			for (unsigned int j=0; j<L3Pulses.size(); j++){
				h_dTlayer->Fill( L3Pulses.at(j).at(1).at(4) - L2Pulses.at(i).at(1).at(4) );
				h_dTlayerL->Fill( L3Pulses.at(j).at(0).at(4) - L2Pulses.at(i).at(0).at(4) );
				if ( L3Pulses.at(j).at(1).at(3) > L2Pulses.at(i).at(1).at(3) ) h_tDiff->Fill(L3Pulses.at(j).at(1).at(3),L3Pulses.at(j).at(1).at(4)-L2Pulses.at(i).at(1).at(4)  );
				if ( L3Pulses.at(j).at(1).at(3) < L2Pulses.at(i).at(1).at(3) ) h_tDiff->Fill(L2Pulses.at(i).at(1).at(3),L3Pulses.at(j).at(1).at(4)-L2Pulses.at(i).at(1).at(4)  );
			}
		   }
	    }
   }

   if ( !muonBarHit ){
	   for (unsigned int i=0; i<L1PulsesNoHit.size(); i++){ 
			h_nPE->Fill(   L1PulsesNoHit.at(i).at(3) );
			h_tDiffSlab->Fill( SlabPulses.at(1).at(4) - L1PulsesNoHit.at(i).at(4) );
			profile->Fill( L1PulsesNoHit.at(i).at(3) , SlabPulses.at(1).at(4) - L1PulsesNoHit.at(i).at(4) );
			h_tDiff->Fill( L1PulsesNoHit.at(i).at(3) , SlabPulses.at(1).at(4) - L1PulsesNoHit.at(i).at(4) );
	   }
	   h_MuHitsSL->Fill( L1PulsesNoHit.size() );
	   for (unsigned int i=0; i<L2PulsesNoHit.size(); i++){ 
			h_nPE->Fill(  L2PulsesNoHit.at(i).at(3) );
			h_tDiffSlab->Fill( SlabPulses.at(3).at(4) - L2PulsesNoHit.at(i).at(4) );
			profile->Fill( L2PulsesNoHit.at(i).at(3) , SlabPulses.at(3).at(4) - L2PulsesNoHit.at(i).at(4) );
			h_tDiff->Fill( L2PulsesNoHit.at(i).at(3) , SlabPulses.at(3).at(4) - L2PulsesNoHit.at(i).at(4) );
	   }
	   h_MuHitsSL->Fill( L2PulsesNoHit.size() );
	   for (unsigned int i=0; i<L3PulsesNoHit.size(); i++){ 
			h_nPE->Fill(   L3PulsesNoHit.at(i).at(3) );
			h_tDiffSlab->Fill( SlabPulses.at(2).at(4) - L3PulsesNoHit.at(i).at(4) );
			profile->Fill( L3PulsesNoHit.at(i).at(3) , SlabPulses.at(2).at(4) - L3PulsesNoHit.at(i).at(4) );
			h_tDiff->Fill( L3PulsesNoHit.at(i).at(3) , SlabPulses.at(2).at(4) - L3PulsesNoHit.at(i).at(4) );
	   }
	   h_MuHitsSL->Fill( L3PulsesNoHit.size() );
//	   for (unsigned int i=0; i<L1Pulses.size(); i++){ 
//			h_dTlayer->Fill( L1Pulses.at(i).at(1).at(4) - L1Pulses.at(i).at(0).at(4) );
//	   }
//	   for (unsigned int i=0; i<L2Pulses.size(); i++){ 
//			h_dTlayer->Fill( L2Pulses.at(i).at(1).at(4) - L2Pulses.at(i).at(2).at(4) );
//	   }
//	   for (unsigned int i=0; i<L3Pulses.size(); i++){ 
//			h_dTlayer->Fill( L3Pulses.at(i).at(0).at(4) - L3Pulses.at(i).at(1).at(4) );
//	   }
   }

//	end of main loop
  }


  if ( muonBarHit ) cout << " probability of 3 bar hit in straight line: " << 1.5*nBarEvts/nSlabEvts << endl;
  if (!muonBarHit ) cout << " probability of no Muon hit in bars: " << 1 - nBarEvts/nSlabEvts << endl;

  outFile->cd();
  profile->Write();
  h_nPE->Write();
  h_tDiffAll->Write();
  h_tDiffSlab->Write();
  h_dTlayer->Write();
  h_dTlayerL->Write();
  h_MuHitsSL->Write();
  h_tDiff->Write();

  outFile->Write();
  outFile->Close();

}
