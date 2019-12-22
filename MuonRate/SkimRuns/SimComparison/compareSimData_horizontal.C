// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/muonAnalysis.cc"
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/milliHists.h"


void compareSimData( TChain *chain, TString output_filename, TString runName, vector<float> channels, bool synch = true ) 
{

  TH1D *hist = new TH1D ("hist","",8,0,8);

  int LrL = 0;
  int RLn = 0;
  int RNnL = 0;
  int LLL = 0;
  int RRR = 0;
  int LNnR = 0;
  int LRn = 0;
  int LrR = 0;

  vector<int> trigCh;
  for (unsigned int i=0; i<channels.size(); i++){
        trigCh.push_back((int)channels[i]);
  }

  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

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

    bool L1 = false;
    bool L2 = false;
    bool L3 = false;
    bool R1 = false;
    bool R2 = false;
    bool R3 = false;


    chain->GetEntry(ia);

//	skip non-synchronised events && require beam On
    bool synchBoards = false;			// only consider events with synched boards
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true; 
    if (!synchBoards && synch) continue;
    if ( !beam ) continue;


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


/////////////////////////////////////////////////////////////////////////////
//	Layer selection
//	Bar selection I, require through going particles in upper, lower channels 
    vector<vector<float>> L1Pulses, L2Pulses, L3Pulses, R1Pulses, R2Pulses, R3Pulses;
    L1Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {0,24,8} );
    L2Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {6,16,12} );
    L3Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {2,22,4} );
    R1Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {1,25,9} );
    R2Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {7,17,13} );
    R3Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {3,23,5} );

    if ( L1Pulses.size() != 0 ) L1 = true;
    if ( L2Pulses.size() != 0 ) L2 = true;
    if ( L3Pulses.size() != 0 ) L3 = true;
    if ( R1Pulses.size() != 0 ) R1 = true;
    if ( R2Pulses.size() != 0 ) R2 = true;
    if ( R3Pulses.size() != 0 ) R3 = true;

    if ( L1 && L2 && !L3 && R1 && !R2 && !R3 ) LrL++;
    if ( L1 && !L2 && !L3 && R1 && R2 && !R3 ) LrR++;
    if ( !L1 && L2 && !L3 && R1 && !R3 ) RLn++;
    if ( L1 && !L3 && !R1 && R2 && !R3 ) LRn++;
    if ( !L1 && L3 && R1 && !R3 ) RNnL++;
    if ( L1 && !L3 && !R1 && R3 ) LNnR++;
    if ( L1 && L2 && L3 && !R1 && !R2 && !R3 ) LLL++;
    if ( !L1 && !L2 && !L3 && R1 && R2 && R3 ) RRR++;


//    if ( L1 && L2 && R1 ) LrL++;
//    if ( L1 && R1 && R2 ) LrR++;
//    if ( L2 && R1 ) RLn++;
//    if ( L1 && R2 ) LRn++;
//    if ( L3 && R1 ) RNnL++;
//    if ( L1 && R3 ) LNnR++;
//    if ( L1 && L2 && L3 ) LLL++;
//    if ( R1 && R2 && R3 ) RRR++;



  }
  float scale = 0.26/0.179;
//  float scale = 1.;


  cout << "Total entries: " << LLL+RRR+LrL+LrR+RLn+LRn+RNnL+LNnR << endl;


  LrL*=scale;
  RLn*=scale;
  RNnL*=scale;
  LLL*=scale;
  RRR*=scale;
  LNnR*=scale;
  LRn*=scale;
  LrR*=scale;

  hist->SetBinContent( 1 , LrL  );
  hist->SetBinContent( 2 , RLn  );
  hist->SetBinContent( 3 , RNnL );
  hist->SetBinContent( 4 , LLL  );
  hist->SetBinContent( 5 , RRR  );
  hist->SetBinContent( 6 , LNnR );
  hist->SetBinContent( 7 , LRn  );
  hist->SetBinContent( 8 , LrR  );

  TCanvas *c = new TCanvas ("canvas","",800,600);
  c->cd();
  hist->Draw("HIST,ERR");
  c->SaveAs("/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/Plots/SimData_hor.png");


  cout << "\n LrL: " << LrL << "\n RLn: " << RLn << "\n RNnL: " << RNnL << "\n LLL: " << LLL << "\n RRR: " << RRR << "\n LNnR: " << LNnR << "\n LRn: " << LRn << "\n LrR: " << LrR << endl;

}
