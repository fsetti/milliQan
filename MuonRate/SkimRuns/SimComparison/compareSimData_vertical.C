// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/muonAnalysis.cc"
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/milliHists.h"


void compareSimData( TChain *chain, TString output_filename, TString runName, vector<float> channels, bool synch = true ) 
{

  float TTT = 0;
  float MMM = 0;
  float LLL = 0;
  float TNnMn = 0;
  float MNnL = 0;
  float TnNmnNl = 0;
  float TmnNlN = 0;
  float LNnNm = 0;
  float MNnT = 0;
  float NlNmnTn = 0;
  float NmlTnN = 0; 


  TH1D *hist = new TH1D ("hist","",11,0,11);

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


    bool T1 = false;
    bool T2 = false;
    bool T3 = false;
    bool M1 = false;
    bool M2 = false;
    bool M3 = false;
    bool L1 = false;
    bool L2 = false;
    bool L3 = false;

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
//	Bar selection II, study vertical angle distributions 
    vector<vector<float>> T1Pulses, T2Pulses, T3Pulses, M1Pulses, M2Pulses, M3Pulses, L1Pulses, L2Pulses, L3Pulses;
    T1Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {0,1} );
    T2Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {6,7} );
    T3Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {2,3} );
    M1Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {24,25} );
    M2Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {16,17} );
    M3Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {22,23} );
    L1Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {8,9} );
    L2Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {12,13} );
    L3Pulses = pulseBarSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {4,5} );

    if ( T1Pulses.size() != 0 ) T1 = true;
    if ( T2Pulses.size() != 0 ) T2 = true;
    if ( T3Pulses.size() != 0 ) T3 = true;
    if ( M1Pulses.size() != 0 ) M1 = true;
    if ( M2Pulses.size() != 0 ) M2 = true;
    if ( M3Pulses.size() != 0 ) M3 = true;
    if ( L1Pulses.size() != 0 ) L1 = true;
    if ( L2Pulses.size() != 0 ) L2 = true;
    if ( L3Pulses.size() != 0 ) L3 = true;

//    if ( T1 && T2 && T3 ) TTT++;
//    if ( M1 && M2 && M3 ) MMM++;
//    if ( L1 && L2 && L3 ) LLL++;
//    if ( T1 && M3 ) TNnMn++;
//    if ( L1 && M3 ) LNnNm++;
//    if ( M1 && L3 ) MNnL++;
//    if ( M1 && T3 ) MNnT++;
//    if ( T1 && M2 && L3 ) TnNmnNl++;
//    if ( L1 && M2 && T3 ) NlNmnTn++;
//    if ( T1 && M1 && L2 ) TmnNlN++;
//    if ( M1 && L1 && T2 ) NmlTnN++;


// set up vetoes boolean
    bool vT123 = false; if ( !T1 && !T2 && !T3 ) vT123 = true;
    bool vM123 = false; if ( !M1 && !M2 && !M3 ) vM123 = true;
    bool vL123 = false; if ( !L1 && !L2 && !L3 ) vL123 = true;
    bool vT23 = false;  if ( !T2 && !T3 ) vT23 = true;
    bool vL23 = false;  if ( !L2 && !L3 ) vL23 = true;
    bool vTM1 = false;  if ( !T1 && !M1 ) vTM1 = true;
    bool vML1 = false;  if ( !M1 && !L1 ) vML1 = true; 
    bool vTM3 = false;  if ( !T3 && !M3 ) vTM3 = true;
    bool vML3 = false;  if ( !M3 && !L3 ) vML3 = true;
    bool vTL1 = false;  if ( !T1 && !L1 ) vTL1 = true;
    bool vTL3 = false;  if ( !T3 && !L3 ) vTL3 = true;



    if ( T1 && M1 && L2 && !T2 && vTM3 ) TmnNlN++;
    if ( T1 && M2 && L3 && !L1 && !T3 ) TnNmnNl++;
    if ( M1 && L3 && vT123 && !L1 && !M3 ) MNnL++;
    if ( T1 && M3 && vML1 && !L2 && !T3 ) TNnMn++;
    if ( L1 && L2 && L3 && vT123 && vM123 ) LLL++;
    if ( M1 && M2 && M3 && vT123 && vL123 ) MMM++;
    if ( T1 && T2 && T3 && vM123 && vL123 ) TTT++;
    if ( L1 && M3 && vTM1 && !T2 && !L3 ) LNnNm++;
    if ( M1 && T3 && vL123 && !T1 && !M2 ) MNnT++;
    if ( L1 && M2 && T3 && !T1 && !L3 ) NlNmnTn++;
    if ( M1 && L1 && T2 && !L2 && vML3 ) NmlTnN++;

  }
  float scale = 0.26/0.179;
//  float scale = 1.;

  cout << "Total entries: " << TTT+MMM+LLL+TNnMn+LNnNm+MNnL+MNnT+TnNmnNl+NlNmnTn+TmnNlN+NmlTnN << endl;

  TmnNlN*=scale;
  TnNmnNl*=scale;
  MNnL*=scale;
  TNnMn*=scale;
  LLL*=scale;
  MMM*=scale;
  TTT*=scale;
  LNnNm*=scale;
  MNnT*=scale;
  NlNmnTn*=scale;
  NmlTnN*=scale;

  hist->SetBinContent( 1 , TmnNlN );
  hist->SetBinContent( 2 , TnNmnNl );
  hist->SetBinContent( 3 , MNnL );
  hist->SetBinContent( 4 , TNnMn );
  hist->SetBinContent( 5 , LLL );
  hist->SetBinContent( 6 , MMM );
  hist->SetBinContent( 7 , TTT );
  hist->SetBinContent( 8 , LNnNm );
  hist->SetBinContent( 9 , MNnT );
  hist->SetBinContent( 10 , NlNmnTn );
  hist->SetBinContent( 11 , NmlTnN );

  TCanvas *c = new TCanvas ("canvas","",800,600);
  c->cd();
  hist->Draw("HIST,ERR");
  c->SaveAs("/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/Plots/SimData_vert.png");


  cout << "\n 1: " << TmnNlN << "\n 2: " << TnNmnNl << "\n 3: " << MNnL << "\n 4: " << TNnMn << "\n 5: " << LLL << "\n 6: " << MMM << "\n 7: " << TTT << "\n 8: " << LNnNm << "\n 9: " << MNnT << "\n 10: " << NlNmnTn << "\n 11: " << NmlTnN << endl;
}
