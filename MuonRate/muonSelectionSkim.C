// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/muonAnalysis.cc"
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/milliHists.h"
#include "TF1.h"


void muonSelectionSkim( TChain *chain, TString output_filename, TString runName, vector<float> channels, bool synch = false ) 
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

  TH1D *h_Time[32], *h_CalTime[32], *h_Height[32], *h_nPE[32], *h_TimeDiff1, *h_TimeDiffThru, *h_nPELargest, *h_nPESmallest, *h_nPEThrough, *h_TimeDiffLayer1, *h_TimeDiffLayer2, *h_TimeDiffLayer3, *h_TimeDiffSecondPulses, *h_TimeDiffSecondPulsesL;  
  TH2D *h_2DdT;
  for (unsigned int c=0; c<nChannels; c++) {
   h_Time[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_Time","Time of pulses in Channel "+ChStr[c],35, tMin, tMax); 
   h_CalTime[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_CalTime","Calibrated time of pulses in Channel "+ChStr[c],35, tMin, tMax); 
   h_Height[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_Height","Height of pulses in Channel "+ChStr[c], 43, 0, 1200); 
   h_nPE[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_nPE","# PEs in Channel "+ChStr[c], 39, 50, 2000 ); 
  }
  h_2DdT = new TH2D ("RunSkim_2Ddt","#Delta_{t} of slabs vs bars", 15, -15, 15, 30, -50, 50);
  h_TimeDiff1 = new TH1D("Run"+runName+"_TimeDiff1","#Delta_{t} calibrated Ch21-18" ,30, -45, 45); 
  h_TimeDiffThru = new TH1D("Run"+runName+"_TimeDiffThru","#Delta_{t} calibrated Ch21-18" ,30, -65, 45); 
  h_TimeDiffSecondPulses = new TH1D("RunSkim_TimeDiffSecondPulses","#Delta_{t} cal, second pulses" ,30, -65, 45); 
  h_TimeDiffSecondPulsesL = new TH1D("RunSkim_TimeDiffSecondPulsesL","#Delta_{t} cal, second pulses" ,30, -65, 45); 
  h_TimeDiffLayer1 = new TH1D("RunSkim_TimeDiffLayer1","#Delta_{t} cal, layer1" ,30, -65, 45); 
  h_TimeDiffLayer2 = new TH1D("RunSkim_TimeDiffLayer2","#Delta_{t} cal, layer2" ,30, -65, 45); 
  h_TimeDiffLayer3 = new TH1D("RunSkim_TimeDiffLayer3","#Delta_{t} cal, layer3" ,30, -65, 45); 
  h_nPELargest = new TH1D("RunSkim_nPEL","nPE largest pulse" ,55, 0, 4e3); 
  h_nPEThrough = new TH1D("RunSkim_nPET","nPE through pulse" ,55, 0, 4e3); 
  h_nPESmallest = new TH1D("RunSkim_nPES","nPE smallest pulse" ,55, 0, 2e3); 
//  h_TimeDiff2 = new TH1D("Run"+runName+"_TimeDiff2","Time diff of pulses in Channel 28-18" ,30, -45, 45); 

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
  if (  synch ) runTime = SumSynchTime;
  double weight = 1./runTime;

//  cout << " SynchTime over nonSynchTime: " << SumSynchTime/SumTime << " , Synch Time: "<< SumSynchTime << ", and non-Synch Time: " << SumTime << " ." << endl;
//  cout << " nSynch over nonSynch: " << (100.*nSynched)/(nSynched+nUnSynched) << endl;


  //Main Event loop
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);

    vector<int> trigChStatus(trigCh.size(),0);  // set no pulses in tag channel at the beginning of every event

//	skip non-synchronised events && require beam On
    bool synchBoards = false;			// only consider events with synched boards
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true; 
    if (!synchBoards && synch) continue;
    if ( !beam ) continue;
    nSelEvts++;

//	apply pre-selection requirements (check muonAnalysis.cc file)
    vector<vector<float>> SlabPulses;
    vector<vector<vector<float>>> BarPulses, totPulses, L1Pulses, L2Pulses, L3Pulses, UpLowPulses;
    totPulses = pulseSlabSelection( chan,height,ptime,nPE,time_module_calibrated, {trigCh} );
    SlabPulses = totPulses.at(0);
    if ( SlabPulses.size() != trigCh.size() ) continue;

//	order pulses by channel
   std::sort(SlabPulses.begin(), SlabPulses.end(), sort_wrt_channel);

//	Fill in the right histograms with earliest occurance of pulses
   for (unsigned int i=0; i<SlabPulses.size(); i++){
	h_Height[(int)SlabPulses.at(i)[0]]->Fill( SlabPulses.at(i)[1], weight);
	h_Time[(int)SlabPulses.at(i)[0]]->Fill( SlabPulses.at(i)[2], weight);
	h_nPE[(int)SlabPulses.at(i)[0]]->Fill( SlabPulses.at(i)[3], weight);
	h_CalTime[(int)SlabPulses.at(i)[0]]->Fill( SlabPulses.at(i)[4], weight);
   }
   float dt_cal = SlabPulses.at(2)[4] - SlabPulses.at(0)[4];
   h_TimeDiff1->Fill( dt_cal );


// require dt_cal between -12 and +12 ns from dt_cal distribution
   if ( fabs( dt_cal ) > 12 ) continue;
    nSlabEvts++;

   L1Pulses = sameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig08,trig09,trig18,trig19} );
   L2Pulses = sameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig612,trig613,trig712,trig713} );
   L3Pulses = sameLayerSelection( chan,height,ptime,nPE,time_module_calibrated, {trig24,trig25,trig34,trig35} );
   UpLowPulses = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig8124,trig9135} );

//  require at least one through track 
   bool throughParticle = false; 
   int nBarHits=0;
   for (unsigned int i=0; i<UpLowPulses.size(); i++){
	if ( UpLowPulses.at(i).size() != 1 && nBarHits == 0){ 
		nBarHits++ ;
		throughParticle = true;
	}
   }
   nBarEvts+=nBarHits; 
   if ( !throughParticle ) continue;


   int nSameLayerHitL1 = 0;
   int nSameLayerHitL2 = 0;
   int nSameLayerHitL3 = 0;
   for (unsigned int i=0; i<L1Pulses.size(); i++){
	if ( L1Pulses.at(i).size() == 1 ) continue;
   	for (unsigned int j=0; j<L2Pulses.size(); j++){
		if ( L2Pulses.at(j).size() == 1 ) continue;
		h_TimeDiffSecondPulses->Fill(  L1Pulses.at(i).at(1)[4] - L2Pulses.at(j).at(1)[4]  );
		h_TimeDiffSecondPulsesL->Fill( L1Pulses.at(i).at(0)[4] - L2Pulses.at(j).at(0)[4] );
		nSameLayerHitL1++; 
	}
   	for (unsigned int j=0; j<L3Pulses.size(); j++){
		if ( L3Pulses.at(j).size() == 1 ) continue;
		h_TimeDiffSecondPulses->Fill(  L1Pulses.at(i).at(1)[4] - L3Pulses.at(j).at(1)[4]  );
		h_TimeDiffSecondPulsesL->Fill( L1Pulses.at(i).at(0)[4] - L3Pulses.at(j).at(0)[4] );
		nSameLayerHitL2++; 
	}
   }
   for (unsigned int i=0; i<L2Pulses.size(); i++){
	if ( L2Pulses.at(i).size() == 1 ) continue;
   	for (unsigned int j=0; j<L3Pulses.size(); j++){
		if ( L3Pulses.at(j).size() == 1 ) continue;
		h_TimeDiffSecondPulses->Fill(  L2Pulses.at(i).at(1)[4] - L3Pulses.at(j).at(1)[4] );
		h_TimeDiffSecondPulsesL->Fill( L2Pulses.at(i).at(0)[4] - L3Pulses.at(j).at(0)[4] );
		nSameLayerHitL3++; 
	}
   }
   if ( nSameLayerHitL1 + nSameLayerHitL2 + nSameLayerHitL3 > 0 ) oneSameLayerEvts++;
   if ( nSameLayerHitL1 + nSameLayerHitL2 + nSameLayerHitL3 > 1 )  twoSameLayerEvts++;
   if ( nSameLayerHitL1 + nSameLayerHitL2 + nSameLayerHitL3 > 2 ) threeSameLayerEvts++;

   for (unsigned int i=0; i<L1Pulses.size(); i++){
	if ( L1Pulses.at(i).size() != 1 ){
		h_TimeDiffLayer1->Fill( L1Pulses.at(i).at(0).at(4) - L1Pulses.at(i).at(1).at(4) );
	}
	if ( L2Pulses.at(i).size() != 1 ){
		h_TimeDiffLayer2->Fill( L2Pulses.at(i).at(0).at(4) - L2Pulses.at(i).at(1).at(4) );
	}
	if ( L3Pulses.at(i).size() != 1 ){
		h_TimeDiffLayer3->Fill( L3Pulses.at(i).at(0).at(4) - L3Pulses.at(i).at(1).at(4) );
	}
    }
   }


  cout << " probability of 3 bar hit in straight line: " << 1.5*nBarEvts/nSlabEvts << endl;
  cout << " ratio of at least 1 hit in same layer over rate of 3 bar hit in straight line: " << oneSameLayerEvts/nBarEvts << endl;
  cout << " ratio of at least 2 hit in same layer over rate of 3 bar hit in straight line: " << twoSameLayerEvts/nBarEvts << endl;
  cout << " ratio of at least 3 hit in same layer over rate of 3 bar hit in straight line: " << threeSameLayerEvts/nBarEvts << endl;
  cout << " probability of multiple hits in at least one layer: " << 1.5*oneSameLayerEvts/nSlabEvts << endl;




// Plot section
//
//
//
  TCanvas *c = new TCanvas("canvas","Run"+runName, 900,700);
  TCanvas *c1 = new TCanvas("canvas1","Run1"+runName, 900,700);
  c1->cd();
  h_TimeDiffThru->SetLineColor(kRed);
  h_TimeDiffThru->SetFillColor(kRed);
  h_TimeDiffThru->SetFillStyle(3001);
//  h_TimeDiffThru->Draw();
//  h_2DdT->SetStats(0);
//  h_2DdT->GetXaxis()->SetTitle("#Delta_{t} slabs [ns]");
//  h_2DdT->GetYaxis()->SetTitle("#Delta_{t} bars [ns]");
//  h_2DdT->Draw("colz");
//  c1->SaveAs("/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/Plots/RunSkim_2DTimeDiff.png");
  c1->Close();

//  c->cd();
//  TF1 *Gauss1 = new TF1("Gauss1","gaus(0)",-12,12);
//  Gauss1->SetParameters(2e3,0,5);
//  TF1 *Gauss2 = new TF1("Gauss2","gaus(3)",-30,-12);
//  Gauss2->SetParameters(150,-25,5);
//  TF1 *total = new TF1 ("total", "gaus(0)+gaus(3)", -30,12);
//  double par[6];
//  h_TimeDiff1->Fit("Gauss1","R");
//  h_TimeDiff1->Fit("Gauss2","R");
//  Gauss1->GetParameters(&par[0]);
//  Gauss2->GetParameters(&par[3]);
//  total->SetParameters(par);
//  total->SetParNames("#alpha_{1}","#mu_{1}","#sigma_{1}","#alpha_{2}","#mu_{2}","#sigma_{2}");
// get # entries between bin 11 and 19
//  int nEntries=0;
//  for (unsigned int i=11; i<20; i++){
//	nEntries+=h_TimeDiff1->GetBinContent(i);
//  }
//  cout << "# entries between been 11 and 19 is: " << nEntries << endl;
//  c->cd();
//  gStyle->SetOptFit(111); 
//  gStyle->SetOptStat("n"); 
//  h_TimeDiff1->SetLineColor(kBlue);
//  h_TimeDiff1->SetFillColor(kBlue);
//  h_TimeDiff1->SetFillStyle(3002);
//  h_TimeDiff1->GetXaxis()->SetTitle("#Delta_{t} [ns]");
//  h_TimeDiff1->GetYaxis()->SetTitle("# of particles");
//  h_TimeDiff1->Fit("total","R");
//  h_TimeDiff1->Draw("HIST");
//  total->SetLineColor(kRed);
//  total->Draw("SAME");
//  c->Update();
//  c->SaveAs("/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/Plots/Run"+runName+"_TimeDiff.png");
//  c->Close();

  outFile->cd();
  for (unsigned int i=0; i<nChannels; i++){
	if ( h_Time[i]->GetEntries() == 0 ) continue;
	h_Time[i]->Write();
	h_CalTime[i]->Write();
	h_Height[i]->Write();
	h_nPE[i]->Write();

  }
  h_nPEThrough->Write();
  h_nPELargest->Write();
  h_nPESmallest->Write();
  h_2DdT->Write();
  h_TimeDiffSecondPulses->Write();
  h_TimeDiffSecondPulsesL->Write();
  h_TimeDiffLayer1->Write();
  h_TimeDiffLayer2->Write();
  h_TimeDiffLayer3->Write();
  outFile->Write();
  outFile->Close();

}
