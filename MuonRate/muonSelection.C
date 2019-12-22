// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/muonAnalysis.cc"
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/milliHists.h"
#include "TF1.h"


void muonSelection( TChain *chain, TString output_filename, TString runName, vector<float> channels, bool synch = false) 
{

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

  TH1D *h_Time[32], *h_CalTime[32], *h_Height[32], *h_nPE[32], *h_TimeDiff1, *h_TimeDiff2; 
  for (unsigned int c=0; c<nChannels; c++) {
   h_Time[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_Time","Time of pulses in Channel "+ChStr[c],35, tMin, tMax); 
   h_CalTime[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_CalTime","Calibrated time of pulses in Channel "+ChStr[c],35, tMin, tMax); 
   h_Height[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_Height","Height of pulses in Channel "+ChStr[c], 43, 0, 1200); 
   h_nPE[c] = new TH1D("Run"+runName+"_CH"+ChStr[c]+"_nPE","# PEs in Channel "+ChStr[c], 39, nPEmin-50, 2000 ); 
  }
  h_TimeDiff1 = new TH1D("Run"+runName+"_TimeDiff1","#Delta_{t} calibrated Ch21-18" ,30, -45, 45); 
//  h_TimeDiff2 = new TH1D("Run"+runName+"_TimeDiff2","Time diff of pulses in Channel 28-18" ,30, -45, 45); 

  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();
  float nSelEvts = 0;

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

  //Main Event loop
  for(int ia = 0; ia<nentries; ia++){

    vector<vector<float>> Pulses;  	// store data for target channel and treat it as TargetPulses pulse
    vector<int> trigChStatus(trigCh.size(),0);  // set no pulses in tag channel at the beginning of every event
    bool synchBoards = false;			// only consider events with synched boards

//	skip non-synchronised events
    chain->GetEntry(ia);
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true; 
    if (!synchBoards && synch) continue;

//	apply pre-selection requirements (check muonAnalysis.cc file)
    vector<vector<vector<float>>> totPulses, L1Pulses, L2Pulses, L3Pulses ;
    totPulses = pulsePreSelection( chan,height,ptime,nPE,time_module_calibrated, {trigCh} );
    Pulses    = totPulses.at(0);

//	apply selection requirements
    vector<vector<float>> selectedPulses;
    selectedPulses = pulseSelection( Pulses, trigCh.size() );
    if ( selectedPulses.size() != trigCh.size() ) continue;


//	order pulses by channel
    std::sort(selectedPulses.begin(), selectedPulses.end(), sort_wrt_channel);
//	Fill in the right histograms with earliest occurance of pulses
    for (unsigned int i=0; i<selectedPulses.size(); i++){
	h_Height[(int)selectedPulses.at(i)[0]]->Fill( selectedPulses.at(i)[1],1./runTime);
	h_Time[(int)selectedPulses.at(i)[0]]->Fill( selectedPulses.at(i)[2],1./runTime);
	h_nPE[(int)selectedPulses.at(i)[0]]->Fill( selectedPulses.at(i)[3],1./runTime);
	h_CalTime[(int)selectedPulses.at(i)[0]]->Fill( selectedPulses.at(i)[4],1./runTime);
    }
    h_TimeDiff1->Fill( selectedPulses.at(2)[4] - selectedPulses.at(0)[4], 1./runTime );
//    h_TimeDiff2->Fill( selectedPulses.at(2)[4] - selectedPulses.at(0)[4] , 1./runTime );
   }


  TCanvas *c_time = new TCanvas("canvas_time","Run"+runName+"_time", 900,700);
  TCanvas *c_height = new TCanvas("canvas_height","Run"+runName+"_height", 900,700);
  TCanvas *c_nPE = new TCanvas("canvas_nPE","Run"+runName+"_nPE", 900,700);
  TCanvas *c = new TCanvas("canvas","Run"+runName, 900,700);
  TCanvas *c1 = new TCanvas("canvas1","Run1"+runName, 900,700);

  TF1 *Gauss1 = new TF1("Gauss1","gaus(0)",-12,12);
  Gauss1->SetParameters(1e-3,0,5);
  TF1 *Gauss2 = new TF1("Gauss2","gaus(3)",-30,-12);
  Gauss2->SetParameters(1e-4,-25,5);
  TF1 *total = new TF1 ("total", "gaus(0)+gaus(3)", -30,12);
  double par[6];
  h_TimeDiff1->Fit("Gauss1","R");
  h_TimeDiff1->Fit("Gauss2","R");
  Gauss1->GetParameters(&par[0]);
  Gauss2->GetParameters(&par[3]);
  total->SetParameters(par);
  total->SetParNames("scale1","#mu_{1}","#sigma_{1}","scale2","#mu_{2}","#sigma_{2}");
  
  c->cd();
  gStyle->SetOptFit(111); 
  gStyle->SetOptStat("n"); 
  h_TimeDiff1->SetLineColor(kBlue);
  h_TimeDiff1->SetFillColor(kBlue);
  h_TimeDiff1->SetFillStyle(3002);
  h_TimeDiff1->GetXaxis()->SetTitle("#Delta_{t} [ns]");
  h_TimeDiff1->GetYaxis()->SetTitle("# of particles");
  h_TimeDiff1->Fit("total","R");
  h_TimeDiff1->Draw("HIST");
  total->Draw("SAME");
  c->Update();
  c->SaveAs("/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/Plots/Run"+runName+"_TimeDiff.png");
  c->Close();
//  c1->cd();
//  h_TimeDiff2->SetLineColor(kRed);
//  h_TimeDiff2->SetFillColor(kRed);
//  h_TimeDiff2->SetFillStyle(3003);
//  h_TimeDiff2->SetStats(0);
//  h_TimeDiff2->Draw("HIST");
//  TLegend *legend1 = new TLegend(.75, .80, .95, .95);
//  legend1->AddEntry(h_TimeDiff2, "Cal");
//  legend1->Draw();
//  c1->SaveAs("/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/Plots/Run"+runName+"_CalTimeDiff.png");
//  c1->Close();
  outFile->cd();
//  h_TimeDiff->Write();
  for (unsigned int i=0; i<nChannels; i++){
	if ( h_Time[i]->GetEntries() == 0 ) continue;
	c_time->cd();
	c_time->SetLogy();
//	h_Time[i]->Draw("HIST");
	c_height->cd();
//	h_Height[i]->Draw("HIST");
	c_nPE->cd();
//	h_nPE[i]->Draw("HIST");
	h_Time[i]->Write();
	h_CalTime[i]->Write();
	h_Height[i]->Write();
	h_nPE[i]->Write();

  }
  c_time->Close();
  c_height->Close();
  c_nPE->Close();
  outFile->Write();
  outFile->Close();

}
