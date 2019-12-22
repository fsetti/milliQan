// Make milliqan validation plots. D. Stuart, Oct. 2017.

#include "milliHists.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "TString.h"
#include "TVector3.h"
using namespace std;


void tree1r(TChain *chain, TString output_filename,TString EventCategory,TString LHCStatus,TString RangeCode, int RunNum)
{
  int debug=0; // Debugging level. 0=none, 100=trace, 300=detailed trace

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";

  // Read the fill information.
  vector<int> fillNum;
  vector<Long64_t> fillStartTime; // in unix epoch time
  vector<Long64_t> fillEndTime; // in unix epoch time
  vector<float> fillLumi; // in pb^-1
  FILE *fillFile;
  fillFile = fopen("/net/cms26/cms26r0/milliqan/milliqanOffline/processedFillList2018.txt","r");
  if (!fillFile) {
    cerr << "Could not open fill file.\n";
  }
  else { // Opened file successfully?  
    int ncol = 0;
    while (ncol >= 0) {
      int intVal; Long64_t dVal; float fVal;
      ncol=fscanf(fillFile,"%d", &intVal);
      fillNum.push_back(intVal);
      ncol=fscanf(fillFile,"%li", &dVal);
      fillStartTime.push_back(dVal);
      ncol=fscanf(fillFile,"%li", &dVal);
      fillEndTime.push_back(dVal);
      ncol=fscanf(fillFile,"%f", &fVal);
      fillLumi.push_back(fVal);
      //cout << "Read fill "<<intVal<<"\n";
      char c;
      c = '0'; while (c != '\n' && ncol>0) ncol=fscanf(fillFile,"%c",&c);
    }
    fclose(fillFile);
  } // Opened file successfully?

  // Record the number of events with specific hit patterns. 
  // These are used to validate the trigger efficiency.
  // In particular, I count events with vertical cosmic like hits, two layer hits, and three layer hits.
  // These are calculated and reported for each event category, which corresponds to different channel thresholds.
  // First define the digitizer 0 sets, 2 hits and then 3 hits
  int nCH0CH8=0; // Within a layer
  int nCH1CH9=0; 
  int nCH0CH9=0;
  int nCH1CH8=0;
  int nCH6CH12=0; 
  int nCH7CH13=0; 
  int nCH6CH13=0; 
  int nCH7CH12=0;
  int nCH2CH4=0; 
  int nCH3CH5=0; 
  int nCH2CH5=0; 
  int nCH3CH4=0; 
  int nCH0CH6=0;   // Between layers
  int nCH1CH7=0;
  int nCH0CH7=0;
  int nCH1CH6=0;
  int nCH8CH12=0;
  int nCH9CH13=0;
  int nCH8CH13=0;
  int nCH9CH12=0;
  int nCH0CH12=0;
  int nCH1CH13=0;
  int nCH0CH13=0;
  int nCH1CH12=0;
  int nCH8CH6=0;
  int nCH9CH7=0;
  int nCH8CH7=0;
  int nCH9CH6=0;
  int nCH6CH2=0;
  int nCH7CH3=0;
  int nCH6CH3=0;
  int nCH7CH2=0;
  int nCH12CH4=0;
  int nCH13CH5=0;
  int nCH12CH5=0;
  int nCH13CH4=0;
  int nCH6CH4=0;
  int nCH7CH5=0;
  int nCH6CH5=0;
  int nCH7CH4=0;
  int nCH12CH2=0;
  int nCH13CH3=0;
  int nCH12CH3=0;
  int nCH13CH2=0;
  int nCH0CH8CH10=0; // 3 hits per layer
  int nCH1CH9CH10=0;
  int nCH0CH9CH10=0;
  int nCH1CH8CH10=0;
  int nCH6CH12CH30=0;
  int nCH7CH13CH30=0;
  int nCH2CH4CH14=0;
  int nCH3CH5CH14=0;
  int nCH2CH5CH14=0;
  int nCH3CH4CH14=0;
  int nCH0CH6CH2=0; // Hits across 3 layers
  int nCH1CH7CH3=0;
  int nCH0CH6CH3=0;
  int nCH1CH7CH2=0;
  int nCH0CH7CH2=0;
  int nCH1CH6CH3=0;
  int nCH8CH12CH4=0;
  int nCH9CH13CH5=0;
  int nCH8CH12CH5=0;
  int nCH9CH13CH4=0;
  int nCH8CH13CH4=0;
  int nCH9CH12CH5=0;
  int nCH0CH6CH4=0;
  int nCH0CH6CH5=0;
  int nCH1CH7CH5=0;
  int nCH1CH7CH4=0;
  // Now define the digitizer 1 sets, 2 hits and then 3 hits
  int nCH11CH8=0; // Two hits in digitizer 1...
  int nCH11CH9=0; 
  int nCH13CH8=0; 
  int nCH13CH9=0;
  int nCH16CH30=0;
  int nCH17CH30=0;
  int nCH16CH19=0;
  int nCH17CH19=0;
  int nCH22CH15=0;
  int nCH23CH15=0;
  int nCH22CH10=0;
  int nCH23CH10=0;
  int nCH27CH24CH29=0; // Three hits within a layer
  int nCH27CH25CH29=0;
  int nCH30CH16CH19=0;
  int nCH30CH17CH19=0;
  int nCH31CH22CH26=0;
  int nCH31CH23CH26=0;
  int nCH24CH16CH22=0; // Hits across 3 layers
  int nCH24CH16CH23=0;
  int nCH24CH17CH23=0;
  int nCH24CH17CH22=0;
  int nCH25CH17CH23=0;
  int nCH25CH16CH23=0;
  int nCH25CH16CH22=0;
  int nCH18CH20=0; // Two hits between slabs
  int nCH20CH28=0;
  int nCH28CH21=0;
  int nCH18CH28=0;
  int nCH20CH21=0;
  int nCH18CH21=0;
  int nCH18CH20CH28=0; // Hits across slab layers
  int nCH20CH28CH21=0;
  int nCH18CH20CH28CH21=0; 
  int nCH18CH20CH28CH21CH0CH6CH2=0; // Thru-going
  int nCH18CH20CH28CH21CH1CH7CH3=0;
  int nCH18CH20CH28CH21CH24CH16CH22=0;
  int nCH18CH20CH28CH21CH25CH17CH23=0;
  int nCH18CH20CH28CH21CH12CH4=0;
  int nCH18CH20CH28CH21CH9CH13CH5=0;
  int nAllSlabsAllLayers=0; // All Slabs and a hit somewhere in each layer
  int nAllSlabsAnyLayer=0; // All Slabs and a hit somewhere in any layer
  int nAllLayersAnySlab=0; // All layers with at least one hit, and at least one slab hit

  // Set voltage and area thresholds for pulses to count as hits.
  // These thresholds differ across channels and differ by event category (0=All, 1=SPE, 2=Small, 3=Cosmic, 4=Thru).
  float VThresholds[32][5]; // Threshold for pulse height for each channel and event category (All,SPE,Small,Cosmic,Thru).
  float DThresholds[32][5]; // Threshold for pulse duration for each channel and event category (All,SPE,Small,Cosmic,Thru).
  float AThresholds[32][5]; // Threshold for pulse area for each channel and event category (All,SPE,Small,Cosmic,Thru).

  // Set the default thresholds
  for (int c=0; c<32; c++) {
    VThresholds[c][0] = 0.; // All, i.e., no cuts beyond pulse finding
    VThresholds[c][1] = 10.; // SPE, upper cut [mV]
    VThresholds[c][2] = 20.; // Small signals, lower cut [mV]
    VThresholds[c][3] = 300.; // Cosmics, lower cut [mV]
    VThresholds[c][4] = 10.; // Any pulse
    DThresholds[c][0] = 0.; // All, i.e., no cuts beyond pulse finding
    DThresholds[c][1] = 15.; // SPE, upper cut [ns]
    DThresholds[c][2] = 5.; // Small signals, lower cut [ns]
    DThresholds[c][3] = 10.; // Cosmics, lower cut [ns]
    DThresholds[c][4] = 10.; // Through going particles, lower cut [ns]
    AThresholds[c][0] = 0.; // All, i.e., no cuts beyond pulse finding
    AThresholds[c][1] = 0.02; // SPE, upper cut [nVs]
    AThresholds[c][2] = 0.1; // Small signals, lower cut [nVs]
    AThresholds[c][3] = 1.; // Cosmics, lower cut [nVs]
    AThresholds[c][4] = 0.01; // Through going particles, lower cut [nVs]
  }
  
  // Special settings
  if (RunNum<31) {
    for (int c=0; c<32; c++) {
      VThresholds[c][3] = 60.; // Cosmics, lower cut [mV]
      VThresholds[c][4] = 1000.; // Through going particles, lower cut [mV]
      AThresholds[c][3] = 5.; // Cosmics, lower cut [nVs]
      AThresholds[c][4] = 20.; // Through going particles, lower cut [nVs]
    }
    VThresholds[8][3] = 10.; // Cosmics, lower cut [mV]
    VThresholds[8][4] = 100.; // Through going particles, lower cut [mV]
    AThresholds[8][3] = 0.4; // Cosmics, lower cut [nVs]
    AThresholds[8][4] = 2.; // Through going particles, lower cut [nVs]
    VThresholds[9][3] = 10.; // Cosmics, lower cut [mV]
    VThresholds[9][4] = 100.; // Through going particles, lower cut [mV]
    AThresholds[9][3] = 1.; // Cosmics, lower cut [nVs]
    AThresholds[9][4] = 2.; // Through going particles, lower cut [nVs]
    VThresholds[10][3] = 10.; // Cosmics, lower cut [mV]
    VThresholds[10][4] = 100.; // Through going particles, lower cut [mV]
    AThresholds[10][3] = 0.4; // Cosmics, lower cut [nVs]
    AThresholds[10][4] = 2.; // Through going particles, lower cut [nVs]
  }

  // Hardcode run dependent cosmic thresholds. These should eventually be read from a file.
  if (RunNum == 24) {
    AThresholds[1][3] = 5.;
    AThresholds[2][3] = 1.5;
    AThresholds[3][3] = 1.;
    AThresholds[4][3] = 2.;
    AThresholds[5][3] = 2.5;
    AThresholds[6][3] = 2.;
    AThresholds[7][3] = 2.;
    AThresholds[8][3] = 0.02;
    AThresholds[9][3] = 0.1;
    AThresholds[10][3] = 0.02;
    AThresholds[11][3] = 1.0;
    AThresholds[12][3] = 3.;
    AThresholds[14][3] = 2.;
    AThresholds[15][3] = 2.;
    VThresholds[1][3] = 100.;
    VThresholds[2][3] = 40.;
    VThresholds[3][3] = 25.;
    VThresholds[4][3] = 35.;
    VThresholds[5][3] = 45.;
    VThresholds[6][3] = 20.;
    VThresholds[7][3] = 50.;
    VThresholds[8][3] = 5.;
    VThresholds[9][3] = 6.;
    VThresholds[10][3] = 6.;
    VThresholds[11][3] = 25.;
    VThresholds[12][3] = 35.;
    VThresholds[14][3] = 35.;
    VThresholds[15][3] = 35.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if (RunNum == 25) {
    AThresholds[1][3] = 11.;
    AThresholds[2][3] = 2.;
    AThresholds[3][3] = 2.;
    AThresholds[4][3] = 5.;
    AThresholds[5][3] = 7.;
    AThresholds[6][3] = 4.;
    AThresholds[7][3] = 4.;
    AThresholds[8][3] = 0.02;
    AThresholds[9][3] = 0.5;
    AThresholds[10][3] = 0.05;
    AThresholds[11][3] = 3.;
    AThresholds[12][3] = 3.;
    AThresholds[14][3] = 4.;
    AThresholds[15][3] = 4.;
    VThresholds[1][3] = 200.;
    VThresholds[2][3] = 80.;
    VThresholds[3][3] = 50.;
    VThresholds[4][3] = 50.;
    VThresholds[5][3] = 50.;
    VThresholds[6][3] = 50.;
    VThresholds[7][3] = 60.;
    VThresholds[8][3] = 3.;
    VThresholds[9][3] = 18.;
    VThresholds[10][3] = 7.;
    VThresholds[11][3] = 60.;
    VThresholds[12][3] = 50.;
    VThresholds[14][3] = 40.;
    VThresholds[15][3] = 50.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if (RunNum == 28) {
    AThresholds[1][3] = 11.;
    AThresholds[2][3] = 2.;
    AThresholds[3][3] = 2.;
    AThresholds[4][3] = 5.;
    AThresholds[5][3] = 7.;
    AThresholds[6][3] = 4.;
    AThresholds[7][3] = 4.;
    AThresholds[8][3] = 0.02;
    AThresholds[9][3] = 0.5;
    AThresholds[10][3] = 0.05;
    AThresholds[11][3] = 3.;
    AThresholds[12][3] = 3.;
    AThresholds[14][3] = 4.;
    AThresholds[15][3] = 4.;
    VThresholds[1][3] = 200.;
    VThresholds[2][3] = 80.;
    VThresholds[3][3] = 50.;
    VThresholds[4][3] = 50.;
    VThresholds[5][3] = 50.;
    VThresholds[6][3] = 50.;
    VThresholds[7][3] = 60.;
    VThresholds[8][3] = 3.;
    VThresholds[9][3] = 18.;
    VThresholds[10][3] = 7.;
    VThresholds[11][3] = 60.;
    VThresholds[12][3] = 50.;
    VThresholds[14][3] = 40.;
    VThresholds[15][3] = 50.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if ((RunNum == 32)||(RunNum == 33)) {
    AThresholds[1][3] = 40.;
    AThresholds[2][3] = 20.;
    AThresholds[3][3] = 10.;
    AThresholds[4][3] = 12.;
    AThresholds[5][3] = 18.;
    AThresholds[6][3] = 12.;
    AThresholds[7][3] = 18.;
    AThresholds[8][3] = 0.7;
    AThresholds[9][3] = 3.;
    AThresholds[10][3] = 1.8;
    AThresholds[11][3] = 10.;
    AThresholds[12][3] = 18.;
    AThresholds[14][3] = 11.;
    AThresholds[15][3] = 12.;
    VThresholds[1][3] = 350.;
    VThresholds[2][3] = 200.;
    VThresholds[3][3] = 350.;
    VThresholds[4][3] = 200.;
    VThresholds[5][3] = 250.;
    VThresholds[6][3] = 220.;
    VThresholds[7][3] = 400.;
    VThresholds[8][3] = 20.;
    VThresholds[9][3] = 80.;
    VThresholds[10][3] = 100.;
    VThresholds[11][3] = 250.;
    VThresholds[12][3] = 250.;
    VThresholds[14][3] = 180.;
    VThresholds[15][3] = 180.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if (RunNum == 48) {
    AThresholds[1][3] = 8.;
    AThresholds[2][3] = 5.;
    AThresholds[3][3] = 1.5;
    AThresholds[4][3] = 8.;
    AThresholds[5][3] = 10.;
    AThresholds[6][3] = 6.;
    AThresholds[7][3] = 4.;
    AThresholds[8][3] = 0.2;
    AThresholds[9][3] = 1.5;
    AThresholds[10][3] = 1.5;
    AThresholds[11][3] = 5.;
    AThresholds[12][3] = 9.;
    AThresholds[14][3] = 8.;
    AThresholds[15][3] = 6.;
    VThresholds[1][3] = 100.;
    VThresholds[2][3] = 100.;
    VThresholds[3][3] = 50.;
    VThresholds[4][3] = 80.;
    VThresholds[5][3] = 80.;
    VThresholds[6][3] = 80.;
    VThresholds[7][3] = 80.;
    VThresholds[8][3] = 8.;
    VThresholds[9][3] = 35.;
    VThresholds[10][3] = 80.;
    VThresholds[11][3] = 80.;
    VThresholds[12][3] = 100.;
    VThresholds[14][3] = 80.;
    VThresholds[15][3] = 80.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if (RunNum == 50) {
    AThresholds[1][3] = 10.;
    AThresholds[2][3] = 8.;
    AThresholds[3][3] = 3.;
    AThresholds[4][3] = 10.;
    AThresholds[5][3] = 15.;
    AThresholds[6][3] = 8.;
    AThresholds[7][3] = 6.;
    AThresholds[8][3] = 0.5;
    AThresholds[9][3] = 2.5;
    AThresholds[10][3] = 0.8;
    AThresholds[11][3] = 9.;
    AThresholds[12][3] = 15.;
    AThresholds[14][3] = 9.;
    AThresholds[15][3] = 10.;
    VThresholds[1][3] = 250.;
    VThresholds[2][3] = 200.;
    VThresholds[3][3] = 70.;
    VThresholds[4][3] = 150.;
    VThresholds[5][3] = 200.;
    VThresholds[6][3] = 120.;
    VThresholds[7][3] = 150.;
    VThresholds[8][3] = 20.;
    VThresholds[9][3] = 60.;
    VThresholds[10][3] = 30.;
    VThresholds[11][3] = 200.;
    VThresholds[12][3] = 200.;
    VThresholds[14][3] = 130.;
    VThresholds[15][3] = 150.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  
  if (RunNum == 61) {
    AThresholds[1][3] = 60.;
    AThresholds[2][3] = 60.;
    AThresholds[3][3] = 60.;
    AThresholds[4][3] = 10.;
    AThresholds[5][3] = 10.;
    AThresholds[6][3] = 10.;
    AThresholds[7][3] = 10.;
    AThresholds[8][3] = 20.;
    AThresholds[9][3] = 50.;
    AThresholds[10][3] = 20.;
    AThresholds[11][3] = 70.;
    AThresholds[12][3] = 10.;
    AThresholds[14][3] = 10.;
    AThresholds[15][3] = 10.;
    VThresholds[1][3] = 600.;
    VThresholds[2][3] = 600.;
    VThresholds[3][3] = 600.;
    VThresholds[4][3] = 100.;
    VThresholds[5][3] = 100.;
    VThresholds[6][3] = 100.;
    VThresholds[7][3] = 600.;
    VThresholds[8][3] = 600.;
    VThresholds[9][3] = 250.;
    VThresholds[10][3] = 600.;
    VThresholds[11][3] = 250.;
    VThresholds[12][3] = 100.;
    VThresholds[14][3] = 100.;
    VThresholds[15][3] = 100.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }

  if (RunNum == 97) {
    AThresholds[1][3] = 2.;
    AThresholds[2][3] = 0.5;
    AThresholds[3][3] = 0.5;
    AThresholds[4][3] = 1.;
    AThresholds[5][3] = 1.;
    AThresholds[6][3] = 1.;
    AThresholds[7][3] = 1.;
    AThresholds[8][3] = 0.01;
    AThresholds[9][3] = 0.05;
    AThresholds[10][3] = 0.01;
    AThresholds[11][3] = 0.5;
    AThresholds[12][3] = 1.;
    AThresholds[14][3] = 1.;
    AThresholds[15][3] = 1.;
    VThresholds[1][3] = 40.;
    VThresholds[2][3] = 20.;
    VThresholds[3][3] = 10.;
    VThresholds[4][3] = 10.;
    VThresholds[5][3] = 10.;
    VThresholds[6][3] = 10.;
    VThresholds[7][3] = 20.;
    VThresholds[8][3] = 2.;
    VThresholds[9][3] = 2.;
    VThresholds[10][3] = 2.;
    VThresholds[11][3] = 10.;
    VThresholds[12][3] = 15.;
    VThresholds[14][3] = 15.;
    VThresholds[15][3] = 15.;
    DThresholds[1][3] = 10.;
    DThresholds[2][3] = 10.;
    DThresholds[3][3] = 10.;
    DThresholds[4][3] = 10.;
    DThresholds[5][3] = 10.;
    DThresholds[6][3] = 10.;
    DThresholds[7][3] = 10.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 10.;
    DThresholds[12][3] = 10.;
    DThresholds[14][3] = 10.;
    DThresholds[15][3] = 10.;
  }
  if (RunNum == 242 || RunNum == 243) {
    AThresholds[0][3] = 50.;
    AThresholds[1][3] = 50.;
    AThresholds[2][3] = 50.;
    AThresholds[3][3] = 50.;
    AThresholds[4][3] = 50.;
    AThresholds[5][3] = 50.;
    AThresholds[6][3] = 50.;
    AThresholds[7][3] = 50.;
    AThresholds[8][3] = 50.;
    AThresholds[9][3] = 50.;
    AThresholds[10][3] = 50.;
    AThresholds[11][3] = 50.;
    AThresholds[12][3] = 50.;
    AThresholds[14][3] = 50.;
    AThresholds[15][3] = 50.;
    VThresholds[0][3] = 500.;
    VThresholds[1][3] = 500.;
    VThresholds[2][3] = 500.;
    VThresholds[3][3] = 500.;
    VThresholds[4][3] = 500.;
    VThresholds[5][3] = 500.;
    VThresholds[6][3] = 500.;
    VThresholds[7][3] = 500.;
    VThresholds[8][3] = 500.;
    VThresholds[9][3] = 500.;
    VThresholds[10][3] = 500.;
    VThresholds[11][3] = 500.;
    VThresholds[12][3] = 500.;
    VThresholds[14][3] = 500.;
    VThresholds[15][3] = 500.;
    DThresholds[0][3] = 70.;
    DThresholds[1][3] = 70.;
    DThresholds[2][3] = 30.;
    DThresholds[3][3] = 40.;
    DThresholds[4][3] = 70.;
    DThresholds[5][3] = 80.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 50.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 30.;
    DThresholds[12][3] = 70.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if (RunNum >= 247 ) {
    AThresholds[0][3] = 2;
    AThresholds[1][3] = 2;
    AThresholds[2][3] = 0.1;
    AThresholds[3][3] = 0.1;
    AThresholds[4][3] = 0.1;
    AThresholds[5][3] = 0.1;
    AThresholds[6][3] = 2;
    AThresholds[7][3] = 2;
    AThresholds[8][3] = 2;
    AThresholds[9][3] = 2;
    AThresholds[10][3] = 2;
    AThresholds[11][3] = 2;
    AThresholds[12][3] = 2;
    AThresholds[13][3] = 0.1;
    AThresholds[14][3] = 2;
    AThresholds[15][3] = 2;
    VThresholds[0][3] = 30.;
    VThresholds[1][3] = 30.;
    VThresholds[2][3] = 7.;
    VThresholds[3][3] = 3.;
    VThresholds[4][3] = 3.;
    VThresholds[5][3] = 3.;
    VThresholds[6][3] = 30.;
    VThresholds[7][3] = 30.;
    VThresholds[8][3] = 30.;
    VThresholds[9][3] = 30.;
    VThresholds[10][3] = 30.;
    VThresholds[11][3] = 30.;
    VThresholds[12][3] = 30.;
    VThresholds[13][3] = 3.;
    VThresholds[14][3] = 30.;
    VThresholds[15][3] = 30.;
    DThresholds[0][3] = 60.;
    DThresholds[1][3] = 60.;
    DThresholds[2][3] = 20.;
    DThresholds[3][3] = 20.;
    DThresholds[4][3] = 20.;
    DThresholds[5][3] = 20.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 60.;
    DThresholds[8][3] = 20.;
    DThresholds[9][3] = 20.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 60.;
    DThresholds[12][3] = 60.;
    DThresholds[13][3] = 20.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if (RunNum >= 253 ) {
    // Needs to be updated.
    AThresholds[0][3] = 7;
    AThresholds[1][3] = 10;
    AThresholds[2][3] = 1.2;
    AThresholds[3][3] = 0.5;
    AThresholds[4][3] = 0.2;
    AThresholds[5][3] = 1.1;
    AThresholds[6][3] = 10;
    AThresholds[7][3] = 12;
    AThresholds[8][3] = 10;
    AThresholds[9][3] = 2;
    AThresholds[10][3] = 5;
    AThresholds[11][3] = 6;
    AThresholds[12][3] = 12;
    AThresholds[13][3] = 0.7;
    AThresholds[14][3] = 10;
    AThresholds[15][3] = 10;
    VThresholds[0][3] = 60.;
    VThresholds[1][3] = 80.;
    VThresholds[2][3] = 30.;
    VThresholds[3][3] = 10.;
    VThresholds[4][3] = 6.;
    VThresholds[5][3] = 20.;
    VThresholds[6][3] = 100.;
    VThresholds[7][3] = 150.;
    VThresholds[8][3] = 120.;
    VThresholds[9][3] = 30.;
    VThresholds[10][3] = 150.;
    VThresholds[11][3] = 70.;
    VThresholds[12][3] = 150.;
    VThresholds[13][3] = 15.;
    VThresholds[14][3] = 150.;
    VThresholds[15][3] = 120.;
    // Duration thresholds not yet optimized
    DThresholds[0][3] = 60.;
    DThresholds[1][3] = 60.;
    DThresholds[2][3] = 20.;
    DThresholds[3][3] = 20.;
    DThresholds[4][3] = 20.;
    DThresholds[5][3] = 20.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 60.;
    DThresholds[8][3] = 20.;
    DThresholds[9][3] = 20.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 60.;
    DThresholds[12][3] = 60.;
    DThresholds[13][3] = 20.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }
  if (RunNum >= 255 ) {
    AThresholds[0][3] = 10;
    AThresholds[1][3] = 10;
    AThresholds[2][3] = 2;
    AThresholds[3][3] = 2;
    AThresholds[4][3] = 0.5;
    AThresholds[5][3] = 2;
    AThresholds[6][3] = 10;
    AThresholds[7][3] = 15;
    AThresholds[8][3] = 10;
    AThresholds[9][3] = 2;
    AThresholds[10][3] = 10;
    AThresholds[11][3] = 10;
    AThresholds[12][3] = 15;
    AThresholds[13][3] = 2;
    AThresholds[14][3] = 15;
    AThresholds[15][3] = 15;
    VThresholds[0][3] = 100.;
    VThresholds[1][3] = 100.;
    VThresholds[2][3] = 60.;
    VThresholds[3][3] = 30.;
    VThresholds[4][3] = 18.;
    VThresholds[5][3] = 50.;
    VThresholds[6][3] = 100.;
    VThresholds[7][3] = 200.;
    VThresholds[8][3] = 200.;
    VThresholds[9][3] = 100.;
    VThresholds[10][3] = 30.;
    VThresholds[11][3] = 30.;
    VThresholds[12][3] = 200.;
    VThresholds[13][3] = 40.;
    VThresholds[14][3] = 200.;
    VThresholds[15][3] = 150.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 60.;
    DThresholds[1][3] = 60.;
    DThresholds[2][3] = 20.;
    DThresholds[3][3] = 20.;
    DThresholds[4][3] = 20.;
    DThresholds[5][3] = 20.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 60.;
    DThresholds[8][3] = 20.;
    DThresholds[9][3] = 20.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 60.;
    DThresholds[12][3] = 60.;
    DThresholds[13][3] = 20.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }

  if (RunNum >= 286 ) {
    AThresholds[0][3] = 1.5;
    AThresholds[1][3] = 2;
    AThresholds[2][3] = 20;
    AThresholds[3][3] = 14;
    AThresholds[4][3] = 4;
    AThresholds[5][3] = 20;
    AThresholds[6][3] = 2;
    AThresholds[7][3] = 2.5;
    AThresholds[8][3] = 2;
    AThresholds[9][3] = 2.5;
    AThresholds[10][3] = 0.5;
    AThresholds[11][3] = 2;
    AThresholds[12][3] = 2.5;
    AThresholds[13][3] = 14;
    AThresholds[14][3] = 2;
    AThresholds[15][3] = 2;
    VThresholds[0][3] = 20.;
    VThresholds[1][3] = 25.;
    VThresholds[2][3] = 600.;
    VThresholds[3][3] = 250.;
    VThresholds[4][3] = 120.;
    VThresholds[5][3] = 400.;
    VThresholds[6][3] = 30.;
    VThresholds[7][3] = 35.;
    VThresholds[8][3] = 30.;
    VThresholds[9][3] = 35.;
    VThresholds[10][3] = 20.;
    VThresholds[11][3] = 25.;
    VThresholds[12][3] = 35.;
    VThresholds[13][3] = 300.;
    VThresholds[14][3] = 30.;
    VThresholds[15][3] = 30.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 60.;
    DThresholds[1][3] = 60.;
    DThresholds[2][3] = 20.;
    DThresholds[3][3] = 20.;
    DThresholds[4][3] = 20.;
    DThresholds[5][3] = 20.;
    DThresholds[6][3] = 60.;
    DThresholds[7][3] = 60.;
    DThresholds[8][3] = 20.;
    DThresholds[9][3] = 20.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 60.;
    DThresholds[12][3] = 60.;
    DThresholds[13][3] = 20.;
    DThresholds[14][3] = 60.;
    DThresholds[15][3] = 60.;
  }


  if (RunNum >= 288 ) { // First run after rechanneling for second digitizer 
    AThresholds[0][3] = 0.1;
    AThresholds[1][3] = 0.1;
    AThresholds[2][3] = 0.1;
    AThresholds[3][3] = 0.1;
    AThresholds[4][3] = 0.1;
    AThresholds[5][3] = 0.1;
    AThresholds[6][3] = 0.1;
    AThresholds[7][3] = 0.1;
    AThresholds[8][3] = 0.1;
    AThresholds[9][3] = 0.1;
    AThresholds[10][3] = 0.1;
    AThresholds[11][3] = 0.1;
    AThresholds[12][3] = 0.1;
    AThresholds[13][3] = 0.1;
    AThresholds[14][3] = 0.1;
    AThresholds[15][3] = 0.1;
    AThresholds[16][3] = 0.1;
    AThresholds[17][3] = 0.1;
    AThresholds[18][3] = 0.1;
    AThresholds[19][3] = 0.1;
    AThresholds[20][3] = 0.1;
    AThresholds[21][3] = 0.1;
    AThresholds[22][3] = 0.1;
    AThresholds[23][3] = 0.1;
    AThresholds[24][3] = 0.1;
    AThresholds[25][3] = 0.1;
    AThresholds[26][3] = 0.1;
    AThresholds[27][3] = 0.1;
    AThresholds[28][3] = 0.1;
    AThresholds[29][3] = 0.1;
    AThresholds[30][3] = 0.1;
    AThresholds[31][3] = 0.1;
    VThresholds[0][3] = 25.;
    VThresholds[1][3] = 25.;
    VThresholds[2][3] = 30.;
    VThresholds[3][3] = 25.;
    VThresholds[4][3] = 30.;
    VThresholds[5][3] = 15.;
    VThresholds[6][3] = 30.;
    VThresholds[7][3] = 30.;
    VThresholds[8][3] = 100.;
    VThresholds[9][3] = 200.;
    VThresholds[10][3] = 6.;
    VThresholds[11][3] = 6.;
    VThresholds[12][3] = 25.;
    VThresholds[13][3] = 25.;
    VThresholds[14][3] = 5.;
    VThresholds[15][3] = 30.;
    VThresholds[16][3] = 30.;
    VThresholds[17][3] = 300.;
    VThresholds[18][3] = 6.;
    VThresholds[19][3] = 5.;
    VThresholds[20][3] = 6.;
    VThresholds[21][3] = 7.;
    VThresholds[22][3] = 20.;
    VThresholds[23][3] = 25.;
    VThresholds[24][3] = 600.;
    VThresholds[25][3] = 200.;
    VThresholds[26][3] = 6.;
    VThresholds[27][3] = 6.;
    VThresholds[28][3] = 5.;
    VThresholds[29][3] = 5.;
    VThresholds[30][3] = 5.;
    VThresholds[31][3] = 5.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 10.;
    DThresholds[1][3] = 10.;
    DThresholds[2][3] = 10.;
    DThresholds[3][3] = 10.;
    DThresholds[4][3] = 10.;
    DThresholds[5][3] = 10.;
    DThresholds[6][3] = 10.;
    DThresholds[7][3] = 10.;
    DThresholds[8][3] = 10.;
    DThresholds[9][3] = 10.;
    DThresholds[10][3] = 10.;
    DThresholds[11][3] = 10.;
    DThresholds[12][3] = 10.;
    DThresholds[13][3] = 10.;
    DThresholds[14][3] = 10.;
    DThresholds[15][3] = 10.;
    DThresholds[16][3] = 10.;
    DThresholds[17][3] = 10.;
    DThresholds[18][3] = 10.;
    DThresholds[19][3] = 10.;
    DThresholds[20][3] = 10.;
    DThresholds[21][3] = 10.;
    DThresholds[22][3] = 10.;
    DThresholds[23][3] = 10.;
    DThresholds[24][3] = 10.;
    DThresholds[25][3] = 10.;
    DThresholds[26][3] = 10.;
    DThresholds[27][3] = 10.;
    DThresholds[28][3] = 10.;
    DThresholds[29][3] = 10.;
    DThresholds[30][3] = 10.;
    DThresholds[31][3] = 10.;
  }

  if (RunNum >= 374 ) { // 1600 on ET and 1450 on all else
    AThresholds[0][3] = 20;
    AThresholds[1][3] = 20;
    AThresholds[2][3] = 20;
    AThresholds[3][3] = 20;
    AThresholds[4][3] = 0.1;
    AThresholds[5][3] = 0.1;
    AThresholds[6][3] = 20;
    AThresholds[7][3] = 20;
    AThresholds[8][3] = 0.1;
    AThresholds[9][3] = 0.1;
    AThresholds[10][3] = 0.1;
    AThresholds[11][3] = 0.1;
    AThresholds[12][3] = 0.1;
    AThresholds[13][3] = 0.1;
    AThresholds[14][3] = 0.1;
    AThresholds[15][3] = 0.1;
    AThresholds[16][3] = 0.1;
    AThresholds[17][3] = 0.1;
    AThresholds[18][3] = 0.1;
    AThresholds[19][3] = 0.1;
    AThresholds[20][3] = 0.1;
    AThresholds[21][3] = 0.1;
    AThresholds[22][3] = 0.1;
    AThresholds[23][3] = 0.1;
    AThresholds[24][3] = 0.1;
    AThresholds[25][3] = 0.1;
    AThresholds[26][3] = 0.1;
    AThresholds[27][3] = 0.1;
    AThresholds[28][3] = 0.1;
    AThresholds[29][3] = 0.1;
    AThresholds[30][3] = 0.1;
    AThresholds[31][3] = 0.1;
    VThresholds[0][3] = 300.;
    VThresholds[1][3] = 300.;
    VThresholds[2][3] = 350.;
    VThresholds[3][3] = 350.;
    VThresholds[4][3] = 30.;
    VThresholds[5][3] = 15.;
    VThresholds[6][3] = 400.;
    VThresholds[7][3] = 400.;
    VThresholds[8][3] = 100.;
    VThresholds[9][3] = 100.;
    VThresholds[10][3] = 30.;
    VThresholds[11][3] = 30.;
    VThresholds[12][3] = 25.;
    VThresholds[13][3] = 25.;
    VThresholds[14][3] = 30.;
    VThresholds[15][3] = 30.;
    VThresholds[16][3] = 30.;
    VThresholds[17][3] = 200.;
    VThresholds[18][3] = 600.;
    VThresholds[19][3] = 30.;
    VThresholds[20][3] = 500.;
    VThresholds[21][3] = 500.;
    VThresholds[22][3] = 20.;
    VThresholds[23][3] = 25.;
    VThresholds[24][3] = 200.;
    VThresholds[25][3] = 200.;
    VThresholds[26][3] = 30.;
    VThresholds[27][3] = 30.;
    VThresholds[28][3] = 500.;
    VThresholds[29][3] = 30.;
    VThresholds[30][3] = 30.;
    VThresholds[31][3] = 30.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 400.; // 878
    DThresholds[1][3] = 400.; // 878
    DThresholds[2][3] = 400.; // 878
    DThresholds[3][3] = 400.; // 878
    DThresholds[4][3] = 400.; // 878
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 400.; // 878
    DThresholds[7][3] = 400.; // 878
    DThresholds[8][3] = 50.; // ET
    DThresholds[9][3] = 50.; // ET
    DThresholds[10][3] = 10.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 400.; // 878
    DThresholds[13][3] = 400.; // 878
    DThresholds[14][3] = 10.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 400.; // 878
    DThresholds[17][3] = 50.; // ET
    DThresholds[18][3] = 400.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 400.; // 878 - L1 slab
    DThresholds[21][3] = 400.; // 878 - Back slab
    DThresholds[22][3] = 40.; // 7725
    DThresholds[23][3] = 400.; // 878
    DThresholds[24][3] = 50.; // ET
    DThresholds[25][3] = 50.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 400.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 10.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 424 ) { // 1600 on ET and 1450 on all else
    AThresholds[0][3] = 3;
    AThresholds[1][3] = 3;
    AThresholds[2][3] = 5;
    AThresholds[3][3] = 5;
    AThresholds[4][3] = 5;
    AThresholds[5][3] = 2;
    AThresholds[6][3] = 4;
    AThresholds[7][3] = 6;
    AThresholds[8][3] = 0.6;
    AThresholds[9][3] = 3;
    AThresholds[10][3] = 2;
    AThresholds[11][3] = 2;
    AThresholds[12][3] = 5;
    AThresholds[13][3] = 5;
    AThresholds[14][3] = 0.1;
    AThresholds[15][3] = 0.1;
    AThresholds[16][3] = 5;
    AThresholds[17][3] = 2;
    AThresholds[18][3] = 60;
    AThresholds[19][3] = 0.1;
    AThresholds[20][3] = 0.5;
    AThresholds[21][3] = 0.1;
    AThresholds[22][3] = 1.5;
    AThresholds[23][3] = 3;
    AThresholds[24][3] = 1.5;
    AThresholds[25][3] = 2;
    AThresholds[26][3] = 0.1;
    AThresholds[27][3] = 1;
    AThresholds[28][3] = 0.2;
    AThresholds[29][3] = 1;
    AThresholds[30][3] = 2;
    AThresholds[31][3] = 1;
    VThresholds[0][3] = 50.;
    VThresholds[1][3] = 60.;
    VThresholds[2][3] = 70.;
    VThresholds[3][3] = 70.;
    VThresholds[4][3] = 70.;
    VThresholds[5][3] = 50.;
    VThresholds[6][3] = 60.;
    VThresholds[7][3] = 80.;
    VThresholds[8][3] = 15.;
    VThresholds[9][3] = 50.;
    VThresholds[10][3] = 40.;
    VThresholds[11][3] = 30.;
    VThresholds[12][3] = 60.;
    VThresholds[13][3] = 60.;
    VThresholds[14][3] = 10.;
    VThresholds[15][3] = 30.;
    VThresholds[16][3] = 80.;
    VThresholds[17][3] = 50.;
    VThresholds[18][3] = 600.;
    VThresholds[19][3] = 10.;
    VThresholds[20][3] = 15.;
    VThresholds[21][3] = 10.;
    VThresholds[22][3] = 50.;
    VThresholds[23][3] = 50.;
    VThresholds[24][3] = 100.;
    VThresholds[25][3] = 30.;
    VThresholds[26][3] = 10.;
    VThresholds[27][3] = 20.;
    VThresholds[28][3] = 10.;
    VThresholds[29][3] = 20.;
    VThresholds[30][3] = 40.;
    VThresholds[31][3] = 30.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 40.; // 878
    DThresholds[1][3] = 60.; // 878
    DThresholds[2][3] = 60.; // 878
    DThresholds[3][3] = 60.; // 878
    DThresholds[4][3] = 60.; // 878
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 60.; // 878
    DThresholds[7][3] = 60.; // 878
    DThresholds[8][3] = 40.; // ET
    DThresholds[9][3] = 40.; // ET
    DThresholds[10][3] = 60.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 60.; // 878
    DThresholds[13][3] = 60.; // 878
    DThresholds[14][3] = 20.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 60.; // 878
    DThresholds[17][3] = 40.; // ET
    DThresholds[18][3] = 400.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 60.; // 878 - L1 slab
    DThresholds[21][3] = 60.; // 878 - Back slab
    DThresholds[22][3] = 20.; // 7725
    DThresholds[23][3] = 40.; // 878
    DThresholds[24][3] = 40.; // ET
    DThresholds[25][3] = 40.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 20.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 50.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 430 ) { // ET=1050, 878=750, 7725=900
    AThresholds[0][3] = 1;
    AThresholds[1][3] = 1;
    AThresholds[2][3] = 2;
    AThresholds[3][3] = 2;
    AThresholds[4][3] = 2;
    AThresholds[5][3] = 5;
    AThresholds[6][3] = 1.5;
    AThresholds[7][3] = 2;
    AThresholds[8][3] = 1;
    AThresholds[9][3] = 5;
    AThresholds[10][3] = 2;
    AThresholds[11][3] = 2;
    AThresholds[12][3] = 2;
    AThresholds[13][3] = 2;
    AThresholds[14][3] = 0.1;
    AThresholds[15][3] = 999; // LHC clock
    AThresholds[16][3] = 2;
    AThresholds[17][3] = 2;
    AThresholds[18][3] = 2;
    AThresholds[19][3] = 0.1;
    AThresholds[20][3] = 0.1;
    AThresholds[21][3] = 0.1;
    AThresholds[22][3] = 2;
    AThresholds[23][3] = 2;
    AThresholds[24][3] = 5;
    AThresholds[25][3] = 3;
    AThresholds[26][3] = 0.1;
    AThresholds[27][3] = 1;
    AThresholds[28][3] = 0.1;
    AThresholds[29][3] = 1;
    AThresholds[30][3] = 2;
    AThresholds[31][3] = 2;
    VThresholds[0][3] = 20.;
    VThresholds[1][3] = 30.;
    VThresholds[2][3] = 30.;
    VThresholds[3][3] = 30.;
    VThresholds[4][3] = 30.;
    VThresholds[5][3] = 30.;
    VThresholds[6][3] = 30.;
    VThresholds[7][3] = 30.;
    VThresholds[8][3] = 10.;
    VThresholds[9][3] = 40.;
    VThresholds[10][3] = 40.;
    VThresholds[11][3] = 30.;
    VThresholds[12][3] = 40.;
    VThresholds[13][3] = 30.;
    VThresholds[14][3] = 20.;
    VThresholds[15][3] = 999.; // LHC clock
    VThresholds[16][3] = 20.;
    VThresholds[17][3] = 30.;
    VThresholds[18][3] = 400.;
    VThresholds[19][3] = 10.;
    VThresholds[20][3] = 10.;
    VThresholds[21][3] = 10.;
    VThresholds[22][3] = 1000.;
    VThresholds[23][3] = 20.;
    VThresholds[24][3] = 50.;
    VThresholds[25][3] = 30.;
    VThresholds[26][3] = 10.;
    VThresholds[27][3] = 20.;
    VThresholds[28][3] = 10.;
    VThresholds[29][3] = 20.;
    VThresholds[30][3] = 40.;
    VThresholds[31][3] = 30.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 40.; // 878
    DThresholds[1][3] = 60.; // 878
    DThresholds[2][3] = 60.; // 878
    DThresholds[3][3] = 60.; // 878
    DThresholds[4][3] = 60.; // 878
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 60.; // 878
    DThresholds[7][3] = 60.; // 878
    DThresholds[8][3] = 40.; // ET
    DThresholds[9][3] = 40.; // ET
    DThresholds[10][3] = 60.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 60.; // 878
    DThresholds[13][3] = 60.; // 878
    DThresholds[14][3] = 20.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 60.; // 878
    DThresholds[17][3] = 40.; // ET
    DThresholds[18][3] = 40.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 60.; // 878 - L1 slab
    DThresholds[21][3] = 60.; // 878 - Back slab
    DThresholds[22][3] = 20.; // 7725
    DThresholds[23][3] = 40.; // 878
    DThresholds[24][3] = 40.; // ET
    DThresholds[25][3] = 40.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 20.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 50.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 485 ) { // ET=1050, 878=750, 7725=900
    AThresholds[0][3] = 2;
    AThresholds[1][3] = 3;
    AThresholds[2][3] = 4;
    AThresholds[3][3] = 4;
    AThresholds[4][3] = 4;
    AThresholds[5][3] = 3;
    AThresholds[6][3] = 3;
    AThresholds[7][3] = 3;
    AThresholds[8][3] = 0.5;
    AThresholds[9][3] = 3;
    AThresholds[10][3] = 2;
    AThresholds[11][3] = 2;
    AThresholds[12][3] = 4;
    AThresholds[13][3] = 4;
    AThresholds[14][3] = 2;
    AThresholds[15][3] = 999; // LHC clock
    AThresholds[16][3] = 4;
    AThresholds[17][3] = 2;
    AThresholds[18][3] = 80;
    AThresholds[19][3] = 2;
    AThresholds[20][3] = 40;
    AThresholds[21][3] = 40;
    AThresholds[22][3] = 20;
    AThresholds[23][3] = 60;
    AThresholds[24][3] = 3;
    AThresholds[25][3] = 1;
    AThresholds[26][3] = 2;
    AThresholds[27][3] = 2;
    AThresholds[28][3] = 20;
    AThresholds[29][3] = 2;
    AThresholds[30][3] = 2;
    AThresholds[31][3] = 2;
    VThresholds[0][3] = 30.;
    VThresholds[1][3] = 40.;
    VThresholds[2][3] = 50.;
    VThresholds[3][3] = 50.;
    VThresholds[4][3] = 30.;
    VThresholds[5][3] = 50.;
    VThresholds[6][3] = 50.;
    VThresholds[7][3] = 50.;
    VThresholds[8][3] = 15.;
    VThresholds[9][3] = 50.;
    VThresholds[10][3] = 40.;
    VThresholds[11][3] = 20.;
    VThresholds[12][3] = 50.;
    VThresholds[13][3] = 50.;
    VThresholds[14][3] = 30.;
    VThresholds[15][3] = 999.; // LHC clock
    VThresholds[16][3] = 50.;
    VThresholds[17][3] = 50.;
    VThresholds[18][3] = 600.;
    VThresholds[19][3] = 50.;
    VThresholds[20][3] = 400.;
    VThresholds[21][3] = 600.;
    VThresholds[22][3] = 1000.;
    VThresholds[23][3] = 400.;
    VThresholds[24][3] = 100.;
    VThresholds[25][3] = 40.;
    VThresholds[26][3] = 30.;
    VThresholds[27][3] = 20.;
    VThresholds[28][3] = 300.;
    VThresholds[29][3] = 20.;
    VThresholds[30][3] = 40.;
    VThresholds[31][3] = 30.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 40.; // 878
    DThresholds[1][3] = 60.; // 878
    DThresholds[2][3] = 60.; // 878
    DThresholds[3][3] = 60.; // 878
    DThresholds[4][3] = 60.; // 878
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 60.; // 878
    DThresholds[7][3] = 60.; // 878
    DThresholds[8][3] = 40.; // ET
    DThresholds[9][3] = 40.; // ET
    DThresholds[10][3] = 60.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 60.; // 878
    DThresholds[13][3] = 60.; // 878
    DThresholds[14][3] = 20.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 60.; // 878
    DThresholds[17][3] = 40.; // ET
    DThresholds[18][3] = 40.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 60.; // 878 - L1 slab
    DThresholds[21][3] = 60.; // 878 - Back slab
    DThresholds[22][3] = 20.; // 7725
    DThresholds[23][3] = 40.; // 878
    DThresholds[24][3] = 40.; // ET
    DThresholds[25][3] = 40.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 20.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 50.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }


  if (RunNum >= 491 ) { // ET=1100, 878=825, 7725=875
    AThresholds[0][3] = 4;
    AThresholds[1][3] = 5;
    AThresholds[2][3] = 5;
    AThresholds[3][3] = 5;
    AThresholds[4][3] = 4;
    AThresholds[5][3] = 3;
    AThresholds[6][3] = 5;
    AThresholds[7][3] = 5;
    AThresholds[8][3] = 1.5;
    AThresholds[9][3] = 9;
    AThresholds[10][3] = 2;
    AThresholds[11][3] = 2;
    AThresholds[12][3] = 5;
    AThresholds[13][3] = 5;
    AThresholds[14][3] = 2;
    AThresholds[15][3] = 999; // LHC clock
    AThresholds[16][3] = 5;
    AThresholds[17][3] = 5;
    AThresholds[18][3] = 40;
    AThresholds[19][3] = 2;
    AThresholds[20][3] = 30;
    AThresholds[21][3] = 30;
    AThresholds[22][3] = 40;
    AThresholds[23][3] = 60;
    AThresholds[24][3] = 8;
    AThresholds[25][3] = 5;
    AThresholds[26][3] = 2;
    AThresholds[27][3] = 1.5;
    AThresholds[28][3] = 15;
    AThresholds[29][3] = 1.5;
    AThresholds[30][3] = 2;
    AThresholds[31][3] = 2;
    VThresholds[0][3] = 50.;
    VThresholds[1][3] = 50.;
    VThresholds[2][3] = 90.;
    VThresholds[3][3] = 60.;
    VThresholds[4][3] = 60.;
    VThresholds[5][3] = 50.;
    VThresholds[6][3] = 70.;
    VThresholds[7][3] = 80.;
    VThresholds[8][3] = 40.;
    VThresholds[9][3] = 100.;
    VThresholds[10][3] = 40.;
    VThresholds[11][3] = 30.;
    VThresholds[12][3] = 60.;
    VThresholds[13][3] = 60.;
    VThresholds[14][3] = 30.;
    VThresholds[15][3] = 999.; // LHC clock
    VThresholds[16][3] = 60.;
    VThresholds[17][3] = 60.;
    VThresholds[18][3] = 600.;
    VThresholds[19][3] = 60.;
    VThresholds[20][3] = 400.;
    VThresholds[21][3] = 400.;
    VThresholds[22][3] = 1200.;
    VThresholds[23][3] = 500.;
    VThresholds[24][3] = 200.;
    VThresholds[25][3] = 100.;
    VThresholds[26][3] = 30.;
    VThresholds[27][3] = 20.;
    VThresholds[28][3] = 300.;
    VThresholds[29][3] = 20.;
    VThresholds[30][3] = 40.;
    VThresholds[31][3] = 30.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 40.; // 878
    DThresholds[1][3] = 40.; // 878
    DThresholds[2][3] = 40.; // 878
    DThresholds[3][3] = 40.; // 878
    DThresholds[4][3] = 4.; // 878  Affected by noise?
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 40.; // 878
    DThresholds[7][3] = 40.; // 878
    DThresholds[8][3] = 40.; // ET
    DThresholds[9][3] = 40.; // ET
    DThresholds[10][3] = 40.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 40.; // 878
    DThresholds[13][3] = 40.; // 878
    DThresholds[14][3] = 20.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 40.; // 878
    DThresholds[17][3] = 40.; // ET
    DThresholds[18][3] = 40.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 40.; // 878 - L1 slab
    DThresholds[21][3] = 40.; // 878 - Back slab
    DThresholds[22][3] = 40.; // 7725
    DThresholds[23][3] = 400.; // 878
    DThresholds[24][3] = 40.; // ET
    DThresholds[25][3] = 40.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 40.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 40.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }


  if (RunNum >= 497 ) { // ET=1100, 878=825, 7725=875, B ramping (ramped) to 3.8T
    AThresholds[0][3] = 2;
    AThresholds[1][3] = 2;
    AThresholds[2][3] = 2;
    AThresholds[3][3] = 2;
    AThresholds[4][3] = 2;
    AThresholds[5][3] = 3;
    AThresholds[6][3] = 2;
    AThresholds[7][3] = 2;
    AThresholds[8][3] = 0.5;
    AThresholds[9][3] = 4;
    AThresholds[10][3] = 0.5;
    AThresholds[11][3] = 0.5;
    AThresholds[12][3] = 2;
    AThresholds[13][3] = 2;
    AThresholds[14][3] = 0.5;
    AThresholds[15][3] = 999; // LHC clock
    AThresholds[16][3] = 2;
    AThresholds[17][3] = 2;
    AThresholds[18][3] = 10;
    AThresholds[19][3] = 2;
    AThresholds[20][3] = 10;
    AThresholds[21][3] = 10;
    AThresholds[22][3] = 20;
    AThresholds[23][3] = 40;
    AThresholds[24][3] = 2;
    AThresholds[25][3] = 2;
    AThresholds[26][3] = 1;
    AThresholds[27][3] = 0.5;
    AThresholds[28][3] = 5;
    AThresholds[29][3] = 0.5;
    AThresholds[30][3] = 1;
    AThresholds[31][3] = 1;
    VThresholds[0][3] = 30.;
    VThresholds[1][3] = 30.;
    VThresholds[2][3] = 30.;
    VThresholds[3][3] = 30.;
    VThresholds[4][3] = 20.;
    VThresholds[5][3] = 40.;
    VThresholds[6][3] = 40.;
    VThresholds[7][3] = 40.;
    VThresholds[8][3] = 10.;
    VThresholds[9][3] = 90.;
    VThresholds[10][3] = 10.;
    VThresholds[11][3] = 20.;
    VThresholds[12][3] = 40.;
    VThresholds[13][3] = 40.;
    VThresholds[14][3] = 10.;
    VThresholds[15][3] = 999.; // LHC clock
    VThresholds[16][3] = 40.;
    VThresholds[17][3] = 40.;
    VThresholds[18][3] = 300.;
    VThresholds[19][3] = 40.;
    VThresholds[20][3] = 300.;
    VThresholds[21][3] = 300.;
    VThresholds[22][3] = 1200.;
    VThresholds[23][3] = 400.;
    VThresholds[24][3] = 100.;
    VThresholds[25][3] = 50.;
    VThresholds[26][3] = 20.;
    VThresholds[27][3] = 5.;
    VThresholds[28][3] = 100.;
    VThresholds[29][3] = 10.;
    VThresholds[30][3] = 20.;
    VThresholds[31][3] = 5.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 40.; // 878
    DThresholds[1][3] = 40.; // 878
    DThresholds[2][3] = 40.; // 878
    DThresholds[3][3] = 40.; // 878
    DThresholds[4][3] = 1.; // 878  Affected by noise?
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 40.; // 878
    DThresholds[7][3] = 40.; // 878
    DThresholds[8][3] = 40.; // ET
    DThresholds[9][3] = 40.; // ET
    DThresholds[10][3] = 40.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 40.; // 878
    DThresholds[13][3] = 40.; // 878
    DThresholds[14][3] = 20.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 40.; // 878
    DThresholds[17][3] = 40.; // ET
    DThresholds[18][3] = 40.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 40.; // 878 - L1 slab
    DThresholds[21][3] = 40.; // 878 - Back slab
    DThresholds[22][3] = 40.; // 7725
    DThresholds[23][3] = 400.; // 878
    DThresholds[24][3] = 40.; // ET
    DThresholds[25][3] = 40.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 40.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 40.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 500 ) { // ET=1100, 878=825, 7725=875, B ramping (ramped) to 3.8T
    AThresholds[0][3] = 2;
    AThresholds[1][3] = 2;
    AThresholds[2][3] = 2;
    AThresholds[3][3] = 2;
    AThresholds[4][3] = 2;
    AThresholds[5][3] = 3;
    AThresholds[6][3] = 2;
    AThresholds[7][3] = 2;
    AThresholds[8][3] = 0.5;
    AThresholds[9][3] = 4;
    AThresholds[10][3] = 0.5;
    AThresholds[11][3] = 0.5;
    AThresholds[12][3] = 2;
    AThresholds[13][3] = 2;
    AThresholds[14][3] = 0.5;
    AThresholds[15][3] = 999; // LHC clock
    AThresholds[16][3] = 2;
    AThresholds[17][3] = 2;
    AThresholds[18][3] = 10;
    AThresholds[19][3] = 2;
    AThresholds[20][3] = 10;
    AThresholds[21][3] = 10;
    AThresholds[22][3] = 20;
    AThresholds[23][3] = 40;
    AThresholds[24][3] = 2;
    AThresholds[25][3] = 2;
    AThresholds[26][3] = 1;
    AThresholds[27][3] = 0.5;
    AThresholds[28][3] = 5;
    AThresholds[29][3] = 0.5;
    AThresholds[30][3] = 1;
    AThresholds[31][3] = 1;
    VThresholds[0][3] = 400.;
    VThresholds[1][3] = 400.;
    VThresholds[2][3] = 400.;
    VThresholds[3][3] = 400.;
    VThresholds[4][3] = 400.;
    VThresholds[5][3] = 400.;
    VThresholds[6][3] = 400.;
    VThresholds[7][3] = 400.;
    VThresholds[8][3] = 400.;
    VThresholds[9][3] = 400.;
    VThresholds[10][3] = 10.;
    VThresholds[11][3] = 20.;
    VThresholds[12][3] = 400.;
    VThresholds[13][3] = 400.;
    VThresholds[14][3] = 10.;
    VThresholds[15][3] = 999.; // LHC clock
    VThresholds[16][3] = 400.;
    VThresholds[17][3] = 400.;
    VThresholds[18][3] = 300.;
    VThresholds[19][3] = 40.;
    VThresholds[20][3] = 300.;
    VThresholds[21][3] = 300.;
    VThresholds[22][3] = 1200.;
    VThresholds[23][3] = 400.;
    VThresholds[24][3] = 400.;
    VThresholds[25][3] = 400.;
    VThresholds[26][3] = 20.;
    VThresholds[27][3] = 5.;
    VThresholds[28][3] = 100.;
    VThresholds[29][3] = 10.;
    VThresholds[30][3] = 20.;
    VThresholds[31][3] = 5.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 40.; // 878
    DThresholds[1][3] = 40.; // 878
    DThresholds[2][3] = 40.; // 878
    DThresholds[3][3] = 40.; // 878
    DThresholds[4][3] = 1.; // 878  Affected by noise?
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 40.; // 878
    DThresholds[7][3] = 40.; // 878
    DThresholds[8][3] = 40.; // ET
    DThresholds[9][3] = 40.; // ET
    DThresholds[10][3] = 40.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 40.; // 878
    DThresholds[13][3] = 40.; // 878
    DThresholds[14][3] = 20.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 40.; // 878
    DThresholds[17][3] = 40.; // ET
    DThresholds[18][3] = 40.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 40.; // 878 - L1 slab
    DThresholds[21][3] = 40.; // 878 - Back slab
    DThresholds[22][3] = 40.; // 7725
    DThresholds[23][3] = 400.; // 878
    DThresholds[24][3] = 40.; // ET
    DThresholds[25][3] = 40.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 40.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 40.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 572 ) { // ET=1100, 878=825, 7725=875, B ramping (ramped) to 3.8T. CONFIRMED VALID
    AThresholds[0][3] = 2.0;
    AThresholds[1][3] = 2.0;
    AThresholds[2][3] = 2.0;
    AThresholds[3][3] = 2.0;
    AThresholds[4][3] = 2.0;
    AThresholds[5][3] = 3.0;
    AThresholds[6][3] = 2.0;
    AThresholds[7][3] = 2.0;
    AThresholds[8][3] = 0.5;
    AThresholds[9][3] = 4.0;
    AThresholds[10][3] = 0.5;
    AThresholds[11][3] = 1.0;
    AThresholds[12][3] = 2.0;
    AThresholds[13][3] = 2.0;
    AThresholds[14][3] = 0.5;
    AThresholds[15][3] = 999; // LHC clock
    AThresholds[16][3] = 4.0;
    AThresholds[17][3] = 4.0;
    AThresholds[18][3] = 20.;
    AThresholds[19][3] = 2.0;
    AThresholds[20][3] = 20.;
    AThresholds[21][3] = 20.;
    AThresholds[22][3] = 2.0;
    AThresholds[23][3] = 40.;
    AThresholds[24][3] = 4.0;
    AThresholds[25][3] = 4.0;
    AThresholds[26][3] = 2.0;
    AThresholds[27][3] = 0.5;
    AThresholds[28][3] = 8.0;
    AThresholds[29][3] = 0.4;
    AThresholds[30][3] = 1.0;
    AThresholds[31][3] = 1.5;
    VThresholds[0][3] = 40.;
    VThresholds[1][3] = 70.;
    VThresholds[2][3] = 80.;
    VThresholds[3][3] = 60.;
    VThresholds[4][3] = 70.;
    VThresholds[5][3] = 70.;
    VThresholds[6][3] = 70.;
    VThresholds[7][3] = 100.;
    VThresholds[8][3] = 15.;
    VThresholds[9][3] = 150.;
    VThresholds[10][3] = 20.;
    VThresholds[11][3] = 25.;
    VThresholds[12][3] = 70.;
    VThresholds[13][3] = 70.;
    VThresholds[14][3] = 10.;
    VThresholds[15][3] = 999.; // LHC clock
    VThresholds[16][3] = 80.;
    VThresholds[17][3] = 100.;
    VThresholds[18][3] = 300.;
    VThresholds[19][3] = 50.;
    VThresholds[20][3] = 300.;
    VThresholds[21][3] = 350.;
    VThresholds[22][3] = 100.;
    VThresholds[23][3] = 400.;
    VThresholds[24][3] = 150.;
    VThresholds[25][3] = 90.;
    VThresholds[26][3] = 30.;
    VThresholds[27][3] = 8.;
    VThresholds[28][3] = 100.;
    VThresholds[29][3] = 8.;
    VThresholds[30][3] = 30.;
    VThresholds[31][3] = 30.;
    // Duration thresholds not updated.
    DThresholds[0][3] = 40.; // 878
    DThresholds[1][3] = 40.; // 878
    DThresholds[2][3] = 40.; // 878
    DThresholds[3][3] = 40.; // 878
    DThresholds[4][3] = 1.; // 878  Affected by noise?
    DThresholds[5][3] = 40.; // 7725
    DThresholds[6][3] = 20.; // 878
    DThresholds[7][3] = 40.; // 878
    DThresholds[8][3] = 40.; // ET
    DThresholds[9][3] = 40.; // ET
    DThresholds[10][3] = 40.; // 878 - L0 top panel
    DThresholds[11][3] = 10.; // 878 - L1 top panel
    DThresholds[12][3] = 40.; // 878
    DThresholds[13][3] = 40.; // 878
    DThresholds[14][3] = 20.; // 878 - L2 top panel
    DThresholds[15][3] = 10.; // LHC clock
    DThresholds[16][3] = 40.; // 878
    DThresholds[17][3] = 40.; // ET
    DThresholds[18][3] = 40.; // 878 - L0 slab
    DThresholds[19][3] = 10.; // 878 - L1 right panel
    DThresholds[20][3] = 40.; // 878 - L1 slab
    DThresholds[21][3] = 40.; // 878 - Back slab
    DThresholds[22][3] = 40.; // 7725
    DThresholds[23][3] = 400.; // 878
    DThresholds[24][3] = 40.; // ET
    DThresholds[25][3] = 40.; // ET
    DThresholds[26][3] = 10.; // 878 - L2 right panel
    DThresholds[27][3] = 10.; // 878 - L0 left panel
    DThresholds[28][3] = 40.; // 878 - L2 slab
    DThresholds[29][3] = 10.; // 878 - L0 right panel
    DThresholds[30][3] = 40.; // 878 - L1 left panel
    DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 582 ) { // ET=1000, 878=775, 7725=850, B=3.8T.
    VThresholds[0][3] = 30.; AThresholds[0][3] = 2.0; DThresholds[0][3] = 40.; // 878bar
    VThresholds[1][3] = 50.; AThresholds[1][3] = 3.5; DThresholds[1][3] = 40.; // 878bar
    VThresholds[2][3] = 50.; AThresholds[2][3] = 4.0; DThresholds[2][3] = 40.; // 878bar
    VThresholds[3][3] = 50.; AThresholds[3][3] = 4.0; DThresholds[3][3] = 40.; // 878bar
    VThresholds[4][3] = 40.; AThresholds[4][3] = 4.0; DThresholds[4][3] = 40.; // 878bar
    VThresholds[5][3] = 50.; AThresholds[5][3] = 2.5; DThresholds[5][3] = 40.; // 7725
    VThresholds[6][3] = 40.; AThresholds[6][3] = 3.0; DThresholds[6][3] = 20.; // 878bar
    VThresholds[7][3] = 60.; AThresholds[7][3] = 5.0; DThresholds[7][3] = 40.; // 878bar
    VThresholds[8][3] = 5.; AThresholds[8][3] = 0.2; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 60.; AThresholds[9][3] = 3.0; DThresholds[9][3] = 40.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 10.; // 878 - L1 left panel
    VThresholds[12][3] = 50.; AThresholds[12][3] = 4.0; DThresholds[12][3] = 40.; // 878bar
    VThresholds[13][3] = 50.; AThresholds[13][3] = 4.0; DThresholds[13][3] = 40.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.4; DThresholds[14][3] = 20.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 50.; AThresholds[16][3] = 4.0; DThresholds[16][3] = 40.; // 878bar
    VThresholds[17][3] = 40.; AThresholds[17][3] = 2.0;  DThresholds[17][3] = 40.; // ET
    VThresholds[18][3] = 200.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 10.; // 878 - L1 right panel
    VThresholds[20][3] = 200.; AThresholds[20][3] = 18.; DThresholds[20][3] = 80.; // 878 - L1 slab
    VThresholds[21][3] = 320.; AThresholds[21][3] = 25.; DThresholds[21][3] = 80.; // 878 - Back slab
    VThresholds[22][3] = 70.; AThresholds[22][3] = 2.3; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 500.; AThresholds[23][3] = 65.; DThresholds[23][3] = 400.; // 878bar
    VThresholds[24][3] = 50.; AThresholds[24][3] = 2.0; DThresholds[24][3] = 40.; // ET
    VThresholds[25][3] = 35.; AThresholds[25][3] = 1.6; DThresholds[25][3] = 40.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 10.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 10.; // 878 - L0 left panel
    VThresholds[28][3] = 100.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 10.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 40.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 616 ) { // ET=1050, 878=750, 7725=900, B=3.8T.
    VThresholds[0][3] = 25.; AThresholds[0][3] = 1.8; DThresholds[0][3] = 40.; // 878bar
    VThresholds[1][3] = 40.; AThresholds[1][3] = 3.0; DThresholds[1][3] = 40.; // 878bar
    VThresholds[2][3] = 45.; AThresholds[2][3] = 3.5; DThresholds[2][3] = 40.; // 878bar
    VThresholds[3][3] = 40.; AThresholds[3][3] = 3.5; DThresholds[3][3] = 40.; // 878bar
    VThresholds[4][3] = 40.; AThresholds[4][3] = 3.5; DThresholds[4][3] = 40.; // 878bar
    VThresholds[5][3] = 90.; AThresholds[5][3] = 3.5; DThresholds[5][3] = 40.; // 7725
    VThresholds[6][3] = 30.; AThresholds[6][3] = 2.5; DThresholds[6][3] = 20.; // 878bar
    VThresholds[7][3] = 55.; AThresholds[7][3] = 4.0; DThresholds[7][3] = 40.; // 878bar
    VThresholds[8][3] = 12.; AThresholds[8][3] = 0.3; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 80.; AThresholds[9][3] = 5.0; DThresholds[9][3] = 40.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 10.; // 878 - L1 left panel
    VThresholds[12][3] = 40.; AThresholds[12][3] = 3.5; DThresholds[12][3] = 40.; // 878bar
    VThresholds[13][3] = 45.; AThresholds[13][3] = 3.5; DThresholds[13][3] = 40.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.4; DThresholds[14][3] = 20.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 40.; AThresholds[16][3] = 3.5; DThresholds[16][3] = 40.; // 878bar
    VThresholds[17][3] = 60.; AThresholds[17][3] = 3.5;  DThresholds[17][3] = 40.; // ET
    VThresholds[18][3] = 200.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 10.; // 878 - L1 right panel
    VThresholds[20][3] = 200.; AThresholds[20][3] = 18.; DThresholds[20][3] = 80.; // 878 - L1 slab
    VThresholds[21][3] = 320.; AThresholds[21][3] = 25.; DThresholds[21][3] = 80.; // 878 - Back slab
    VThresholds[22][3] = 140.; AThresholds[22][3] = 4.0; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 500.; AThresholds[23][3] = 70.; DThresholds[23][3] = 400.; // 878bar
    VThresholds[24][3] = 80.; AThresholds[24][3] = 3.0; DThresholds[24][3] = 40.; // ET
    VThresholds[25][3] = 50.; AThresholds[25][3] = 2.7; DThresholds[25][3] = 40.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 10.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 10.; // 878 - L0 left panel
    VThresholds[28][3] = 100.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 10.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 40.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 651) { // ET=1150, 878=800, 7725=800, B=3.8T.
    VThresholds[0][3] = 40.; AThresholds[0][3] = 1.8; DThresholds[0][3] = 40.; // 878bar
    VThresholds[1][3] = 60.; AThresholds[1][3] = 3.0; DThresholds[1][3] = 40.; // 878bar
    VThresholds[2][3] = 45.; AThresholds[2][3] = 3.5; DThresholds[2][3] = 40.; // 878bar
    VThresholds[3][3] = 40.; AThresholds[3][3] = 3.5; DThresholds[3][3] = 40.; // 878bar
    VThresholds[4][3] = 40.; AThresholds[4][3] = 3.5; DThresholds[4][3] = 40.; // 878bar
    VThresholds[5][3] = 40.; AThresholds[5][3] = 1.; DThresholds[5][3] = 40.; // 7725
    VThresholds[6][3] = 30.; AThresholds[6][3] = 2.5; DThresholds[6][3] = 20.; // 878bar
    VThresholds[7][3] = 55.; AThresholds[7][3] = 4.0; DThresholds[7][3] = 40.; // 878bar
    VThresholds[8][3] = 8.; AThresholds[8][3] = 0.3; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 80.; AThresholds[9][3] = 5.0; DThresholds[9][3] = 40.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 10.; // 878 - L1 left panel
    VThresholds[12][3] = 40.; AThresholds[12][3] = 3.5; DThresholds[12][3] = 40.; // 878bar
    VThresholds[13][3] = 45.; AThresholds[13][3] = 3.5; DThresholds[13][3] = 40.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.4; DThresholds[14][3] = 20.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 40.; AThresholds[16][3] = 3.5; DThresholds[16][3] = 40.; // 878bar
    VThresholds[17][3] = 60.; AThresholds[17][3] = 3.5;  DThresholds[17][3] = 40.; // ET
    VThresholds[18][3] = 200.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 10.; // 878 - L1 right panel
    VThresholds[20][3] = 200.; AThresholds[20][3] = 18.; DThresholds[20][3] = 80.; // 878 - L1 slab
    VThresholds[21][3] = 320.; AThresholds[21][3] = 25.; DThresholds[21][3] = 80.; // 878 - Back slab
    VThresholds[22][3] = 40.; AThresholds[22][3] = 1.0; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 500.; AThresholds[23][3] = 70.; DThresholds[23][3] = 400.; // 878bar
    VThresholds[24][3] = 80.; AThresholds[24][3] = 3.0; DThresholds[24][3] = 40.; // ET
    VThresholds[25][3] = 50.; AThresholds[25][3] = 2.5; DThresholds[25][3] = 40.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 10.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 10.; // 878 - L0 left panel
    VThresholds[28][3] = 100.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 10.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 40.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  // Not confirmed
  if (RunNum >= 670 ) { // ET=1200, 878=700, 7725=950, B=3.8T.
    VThresholds[0][3] = 8.; AThresholds[0][3] = 1.8; DThresholds[0][3] = 40.; // 878bar
    VThresholds[1][3] = 10.; AThresholds[1][3] = 3.0; DThresholds[1][3] = 40.; // 878bar
    VThresholds[2][3] = 10.; AThresholds[2][3] = 3.5; DThresholds[2][3] = 40.; // 878bar
    VThresholds[3][3] = 10.; AThresholds[3][3] = 3.5; DThresholds[3][3] = 40.; // 878bar
    VThresholds[4][3] = 10.; AThresholds[4][3] = 3.5; DThresholds[4][3] = 40.; // 878bar
    VThresholds[5][3] = 100.; AThresholds[5][3] = 3.5; DThresholds[5][3] = 40.; // 7725
    VThresholds[6][3] = 10.; AThresholds[6][3] = 2.5; DThresholds[6][3] = 20.; // 878bar
    VThresholds[7][3] = 10.; AThresholds[7][3] = 4.0; DThresholds[7][3] = 40.; // 878bar
    VThresholds[8][3] = 15.; AThresholds[8][3] = 0.3; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 300.; AThresholds[9][3] = 5.0; DThresholds[9][3] = 40.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 10.; // 878 - L1 left panel
    VThresholds[12][3] = 10.; AThresholds[12][3] = 3.5; DThresholds[12][3] = 40.; // 878bar
    VThresholds[13][3] = 10.; AThresholds[13][3] = 3.5; DThresholds[13][3] = 40.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.4; DThresholds[14][3] = 20.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 20.; AThresholds[16][3] = 3.5; DThresholds[16][3] = 40.; // 878bar
    VThresholds[17][3] = 120.; AThresholds[17][3] = 3.5;  DThresholds[17][3] = 40.; // ET
    VThresholds[18][3] = 200.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 10.; // 878 - L1 right panel
    VThresholds[20][3] = 200.; AThresholds[20][3] = 18.; DThresholds[20][3] = 80.; // 878 - L1 slab
    VThresholds[21][3] = 320.; AThresholds[21][3] = 25.; DThresholds[21][3] = 80.; // 878 - Back slab
    VThresholds[22][3] = 120.; AThresholds[22][3] = 4.0; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 500.; AThresholds[23][3] = 70.; DThresholds[23][3] = 400.; // 878bar
    VThresholds[24][3] = 300.; AThresholds[24][3] = 3.0; DThresholds[24][3] = 40.; // ET
    VThresholds[25][3] = 150.; AThresholds[25][3] = 2.5; DThresholds[25][3] = 40.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 10.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 10.; // 878 - L0 left panel
    VThresholds[28][3] = 100.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 10.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 40.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 725) { // ET=1150, 878=800, 7725=800, B=3.8T.
    VThresholds[0][3] = 40.; AThresholds[0][3] = 1.8; DThresholds[0][3] = 40.; // 878bar
    VThresholds[1][3] = 60.; AThresholds[1][3] = 3.0; DThresholds[1][3] = 40.; // 878bar
    VThresholds[2][3] = 45.; AThresholds[2][3] = 3.5; DThresholds[2][3] = 40.; // 878bar
    VThresholds[3][3] = 40.; AThresholds[3][3] = 3.5; DThresholds[3][3] = 40.; // 878bar
    VThresholds[4][3] = 40.; AThresholds[4][3] = 3.5; DThresholds[4][3] = 40.; // 878bar
    VThresholds[5][3] = 40.; AThresholds[5][3] = 1.; DThresholds[5][3] = 40.; // 7725
    VThresholds[6][3] = 30.; AThresholds[6][3] = 2.5; DThresholds[6][3] = 20.; // 878bar
    VThresholds[7][3] = 55.; AThresholds[7][3] = 4.0; DThresholds[7][3] = 40.; // 878bar
    VThresholds[8][3] = 8.; AThresholds[8][3] = 0.3; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 80.; AThresholds[9][3] = 5.0; DThresholds[9][3] = 40.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 10.; // 878 - L1 left panel
    VThresholds[12][3] = 40.; AThresholds[12][3] = 3.5; DThresholds[12][3] = 40.; // 878bar
    VThresholds[13][3] = 45.; AThresholds[13][3] = 3.5; DThresholds[13][3] = 40.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.4; DThresholds[14][3] = 20.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 40.; AThresholds[16][3] = 3.5; DThresholds[16][3] = 40.; // 878bar
    VThresholds[17][3] = 60.; AThresholds[17][3] = 3.5;  DThresholds[17][3] = 40.; // ET
    VThresholds[18][3] = 200.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 10.; // 878 - L1 right panel
    VThresholds[20][3] = 200.; AThresholds[20][3] = 18.; DThresholds[20][3] = 80.; // 878 - L1 slab
    VThresholds[21][3] = 320.; AThresholds[21][3] = 25.; DThresholds[21][3] = 80.; // 878 - Back slab
    VThresholds[22][3] = 40.; AThresholds[22][3] = 1.0; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 500.; AThresholds[23][3] = 70.; DThresholds[23][3] = 400.; // 878bar
    VThresholds[24][3] = 80.; AThresholds[24][3] = 3.0; DThresholds[24][3] = 40.; // ET
    VThresholds[25][3] = 50.; AThresholds[25][3] = 2.5; DThresholds[25][3] = 40.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 10.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 10.; // 878 - L0 left panel
    VThresholds[28][3] = 100.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 10.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 40.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 10.; // 878 - L2 left panel
  }

  if (RunNum >= 728) { // ET=1300, 878box1=1100V, L3TL_878=1100V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1450V, B=3.8T.
    VThresholds[0][3] = 270.; AThresholds[0][3] = 18.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 250.; AThresholds[1][3] = 27.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 320.; AThresholds[2][3] = 28.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 350.; AThresholds[3][3] = 30.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 260.; AThresholds[4][3] = 28.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 300.; AThresholds[6][3] = 22.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 380.; AThresholds[7][3] = 35.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 95.; AThresholds[8][3] = 4.; DThresholds[8][3] = 50.; // ET
    VThresholds[9][3] = 850.; AThresholds[9][3] = 40.; DThresholds[9][3] = 100.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 40.; // 878 - L1 left panel
    VThresholds[12][3] = 320.; AThresholds[12][3] = 26.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 300.; AThresholds[13][3] = 30.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.3; DThresholds[14][3] = 40.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 420.; AThresholds[16][3] = 30.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 620.; AThresholds[17][3] = 30.;  DThresholds[17][3] = 100.; // ET
    VThresholds[18][3] = 250.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 70.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 380.; AThresholds[21][3] = 25.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 14.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 70.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 780.; AThresholds[24][3] = 28.; DThresholds[24][3] = 100.; // ET
    VThresholds[25][3] = 500.; AThresholds[25][3] = 22.; DThresholds[25][3] = 100.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 60.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 30.; // 878 - L0 left panel
    VThresholds[28][3] = 90.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 50.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 60.; // 878 - L2 left panel
  }

  // Some changes in ET HV but no cosmic runs

  if (RunNum >= 802) { // ET=1150, 878box1=1200V, L3TL_878=1200V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1450V, B=3.8T.
    VThresholds[0][3] = 370.; AThresholds[0][3] = 20.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 300.; AThresholds[1][3] = 30.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 320.; AThresholds[2][3] = 30.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 350.; AThresholds[3][3] = 30.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 260.; AThresholds[4][3] = 28.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 300.; AThresholds[6][3] = 30.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 380.; AThresholds[7][3] = 50.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 15.; AThresholds[8][3] = 0.5; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 200.; AThresholds[9][3] = 4.; DThresholds[9][3] = 80.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 40.; // 878 - L1 left panel
    VThresholds[12][3] = 320.; AThresholds[12][3] = 26.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 300.; AThresholds[13][3] = 30.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.3; DThresholds[14][3] = 40.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 420.; AThresholds[16][3] = 30.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 150.; AThresholds[17][3] = 5.;  DThresholds[17][3] = 50.; // ET
    VThresholds[18][3] = 250.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 50.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 380.; AThresholds[21][3] = 25.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 10.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 50.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 200.; AThresholds[24][3] = 5.; DThresholds[24][3] = 60.; // ET
    VThresholds[25][3] = 100.; AThresholds[25][3] = 5.; DThresholds[25][3] = 60.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 60.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 30.; // 878 - L0 left panel
    VThresholds[28][3] = 90.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 50.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 60.; // 878 - L2 left panel
  }

  if (RunNum >= 804) { // ET=1100, 878box1=1250V, L3TL_878=1250V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1450V, B=3.8T.
    VThresholds[0][3] = 370.; AThresholds[0][3] = 20.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 350.; AThresholds[1][3] = 30.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 400.; AThresholds[2][3] = 30.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 400.; AThresholds[3][3] = 30.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 360.; AThresholds[4][3] = 28.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 380.; AThresholds[6][3] = 30.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 500.; AThresholds[7][3] = 50.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 10.; AThresholds[8][3] = 0.8; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 140.; AThresholds[9][3] = 6.; DThresholds[9][3] = 80.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 40.; // 878 - L1 left panel
    VThresholds[12][3] = 320.; AThresholds[12][3] = 26.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 300.; AThresholds[13][3] = 30.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.3; DThresholds[14][3] = 40.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 600.; AThresholds[16][3] = 30.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 100.; AThresholds[17][3] = 2.;  DThresholds[17][3] = 50.; // ET
    VThresholds[18][3] = 250.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 50.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 380.; AThresholds[21][3] = 25.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 10.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 50.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 90.; AThresholds[24][3] = 4.; DThresholds[24][3] = 60.; // ET
    VThresholds[25][3] = 80.; AThresholds[25][3] = 4.; DThresholds[25][3] = 60.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 60.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 30.; // 878 - L0 left panel
    VThresholds[28][3] = 90.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 50.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 60.; // 878 - L2 left panel
    VThresholds[0][4] = 370.; AThresholds[0][4] = 20.; DThresholds[0][4] = 80.; // 878bar
    VThresholds[1][4] = 350.; AThresholds[1][4] = 30.; DThresholds[1][4] = 100.; // 878bar
    VThresholds[2][4] = 400.; AThresholds[2][4] = 30.; DThresholds[2][4] = 100.; // 878bar
    VThresholds[3][4] = 400.; AThresholds[3][4] = 30.; DThresholds[3][4] = 100.; // 878bar
    VThresholds[4][4] = 360.; AThresholds[4][4] = 28.; DThresholds[4][4] = 100.; // 878bar
    VThresholds[5][4] = 700.; AThresholds[5][4] = 28.; DThresholds[5][4] = 80.; // 7725
    VThresholds[6][4] = 380.; AThresholds[6][4] = 30.; DThresholds[6][4] = 100.; // 878bar
    VThresholds[7][4] = 500.; AThresholds[7][4] = 50.; DThresholds[7][4] = 100.; // 878bar
    VThresholds[8][4] = 10.; AThresholds[8][4] = 0.8; DThresholds[8][4] = 40.; // ET
    VThresholds[9][4] = 140.; AThresholds[9][4] = 6.; DThresholds[9][4] = 80.; // ET
    VThresholds[10][4] = 10.; AThresholds[10][4] = 0.6; DThresholds[10][4] = 40.; // 878 - L0 top panel
    VThresholds[11][4] = 25.; AThresholds[11][4] = 1.5; DThresholds[11][4] = 40.; // 878 - L1 left panel
    VThresholds[12][4] = 320.; AThresholds[12][4] = 26.; DThresholds[12][4] = 100.; // 878bar
    VThresholds[13][4] = 300.; AThresholds[13][4] = 30.; DThresholds[13][4] = 100.; // 878bar
    VThresholds[14][4] = 10.; AThresholds[14][4] = 0.3; DThresholds[14][4] = 40.; // 878 - L2 top panel
    VThresholds[15][4] = 999.; AThresholds[15][4] = 999; DThresholds[15][4] = 10.; // LHC clock
    VThresholds[16][4] = 600.; AThresholds[16][4] = 30.; DThresholds[16][4] = 100.; // 878bar
    VThresholds[17][4] = 100.; AThresholds[17][4] = 2.;  DThresholds[17][4] = 50.; // ET
    VThresholds[18][4] = 250.; AThresholds[18][4] = 20.; DThresholds[18][4] = 80.; // 878 - L0 slab
    VThresholds[19][4] = 40.; AThresholds[19][4] = 3.0; DThresholds[19][4] = 50.; // 878 - L1 right panel
    VThresholds[20][4] = 250.; AThresholds[20][4] = 18.; DThresholds[20][4] = 100.; // 878 - L1 slab
    VThresholds[21][4] = 380.; AThresholds[21][4] = 25.; DThresholds[21][4] = 100.; // 878 - Back slab
    VThresholds[22][4] = 500.; AThresholds[22][4] = 10.; DThresholds[22][4] = 40.; // 7725
    VThresholds[23][4] = 580.; AThresholds[23][4] = 50.; DThresholds[23][4] = 200.; // 878bar
    VThresholds[24][4] = 90.; AThresholds[24][4] = 4.; DThresholds[24][4] = 60.; // ET
    VThresholds[25][4] = 80.; AThresholds[25][4] = 4.; DThresholds[25][4] = 60.; // ET
    VThresholds[26][4] = 30.; AThresholds[26][4] = 1.0; DThresholds[26][4] = 60.; // 878 - L2 right panel
    VThresholds[27][4] = 8.; AThresholds[27][4] = 0.4; DThresholds[27][4] = 30.; // 878 - L0 left panel
    VThresholds[28][4] = 90.; AThresholds[28][4] = 7.0; DThresholds[28][4] = 80.; // 878 - L2 slab
    VThresholds[29][4] = 6.; AThresholds[29][4] = 0.4; DThresholds[29][4] = 20.; // 878 - L0 right panel
    VThresholds[30][4] = 30.; AThresholds[30][4] = 1.0; DThresholds[30][4] = 50.; // 878 - L1 top panel
    VThresholds[31][4] = 25.; AThresholds[31][4] = 1.5; DThresholds[31][4] = 60.; // 878 - L2 left panel
  }

  if (RunNum >= 805) { // ET=1100, 878box1=1300V, L3TL_878=1300V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1450V, B=3.8T.
    VThresholds[0][3] = 450.; AThresholds[0][3] = 40.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 380.; AThresholds[1][3] = 40.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 500.; AThresholds[2][3] = 50.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 520.; AThresholds[3][3] = 50.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 400.; AThresholds[4][3] = 45.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 450.; AThresholds[6][3] = 45.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 580.; AThresholds[7][3] = 65.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 10.; AThresholds[8][3] = 0.8; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 140.; AThresholds[9][3] = 6.; DThresholds[9][3] = 80.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 40.; // 878 - L1 left panel
    VThresholds[12][3] = 550.; AThresholds[12][3] = 55.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 480.; AThresholds[13][3] = 55.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.3; DThresholds[14][3] = 40.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 650.; AThresholds[16][3] = 60.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 100.; AThresholds[17][3] = 2.;  DThresholds[17][3] = 50.; // ET
    VThresholds[18][3] = 250.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 50.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 380.; AThresholds[21][3] = 25.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 10.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 60.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 90.; AThresholds[24][3] = 4.; DThresholds[24][3] = 60.; // ET
    VThresholds[25][3] = 80.; AThresholds[25][3] = 4.; DThresholds[25][3] = 60.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 60.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 30.; // 878 - L0 left panel
    VThresholds[28][3] = 90.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 50.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 60.; // 878 - L2 left panel
    VThresholds[0][4] = 450.; AThresholds[0][4] = 40.; DThresholds[0][4] = 80.; // 878bar
    VThresholds[1][4] = 380.; AThresholds[1][4] = 40.; DThresholds[1][4] = 100.; // 878bar
    VThresholds[2][4] = 500.; AThresholds[2][4] = 50.; DThresholds[2][4] = 100.; // 878bar
    VThresholds[3][4] = 520.; AThresholds[3][4] = 50.; DThresholds[3][4] = 100.; // 878bar
    VThresholds[4][4] = 400.; AThresholds[4][4] = 28.; DThresholds[4][4] = 100.; // 878bar
    VThresholds[5][4] = 700.; AThresholds[5][4] = 28.; DThresholds[5][4] = 80.; // 7725
    VThresholds[6][4] = 450.; AThresholds[6][4] = 45.; DThresholds[6][4] = 100.; // 878bar
    VThresholds[7][4] = 580.; AThresholds[7][4] = 65.; DThresholds[7][4] = 100.; // 878bar
    VThresholds[8][4] = 10.; AThresholds[8][4] = 0.8; DThresholds[8][4] = 40.; // ET
    VThresholds[9][4] = 140.; AThresholds[9][4] = 6.; DThresholds[9][4] = 80.; // ET
    VThresholds[10][4] = 10.; AThresholds[10][4] = 0.6; DThresholds[10][4] = 40.; // 878 - L0 top panel
    VThresholds[11][4] = 25.; AThresholds[11][4] = 1.5; DThresholds[11][4] = 40.; // 878 - L1 left panel
    VThresholds[12][4] = 550.; AThresholds[12][4] = 55.; DThresholds[12][4] = 100.; // 878bar
    VThresholds[13][4] = 480.; AThresholds[13][4] = 55.; DThresholds[13][4] = 100.; // 878bar
    VThresholds[14][4] = 10.; AThresholds[14][4] = 0.3; DThresholds[14][4] = 40.; // 878 - L2 top panel
    VThresholds[15][4] = 999.; AThresholds[15][4] = 999; DThresholds[15][4] = 10.; // LHC clock
    VThresholds[16][4] = 650.; AThresholds[16][4] = 60.; DThresholds[16][4] = 100.; // 878bar
    VThresholds[17][4] = 100.; AThresholds[17][4] = 2.;  DThresholds[17][4] = 50.; // ET
    VThresholds[18][4] = 250.; AThresholds[18][4] = 20.; DThresholds[18][4] = 80.; // 878 - L0 slab
    VThresholds[19][4] = 40.; AThresholds[19][4] = 3.0; DThresholds[19][4] = 50.; // 878 - L1 right panel
    VThresholds[20][4] = 250.; AThresholds[20][4] = 18.; DThresholds[20][4] = 100.; // 878 - L1 slab
    VThresholds[21][4] = 380.; AThresholds[21][4] = 25.; DThresholds[21][4] = 100.; // 878 - Back slab
    VThresholds[22][4] = 500.; AThresholds[22][4] = 10.; DThresholds[22][4] = 40.; // 7725
    VThresholds[23][4] = 580.; AThresholds[23][4] = 60.; DThresholds[23][4] = 200.; // 878bar
    VThresholds[24][4] = 90.; AThresholds[24][4] = 4.; DThresholds[24][4] = 60.; // ET
    VThresholds[25][4] = 80.; AThresholds[25][4] = 4.; DThresholds[25][4] = 60.; // ET
    VThresholds[26][4] = 30.; AThresholds[26][4] = 1.0; DThresholds[26][4] = 60.; // 878 - L2 right panel
    VThresholds[27][4] = 8.; AThresholds[27][4] = 0.4; DThresholds[27][4] = 30.; // 878 - L0 left panel
    VThresholds[28][4] = 90.; AThresholds[28][4] = 7.0; DThresholds[28][4] = 80.; // 878 - L2 slab
    VThresholds[29][4] = 6.; AThresholds[29][4] = 0.4; DThresholds[29][4] = 20.; // 878 - L0 right panel
    VThresholds[30][4] = 30.; AThresholds[30][4] = 1.0; DThresholds[30][4] = 50.; // 878 - L1 top panel
    VThresholds[31][4] = 25.; AThresholds[31][4] = 1.5; DThresholds[31][4] = 60.; // 878 - L2 left panel
  }

  if (RunNum >= 806) { // ET=1100, 878box1=1350V, L3TL_878=1350V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1450V, B=3.8T.
    VThresholds[0][3] = 470.; AThresholds[0][3] = 40.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 400.; AThresholds[1][3] = 40.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 550.; AThresholds[2][3] = 50.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 580.; AThresholds[3][3] = 50.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 450.; AThresholds[4][3] = 45.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 500.; AThresholds[6][3] = 45.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 620.; AThresholds[7][3] = 65.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 10.; AThresholds[8][3] = 0.8; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 140.; AThresholds[9][3] = 6.; DThresholds[9][3] = 80.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 40.; // 878 - L1 left panel
    VThresholds[12][3] = 600.; AThresholds[12][3] = 55.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 550.; AThresholds[13][3] = 55.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.3; DThresholds[14][3] = 40.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 720.; AThresholds[16][3] = 60.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 100.; AThresholds[17][3] = 2.;  DThresholds[17][3] = 50.; // ET
    VThresholds[18][3] = 250.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 50.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 380.; AThresholds[21][3] = 25.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 10.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 60.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 90.; AThresholds[24][3] = 4.; DThresholds[24][3] = 60.; // ET
    VThresholds[25][3] = 80.; AThresholds[25][3] = 4.; DThresholds[25][3] = 60.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 60.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 30.; // 878 - L0 left panel
    VThresholds[28][3] = 90.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 50.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 60.; // 878 - L2 left panel
    VThresholds[0][4] = 470.; AThresholds[0][4] = 40.; DThresholds[0][4] = 80.; // 878bar
    VThresholds[1][4] = 400.; AThresholds[1][4] = 40.; DThresholds[1][4] = 100.; // 878bar
    VThresholds[2][4] = 550.; AThresholds[2][4] = 50.; DThresholds[2][4] = 100.; // 878bar
    VThresholds[3][4] = 580.; AThresholds[3][4] = 50.; DThresholds[3][4] = 100.; // 878bar
    VThresholds[4][4] = 450.; AThresholds[4][4] = 28.; DThresholds[4][4] = 100.; // 878bar
    VThresholds[5][4] = 700.; AThresholds[5][4] = 28.; DThresholds[5][4] = 80.; // 7725
    VThresholds[6][4] = 500.; AThresholds[6][4] = 45.; DThresholds[6][4] = 100.; // 878bar
    VThresholds[7][4] = 620.; AThresholds[7][4] = 65.; DThresholds[7][4] = 100.; // 878bar
    VThresholds[8][4] = 10.; AThresholds[8][4] = 0.8; DThresholds[8][4] = 40.; // ET
    VThresholds[9][4] = 140.; AThresholds[9][4] = 6.; DThresholds[9][4] = 80.; // ET
    VThresholds[10][4] = 10.; AThresholds[10][4] = 0.6; DThresholds[10][4] = 40.; // 878 - L0 top panel
    VThresholds[11][4] = 25.; AThresholds[11][4] = 1.5; DThresholds[11][4] = 40.; // 878 - L1 left panel
    VThresholds[12][4] = 600.; AThresholds[12][4] = 55.; DThresholds[12][4] = 100.; // 878bar
    VThresholds[13][4] = 550.; AThresholds[13][4] = 55.; DThresholds[13][4] = 100.; // 878bar
    VThresholds[14][4] = 10.; AThresholds[14][4] = 0.3; DThresholds[14][4] = 40.; // 878 - L2 top panel
    VThresholds[15][4] = 999.; AThresholds[15][4] = 999; DThresholds[15][4] = 10.; // LHC clock
    VThresholds[16][4] = 720.; AThresholds[16][4] = 60.; DThresholds[16][4] = 100.; // 878bar
    VThresholds[17][4] = 100.; AThresholds[17][4] = 2.;  DThresholds[17][4] = 50.; // ET
    VThresholds[18][4] = 250.; AThresholds[18][4] = 20.; DThresholds[18][4] = 80.; // 878 - L0 slab
    VThresholds[19][4] = 40.; AThresholds[19][4] = 3.0; DThresholds[19][4] = 50.; // 878 - L1 right panel
    VThresholds[20][4] = 250.; AThresholds[20][4] = 18.; DThresholds[20][4] = 100.; // 878 - L1 slab
    VThresholds[21][4] = 380.; AThresholds[21][4] = 25.; DThresholds[21][4] = 100.; // 878 - Back slab
    VThresholds[22][4] = 500.; AThresholds[22][4] = 10.; DThresholds[22][4] = 40.; // 7725
    VThresholds[23][4] = 580.; AThresholds[23][4] = 60.; DThresholds[23][4] = 200.; // 878bar
    VThresholds[24][4] = 90.; AThresholds[24][4] = 4.; DThresholds[24][4] = 60.; // ET
    VThresholds[25][4] = 80.; AThresholds[25][4] = 4.; DThresholds[25][4] = 60.; // ET
    VThresholds[26][4] = 30.; AThresholds[26][4] = 1.0; DThresholds[26][4] = 60.; // 878 - L2 right panel
    VThresholds[27][4] = 8.; AThresholds[27][4] = 0.4; DThresholds[27][4] = 30.; // 878 - L0 left panel
    VThresholds[28][4] = 90.; AThresholds[28][4] = 7.0; DThresholds[28][4] = 80.; // 878 - L2 slab
    VThresholds[29][4] = 6.; AThresholds[29][4] = 0.4; DThresholds[29][4] = 20.; // 878 - L0 right panel
    VThresholds[30][4] = 30.; AThresholds[30][4] = 1.0; DThresholds[30][4] = 50.; // 878 - L1 top panel
    VThresholds[31][4] = 25.; AThresholds[31][4] = 1.5; DThresholds[31][4] = 60.; // 878 - L2 left panel
  }

  if (RunNum >= 815) { // ET=1100, 878box1=1400V, L3TL_878=1400V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1450V, B=3.8T.
    VThresholds[0][3] = 520.; AThresholds[0][3] = 50.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 450.; AThresholds[1][3] = 60.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 600.; AThresholds[2][3] = 60.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 620.; AThresholds[3][3] = 70.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 500.; AThresholds[4][3] = 60.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 550.; AThresholds[6][3] = 60.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 680.; AThresholds[7][3] = 80.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 10.; AThresholds[8][3] = 0.8; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 140.; AThresholds[9][3] = 6.; DThresholds[9][3] = 80.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 40.; // 878 - L1 left panel
    VThresholds[12][3] = 650.; AThresholds[12][3] = 70.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 580.; AThresholds[13][3] = 70.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.3; DThresholds[14][3] = 40.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 800.; AThresholds[16][3] = 80.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 100.; AThresholds[17][3] = 2.;  DThresholds[17][3] = 50.; // ET
    VThresholds[18][3] = 250.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 50.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 380.; AThresholds[21][3] = 25.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 10.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 60.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 90.; AThresholds[24][3] = 4.; DThresholds[24][3] = 60.; // ET
    VThresholds[25][3] = 80.; AThresholds[25][3] = 4.; DThresholds[25][3] = 60.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 60.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 30.; // 878 - L0 left panel
    VThresholds[28][3] = 90.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 50.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 60.; // 878 - L2 left panel
    VThresholds[0][4] = 520.; AThresholds[0][4] = 50.; DThresholds[0][4] = 80.; // 878bar
    VThresholds[1][4] = 450.; AThresholds[1][4] = 60.; DThresholds[1][4] = 100.; // 878bar
    VThresholds[2][4] = 600.; AThresholds[2][4] = 60.; DThresholds[2][4] = 100.; // 878bar
    VThresholds[3][4] = 620.; AThresholds[3][4] = 70.; DThresholds[3][4] = 100.; // 878bar
    VThresholds[4][4] = 500.; AThresholds[4][4] = 60.; DThresholds[4][4] = 100.; // 878bar
    VThresholds[5][4] = 700.; AThresholds[5][4] = 28.; DThresholds[5][4] = 80.; // 7725
    VThresholds[6][4] = 550.; AThresholds[6][4] = 60.; DThresholds[6][4] = 100.; // 878bar
    VThresholds[7][4] = 680.; AThresholds[7][4] = 80.; DThresholds[7][4] = 100.; // 878bar
    VThresholds[8][4] = 10.; AThresholds[8][4] = 0.8; DThresholds[8][4] = 40.; // ET
    VThresholds[9][4] = 140.; AThresholds[9][4] = 6.; DThresholds[9][4] = 80.; // ET
    VThresholds[10][4] = 10.; AThresholds[10][4] = 0.6; DThresholds[10][4] = 40.; // 878 - L0 top panel
    VThresholds[11][4] = 25.; AThresholds[11][4] = 1.5; DThresholds[11][4] = 40.; // 878 - L1 left panel
    VThresholds[12][4] = 650.; AThresholds[12][4] = 70.; DThresholds[12][4] = 100.; // 878bar
    VThresholds[13][4] = 580.; AThresholds[13][4] = 70.; DThresholds[13][4] = 100.; // 878bar
    VThresholds[14][4] = 10.; AThresholds[14][4] = 0.3; DThresholds[14][4] = 40.; // 878 - L2 top panel
    VThresholds[15][4] = 999.; AThresholds[15][4] = 999; DThresholds[15][4] = 10.; // LHC clock
    VThresholds[16][4] = 800.; AThresholds[16][4] = 80.; DThresholds[16][4] = 100.; // 878bar
    VThresholds[17][4] = 100.; AThresholds[17][4] = 2.;  DThresholds[17][4] = 50.; // ET
    VThresholds[18][4] = 250.; AThresholds[18][4] = 20.; DThresholds[18][4] = 80.; // 878 - L0 slab
    VThresholds[19][4] = 40.; AThresholds[19][4] = 3.0; DThresholds[19][4] = 50.; // 878 - L1 right panel
    VThresholds[20][4] = 250.; AThresholds[20][4] = 18.; DThresholds[20][4] = 100.; // 878 - L1 slab
    VThresholds[21][4] = 380.; AThresholds[21][4] = 25.; DThresholds[21][4] = 100.; // 878 - Back slab
    VThresholds[22][4] = 500.; AThresholds[22][4] = 10.; DThresholds[22][4] = 40.; // 7725
    VThresholds[23][4] = 580.; AThresholds[23][4] = 60.; DThresholds[23][4] = 200.; // 878bar
    VThresholds[24][4] = 90.; AThresholds[24][4] = 4.; DThresholds[24][4] = 60.; // ET
    VThresholds[25][4] = 80.; AThresholds[25][4] = 4.; DThresholds[25][4] = 60.; // ET
    VThresholds[26][4] = 30.; AThresholds[26][4] = 1.0; DThresholds[26][4] = 60.; // 878 - L2 right panel
    VThresholds[27][4] = 8.; AThresholds[27][4] = 0.4; DThresholds[27][4] = 30.; // 878 - L0 left panel
    VThresholds[28][4] = 90.; AThresholds[28][4] = 7.0; DThresholds[28][4] = 80.; // 878 - L2 slab
    VThresholds[29][4] = 6.; AThresholds[29][4] = 0.4; DThresholds[29][4] = 20.; // 878 - L0 right panel
    VThresholds[30][4] = 30.; AThresholds[30][4] = 1.0; DThresholds[30][4] = 50.; // 878 - L1 top panel
    VThresholds[31][4] = 25.; AThresholds[31][4] = 1.5; DThresholds[31][4] = 60.; // 878 - L2 left panel
  }

  if (RunNum >= 816) { // ET=1100, 878box1=1450V, L3TL_878=1450V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1450V, B=3.8T.
    VThresholds[0][3] = 600.; AThresholds[0][3] = 60.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 500.; AThresholds[1][3] = 65.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 680.; AThresholds[2][3] = 70.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 700.; AThresholds[3][3] = 80.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 550.; AThresholds[4][3] = 75.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 580.; AThresholds[6][3] = 65.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 720.; AThresholds[7][3] = 90.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 10.; AThresholds[8][3] = 0.8; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 140.; AThresholds[9][3] = 6.; DThresholds[9][3] = 80.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 25.; AThresholds[11][3] = 1.5; DThresholds[11][3] = 40.; // 878 - L1 left panel
    VThresholds[12][3] = 700.; AThresholds[12][3] = 80.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 650.; AThresholds[13][3] = 80.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.3; DThresholds[14][3] = 40.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 880.; AThresholds[16][3] = 90.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 100.; AThresholds[17][3] = 2.;  DThresholds[17][3] = 50.; // ET
    VThresholds[18][3] = 250.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 40.; AThresholds[19][3] = 3.0; DThresholds[19][3] = 50.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 380.; AThresholds[21][3] = 25.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 10.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 60.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 90.; AThresholds[24][3] = 4.; DThresholds[24][3] = 60.; // ET
    VThresholds[25][3] = 80.; AThresholds[25][3] = 4.; DThresholds[25][3] = 60.; // ET
    VThresholds[26][3] = 30.; AThresholds[26][3] = 1.0; DThresholds[26][3] = 60.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.4; DThresholds[27][3] = 30.; // 878 - L0 left panel
    VThresholds[28][3] = 90.; AThresholds[28][3] = 7.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.4; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 30.; AThresholds[30][3] = 1.0; DThresholds[30][3] = 50.; // 878 - L1 top panel
    VThresholds[31][3] = 25.; AThresholds[31][3] = 1.5; DThresholds[31][3] = 60.; // 878 - L2 left panel
    VThresholds[0][4] = 600.; AThresholds[0][4] = 60.; DThresholds[0][4] = 80.; // 878bar
    VThresholds[1][4] = 500.; AThresholds[1][4] = 65.; DThresholds[1][4] = 100.; // 878bar
    VThresholds[2][4] = 680.; AThresholds[2][4] = 70.; DThresholds[2][4] = 100.; // 878bar
    VThresholds[3][4] = 700.; AThresholds[3][4] = 80.; DThresholds[3][4] = 100.; // 878bar
    VThresholds[4][4] = 550.; AThresholds[4][4] = 75.; DThresholds[4][4] = 100.; // 878bar
    VThresholds[5][4] = 700.; AThresholds[5][4] = 28.; DThresholds[5][4] = 80.; // 7725
    VThresholds[6][4] = 580.; AThresholds[6][4] = 65.; DThresholds[6][4] = 100.; // 878bar
    VThresholds[7][4] = 720.; AThresholds[7][4] = 90.; DThresholds[7][4] = 100.; // 878bar
    VThresholds[8][4] = 10.; AThresholds[8][4] = 0.8; DThresholds[8][4] = 40.; // ET
    VThresholds[9][4] = 140.; AThresholds[9][4] = 6.; DThresholds[9][4] = 80.; // ET
    VThresholds[10][4] = 10.; AThresholds[10][4] = 0.6; DThresholds[10][4] = 40.; // 878 - L0 top panel
    VThresholds[11][4] = 25.; AThresholds[11][4] = 1.5; DThresholds[11][4] = 40.; // 878 - L1 left panel
    VThresholds[12][4] = 700.; AThresholds[12][4] = 80.; DThresholds[12][4] = 100.; // 878bar
    VThresholds[13][4] = 650.; AThresholds[13][4] = 80.; DThresholds[13][4] = 100.; // 878bar
    VThresholds[14][4] = 10.; AThresholds[14][4] = 0.3; DThresholds[14][4] = 40.; // 878 - L2 top panel
    VThresholds[15][4] = 999.; AThresholds[15][4] = 999; DThresholds[15][4] = 10.; // LHC clock
    VThresholds[16][4] = 880.; AThresholds[16][4] = 90.; DThresholds[16][4] = 100.; // 878bar
    VThresholds[17][4] = 100.; AThresholds[17][4] = 2.;  DThresholds[17][4] = 50.; // ET
    VThresholds[18][4] = 250.; AThresholds[18][4] = 20.; DThresholds[18][4] = 80.; // 878 - L0 slab
    VThresholds[19][4] = 40.; AThresholds[19][4] = 3.0; DThresholds[19][4] = 50.; // 878 - L1 right panel
    VThresholds[20][4] = 250.; AThresholds[20][4] = 18.; DThresholds[20][4] = 100.; // 878 - L1 slab
    VThresholds[21][4] = 380.; AThresholds[21][4] = 25.; DThresholds[21][4] = 100.; // 878 - Back slab
    VThresholds[22][4] = 500.; AThresholds[22][4] = 10.; DThresholds[22][4] = 40.; // 7725
    VThresholds[23][4] = 580.; AThresholds[23][4] = 60.; DThresholds[23][4] = 200.; // 878bar
    VThresholds[24][4] = 90.; AThresholds[24][4] = 4.; DThresholds[24][4] = 60.; // ET
    VThresholds[25][4] = 80.; AThresholds[25][4] = 4.; DThresholds[25][4] = 60.; // ET
    VThresholds[26][4] = 30.; AThresholds[26][4] = 1.0; DThresholds[26][4] = 60.; // 878 - L2 right panel
    VThresholds[27][4] = 8.; AThresholds[27][4] = 0.4; DThresholds[27][4] = 30.; // 878 - L0 left panel
    VThresholds[28][4] = 90.; AThresholds[28][4] = 7.0; DThresholds[28][4] = 80.; // 878 - L2 slab
    VThresholds[29][4] = 6.; AThresholds[29][4] = 0.4; DThresholds[29][4] = 20.; // 878 - L0 right panel
    VThresholds[30][4] = 30.; AThresholds[30][4] = 1.0; DThresholds[30][4] = 50.; // 878 - L1 top panel
    VThresholds[31][4] = 25.; AThresholds[31][4] = 1.5; DThresholds[31][4] = 60.; // 878 - L2 left panel
  }

  if (RunNum >= 816) { // ET=1100, 878box1=1450V, L3TL_878=1450V, L3BR_7725=1120V, L3ML_7725=1030V, Bars&Slabs=1400V, B=3.8T.
    VThresholds[0][3] = 600.; AThresholds[0][3] = 60.; DThresholds[0][3] = 80.; // 878bar
    VThresholds[1][3] = 500.; AThresholds[1][3] = 65.; DThresholds[1][3] = 100.; // 878bar
    VThresholds[2][3] = 680.; AThresholds[2][3] = 70.; DThresholds[2][3] = 100.; // 878bar
    VThresholds[3][3] = 700.; AThresholds[3][3] = 80.; DThresholds[3][3] = 100.; // 878bar
    VThresholds[4][3] = 550.; AThresholds[4][3] = 75.; DThresholds[4][3] = 100.; // 878bar
    VThresholds[5][3] = 700.; AThresholds[5][3] = 28.; DThresholds[5][3] = 80.; // 7725
    VThresholds[6][3] = 580.; AThresholds[6][3] = 65.; DThresholds[6][3] = 100.; // 878bar
    VThresholds[7][3] = 720.; AThresholds[7][3] = 90.; DThresholds[7][3] = 100.; // 878bar
    VThresholds[8][3] = 10.; AThresholds[8][3] = 0.8; DThresholds[8][3] = 40.; // ET
    VThresholds[9][3] = 140.; AThresholds[9][3] = 6.; DThresholds[9][3] = 80.; // ET
    VThresholds[10][3] = 8.; AThresholds[10][3] = 0.1; DThresholds[10][3] = 20.; // 878 - L0 top panel
    VThresholds[11][3] = 8.; AThresholds[11][3] = 0.1; DThresholds[11][3] = 20.; // 878 - L1 left panel
    VThresholds[12][3] = 700.; AThresholds[12][3] = 80.; DThresholds[12][3] = 100.; // 878bar
    VThresholds[13][3] = 650.; AThresholds[13][3] = 80.; DThresholds[13][3] = 100.; // 878bar
    VThresholds[14][3] = 8.; AThresholds[14][3] = 0.1; DThresholds[14][3] = 20.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 880.; AThresholds[16][3] = 90.; DThresholds[16][3] = 100.; // 878bar
    VThresholds[17][3] = 100.; AThresholds[17][3] = 2.;  DThresholds[17][3] = 50.; // ET
    VThresholds[18][3] = 200.; AThresholds[18][3] = 10.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 8.; AThresholds[19][3] = 0.1; DThresholds[19][3] = 20.; // 878 - L1 right panel
    VThresholds[20][3] = 200.; AThresholds[20][3] = 10.; DThresholds[20][3] = 100.; // 878 - L1 slab
    VThresholds[21][3] = 300.; AThresholds[21][3] = 10.; DThresholds[21][3] = 100.; // 878 - Back slab
    VThresholds[22][3] = 500.; AThresholds[22][3] = 10.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 580.; AThresholds[23][3] = 60.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 90.; AThresholds[24][3] = 4.; DThresholds[24][3] = 60.; // ET
    VThresholds[25][3] = 80.; AThresholds[25][3] = 4.; DThresholds[25][3] = 60.; // ET
    VThresholds[26][3] = 8.; AThresholds[26][3] = 0.1; DThresholds[26][3] = 20.; // 878 - L2 right panel
    VThresholds[27][3] = 8.; AThresholds[27][3] = 0.1; DThresholds[27][3] = 20.; // 878 - L0 left panel
    VThresholds[28][3] = 60.; AThresholds[28][3] = 4.0; DThresholds[28][3] = 80.; // 878 - L2 slab
    VThresholds[29][3] = 6.; AThresholds[29][3] = 0.1; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 8.; AThresholds[30][3] = 0.1; DThresholds[30][3] = 20.; // 878 - L1 top panel
    VThresholds[31][3] = 8.; AThresholds[31][3] = 0.1; DThresholds[31][3] = 20.; // 878 - L2 left panel
    VThresholds[0][4] = 600.; AThresholds[0][4] = 60.; DThresholds[0][4] = 80.; // 878bar
    VThresholds[1][4] = 500.; AThresholds[1][4] = 65.; DThresholds[1][4] = 100.; // 878bar
    VThresholds[2][4] = 680.; AThresholds[2][4] = 70.; DThresholds[2][4] = 100.; // 878bar
    VThresholds[3][4] = 700.; AThresholds[3][4] = 80.; DThresholds[3][4] = 100.; // 878bar
    VThresholds[4][4] = 550.; AThresholds[4][4] = 75.; DThresholds[4][4] = 100.; // 878bar
    VThresholds[5][4] = 700.; AThresholds[5][4] = 28.; DThresholds[5][4] = 80.; // 7725
    VThresholds[6][4] = 580.; AThresholds[6][4] = 65.; DThresholds[6][4] = 100.; // 878bar
    VThresholds[7][4] = 720.; AThresholds[7][4] = 90.; DThresholds[7][4] = 100.; // 878bar
    VThresholds[8][4] = 10.; AThresholds[8][4] = 0.8; DThresholds[8][4] = 40.; // ET
    VThresholds[9][4] = 140.; AThresholds[9][4] = 6.; DThresholds[9][4] = 80.; // ET
    VThresholds[10][4] = 8.; AThresholds[10][4] = 0.1; DThresholds[10][4] = 20.; // 878 - L0 top panel
    VThresholds[11][4] = 8.; AThresholds[11][4] = 0.1; DThresholds[11][4] = 20.; // 878 - L1 left panel
    VThresholds[12][4] = 700.; AThresholds[12][4] = 80.; DThresholds[12][4] = 100.; // 878bar
    VThresholds[13][4] = 650.; AThresholds[13][4] = 80.; DThresholds[13][4] = 100.; // 878bar
    VThresholds[14][4] = 8.; AThresholds[14][4] = 0.1; DThresholds[14][4] = 20.; // 878 - L2 top panel
    VThresholds[15][4] = 999.; AThresholds[15][4] = 999; DThresholds[15][4] = 10.; // LHC clock
    VThresholds[16][4] = 880.; AThresholds[16][4] = 90.; DThresholds[16][4] = 100.; // 878bar
    VThresholds[17][4] = 100.; AThresholds[17][4] = 2.;  DThresholds[17][4] = 50.; // ET
    VThresholds[18][4] = 200.; AThresholds[18][4] = 10.; DThresholds[18][4] = 80.; // 878 - L0 slab
    VThresholds[19][4] = 8.; AThresholds[19][4] = 0.1; DThresholds[19][4] = 20.; // 878 - L1 right panel
    VThresholds[20][4] = 200.; AThresholds[20][4] = 10.; DThresholds[20][4] = 100.; // 878 - L1 slab
    VThresholds[21][4] = 300.; AThresholds[21][4] = 10.; DThresholds[21][4] = 100.; // 878 - Back slab
    VThresholds[22][4] = 500.; AThresholds[22][4] = 10.; DThresholds[22][4] = 40.; // 7725
    VThresholds[23][4] = 580.; AThresholds[23][4] = 60.; DThresholds[23][4] = 200.; // 878bar
    VThresholds[24][4] = 90.; AThresholds[24][4] = 4.; DThresholds[24][4] = 60.; // ET
    VThresholds[25][4] = 80.; AThresholds[25][4] = 4.; DThresholds[25][4] = 60.; // ET
    VThresholds[26][4] = 8.; AThresholds[26][4] = 0.1; DThresholds[26][4] = 20.; // 878 - L2 right panel
    VThresholds[27][4] = 8.; AThresholds[27][4] = 0.1; DThresholds[27][4] = 20.; // 878 - L0 left panel
    VThresholds[28][4] = 60.; AThresholds[28][4] = 4.0; DThresholds[28][4] = 80.; // 878 - L2 slab
    VThresholds[29][4] = 6.; AThresholds[29][4] = 0.1; DThresholds[29][4] = 20.; // 878 - L0 right panel
    VThresholds[30][4] = 8.; AThresholds[30][4] = 0.1; DThresholds[30][4] = 20.; // 878 - L1 top panel
    VThresholds[31][4] = 8.; AThresholds[31][4] = 1.5; DThresholds[31][4] = 20.; // 878 - L2 left panel
  }

  if (RunNum >= 931) { // ET=1450, All878s=1450V, L3BR_7725=1220V, L3ML_7725=1130V, B=3.8T.
    VThresholds[0][3] = 500.; AThresholds[0][3] = 70.; DThresholds[0][3] = 200.; // 878bar
    VThresholds[1][3] = 600.; AThresholds[1][3] = 60.; DThresholds[1][3] = 200.; // 878bar
    VThresholds[2][3] = 670.; AThresholds[2][3] = 75.; DThresholds[2][3] = 200.; // 878bar
    VThresholds[3][3] = 700.; AThresholds[3][3] = 82.; DThresholds[3][3] = 200.; // 878bar
    VThresholds[4][3] = 570.; AThresholds[4][3] = 80.; DThresholds[4][3] = 200.; // 878bar
    VThresholds[5][3] = 1100.; AThresholds[5][3] = 50.; DThresholds[5][3] = 100.; // 7725
    VThresholds[6][3] = 600.; AThresholds[6][3] = 65.; DThresholds[6][3] = 200.; // 878bar
    VThresholds[7][3] = 750.; AThresholds[7][3] = 95.; DThresholds[7][3] = 200.; // 878bar
    VThresholds[8][3] = 750.; AThresholds[8][3] = 60.; DThresholds[8][3] = 200.; // R878bar (new)
    VThresholds[9][3] = 1100.; AThresholds[9][3] = 85.; DThresholds[9][3] = 200.; // ET
    VThresholds[10][3] = 10.; AThresholds[10][3] = 0.6; DThresholds[10][3] = 40.; // 878 - L0 top panel
    VThresholds[11][3] = 10.; AThresholds[11][3] = 1.0; DThresholds[11][3] = 40.; // 878 - L1 top panel
    VThresholds[12][3] = 720.; AThresholds[12][3] = 80.; DThresholds[12][3] = 200.; // 878bar
    VThresholds[13][3] = 650.; AThresholds[13][3] = 80.; DThresholds[13][3] = 200.; // 878bar
    VThresholds[14][3] = 10.; AThresholds[14][3] = 0.5; DThresholds[14][3] = 30.; // 878 - L2 top panel
    VThresholds[15][3] = 999.; AThresholds[15][3] = 999; DThresholds[15][3] = 10.; // LHC clock
    VThresholds[16][3] = 850.; AThresholds[16][3] = 90.; DThresholds[16][3] = 200.; // 878bar
    VThresholds[17][3] = 1100.; AThresholds[17][3] = 70.;  DThresholds[17][3] = 150.; // ET
    VThresholds[18][3] = 200.; AThresholds[18][3] = 20.; DThresholds[18][3] = 80.; // 878 - L0 slab
    VThresholds[19][3] = 20.; AThresholds[19][3] = 1.0; DThresholds[19][3] = 20.; // 878 - L1 right panel
    VThresholds[20][3] = 250.; AThresholds[20][3] = 18.; DThresholds[20][3] = 80.; // 878 - L1 slab
    VThresholds[21][3] = 300.; AThresholds[21][3] = 15.; DThresholds[21][3] = 80.; // 878 - Back slab
    VThresholds[22][3] = 1100.; AThresholds[22][3] = 25.; DThresholds[22][3] = 40.; // 7725
    VThresholds[23][3] = 560.; AThresholds[23][3] = 65.; DThresholds[23][3] = 200.; // 878bar
    VThresholds[24][3] = 1100.; AThresholds[24][3] = 70.; DThresholds[24][3] = 20.; // ET
    VThresholds[25][3] = 1100.; AThresholds[25][3] = 60.; DThresholds[25][3] = 20.; // ET
    VThresholds[26][3] = 20.; AThresholds[26][3] = 1.; DThresholds[26][3] = 20.; // 878 - L2 right panel
    VThresholds[27][3] = 20.; AThresholds[27][3] = 1.; DThresholds[27][3] = 20.; // 878 - L0 left panel
    VThresholds[28][3] = 140.; AThresholds[28][3] = 8.0; DThresholds[28][3] = 40.; // 878 - L2 slab
    VThresholds[29][3] = 8.; AThresholds[29][3] = 0.2; DThresholds[29][3] = 20.; // 878 - L0 right panel
    VThresholds[30][3] = 20.; AThresholds[30][3] = 1.; DThresholds[30][3] = 20.; // 878 - L1 left panel
    VThresholds[31][3] = 20.; AThresholds[31][3] = 1.; DThresholds[31][3] = 20.; // 878 - L2 left panel
    for(int c=0; c<32; c++) {
      VThresholds[c][4] = VThresholds[c][3];
      AThresholds[c][4] = AThresholds[c][3];
      DThresholds[c][4] = DThresholds[c][3];
    }
  }

  // // Read the thresholds and other settings from the data base file
  // string DBfile = getenv("MILLIQANDB");
  // FILE* f;
  // f = fopen(DBfile.c_str(),"r");
  // if ( !f ) {
  //   cerr << "I could not open the database file "<<DBfile<<" to read settings.\n";
  //   return;
  // }
  // char word[400];  // The current line of input in character format
  // std::string Word; // The current line of input in string format
  // int ncol = 1;
  // while (ncol>=0) { // Until end of file
  //   ncol = fscanf(f,"%s",word);
  //   Word = word;
  //   if (debug>0) cout << ":"<<Word<<":\n";
  //   tmps = Word.erase(14);
  //   if (tmps != "LECROYWR104MXi") {
  // } // Until end of file
    
  if (debug>100) cout << "Initializing Chain for run number " << RunNum<<".\n";
  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();
  cout<<"The number of entries is: "<<nentries<<endl;

  TFile *outFile = new TFile(output_filename.Data(),"UPDATE");
  outFile->cd();

  //Histogram Declarations
  if (debug>100) cout << "Booking histograms.\n";

  // Bins and ranges
  int nBins; float minX, maxX;
  
  // Histogram names are Cuts_LHCStatus_Variable_Chan_Range
  
  // Plots of the max sample height, with several different range choices
  if (RangeCode == "A") { nBins = 140; minX = 0.; maxX = 1400.; }
  if (RangeCode == "B") { nBins = 50; minX = 0.; maxX = 500.; }
  if (RangeCode == "C") { nBins = 50; minX = 0.; maxX = 100.; }
  if (RangeCode == "D") { nBins = 40; minX = 0.; maxX = 40.; }
  TH1D *h_Max[32];
  for (int c=0; c<32; c++)
   h_Max[c] = new TH1D(EventCategory+"_"+LHCStatus+"_Max_"+ChStr[c]+"_"+RangeCode,"Max sample in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the after pulse area (1/1000) per channel, with several different range choices
  if (RangeCode == "A") { nBins = 200; minX = 0.; maxX = 100.; }
  if (RangeCode == "B") { nBins = 200; minX = 0.; maxX = 50.; }
  if (RangeCode == "C") { nBins = 50; minX = 0.; maxX = 2.; }
  if (RangeCode == "D") { nBins = 50; minX = 0.; maxX = 0.2; }
  TH1D *h_AfterPulseArea[32];
  for (int c=0; c<32; c++)
    h_AfterPulseArea[c] = new TH1D(EventCategory+"_"+LHCStatus+"_AfterPulseArea_"+ChStr[c]+"_"+RangeCode,"Area of afterpulses in Channel "+ChStr[c],nBins,minX,maxX);
 
  // Plots of the after pulse height per channel, with several different range choices
  if (RangeCode == "A") { nBins = 60; minX = 0.; maxX = 240.; }
  if (RangeCode == "B") { nBins = 60; minX = 0.; maxX = 60.; }
  if (RangeCode == "C") { nBins = 40; minX = 0.; maxX = 20.; }
  if (RangeCode == "D") { nBins = 45; minX = 0.; maxX = 15.; }
  TH1D *h_AfterPulseHeight[32];
  for (int c=0; c<32; c++)
    h_AfterPulseHeight[c] = new TH1D(EventCategory+"_"+LHCStatus+"_AfterPulseHeight_"+ChStr[c]+"_"+RangeCode,"Height of afterpulses in Channel "+ChStr[c],nBins,minX,maxX);
   
  // Plots of the after pulse duration per channel, with several different range choices
  if (RangeCode == "A") { nBins = 128; minX = 0.; maxX = 640.; }
  if (RangeCode == "B") { nBins = 64; minX = 0.; maxX = 160.; }
  if (RangeCode == "C") { nBins = 32; minX = 0.; maxX = 40.; }
  if (RangeCode == "D") { nBins = 32; minX = 0.; maxX = 20.; }
  TH1D *h_AfterPulseDuration[32];
  for (int c=0; c<32; c++)
    h_AfterPulseDuration[c] = new TH1D(EventCategory+"_"+LHCStatus+"_AfterPulseDuration_"+ChStr[c]+"_"+RangeCode,"Duration of afterpulses in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the first pulse area (1/1000) per channel, with several different range choices
  if (RangeCode == "A") { nBins = 400; minX = 0.; maxX = 200.; }
  if (RangeCode == "B") { nBins = 400; minX = 0.; maxX = 100.; }
  if (RangeCode == "C") { nBins = 400; minX = 0.; maxX = 25.; }
  if (RangeCode == "D") { nBins = 200; minX = 0.; maxX = 5.; }
  TH1D *h_FirstPulseArea[32];
  for (int c=0; c<32; c++)
    h_FirstPulseArea[c] = new TH1D(EventCategory+"_"+LHCStatus+"_FirstPulseArea_"+ChStr[c]+"_"+RangeCode,"Area of first pulse in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the afterpulse time, wrt first pulse, per channel, with several different range choices
  if (RangeCode == "A") { nBins = 1024; minX = 0.; maxX = 2*1024*0.625; }
  if (RangeCode == "B") { nBins = 512; minX = 0.; maxX = 2*1024*0.625; }
  if (RangeCode == "C") { nBins = 128; minX = 0.; maxX = 1024*0.625; }
  if (RangeCode == "D") { nBins = 100; minX = 100.; maxX = 300.; }
  TH1D *h_AfterPulseTime[32];
  for (int c=0; c<32; c++)
    h_AfterPulseTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_AfterPulseTime_"+ChStr[c]+"_"+RangeCode,"Time delay of afterpulses in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the first pulse height per channel, with several different range choices
  if (RangeCode == "A") { nBins = 260; minX = 0.; maxX = 1300.; }
  if (RangeCode == "B") { nBins = 200; minX = 0.; maxX = 500.; }
  if (RangeCode == "C") { nBins = 150; minX = 0.; maxX = 150.; }
  if (RangeCode == "D") { nBins = 120; minX = 0.; maxX = 40.; }
  TH1D *h_FirstPulseHeight[32];
  for (int c=0; c<32; c++)
   h_FirstPulseHeight[c] = new TH1D(EventCategory+"_"+LHCStatus+"_FirstPulseHeight_"+ChStr[c]+"_"+RangeCode,"Height of first pulse in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the first pulse duration per channel, with several different range choices
  if (RangeCode == "A") { nBins = 128; minX = 0.; maxX = 640.; }
  if (RangeCode == "B") { nBins = 64; minX = 0.; maxX = 160.; }
  if (RangeCode == "C") { nBins = 32; minX = 0.; maxX = 40.; }
  if (RangeCode == "D") { nBins = 32; minX = 0.; maxX = 20.; }
  TH1D *h_FirstDuration[32];
  for (int c=0; c<32; c++)
    h_FirstDuration[c] = new TH1D(EventCategory+"_"+LHCStatus+"_FirstPulseDuration_"+ChStr[c]+"_"+RangeCode,"Duration of first pulse in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the first pulse time per channel, with several different range choices
  if (RangeCode == "A") { nBins = 512; minX = 0.; maxX = 1024*0.625; }
  if (RangeCode == "B") { nBins = 50; minX = 0.; maxX = 200.; }
  if (RangeCode == "C") { nBins = 50; minX = 150.; maxX = 300.; }
  if (RangeCode == "D") { nBins = 50; minX = 300.; maxX = 450.; }
  TH1D *h_FirstTime[32];
  for (int c=0; c<32; c++)
    h_FirstTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_FirstPulseTime_"+ChStr[c]+"_"+RangeCode,"Time of first pulse in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the time difference wrt the closest LHC clock pulse
  if (RangeCode == "A") { nBins = 50; minX = -25.; maxX = 25.; }
  if (RangeCode == "B") { nBins = 30; minX = -15.; maxX = 15.; }
  if (RangeCode == "C") { nBins = 15; minX = -15.; maxX = 15.; }
  if (RangeCode == "D") { nBins = 10; minX = -15.; maxX = 15.; }
  TH1D *h_LHCTime[32];
  TH1D *h_LHCTimeCalib[32];
  for (int c=0; c<32; c++) {
    h_LHCTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_LHCTime_"+ChStr[c]+"_"+RangeCode,"Time wrt LHC clock of largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_LHCTimeCalib[c] = new TH1D(EventCategory+"_"+LHCStatus+"_LHCTimeCalib_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt LHC clock of largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
  }

  // Plots of the time difference between each channel's hit and a hit in specific other channels.
  // I compare to each of the slabs (18,20,28,21) and to the top right corner bar in each layer (1,7,3)
  // Compare only the largest pulse in each channel.
  if (RangeCode == "A") { nBins = 200; minX = -200.; maxX = 200.; }
  if (RangeCode == "B") { nBins = 50; minX = -200.; maxX = 200.; }
  if (RangeCode == "C") { nBins = 50; minX = -50.; maxX = 50.; }
  if (RangeCode == "D") { nBins = 30; minX = -30.; maxX = 30.; }
  TH1D *h_TimeWRTCH18[32];
  TH1D *h_TimeCalibWRTCH18[32];
  TH1D *h_TimeWRTCH20[32];
  TH1D *h_TimeCalibWRTCH20[32];
  TH1D *h_TimeWRTCH28[32];
  TH1D *h_TimeCalibWRTCH28[32];
  TH1D *h_TimeWRTCH21[32];
  TH1D *h_TimeCalibWRTCH21[32];
  TH1D *h_TimeWRTCH7[32];
  TH1D *h_TimeCalibWRTCH7[32];
  TH1D *h_TimeWRTCH0[32];
  TH1D *h_TimeCalibWRTCH0[32];
  TH1D *h_TimeWRTCH1[32];
  TH1D *h_TimeCalibWRTCH1[32];
  TH1D *h_TimeWRTCH2[32];
  TH1D *h_TimeCalibWRTCH2[32];
  TH1D *h_TimeWRTCH3[32];
  TH1D *h_TimeCalibWRTCH3[32];
  TH1D *h_TimeWRTCH4[32];
  TH1D *h_TimeCalibWRTCH4[32];
  for (int c=0; c<32; c++) {
    h_TimeWRTCH18[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH18_"+ChStr[c]+"_"+RangeCode,"Time wrt CH18 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH18[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH18_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH18 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH20[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH20_"+ChStr[c]+"_"+RangeCode,"Time wrt CH20 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH20[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH20_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH20 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH28[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH28_"+ChStr[c]+"_"+RangeCode,"Time wrt CH28 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH28[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH28_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH28 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH21[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH21_"+ChStr[c]+"_"+RangeCode,"Time wrt CH21 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH21[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH21_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH21 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH7[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH7_"+ChStr[c]+"_"+RangeCode,"Time wrt CH7 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH7[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH7_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH7 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH0[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH0_"+ChStr[c]+"_"+RangeCode,"Time wrt CH0 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH0[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH0_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH0 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH1[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH1_"+ChStr[c]+"_"+RangeCode,"Time wrt CH1 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH1[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH1_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH1 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH2[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH2_"+ChStr[c]+"_"+RangeCode,"Time wrt CH2 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH2[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH2_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH2 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH3[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH3_"+ChStr[c]+"_"+RangeCode,"Time wrt CH3 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH3[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH3_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH3 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeWRTCH4[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeWRTCH4_"+ChStr[c]+"_"+RangeCode,"Time wrt CH4 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
    h_TimeCalibWRTCH4[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TimeCalibWRTCH4_"+ChStr[c]+"_"+RangeCode,"Calibrated time wrt CH4 for largest pulse in Channel "+ChStr[c],nBins,minX,maxX);
  }

  // // Plots of the first pulse delta time (wrt earliest hit) per channel, with several different range choices
  // if (RangeCode == "A") { nBins = 480; minX = 0.; maxX = 300.; }
  // if (RangeCode == "B") { nBins = 240; minX = 0.; maxX = 150.; }
  // if (RangeCode == "C") { nBins = 80; minX = 0.; maxX = 50.; }
  // if (RangeCode == "D") { nBins = 40; minX = 0.; maxX = 25.; }
  // TH1D *h_FirstDeltaTime[32];
  // for (int c=0; c<32; c++)
  //   h_FirstDeltaTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_FirstPulseDeltaTime_"+ChStr[c]+"_"+RangeCode,"Time after trigger of first pulse in Channel "+ChStr[c],nBins,minX,maxX);
  
  // // Plots of the time wrt 300 ns of the closest hit. This tests if there is a hit near the trigger time, as should be required.  The histogram gets filled for the channel that contains that nearest hit
  // if (RangeCode == "A") { nBins = 40; minX = -200.; maxX = 200.; }
  // if (RangeCode == "B") { nBins = 40; minX = -50.; maxX = 50.; }
  // if (RangeCode == "C") { nBins = 32; minX = -20.; maxX = 20.; }
  // if (RangeCode == "D") { nBins = 16; minX = -20.; maxX = 20.; }
  // TH1D *h_Nearest1DeltaTime[32];
  // for (int c=0; c<32; c++)
  //   h_Nearest1DeltaTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_Nearest1DeltaTime_"+ChStr[c]+"_"+RangeCode,"Time wrt 300ns of nearest pulse when it is in Channel "+ChStr[c],nBins,minX,maxX);
  
  // // Plots of the time wrt 300 ns of the second closest hit *before* 300ns. This tests if there is a second hit within the trigger window, as should be required.  The histogram gets filled for the channel that contains that 2nd nearest hit
  // if (RangeCode == "A") { nBins = 40; minX = -200.; maxX = 0.; }
  // if (RangeCode == "B") { nBins = 40; minX = -50.; maxX = 0.; }
  // if (RangeCode == "C") { nBins = 32; minX = -20.; maxX = 0.; }
  // if (RangeCode == "D") { nBins = 16; minX = -20.; maxX = 0.; }
  // TH1D *h_Nearest2DeltaTime[32];
  // for (int c=0; c<32; c++)
  //  h_Nearest2DeltaTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_Nearest2DeltaTime_"+ChStr[c]+"_"+RangeCode,"Time wrt 300ns of 2nd nearest pulse when it is in Channel "+ChStr[c],nBins,minX,maxX);

  // // Plots of the time difference between the nearest and second nearest hit to 300ns. This tests if there is a second hit within the trigger window, as should be required.  The histogram gets filled for the channel that contains that 2nd nearest hit
  // if (RangeCode == "A") { nBins = 40; minX = -200.; maxX = 50.; }
  // if (RangeCode == "B") { nBins = 40; minX = -50.; maxX = 50.; }
  // if (RangeCode == "C") { nBins = 32; minX = -20.; maxX = 10.; }
  // if (RangeCode == "D") { nBins = 16; minX = -20.; maxX = 10.; }
  // TH1D *h_Nearest12DeltaTime[32];
  // for (int c=0; c<32; c++)
  //  h_Nearest12DeltaTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_Nearest12DeltaTime_"+ChStr[c]+"_"+RangeCode,"Time wrt trigger hit of 2nd nearest pulse when it is in Channel "+ChStr[c],nBins,minX,maxX);
  
  // Plots of the total pulse area (1/1000) per channel, with several different range choices
  if (RangeCode == "A") { nBins = 100; minX = 0.; maxX = 300.; }
  if (RangeCode == "B") { nBins = 100; minX = 0.; maxX = 100.; }
  if (RangeCode == "C") { nBins = 100; minX = 0.; maxX = 25.; }
  if (RangeCode == "D") { nBins = 100; minX = 0.; maxX = 5.; }
  TH1D *h_TotalPulseArea[32];
  for (int c=0; c<32; c++) 
   h_TotalPulseArea[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TotalPulseArea_"+ChStr[c]+"_"+RangeCode,"Total area of pulses in Channel "+ChStr[c],nBins,minX,maxX);

  // // Plots of the time difference between a pulse in one bar and a pulse in its upper or lower neighbor (if there is one).  For the sidecars, this is never true.
  // if (RangeCode == "A") { nBins = 128; minX = -40.; maxX = 40.; }
  // if (RangeCode == "B") { nBins = 64; minX = -20.; maxX = 20.; }
  // if (RangeCode == "C") { nBins = 32; minX = -20.; maxX = 20.; }
  // if (RangeCode == "D") { nBins = 32; minX = -10.; maxX = 10.; }
  // TH1D *h_VertNeighborDeltaTime[32];
  // for (int c=0; c<32; c++) 
  //  h_VertNeighborDeltaTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_VertNeighborDeltaTime_"+ChStr[c]+"_"+RangeCode,"Time difference to vertical neighbor hit in Channel "+ChStr[c],nBins,minX,maxX);

  // // Plots of the time difference between a pulse in one bar and a pulse in its angled upper or lower neighbor (if there is one).  For the sidecars, only fill in their channel otherwise fill in both.
  // if (RangeCode == "A") { nBins = 128; minX = -40.; maxX = 40.; }
  // if (RangeCode == "B") { nBins = 64; minX = -20.; maxX = 20.; }
  // if (RangeCode == "C") { nBins = 64; minX = -20.; maxX = 20.; }
  // if (RangeCode == "D") { nBins = 32; minX = -10.; maxX = 10.; }
  // TH1D *h_AngleNeighborDeltaTime[32];
  // for (int c=0; c<32; c++) 
  //   h_AngleNeighborDeltaTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_AngleNeighborDeltaTime_"+ChStr[c]+"_"+RangeCode,"Time difference to angled neighbor hit in Channel "+ChStr[c],nBins,minX,maxX);

  // // Plots of the time difference between a pulse in one bar and a pulse in its longitudinal neighbors (between layers).
  // if (RangeCode == "A") { nBins = 128; minX = -40.; maxX = 40.; }
  // if (RangeCode == "B") { nBins = 64; minX = -20.; maxX = 20.; }
  // if (RangeCode == "C") { nBins = 64; minX = -20.; maxX = 20.; }
  // if (RangeCode == "D") { nBins = 32; minX = -10.; maxX = 10.; }
  // TH1D *h_LayerNeighborDeltaTime[32];
  // for (int c=0; c<32; c++) 
  //   h_LayerNeighborDeltaTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_LayerNeighborDeltaTime_"+ChStr[c]+"_"+RangeCode,"Time difference to layer neighbor hit in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the number of pulses per channel, with several different range choices
  if (RangeCode == "A") { nBins = 20; minX = 0.; maxX = 20.; }
  if (RangeCode == "B") { nBins = 10; minX = 0.; maxX = 10.; }
  if (RangeCode == "C") { nBins = 5; minX = 0.; maxX = 5.; }
  if (RangeCode == "D") { nBins = 3; minX = 0.; maxX = 3.; }
  TH1D *h_NPulses[32];
  for (int c=0; c<32; c++)
   h_NPulses[c] = new TH1D(EventCategory+"_"+LHCStatus+"_NPulses_"+ChStr[c]+"_"+RangeCode,"Number of pulses in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of paired hits in other channels when this channel is hit
  if (RangeCode == "A") { nBins = 34; minX = -0.5; maxX = 33.5; }
  if (RangeCode == "B") { nBins = 18; minX = -0.5; maxX = 17.5; }
  if (RangeCode == "C") { nBins = 18; minX = 15.5; maxX = 33.5; }
  if (RangeCode == "D") { nBins = 36; minX = -0.5; maxX = 35.5; }
  TH1D *h_NPairedHits[32];
  for (int c=0; c<32; c++)
   h_NPairedHits[c] = new TH1D(EventCategory+"_"+LHCStatus+"_NPairedHits_"+ChStr[c]+"_"+RangeCode,"Number of coincident hits in other channels for Channel "+ChStr[c],nBins,minX,maxX);

  ///////////// Plots that are not channel dependent 

  // Plots of hit triplets when no other channel is hit
  // Bins are as follows:
  //  1 = CH0-6-2
  //  2 = CH1-7-3
  //  3 = CH24-16-22
  //  4 = CH25-17-23
  //  5 = CH8-12-4
  //  6 = CH9-13-5
  //  7 = CH10-30-14
  //  8 = CH29-19-26
  //  9 = CH27-11-31
  // 10 = CH18-20-28
  // 11 = CH18-21-28
  if (RangeCode == "A") { nBins = 15; minX = -0.5; maxX = 14.5; }
  if (RangeCode == "B") { nBins = 30; minX = -0.5; maxX = 14.5; }
  if (RangeCode == "C") { nBins = 15; minX = -0.5; maxX = 14.5; }
  if (RangeCode == "D") { nBins = 15; minX = -0.5; maxX = 14.5; }
  TH1D *h_NTripletHits;
  TH1D *h_NTripletHits3;
  TH1D *h_NTripletHits4;
  h_NTripletHits = new TH1D(EventCategory+"_"+LHCStatus+"_NTripletHits_"+RangeCode,"Number of triplet hits ",nBins,minX,maxX);
  h_NTripletHits3 = new TH1D(EventCategory+"_"+LHCStatus+"_NTripletHits3_"+RangeCode,"Number of triplet hits in 3 hit events ",nBins,minX,maxX);
  h_NTripletHits4 = new TH1D(EventCategory+"_"+LHCStatus+"_NTripletHits4_"+RangeCode,"Number of triplet hits in 4 hit events",nBins,minX,maxX);

  // Plots of number of hit bars
  if (RangeCode == "A") { nBins = 34; minX = -0.5; maxX = 33.5; }
  if (RangeCode == "B") { nBins = 18; minX = -0.5; maxX = 17.5; }
  if (RangeCode == "C") { nBins = 18; minX = 15.5; maxX = 33.5; }
  if (RangeCode == "D") { nBins = 36; minX = -0.5; maxX = 35.5; }
  TH1D *h_NHitBars;
  h_NHitBars = new TH1D(EventCategory+"_"+LHCStatus+"_NHitBars_"+RangeCode,"Number of hit bars",nBins,minX,maxX);

  // Plots of number of hit layers (bars only, not slabs or panels)
  if (RangeCode == "A") { nBins = 4; minX = -0.5; maxX = 3.5; }
  if (RangeCode == "B") { nBins = 4; minX = -0.5; maxX = 3.5; }
  if (RangeCode == "C") { nBins = 4; minX = -0.5; maxX = 3.5; }
  if (RangeCode == "D") { nBins = 4; minX = -0.5; maxX = 3.5; }
  TH1D *h_NHitLayers;
  h_NHitLayers = new TH1D(EventCategory+"_"+LHCStatus+"_NHitLayers_"+RangeCode,"Number of layers with hit bars",nBins,minX,maxX);

  // Plots of the total pulse area (1/1000) across the whole detector
  if (RangeCode == "A") { nBins = 100; minX = 0.; maxX = 2000.; }
  if (RangeCode == "B") { nBins = 100; minX = 0.; maxX = 500.; }
  if (RangeCode == "C") { nBins = 100; minX = 0.; maxX = 200.; }
  if (RangeCode == "D") { nBins = 100; minX = 0.; maxX = 50.; }
  TH1D *h_TotalGlobalPulseArea;
  h_TotalGlobalPulseArea = new TH1D(EventCategory+"_"+LHCStatus+"_TotalGlobalPulseArea_"+RangeCode,"Total area of all pulses across detector",nBins,minX,maxX);

  if (RangeCode == "A") { nBins = nentries/10; minX = 0.; maxX = 1.*(nentries); }
  if (RangeCode == "B") { nBins = nentries/50; minX = 0.; maxX = 1.*(nentries); }
  if (RangeCode == "C") { nBins = nentries/100; minX = 0.; maxX = 1.*(nentries); }
  if (RangeCode == "D") { nBins = nentries/1000; minX = 0.; maxX = 1.*(nentries); }
  TH1D *h_EventNumberUnsynched;
  h_EventNumberUnsynched = new TH1D(EventCategory+"_"+LHCStatus+"_EventNumberUnsynched_"+RangeCode,"Unsynched event numbers",nBins,minX,maxX);

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
        h_EventNumberUnsynched->Fill(float(ia),1.);
      }
      if (prevTime>0 && ((filenum==1 && event>3) || (filenum>1 && event>1))) {
        SumTime += event_time_fromTDC-prevTime;
        if (synched) SumSynchTime += event_time_fromTDC-prevTime;
      }
      prevTime = event_time_fromTDC;
     if (event_time_fromTDC < minTime) { // New first event
      minTime = event_time_fromTDC;
      firstEvt = ia;
      if (RunNum>24) FirstTime = *event_t_string;
     } // New first event
     if (event_time_fromTDC > maxTime) { // New last event
      maxTime = event_time_fromTDC;
      lastEvt = ia;
      if (RunNum>24) LastTime = *event_t_string;
     } // New last event
   } // Ignore the first events, which may be buggy
  } // Loop over all events
  double minHour = minTime/3600; double maxHour = maxTime/3600;
  double runTime = SumTime;

  // if (EventCategory=="All" && RangeCode=="A") {
  //  prevTime = minTime;
  //  for(int ia = 0; ia<nentries; ia++){ // Loop over all events
  //   chain->GetEntry(ia);
  //   cout << "TIMEANT: "<<filenum<<" "<<event<<" "<<10000*(filenum-1)+event<<" "<<Long64_t(event_time_fromTDC)<<" "<<event_time_fromTDC-minTime<<" "<<1000*(event_time_fromTDC-prevTime)<<"\n";
  //   prevTime = event_time_fromTDC;
  //  } // Loop over all events
  // }

  // Determine the LHC fills present during this run. 
  // For now, only include fully contained fills...
  float totalLumi = 0.;
  cout << "Epoch time range = "<<Long64_t(minTime)<<" - "<<Long64_t(maxTime)<<"\n";
  for (unsigned int iFill=0; iFill<fillNum.size(); iFill++) {
    if (fillStartTime[iFill]>minTime && fillStartTime[iFill]<maxTime && fillEndTime[iFill]<maxTime) {
      // Make sure that we have not already used this fill.
      bool newFill = true;
      for (unsigned int iFill0 = 0; iFill0<iFill; iFill0++)
      	if (fillNum[iFill0]==fillNum[iFill])
      	  newFill = false;
      if (newFill) {
       cout << "Including lumi from fill "<<fillNum[iFill]<<" "<<fillStartTime[iFill]<<" "<<fillEndTime[iFill]<<" "<<fillLumi[iFill]<<"\n";
       totalLumi += fillLumi[iFill];
      }
    }  vector<int> fillNum;
  }  
  // Print a summary of the run time and lumi
  cout << "Summary: Run "<<RunNum<<" "<<nentries<<" events over "<<setprecision(3)<<(maxTime-minTime)/3600.<<" hours between "<<FirstTime<<" and "<<LastTime<<"; Trigger rate = "<<nentries/(maxTime-minTime)<<" Hz. Livetime="<<int(runTime)<<" sec = "<<setprecision(4)<<runTime/3600.<<" hr. Lumi="<<totalLumi<<"/pb. Scale to /hour = "<<setprecision(6)<<3600./(runTime)<<"\n";
  cout << "Synch: &epsilon;<sub>sync</sub> = " <<setprecision(3)<< (100.*nSynched)/(nSynched+nUnSynched) << " &plusmn; " << 100.*sqrt(nSynched)/(nSynched+nUnSynched) << " %";
  cout << "  T<sub>sync</sub> = "<<setprecision(3)<<SumSynchTime/3600.<<" hr\n";
  if (RangeCode == "A") { nBins = 1000; minX = minHour - 0.5; maxX = maxHour + 2.; }
  if (RangeCode == "B") { nBins = 200; minX = minHour - 0.5; maxX = maxHour + 2.; }
  if (RangeCode == "C") { nBins = 60; minX = minHour - 0.5; maxX = maxHour + 2.; }
  if (RangeCode == "D") { nBins = 10; minX = minHour - 0.5; maxX = maxHour + 2.; }
  TH1D *h_EventTime[32];
  for (int c=0; c<32; c++)
   h_EventTime[c] = new TH1D(EventCategory+"_"+LHCStatus+"_EventTime_"+ChStr[c]+"_"+RangeCode,"Time of events with a hit in Channel "+ChStr[c],nBins,minX,maxX);

  if (RangeCode == "A") { nBins = 500; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "B") { nBins = 200; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "C") { nBins = 60; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "D") { nBins = 10; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  TH1D *h_EventTimeSeconds[32];
  for (int c=0; c<32; c++)
   h_EventTimeSeconds[c] = new TH1D(EventCategory+"_"+LHCStatus+"_EventTimeSeconds_"+ChStr[c]+"_"+RangeCode,"Seconds into run for events with a hit in Channel"+ChStr[c],nBins,minX,maxX);

  if (RangeCode == "A") { nBins = 500; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "B") { nBins = 200; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "C") { nBins = 60; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "D") { nBins = 10; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  TH1D *h_EventTimeSynchedSeconds[32];
  for (int c=0; c<32; c++)
   h_EventTimeSynchedSeconds[c] = new TH1D(EventCategory+"_"+LHCStatus+"_EventTimeSynchedSeconds_"+ChStr[c]+"_"+RangeCode,"Seconds into run for synched events with a hit in Channel"+ChStr[c],nBins,minX,maxX);

  if (RangeCode == "A") { nBins = 500; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "B") { nBins = 200; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "C") { nBins = 60; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  if (RangeCode == "D") { nBins = 10; minX = 0.; maxX = 1.2*(maxTime-minTime); }
  TH1D *h_EventTimeUnsynchedSeconds[32];
  for (int c=0; c<32; c++)
   h_EventTimeUnsynchedSeconds[c] = new TH1D(EventCategory+"_"+LHCStatus+"_EventTimeUnsynchedSeconds_"+ChStr[c]+"_"+RangeCode,"Seconds into run for unsynched events with a hit in Channel"+ChStr[c],nBins,minX,maxX);

  if (RangeCode == "A") { nBins = maxFile+5; minX = - 0.5; maxX = maxFile + 4.5; }
  if (RangeCode == "B") { nBins = maxFile+5; minX = - 0.5; maxX = maxFile + 4.5; }
  if (RangeCode == "C") { nBins = maxFile+5; minX = - 0.5; maxX = maxFile + 4.5; }
  if (RangeCode == "D") { nBins = maxFile+5; minX = - 0.5; maxX = maxFile + 4.5; }
  TH1D *h_FileNum[32];
  for (int c=0; c<32; c++)
   h_FileNum[c] = new TH1D(EventCategory+"_"+LHCStatus+"_FileNum_"+ChStr[c]+"_"+RangeCode,"Filenumber of events with a hit in Channel "+ChStr[c],nBins,minX,maxX);

  // Histograms of the sideband mean
  if (RangeCode == "A") { nBins = 1000; minX = -20.; maxX = 20.; }
  if (RangeCode == "B") { nBins = 200; minX = -10.; maxX = 10.; }
  if (RangeCode == "C") { nBins = 100; minX = -5.; maxX = 5.; }
  if (RangeCode == "D") { nBins = 80; minX = -2.; maxX = 2.; }
  TH1D *h_SidebandMean[32];
  for (int c=0; c<32; c++)
   h_SidebandMean[c] = new TH1D(EventCategory+"_"+LHCStatus+"_SidebandMean_"+ChStr[c]+"_"+RangeCode,"Sideband mean in Channel "+ChStr[c],nBins,minX,maxX);

  // Histograms of the sideband RMS
  if (RangeCode == "A") { nBins = 1000; minX = 0.; maxX = 20.; }
  if (RangeCode == "B") { nBins = 200; minX = 0.; maxX = 10.; }
  if (RangeCode == "C") { nBins = 100; minX = 0.; maxX = 5.; }
  if (RangeCode == "D") { nBins = 80; minX = 0.; maxX = 2.; }
  TH1D *h_SidebandRMS[32];
  for (int c=0; c<32; c++)
   h_SidebandRMS[c] = new TH1D(EventCategory+"_"+LHCStatus+"_SidebandRMS_"+ChStr[c]+"_"+RangeCode,"Sideband RMS in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the max sample height for tagged channels, with several different range choices
  if (RangeCode == "A") { nBins = 70; minX = 0.; maxX = 1400.; }
  if (RangeCode == "B") { nBins = 50; minX = 0.; maxX = 500.; }
  if (RangeCode == "C") { nBins = 25; minX = 0.; maxX = 100.; }
  if (RangeCode == "D") { nBins = 20; minX = 0.; maxX = 40.; }
  TH1D *h_TaggedMax[32];
  for (int c=0; c<32; c++)
   h_TaggedMax[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TaggedMax_"+ChStr[c]+"_"+RangeCode,"Max sample in tagged Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the tagged pulse area (1/1000) per channel, with several different range choices
  if (RangeCode == "A") { nBins = 100; minX = 0.; maxX = 200.; }
  if (RangeCode == "B") { nBins = 100; minX = 0.; maxX = 100.; }
  if (RangeCode == "C") { nBins = 50; minX = 0.; maxX = 25.; }
  if (RangeCode == "D") { nBins = 50; minX = 0.; maxX = 5.; }
  TH1D *h_TaggedPulseArea[32];
  for (int c=0; c<32; c++)
    h_TaggedPulseArea[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TaggedPulseArea_"+ChStr[c]+"_"+RangeCode,"Area of neighbor tagged pulse in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the tagged total pulse area (1/1000) per channel, with several different range choices
  if (RangeCode == "A") { nBins = 100; minX = 0.; maxX = 200.; }
  if (RangeCode == "B") { nBins = 100; minX = 0.; maxX = 100.; }
  if (RangeCode == "C") { nBins = 100; minX = 0.; maxX = 25.; }
  if (RangeCode == "D") { nBins = 100; minX = 0.; maxX = 5.; }
  TH1D *h_TaggedTotalPulseArea[32];
  for (int c=0; c<32; c++)
    h_TaggedTotalPulseArea[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TaggedTotalPulseArea_"+ChStr[c]+"_"+RangeCode,"Total area of neighbor tagged pulses in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the tagged pulse height per channel, with several different range choices
  if (RangeCode == "A") { nBins = 130; minX = 0.; maxX = 1300.; }
  if (RangeCode == "B") { nBins = 100; minX = 0.; maxX = 500.; }
  if (RangeCode == "C") { nBins = 50; minX = 0.; maxX = 150.; }
  if (RangeCode == "D") { nBins = 40; minX = 0.; maxX = 40.; }
  TH1D *h_TaggedPulseHeight[32];
  for (int c=0; c<32; c++)
    h_TaggedPulseHeight[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TaggedPulseHeight_"+ChStr[c]+"_"+RangeCode,"Height of neighbor tagged pulse in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the tagged pulse duration per channel, with several different range choices
  if (RangeCode == "A") { nBins = 128; minX = 0.; maxX = 640.; }
  if (RangeCode == "B") { nBins = 64; minX = 0.; maxX = 160.; }
  if (RangeCode == "C") { nBins = 32; minX = 0.; maxX = 40.; }
  if (RangeCode == "D") { nBins = 32; minX = 0.; maxX = 20.; }
  TH1D *h_TaggedDuration[32];
  for (int c=0; c<32; c++)
    h_TaggedDuration[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TaggedPulseDuration_"+ChStr[c]+"_"+RangeCode,"Duration of neighbor tagged pulse in Channel "+ChStr[c],nBins,minX,maxX);

  // Plots of the pulse height in each channel that is part of a three-fold coincidence with its inter-layer partners
  if (RangeCode == "A") { nBins = 140; minX = 0.; maxX = 1400.; }
  if (RangeCode == "B") { nBins = 50; minX = 0.; maxX = 500.; }
  if (RangeCode == "C") { nBins = 50; minX = 0.; maxX = 100.; }
  if (RangeCode == "D") { nBins = 40; minX = 0.; maxX = 40.; }
  TString partnerThresholdStr = " 10mVPartnerThreshold";
  if (EventCategory == "SPE") partnerThresholdStr = " 20mVPartnerThreshold";
  if (EventCategory == "Small" || EventCategory == "SmallTag") partnerThresholdStr = " 30mVPartnerThreshold";
  if (EventCategory == "Cosmic" || EventCategory == "CosmicTag") partnerThresholdStr = " 40mVPartnerThreshold";
  if (EventCategory == "Thru" || EventCategory == "ThruTag") partnerThresholdStr = " 50mVPartnerThreshold";
  TH1D *h_TripleLayerMax[32];
  for (int c=0; c<32; c++)
    h_TripleLayerMax[c] = new TH1D(EventCategory+"_"+LHCStatus+"_TripleLayerMax_"+ChStr[c]+"_"+RangeCode,"Max sample in TripleLayer Channel "+ChStr[c]+partnerThresholdStr,nBins,minX,maxX);

  // Plots of the pulse height in each channel that is part of a two-fold coincidence with its inter-layer partners
  if (RangeCode == "A") { nBins = 140; minX = 0.; maxX = 1400.; }
  if (RangeCode == "B") { nBins = 50; minX = 0.; maxX = 500.; }
  if (RangeCode == "C") { nBins = 50; minX = 0.; maxX = 100.; }
  if (RangeCode == "D") { nBins = 40; minX = 0.; maxX = 40.; }
  TH1D *h_DoubleLayerMax[32];
  for (int c=0; c<32; c++)
    h_DoubleLayerMax[c] = new TH1D(EventCategory+"_"+LHCStatus+"_DoubleLayerMax_"+ChStr[c]+"_"+RangeCode,"Max sample in DoubleLayer Channel "+ChStr[c]+partnerThresholdStr,nBins,minX,maxX);

//  TH1D *h_OccupancyMap_0 = new TH1D(EventCategory+"_"+LHCStatus+"_OccupancyMap_"+RangeCode,"Occupancy Map",16,-0.5,15.5);
  // Save information about the trigger rate for different threshold options.
  // The goal is to produce an ascii ntuple with the rate (events/hour) with different choices
  // of trigger thresholds for each channel. I measure this brute force by looping over all 16 threshold options
  // and count the number of events that would pass either a single channel or double coincidence trigger
  // using these thresholds.  There is no timing information used, just max>threshold.
  int nSingleChannel[32][600]; // 32 channels and thresholds from 0 to 600 mV
  int nDoubleCoincidence[32][600]; // 32 channels and thresholds from 0 to 600 mV
  // Initialize the trigger emulation counters
  for (int c=0; c<32; c++)
    for (int iThr=0; iThr<600; iThr++) {
      nSingleChannel[c][iThr] = 0;
      nDoubleCoincidence[c][iThr] = 0;
    }

  prevTime = 0.; // Initialize the previout event time; used to determine time between events
  long timeBetweenEvents = 0; // microseconds
  long prev_event_trigger_time_tag_b0 = 0;
  long prev_event_trigger_time_tag_b1 = 0;
  long prev_gTDC_b0 = 0;
  long prev_gTDC_b1 = 0;

  //Main Event loop
  for(int ia = 0; ia<nentries; ia++){
    if (debug>100 || (ia%100000 == 0)) {cout << "Entry #"<<ia<<"\n"; cout.flush();}

    timeBetweenEvents = long(1000000*(event_time_fromTDC-prevTime));
    prevTime = event_time_fromTDC;

    chain->GetEntry(ia);

    bool GoodEvent = true;

    synched = false;
    if (event_trigger_time_tag_b1==event_trigger_time_tag_b0)
      synched = true;

    if (run==24 && event<20) GoodEvent = false;
    if (event==0 && filenum==1) GoodEvent = false; // The first event is sometimes buggy

    int evtCategoryCode = 0;
    if (EventCategory == "All") evtCategoryCode = 0;
    if (EventCategory == "SPE") evtCategoryCode = 1;
    if ((EventCategory == "Small") || (EventCategory == "SmallTag")) evtCategoryCode = 2;
    if ((EventCategory == "Cosmic") || (EventCategory == "CosmicTag")) evtCategoryCode = 3;
    if ((EventCategory == "Thru") || (EventCategory == "ThruTag")) evtCategoryCode = 4;

    // Calculate total area of pulses and save a more convenient maxSample per channel
    float TotalArea[32]; // in nVs
    float TotalGlobalArea = 0.;
    float maxSample[32]; // in mV
    for (int c=0; c<32; c++)
      TotalArea[c] = 0.; // Sum of area of pulses in the channel
    for (int i=0; i<int(chan->size()); i++) {
      TotalArea[int(chan->at(i))] += area->at(i)/1000.;
      TotalGlobalArea += area->at(i)/1000.;
    }
    for (int c=0; c<32; c++)
     maxSample[c] = max_sample->at(c);

    // Update the trigger emulation counters, but only for the first category of analysis
//    if (EventCategory == "All" && LHCStatus == "All" && RangeCode == "A") { // Update trigger emulation counters?
//      // Loop through every trigger threshold option for every channel
//      int TrigForChan[32];
//      for (TrigForChan[0]=0; TrigForChan[0]<600; TrigForChan[0]++)
//      for (TrigForChan[1]=0; TrigForChan[1]<600; TrigForChan[1]++)
//      for (TrigForChan[2]=0; TrigForChan[2]<600; TrigForChan[2]++)
//      for (TrigForChan[3]=0; TrigForChan[3]<600; TrigForChan[3]++)
//      for (TrigForChan[4]=0; TrigForChan[4]<600; TrigForChan[4]++)
//      for (TrigForChan[5]=0; TrigForChan[5]<600; TrigForChan[5]++)
//      for (TrigForChan[6]=0; TrigForChan[6]<600; TrigForChan[6]++)
//      for (TrigForChan[7]=0; TrigForChan[7]<600; TrigForChan[7]++)
//      for (TrigForChan[8]=0; TrigForChan[8]<600; TrigForChan[8]++)
//      for (TrigForChan[9]=0; TrigForChan[9]<600; TrigForChan[9]++)
//      for (TrigForChan[10]=0; TrigForChan[10]<600; TrigForChan[10]++)
//      for (TrigForChan[11]=0; TrigForChan[11]<600; TrigForChan[11]++)
//      for (TrigForChan[12]=0; TrigForChan[12]<600; TrigForChan[12]++)
//      for (TrigForChan[13]=0; TrigForChan[13]<600; TrigForChan[13]++)
//      for (TrigForChan[14]=0; TrigForChan[14]<600; TrigForChan[14]++)
//      for (TrigForChan[15]=0; TrigForChan[15]<600; TrigForChan[15]++) {
//        // Determine the number of groups above threshold for this option.
//        int groupTrig[8];
//        for (int c=0; c<8; c++) groupTrig[c] = false;
//        if (maxSample[0] > float(TrigForChan[0]) || 
//            maxSample[1] > float(TrigForChan[1]) ) groupTrig[0] = true; //0=group1
//        if (maxSample[2] > float(TrigForChan[2]) || 
//            maxSample[3] > float(TrigForChan[3]) ) groupTrig[1] = true; //1=group2
//        if (maxSample[4] > float(TrigForChan[4]) || 
//            maxSample[5] > float(TrigForChan[5]) ) groupTrig[2] = true; //2=group3
//        if (maxSample[6] > float(TrigForChan[6]) || 
//            maxSample[7] > float(TrigForChan[7]) ) groupTrig[3] = true; //3=group4
//        if (maxSample[8] > float(TrigForChan[8]) || 
//            maxSample[9] > float(TrigForChan[9]) ) groupTrig[4] = true; //4=group5
//        if (maxSample[10] > float(TrigForChan[10]) || 
//            maxSample[11] > float(TrigForChan[11]) ) groupTrig[5] = true; //5=group6
//        if (maxSample[12] > float(TrigForChan[12]) || 
//            maxSample[13] > float(TrigForChan[13]) ) groupTrig[6] = true; //6=group7
//        if (maxSample[14] > float(TrigForChan[14]) || 
//            maxSample[15] > float(TrigForChan[15]) ) groupTrig[7] = true; //7=group8
//          nSingleChannel[c][iThr] = 0;
//          nDoubleCoincidence[c][iThr] = 0;
//      }
//    } // Update trigger emulation counters?


    // Build a map of hits above thresholds appropriate for the passed event category, and record the pulse number; consider only the first pulse in each channel.
    bool HitChan[32]; // Simple flag showing existence of a hit
    int ChanHitIndex[32]; // Index to the pulse with the hit
    int nHits=0; // Number of channels with hits
    for (int c=0; c<32; c++) {
      ChanHitIndex[c] = -1; // No hit
      HitChan[c] = false;
    }
    for (int c=0; c<32; c++)
     for (int i=0; i<int(chan->size()); i++) // Pulse loop 
	  if ( (chan->at(i) == c) &&
	     //(ipulse->at(i) == 0) &&
	     (duration->at(i) > DThresholds[c][evtCategoryCode]) &&
	     // (height->at(i) > VThresholds[c][evtCategoryCode]) &&
	     // (area->at(i)/1000. > AThresholds[c][evtCategoryCode]) &&
	     (maxSample[c] > VThresholds[c][evtCategoryCode]) &&
	     (TotalArea[c] > AThresholds[c][evtCategoryCode])
	     ) { // Appropriate hit in this channel?
	        HitChan[c] = true;
	        ChanHitIndex[c] = i;
	       } // Appropriate hit in this channel?

    for (int c=0; c<32; c++)
      if (HitChan[c]) nHits++;

    // Plot of the number of hit bars
    h_NHitBars->Fill(1.*nHits,1.);

    // Plot of the global total pulse area (across all bars)
    h_TotalGlobalPulseArea->Fill(TotalGlobalArea,1.);
    
    // Check if this event passes the requirements for the selected event category
    if ((evtCategoryCode == 0) && (nHits<1)) GoodEvent = false; // All category
    if ((evtCategoryCode == 1) && (nHits>0)) GoodEvent = false; // SPE category; No hits above the threshold
    if ((EventCategory == "Small") && (nHits<2)) GoodEvent = false; // Small category
    if ((EventCategory == "SmallTag") && (nHits<1)) GoodEvent = false; // SmallTag category
    if ((EventCategory == "Cosmic") && (nHits<2)) GoodEvent = false; // Cosmic category
    if ((EventCategory == "CosmicTag") && (nHits<2)) GoodEvent = false; // CosmicTag category
    if ((evtCategoryCode == 4) && !(HitChan[18] && HitChan[20] && HitChan[28] && HitChan[21])) GoodEvent = false; // Through going particle category

    // Veto events based on the LHCStatus string
    if (LHCStatus=="Beam" && !beam) GoodEvent = false;
    if (LHCStatus=="NoBeam" && beam) GoodEvent = false;

    // // Determine if the hit channel pattern is consistent with neighboring hits for cosmics.
    // bool cosmicneighbor = false;
    // // Hit patterns within layer0
    // if (HitChan[2] && HitChan[0]) cosmicneighbor = true;
    // if (HitChan[3] && HitChan[1]) cosmicneighbor = true;
    // if (HitChan[2] && HitChan[1]) cosmicneighbor = true;
    // if (HitChan[3] && HitChan[0]) cosmicneighbor = true;
    // if (HitChan[3] && HitChan[4]) cosmicneighbor = true;
    // if (HitChan[3] && HitChan[5]) cosmicneighbor = true;
    // if (HitChan[2] && HitChan[5]) cosmicneighbor = true;
    // // Hit patterns within layer1
    // if (HitChan[6] && HitChan[4]) cosmicneighbor = true;
    // if (HitChan[7] && HitChan[5]) cosmicneighbor = true;
    // if (HitChan[6] && HitChan[5]) cosmicneighbor = true;
    // if (HitChan[7] && HitChan[4]) cosmicneighbor = true;
    // if (HitChan[7] && HitChan[12]) cosmicneighbor = true;
    // // Hit patterns within layer2
    // if (HitChan[10] && HitChan[8]) cosmicneighbor = true;
    // if (HitChan[11] && HitChan[9]) cosmicneighbor = true;
    // if (HitChan[10] && HitChan[9]) cosmicneighbor = true;
    // if (HitChan[8] && HitChan[11]) cosmicneighbor = true;
    // if (HitChan[11] && HitChan[14]) cosmicneighbor = true;
    // // Hit patterns between layer0 and layer1
    // if (HitChan[2] && HitChan[4]) cosmicneighbor = true;
    // if (HitChan[0] && HitChan[6]) cosmicneighbor = true;
    // if (HitChan[3] && HitChan[5]) cosmicneighbor = true;
    // if (HitChan[1] && HitChan[7]) cosmicneighbor = true;
    // if (HitChan[15] && HitChan[12]) cosmicneighbor = true;
    // // Hit patterns between layer1 and layer2
    // if (HitChan[6] && HitChan[8]) cosmicneighbor = true;
    // if (HitChan[7] && HitChan[9]) cosmicneighbor = true;
    // if (HitChan[4] && HitChan[10]) cosmicneighbor = true;
    // if (HitChan[5] && HitChan[11]) cosmicneighbor = true;
    // if (HitChan[12] && HitChan[14]) cosmicneighbor = true;
    // if ((EventCategory == "CosmicTag") && (!cosmicneighbor)) GoodEvent = false;

    // Temporary code for veto panel stack test
    // Count the number of channels that have more than 7 mV in them. Mark the event as good for CosmicTag if there are at least 6 such channels.
    if ((EventCategory=="CosmicTag") && (RunNum>231 && RunNum<246)) {
     int nStackHits = 0;
     for (int c=0; c<32; c++)
       if (maxSample[c]>7.)
      	 nStackHits++;
     if (nStackHits>=6) {
       GoodEvent = true;
       //cout << "CosmicTag ia= "<<ia<<" event= "<<event<<" nStackHits= "<<nStackHits<<"\n";
     }
     else
       GoodEvent = false;
    }


    if (GoodEvent) { // Good event?

     // Plots of the raw MAX pulse height in each channel
     for (int c=0; c<32; c++)
      h_Max[c]->Fill(maxSample[c],1.);

     // Plots of the sideband mean
     for (unsigned int i = 0; i<sideband_mean->size(); i++)
	  if (ipulse->at(i) == 0)
	    h_SidebandMean[int(chan->at(i))]->Fill(sideband_mean->at(i),1.);

     // Plots of the sideband RMS
     for (unsigned int i = 0; i<sideband_mean->size(); i++)
	  if (ipulse->at(i) == 0 ) 
	   h_SidebandRMS[int(chan->at(i))]->Fill(sideband_RMS->at(i),1.);
	        
     // Plots of event time and filenum for channels with hits
     for (int c=0; c<32; c++)
      if (HitChan[c]) {
       h_FileNum[c]->Fill(1.*filenum,1.);
       h_EventTime[c]->Fill(event_time_fromTDC/3600.,1.);
	   h_EventTimeSeconds[c]->Fill(event_time_fromTDC-minTime,1.);
	   if (synched) 
	   	h_EventTimeSynchedSeconds[c]->Fill(event_time_fromTDC-minTime,1.);
       else
	   	h_EventTimeUnsynchedSeconds[c]->Fill(event_time_fromTDC-minTime,1.);
      }

  //    // Make plots of the time wrt 300 ns of the nearest and second nearest hit, where the second nearest hit is required to be *before* the closest hit.
  //    int i1 = -1; // Index to nearest hit
  //    int i2 = -1; // Index to 2nd nearest hit
  //    // Find closest
  //    for (int i=0; i<int(ptime->size()); i++) // Pulse loop 
	 //  if (i1<0 || fabs(ptime->at(i)-300.)<fabs(ptime->at(i1)-300.) ) // Nearer?
	 //   i1 = i;
  //    if (i1>=0) { // Found a nearest hit?
	 //  h_Nearest1DeltaTime[int(chan->at(i1))]->Fill(float(ptime->at(i1)-300.),1.);
	 //  // Find 2nd closest
	 //  for (int i=0; i<int(ptime->size()); i++) // Pulse loop 
	 //   if (ptime->at(i)<ptime->at(i1)) // Before the trigger hit?
	 //    if (i2<0 || ptime->at(i1)-ptime->at(i)<ptime->at(i1)-ptime->at(i2) ) // Nearer?
	 //     i2 = i;
	 //  if (i2>=0) { // Found a 2nd nearest hit?
	 //  // Fill histograms of the time of this 2nd nearest hit.
	 //  h_Nearest2DeltaTime[int(chan->at(i2))]->Fill(float(ptime->at(i2)-300.),1.);
	 //  // Fill histograms of the time of this 2nd nearest hit minus the nearest.
	 //  h_Nearest12DeltaTime[int(chan->at(i2))]->Fill(float(ptime->at(i2)-ptime->at(i1)),1.);
	 // } // Found a 2nd nearest hit?
  //   } // Found a nearest hit?
      
    // Figure out the event time as the time of the earliest pulse in any channel.
    float t0 = 9.E9;  
    for (int i=0; i<int(ptime->size()); i++) // Pulse loop 
	 if (ipulse->at(i) == 0 ) // First pulse in this channel
	  if (ptime->at(i)<t0) t0 = ptime->at(i);
      
    // Plots of the number of pulses, first pulse time and area in each channel
    for (int i=0; i<int(ptime->size()); i++) { // Pulse loop 
	   if (ipulse->at(i) == 0 ) { // First pulse in this channel
	    h_FirstTime[int(chan->at(i))]->Fill(float(ptime->at(i)),1.);
	    // h_FirstDeltaTime[int(chan->at(i))]->Fill(float(ptime->at(i))-t0,1.);
	    h_FirstPulseArea[int(chan->at(i))]->Fill(float(area->at(i))/1000.,1.);
	    h_FirstPulseHeight[int(chan->at(i))]->Fill(float(height->at(i)),1.);
	    h_FirstDuration[int(chan->at(i))]->Fill(float(duration->at(i)),1.);
	    h_NPulses[int(chan->at(i))]->Fill(float(npulses->at(i)),1.);
	   } // First pulse in this channel
    } // Pulse loop

    // Plot the pulse time relative to the closest LHC clock edge
    for (int c=0; c<32; c++)
      if (c!=15) { // Channel loop
        // First find the index of the largest pulse from this channel in the event
        int iLargest = -1;
        for (int i=0; i<int(ptime->size()); i++) 
          if (int(chan->at(i))==c) {
            if (iLargest<0)
              iLargest = i;
            else if (area->at(i)>area->at(iLargest))
              iLargest = i;
          }
        // Now find the nearest LHC clock edge (pulse)
        if (iLargest>=0 && area->at(iLargest)>AThresholds[c][evtCategoryCode] && 
          maxSample[c]>VThresholds[c][evtCategoryCode]) { // A pulse found?
         float dLHCtime = 999.;
         for (int iL=0; iL<int(ptime->size()); iL++) 
          if (int(chan->at(iL)) == 15) { // LHC pulse loop
           if (fabs(ptime->at(iLargest)-ptime->at(iL))<fabs(dLHCtime))
            dLHCtime = ptime->at(iLargest)-ptime->at(iL);
          } // LHC pulse loop
          h_LHCTime[c]->Fill(dLHCtime);
          // Repeat using calibrated time
          dLHCtime = 999.;
          for (int iL=0; iL<int(ptime->size()); iL++) 
           if (int(chan->at(iL)) == 15) { // LHC pulse loop
            if (fabs(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iL))<fabs(dLHCtime))
             dLHCtime = time_module_calibrated->at(iLargest)-time_module_calibrated->at(iL);
          } // LHC pulse loop
          h_LHCTimeCalib[c]->Fill(dLHCtime);
        } // A pulse found?
    } // Channel loop

    // Plot the pulse time relative to the largest pulse in various channels
    for (int iprobe=0; iprobe<10; iprobe++) { // Loop through 1,7,3,18,20,28,21,0,2,4
      int prbCH = 1;
      if (iprobe==0) prbCH=1;
      if (iprobe==1) prbCH=7;
      if (iprobe==2) prbCH=3;
      if (iprobe==3) prbCH=18;
      if (iprobe==4) prbCH=20;
      if (iprobe==5) prbCH=28;
      if (iprobe==6) prbCH=21;
      if (iprobe==7) prbCH=0;
      if (iprobe==8) prbCH=2;
      if (iprobe==9) prbCH=4;
      // First find the index to the largest pulse in the probe channel
      int iCH = -1;
      for (int i=0; i<int(ptime->size()); i++) 
        if (int(chan->at(i))==prbCH) {
          if (iCH<0 || area->at(i)>area->at(iCH))
            iCH = i;
      }
      // Now loop over all channels
      if (iCH>=0 && area->at(iCH)>AThresholds[prbCH][evtCategoryCode] && 
          maxSample[prbCH]>VThresholds[prbCH][evtCategoryCode]) { // Found a pulse in probe channel?
        for (int c=0; c<32; c++)
          if (c!=15 && c!=prbCH) { // Channel loop
            // First find the index to the largest pulse from this channel in the event
            int iLargest = -1;
            for (int i=0; i<int(ptime->size()); i++) 
              if (int(chan->at(i))==c) {
                if (iLargest<0 || area->at(i)>area->at(iLargest))
                  iLargest = i;
            }
            if (iLargest>=0 && area->at(iLargest)>AThresholds[c][evtCategoryCode] && maxSample[c]>VThresholds[c][evtCategoryCode]) {
              if (iprobe==0) {
               h_TimeWRTCH1[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH1[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==7) {
               h_TimeWRTCH0[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH0[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==8) {
               h_TimeWRTCH2[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH2[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==9) {
               h_TimeWRTCH4[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH4[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==1) {
               h_TimeWRTCH7[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH7[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==2) {
               h_TimeWRTCH3[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH3[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==3) {
               h_TimeWRTCH18[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH18[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==4) {
               h_TimeWRTCH20[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH20[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==5) {
               h_TimeWRTCH28[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH28[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
              if (iprobe==6) {
               h_TimeWRTCH21[c]->Fill(ptime->at(iLargest)-ptime->at(iCH));
               h_TimeCalibWRTCH21[c]->Fill(time_module_calibrated->at(iLargest)-time_module_calibrated->at(iCH));
              }
            }
        } // Channel loop
      } // Found a pulse in probe channel?
    } // Loop through 1,7,3,18,20,28,21,0,2,4

    // Plots of afterpulses.
    // Only consider clean after pulses with at least 10 ns between it and the end of the previous pulse. 
    float firstTime[32]; // Time of first pulse in channel
    for (int i=0; i<int(ptime->size()); i++)
      if (ipulse->at(i)==0)
        firstTime[int(chan->at(i))] = ptime->at(i);
    for (int i=0; i<int(ptime->size()); i++) { // Pulse loop 
	   if (ipulse->at(i) != 0 ) { // Not first pulse in this channel?
	    if (ptime->at(i)>ptime->at(i-1)+duration->at(i-1)+20. &&
	      duration->at(i)>7.) { // Clean afterpulse
	      h_AfterPulseArea[int(chan->at(i))]->Fill(float(area->at(i))/1000.,1.);
	      h_AfterPulseHeight[int(chan->at(i))]->Fill(float(height->at(i)),1.);
        h_AfterPulseDuration[int(chan->at(i))]->Fill(float(duration->at(i)),1.);
        h_AfterPulseTime[int(chan->at(i))]->Fill(float(ptime->at(i))-firstTime[int(chan->at(i))],1.);
	    } // Clean afterpulse?
	   } // Not first pulse in this channel?
    } // Pulse loop

    // Plot charge, height, and duration of tagged channels.
    // Tagged channels are those with a cosmic partner containing a hit.

    bool NeighborHit[32]; // True for each channel if there is a hit in its neighbor
    for (int c=0; c<32; c++)
	  NeighborHit[c] = false;

    // The channel mapping changed over time, e.g., in TS1, so the code below is run dependent.

    if (RunNum <= 863) { // Pre TS1
     // Cosmic tags; loose selection for Cosmics and tighter selection later for CosmicTag
     if (EventCategory == "CosmicTag") { // Tighter selection for CosmicTag category
      // Here we use tight requirements of 3 hits around the channel (except slabs, which require 2) to obtain a pure sample.
      if ((HitChan[24] && HitChan[8] && HitChan[10]) ||
          (run>817 && HitChan[24] && HitChan[10])) // CH8 dies after run 817
	      NeighborHit[0] = true;
      if (HitChan[25] && HitChan[9] && HitChan[10])
 	      NeighborHit[1] = true;
      if (HitChan[22] && HitChan[4] && HitChan[14]) 
	      NeighborHit[2] = true ;
      if (HitChan[23] && HitChan[5] && HitChan[14])
	      NeighborHit[3] = true;
      if (HitChan[2] && HitChan[22] && HitChan[14])
	      NeighborHit[4] = true;
      if (HitChan[3] && HitChan[23] && HitChan[14])
	      NeighborHit[5] = true;
      if (HitChan[16] && HitChan[12] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)))
	      NeighborHit[6] = true;
      if (HitChan[17] && HitChan[13] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)))
	      NeighborHit[7] = true;
      if (HitChan[0] && HitChan[24] && HitChan[10])
	      NeighborHit[8] = true;
      if (HitChan[1] && HitChan[25] && HitChan[10])
        NeighborHit[9] = true;
      if (HitChan[6] && HitChan[16] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)))
        NeighborHit[12] = true;
      if (HitChan[7] && HitChan[17] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)))
        NeighborHit[13] = true;
      if (HitChan[6] && HitChan[12] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)))
        NeighborHit[16] = true;
      if (HitChan[7] && HitChan[13] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)))
        NeighborHit[17] = true;
      if (HitChan[2] && HitChan[4] && HitChan[14])
        NeighborHit[22] = true;
      if (HitChan[3] && HitChan[5] && HitChan[14])
        NeighborHit[23] = true;
      if ((HitChan[0] && HitChan[8] && HitChan[10]) ||
          (run>817 && HitChan[0] && HitChan[10])) // CH8 dies after run 817
        NeighborHit[24] = true;
      if (HitChan[1] && HitChan[9] && HitChan[10])
        NeighborHit[25] = true;

      // Top veto panels
      if ((HitChan[0]&&HitChan[24]&&HitChan[8]) || (HitChan[1]&&HitChan[25]&&HitChan[9]))
        NeighborHit[10] = true; // Panel L0 top 
      if ((HitChan[6]&&HitChan[16]&&HitChan[12]) || (HitChan[7]&&HitChan[17]&&HitChan[13]))
         NeighborHit[30] = true; // Panel L1 top
      if ((HitChan[2]&&HitChan[22]&&HitChan[4]) || (HitChan[3]&&HitChan[23]&&HitChan[5]))
         NeighborHit[14] = true; // Panel L2 top

      // Side veto panels
      if (HitChan[29] && 
          (HitChan[0]||HitChan[24]||HitChan[8]) && 
          (HitChan[1]||HitChan[25]||HitChan[9]))
        NeighborHit[27] = true; // Panel L0 left
      if (HitChan[27] && 
          (HitChan[0]||HitChan[24]||HitChan[8]) && 
          (HitChan[1]||HitChan[25]||HitChan[9]))
        NeighborHit[29] = true; // Panel L0 right
      if (HitChan[19] && 
          (HitChan[6]||HitChan[16]||HitChan[12]) && 
          (HitChan[7]||HitChan[17]||HitChan[13]))
         NeighborHit[11] = true; // Panel L1 left
      if (((HitChan[30]&&RunNum>900)||(HitChan[11]&&RunNum<900)) && 
          (HitChan[6]||HitChan[16]||HitChan[12]) && 
          (HitChan[7]||HitChan[17]||HitChan[13]))
        NeighborHit[19] = true; // Panel L1 right
      if ((HitChan[31] && 
          (HitChan[2]||HitChan[22]||HitChan[4]) && 
          (HitChan[3]||HitChan[23]||HitChan[5])) ||
        	(run>833 && (HitChan[2]||HitChan[22]||HitChan[4]) && 
          (HitChan[3]||HitChan[23]||HitChan[5])))
       NeighborHit[26] = true; // Panel L2 right
      if (HitChan[26] && 
          (HitChan[2]||HitChan[22]||HitChan[4]) && 
          (HitChan[3]||HitChan[23]||HitChan[5]))
        NeighborHit[31] = true; // Pannel L2 left

      // Slabs
      if ( (HitChan[24] && HitChan[8]) || (HitChan[25] && HitChan[9]) || 
           (HitChan[0] && HitChan[24]) || (HitChan[1] && HitChan[25]) ) {
       NeighborHit[18] = true; // Slab1 
       NeighborHit[20] = true; // Slab0
      }
      if ( (HitChan[16] && HitChan[12]) || (HitChan[17] && HitChan[13]) || 
           (HitChan[6] && HitChan[16]) || (HitChan[7] && HitChan[17]) ) {
       NeighborHit[20] = true; // Slab0
       NeighborHit[28] = true; // Slab2
      }
      if ( (HitChan[22] && HitChan[4]) || (HitChan[23] && HitChan[5]) || 
           (HitChan[2] && HitChan[22]) || (HitChan[3] && HitChan[23]) ) {
       NeighborHit[28] = true; // Slab2
       NeighborHit[21] = true; // Slab3
      }
     } // Tighter selection for CosmicTag category

     else { // Looser selection for categories other than CosmicTag
      // Here we require ANY TWO hits in line with the probe channel
      // Used to obtain an efficient sample, and in particular one that tests trigger validity
      if ((HitChan[24] && HitChan[8]) || (HitChan[24] && HitChan[10]) || (HitChan[8] && HitChan[10]))
	      NeighborHit[0] = true;
      if ((HitChan[25] && HitChan[9]) || (HitChan[25] && HitChan[10]) || (HitChan[9] && HitChan[10]))
	      NeighborHit[1] = true;
      if ((HitChan[22] && HitChan[4]) || (HitChan[22] && HitChan[14]) || (HitChan[4] && HitChan[14]))
	      NeighborHit[2] = true;
      if ((HitChan[23] && HitChan[5]) || (HitChan[23] && HitChan[14]) || (HitChan[5] && HitChan[14]))
	      NeighborHit[3] = true;
      if ((HitChan[2] && HitChan[22]) || (HitChan[22] && HitChan[14]) || (HitChan[2] && HitChan[14]))
	      NeighborHit[4] = true;
      if ((HitChan[3] && HitChan[23]) || (HitChan[23] && HitChan[14]) || (HitChan[3] && HitChan[14]))
	      NeighborHit[5] = true;
      if ((HitChan[16] && HitChan[12]) || 
         (HitChan[16] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900))) || 
         (HitChan[12] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900))))
	      NeighborHit[6] = true;
      if ((HitChan[17] && HitChan[13]) || 
         (HitChan[17] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900))) || 
         (HitChan[13] && ((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900))))
	      NeighborHit[7] = true;
      if ((HitChan[0] && HitChan[24]) || (HitChan[10] && HitChan[24]) || (HitChan[10] && HitChan[0]))
	      NeighborHit[8] = true;
      if ((HitChan[1] && HitChan[25]) || (HitChan[10] && HitChan[25]) || (HitChan[10] && HitChan[1]))
        NeighborHit[9] = true;
      if ((HitChan[6] && HitChan[16]) || 
         (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[16]) || 
         (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[6]))
        NeighborHit[12] = true;
      if ((HitChan[7] && HitChan[17]) || 
          (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[17]) || 
          (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[7]))
        NeighborHit[13] = true;
      if ((HitChan[6] && HitChan[12]) || 
          (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[12]) || 
          (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[6]))
        NeighborHit[16] = true;
      if ((HitChan[7] && HitChan[13]) || 
          (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[13]) || 
          (((HitChan[30]&&RunNum<900)||(HitChan[11]&&RunNum>900)) && HitChan[7]))
        NeighborHit[17] = true;
      if ((HitChan[2] && HitChan[4]) || (HitChan[14] && HitChan[4]) || (HitChan[14] && HitChan[2]))
        NeighborHit[22] = true;
      if ((HitChan[3] && HitChan[5]) || (HitChan[14] && HitChan[5]) || (HitChan[14] && HitChan[3]))
        NeighborHit[23] = true;
      if ((HitChan[0] && HitChan[8]) || (HitChan[10] && HitChan[8]) || (HitChan[10] && HitChan[0]))
        NeighborHit[24] = true;
      if ((HitChan[1] && HitChan[9]) || (HitChan[10] && HitChan[9]) || (HitChan[10] && HitChan[1]))
        NeighborHit[25] = true;

      // Veto panels
      if ((HitChan[0]&&HitChan[24]) || (HitChan[0]&&HitChan[8])  || 
          (HitChan[24]&&HitChan[8]) || (HitChan[1]&&HitChan[25]) || 
          (HitChan[1]&&HitChan[9]) || (HitChan[25]&&HitChan[9]))
        NeighborHit[10] = true; // Panel L0 top
      if ((HitChan[6]&&HitChan[16]) || (HitChan[6]&&HitChan[12])  || 
          (HitChan[16]&&HitChan[12]) || (HitChan[7]&&HitChan[17]) || 
          (HitChan[7]&&HitChan[13]) || (HitChan[17]&&HitChan[13])) {
       if (RunNum<900) 
         NeighborHit[30] = true; // Panel L1 top
       else
         NeighborHit[11] = true; // Panel L1 top
      }
      if ((HitChan[2]&&HitChan[22]) || (HitChan[2]&&HitChan[4])  || 
          (HitChan[22]&&HitChan[4]) || (HitChan[3]&&HitChan[23]) || 
          (HitChan[3]&&HitChan[5]) || (HitChan[23]&&HitChan[5]))
        NeighborHit[14] = true; // Panel L2 top

      // Side veto panels
      if ((HitChan[0]||HitChan[24]||HitChan[8]) && (HitChan[1]||HitChan[25]||HitChan[9])) {
        NeighborHit[27] = true; // Panel L0 left
        NeighborHit[29] = true; // Panel L0 right
      }
      if ((HitChan[6]||HitChan[16]||HitChan[12]) && (HitChan[7]||HitChan[17]||HitChan[13])) {
        NeighborHit[19] = true; // Panel L1 right
        if (RunNum<900) 
         NeighborHit[11] = true; // Panel L1 left
        else
         NeighborHit[30] = true; // Panel L1 left
      }
      if ((HitChan[2]||HitChan[22]||HitChan[4]) && (HitChan[3]||HitChan[23]||HitChan[5])) {
        NeighborHit[31] = true; // Panel L2 left
        NeighborHit[26] = true; // Panel L2 right
      }

      // Slabs
      if ( (HitChan[24] && HitChan[8]) || (HitChan[25] && HitChan[9]) || 
           (HitChan[0] && HitChan[24]) || (HitChan[1] && HitChan[25]) ) {
        NeighborHit[18] = true; // Slab 0
        NeighborHit[20] = true; // Slab 1
      }
      if ( (HitChan[16] && HitChan[12]) || (HitChan[17] && HitChan[13]) || 
           (HitChan[6] && HitChan[16]) || (HitChan[7] && HitChan[17]) ) {
        NeighborHit[20] = true; // Slab 1
        NeighborHit[28] = true; // Slab 2
      }
      if ( (HitChan[22] && HitChan[4]) || (HitChan[23] && HitChan[5]) || 
          (HitChan[2] && HitChan[22]) || (HitChan[3] && HitChan[23]) ) {
        NeighborHit[28] = true; // Slab 2
        NeighborHit[21] = true; // Slab 3
      }
     } // Looser selection for categories other than CosmicTag
    } // Run <= 863 (pre TS1)

    else if (RunNum>=864 && RunNum<911) { // Special selection for TS1 runs
      // Here we require ANY two hits in the L1 channels placed below L2
      if ((HitChan[10] && HitChan[27]) || 
          ((HitChan[10]||HitChan[27]) && (HitChan[0]||HitChan[1]||HitChan[24]||HitChan[25]))) {
        NeighborHit[0] = true;
        NeighborHit[1] = true;
        NeighborHit[2] = true;
        NeighborHit[3] = true;
        NeighborHit[4] = true;
        NeighborHit[5] = true;
        NeighborHit[6] = true;
        NeighborHit[7] = true;
        NeighborHit[8] = true;
        NeighborHit[9] = true;
        NeighborHit[10] = true;
        NeighborHit[11] = true;
        NeighborHit[12] = true;
        NeighborHit[13] = true;
        NeighborHit[14] = true;
        NeighborHit[16] = true;
        NeighborHit[17] = true;
        NeighborHit[18] = true;
        NeighborHit[19] = true;
        NeighborHit[20] = true;
        NeighborHit[21] = true;
        NeighborHit[22] = true;
        NeighborHit[23] = true;
        NeighborHit[24] = true;
        NeighborHit[25] = true;
        NeighborHit[26] = true;
        NeighborHit[27] = true;
        NeighborHit[28] = true;
        NeighborHit[29] = true;
        NeighborHit[30] = true;
        NeighborHit[31] = true;
     }
    } // Special selection for TS1 runs

    else { // Post TS1 runs?
     // Cosmic tags; loose selection for Cosmics and tighter selection later for CosmicTag
     if (EventCategory == "CosmicTag") { // Tighter selection for CosmicTag category
      // Here we use tight requirements of 3 hits around the channel (except slabs, which require 2) to obtain a pure sample.
      if (HitChan[24] && HitChan[8] && HitChan[10])
        NeighborHit[0] = true;
      if (HitChan[25] && HitChan[9] && HitChan[10])
        NeighborHit[1] = true;
      if (HitChan[22] && HitChan[4] && HitChan[14]) 
        NeighborHit[2] = true ;
      if (HitChan[23] && HitChan[5] && HitChan[14])
        NeighborHit[3] = true;
      if (HitChan[2] && HitChan[22] && HitChan[14])
        NeighborHit[4] = true;
      if (HitChan[3] && HitChan[23] && HitChan[14])
        NeighborHit[5] = true;
      if (HitChan[16] && HitChan[12] && HitChan[11])
        NeighborHit[6] = true;
      if (HitChan[17] && HitChan[13] && HitChan[11])
        NeighborHit[7] = true;
      if (HitChan[0] && HitChan[24] && HitChan[10])
        NeighborHit[8] = true;
      if (HitChan[1] && HitChan[25] && HitChan[10])
        NeighborHit[9] = true;
      if (HitChan[6] && HitChan[16] && HitChan[11])
        NeighborHit[12] = true;
      if (HitChan[7] && HitChan[17] && HitChan[11])
        NeighborHit[13] = true;
      if (HitChan[6] && HitChan[12] && HitChan[11])
        NeighborHit[16] = true;
      if (HitChan[7] && HitChan[13] && HitChan[11])
        NeighborHit[17] = true;
      if (HitChan[2] && HitChan[4] && HitChan[14])
        NeighborHit[22] = true;
      if (HitChan[3] && HitChan[5] && HitChan[14])
        NeighborHit[23] = true;
      if (HitChan[0] && HitChan[8] && HitChan[10])
        NeighborHit[24] = true;
      if (HitChan[1] && HitChan[9] && HitChan[10])
        NeighborHit[25] = true;

      // Top veto panels
      if ((HitChan[0]&&HitChan[24]&&HitChan[8]) || (HitChan[1]&&HitChan[25]&&HitChan[9]))
        NeighborHit[10] = true; // Panel L0 top 
      if ((HitChan[6]&&HitChan[16]&&HitChan[12]) || (HitChan[7]&&HitChan[17]&&HitChan[13]))
         NeighborHit[11] = true; // Panel L1 top
      if ((HitChan[2]&&HitChan[22]&&HitChan[4]) || (HitChan[3]&&HitChan[23]&&HitChan[5]))
         NeighborHit[14] = true; // Panel L2 top

      // Side veto panels
      if (HitChan[29] && 
          (HitChan[0]||HitChan[24]||HitChan[8]) && 
          (HitChan[1]||HitChan[25]||HitChan[9]))
        NeighborHit[27] = true; // Panel L0 left
      if (HitChan[27] && 
          (HitChan[0]||HitChan[24]||HitChan[8]) && 
          (HitChan[1]||HitChan[25]||HitChan[9]))
        NeighborHit[29] = true; // Panel L0 right
      if (HitChan[19] && 
          (HitChan[6]||HitChan[16]||HitChan[12]) && 
          (HitChan[7]||HitChan[17]||HitChan[13]))
        NeighborHit[30] = true; // Panel L1 left
      if (HitChan[30] && 
          (HitChan[6]||HitChan[16]||HitChan[12]) && 
          (HitChan[7]||HitChan[17]||HitChan[13]))
        NeighborHit[19] = true; // Panel L1 right
      if ((HitChan[31] && 
          (HitChan[2]||HitChan[22]||HitChan[4]) && 
          (HitChan[3]||HitChan[23]||HitChan[5])) ||
          (run>833 && (HitChan[2]||HitChan[22]||HitChan[4]) && 
          (HitChan[3]||HitChan[23]||HitChan[5])))
        NeighborHit[26] = true; // Panel L2 right
      if (HitChan[26] && 
          (HitChan[2]||HitChan[22]||HitChan[4]) && 
          (HitChan[3]||HitChan[23]||HitChan[5]))
        NeighborHit[31] = true; // Pannel L2 left

      // Slabs
      if ( (HitChan[24] && HitChan[8] && !HitChan[0] && !HitChan[1] && !HitChan[10]) || 
           (HitChan[25] && HitChan[9] && !HitChan[0] && !HitChan[1] && !HitChan[10]) || 
           (HitChan[10] && HitChan[0] && HitChan[24] && !HitChan[8]) || 
           (HitChan[10] && HitChan[1] && HitChan[25] && !HitChan[9]) ) {
        NeighborHit[18] = true; // Slab1 
        NeighborHit[20] = true; // Slab0
      }
      if ( (HitChan[16] && HitChan[12] && !HitChan[6] && !HitChan[11]) || 
           (HitChan[17] && HitChan[13] && !HitChan[7] && !HitChan[11]) || 
           (HitChan[6] && HitChan[16] && HitChan[11] && !HitChan[12]) || 
           (HitChan[7] && HitChan[17] && HitChan[11] && !HitChan[13]) ) {
        NeighborHit[20] = true; // Slab0
        NeighborHit[28] = true; // Slab2
      }
      if ( (HitChan[22] && HitChan[4] && !HitChan[2] && !HitChan[14]) || 
           (HitChan[23] && HitChan[5] && !HitChan[3] && !HitChan[14]) || 
           (HitChan[2] && HitChan[22] && HitChan[14]) || 
           (HitChan[3] && HitChan[23] && HitChan[14]) ) {
       NeighborHit[28] = true; // Slab2
       NeighborHit[21] = true; // Slab3
      }
     } // Tighter selection for CosmicTag category

     else { // Looser selection for categories other than CosmicTag
      // Here we require ANY TWO hits in line with the probe channel
      // Used to obtain an efficient sample, and in particular one that tests trigger validity
      if ((HitChan[24] && HitChan[8]) || (HitChan[24] && HitChan[10]) || (HitChan[8] && HitChan[10]))
        NeighborHit[0] = true;
      if ((HitChan[25] && HitChan[9]) || (HitChan[25] && HitChan[10]) || (HitChan[9] && HitChan[10]))
        NeighborHit[1] = true;
      if ((HitChan[22] && HitChan[4]) || (HitChan[22] && HitChan[14]) || (HitChan[4] && HitChan[14]))
        NeighborHit[2] = true;
      if ((HitChan[23] && HitChan[5]) || (HitChan[23] && HitChan[14]) || (HitChan[5] && HitChan[14]))
        NeighborHit[3] = true;
      if ((HitChan[2] && HitChan[22]) || (HitChan[22] && HitChan[14]) || (HitChan[2] && HitChan[14]))
        NeighborHit[4] = true;
      if ((HitChan[3] && HitChan[23]) || (HitChan[23] && HitChan[14]) || (HitChan[3] && HitChan[14]))
        NeighborHit[5] = true;
      if ((HitChan[16] && HitChan[12]) || 
          (HitChan[16] && HitChan[11]) || 
          (HitChan[12] && HitChan[11]))
        NeighborHit[6] = true;
      if ((HitChan[17] && HitChan[13]) || 
          (HitChan[17] && HitChan[11]) || 
          (HitChan[13] && HitChan[11]))
        NeighborHit[7] = true;
      if ((HitChan[0] && HitChan[24]) || 
          (HitChan[10] && HitChan[24]) || 
          (HitChan[10] && HitChan[0]))
        NeighborHit[8] = true;
      if ((HitChan[1] && HitChan[25]) || 
          (HitChan[10] && HitChan[25]) || 
          (HitChan[10] && HitChan[1]))
        NeighborHit[9] = true;
      if ((HitChan[6] && HitChan[16]) || 
          (HitChan[11] && HitChan[16]) || 
          (HitChan[11] && HitChan[6]))
        NeighborHit[12] = true;
      if ((HitChan[7] && HitChan[17]) || 
          (HitChan[11] && HitChan[17]) || 
          (HitChan[11] && HitChan[7]))
        NeighborHit[13] = true;
      if ((HitChan[6] && HitChan[12]) || 
          (HitChan[11] && HitChan[12]) || 
          (HitChan[11] && HitChan[6]))
        NeighborHit[16] = true;
      if ((HitChan[7] && HitChan[13]) || 
          (HitChan[11] && HitChan[13]) || 
          (HitChan[11] && HitChan[7]))
        NeighborHit[17] = true;
      if ((HitChan[2] && HitChan[4]) || (HitChan[14] && HitChan[4]) || (HitChan[14] && HitChan[2]))
        NeighborHit[22] = true;
      if ((HitChan[3] && HitChan[5]) || (HitChan[14] && HitChan[5]) || (HitChan[14] && HitChan[3]))
        NeighborHit[23] = true;
      if ((HitChan[0] && HitChan[8]) || (HitChan[10] && HitChan[8]) || (HitChan[10] && HitChan[0]))
        NeighborHit[24] = true;
      if ((HitChan[1] && HitChan[9]) || (HitChan[10] && HitChan[9]) || (HitChan[10] && HitChan[1]))
        NeighborHit[25] = true;

      // Veto panels
      if ((HitChan[0]&&HitChan[24]) || (HitChan[0]&&HitChan[8])  || 
          (HitChan[24]&&HitChan[8]) || (HitChan[1]&&HitChan[25]) || 
          (HitChan[1]&&HitChan[9]) || (HitChan[25]&&HitChan[9]))
        NeighborHit[10] = true; // Panel L0 top
      if ((HitChan[6]&&HitChan[16]) || (HitChan[6]&&HitChan[12])  || 
          (HitChan[16]&&HitChan[12]) || (HitChan[7]&&HitChan[17]) || 
          (HitChan[7]&&HitChan[13]) || (HitChan[17]&&HitChan[13])) {
        NeighborHit[11] = true; // Panel L1 top
      }
      if ((HitChan[2]&&HitChan[22]) || (HitChan[2]&&HitChan[4])  || 
          (HitChan[22]&&HitChan[4]) || (HitChan[3]&&HitChan[23]) || 
          (HitChan[3]&&HitChan[5]) || (HitChan[23]&&HitChan[5]))
        NeighborHit[14] = true; // Panel L2 top

      // Side veto panels
      if ((HitChan[0]||HitChan[24]||HitChan[8]) && (HitChan[1]||HitChan[25]||HitChan[9])) {
        NeighborHit[27] = true; // Panel L0 left
        NeighborHit[29] = true; // Panel L0 right
      }
      if ((HitChan[6]||HitChan[16]||HitChan[12]) && (HitChan[7]||HitChan[17]||HitChan[13])) {
        NeighborHit[19] = true; // Panel L1 right
        NeighborHit[30] = true; // Panel L1 left
      }
      if ((HitChan[2]||HitChan[22]||HitChan[4]) && (HitChan[3]||HitChan[23]||HitChan[5])) {
        NeighborHit[31] = true; // Panel L2 left
        NeighborHit[26] = true; // Panel L2 right
      }

      // Slabs
      if ( (HitChan[24] && HitChan[8]) || (HitChan[25] && HitChan[9]) || 
           (HitChan[0] && HitChan[24]) || (HitChan[1] && HitChan[25]) ) {
        NeighborHit[18] = true; // Slab 0
        NeighborHit[20] = true; // Slab 1
      }
      if ( (HitChan[16] && HitChan[12]) || (HitChan[17] && HitChan[13]) || 
           (HitChan[6] && HitChan[16]) || (HitChan[7] && HitChan[17]) ) {
        NeighborHit[20] = true; // Slab 1
        NeighborHit[28] = true; // Slab 2
      }
      if ( (HitChan[22] && HitChan[4]) || (HitChan[23] && HitChan[5]) || 
          (HitChan[2] && HitChan[22]) || (HitChan[3] && HitChan[23]) ) {
        NeighborHit[28] = true; // Slab 2
        NeighborHit[21] = true; // Slab 3
      }
     } // Looser selection for categories other than CosmicTag
    } // Post TS1 runs?

 	  // Determine which layers and slabs are hit
    bool layHit[3]; bool slabHit[4]; 
    layHit[0] = false; layHit[1] = false; layHit[2] = false;
    slabHit[0] = false; slabHit[1] = false; slabHit[2] = false; slabHit[3] = false; 
    if (HitChan[18]) slabHit[0] = true;
    if (HitChan[20]) slabHit[1] = true;
    if (HitChan[28]) slabHit[2] = true;
    if (HitChan[21]) slabHit[3] = true;
    if (HitChan[0] || HitChan[1] || HitChan[24] || HitChan[25] || HitChan[8] || HitChan[9]) layHit[0] = true;
    if (HitChan[6] || HitChan[7] || HitChan[16] || HitChan[17] || HitChan[12] || HitChan[13]) layHit[1] = true;
    if (HitChan[2] || HitChan[3] || HitChan[22] || HitChan[23] || HitChan[4] || HitChan[5]) layHit[2] = true;
    int nSlabHits = 0;
    for (int i=0; i<4; i++)
    	if (slabHit[i])
    		nSlabHits++;
    	int nLayHits = 0;
    	for (int i=0; i<3; i++)
    		if (layHit[i])
    			nLayHits++;

    // Plot number of layers hit
    h_NHitLayers->Fill(1.*nLayHits,1.);

    // Dump some crude event displays for cosmics
    bool showDisplay = true;
    if (showDisplay && evtCategoryCode==3 && (HitChan[10] || HitChan[11] || HitChan[14])) { // Cosmic step
      if (nHits>=3) { // At least 3 hits?
        cout << "\n";
        cout << "Display: file="<<filenum<<" event="<<event<<"\n";
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L2: ";
        cout << "     ["<<HitChan[21] <<"]\n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L2: ";
        cout << "  --" << HitChan[14]<< HitChan[14]<< HitChan[14]<< HitChan[14]<< HitChan[14]<<"--\n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L2: ";
        cout << "  /"<<HitChan[31] << "/" << HitChan[2] <<" "<< HitChan[3] <<"\\"<<HitChan[26]<<"\\ \n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L2: ";
        cout << " /"<<HitChan[31] << "/ " << HitChan[22] <<" "<< HitChan[23] <<" \\"<<HitChan[26]<<"\\ \n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L2: ";
        cout << "/"<<HitChan[31] << "/  " << HitChan[4] <<" "<< HitChan[5] <<"  \\"<<HitChan[26]<<"\\ \n"; 

        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L1: ";
        cout << "     ["<<HitChan[28] <<"]\n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L1: ";
        cout << "  --" << HitChan[11]<< HitChan[11]<< HitChan[11]<< HitChan[11]<< HitChan[11]<<"--\n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L1: ";
        cout << "  /"<<HitChan[30] << "/" << HitChan[6] <<" "<< HitChan[7] <<"\\"<<HitChan[19]<<"\\ \n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L1: ";
        cout << " /"<<HitChan[30] << "/ " << HitChan[16] <<" "<< HitChan[17] <<" \\"<<HitChan[19]<<"\\ \n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L1: ";
        cout << "/"<<HitChan[30] << "/  " << HitChan[12] <<" "<< HitChan[13] <<"  \\"<<HitChan[19]<<"\\ \n"; 

        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L0: ";
        cout << "     ["<<HitChan[20] <<"]\n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L0: ";
        cout << "  --" << HitChan[10]<< HitChan[10]<< HitChan[10]<< HitChan[10]<< HitChan[10]<<"--\n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L0: ";
        cout << "  /"<<HitChan[27] << "/" << HitChan[0] <<" "<< HitChan[1] <<"\\"<<HitChan[29]<<"\\ \n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L0: ";
        cout << " /"<<HitChan[27] << "/ " << HitChan[24] <<" "<< HitChan[25] <<" \\"<<HitChan[29]<<"\\ \n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L0: ";
        cout << "/"<<HitChan[27] << "/  " << HitChan[8] <<" "<< HitChan[9] <<"  \\"<<HitChan[29]<<"\\ \n"; 
        cout << "Display: nS="<<nSlabHits<<" nL="<<nLayHits<<" "<<layHit[0]<<layHit[1]<<layHit[2]<<" "<<filenum<<" "<<event<<" L0: ";
        cout << "     ["<<HitChan[18] <<"]\n"; 
      } // nHits>=4
    } // Cosmic step

    // Plots of max sample for tagged channels
    for (int c=0; c<32; c++)
     if (NeighborHit[c] && maxSample[c]>0.5*VThresholds[c][evtCategoryCode]) h_TaggedMax[c]->Fill(maxSample[c],1.);

    // Now make plots of pulse area, height, and duration for tagged channels.
    float TotalArea[32];
    for (int c=0; c<32; c++)
	 TotalArea[c] = 0.; // Sum of area of pulses in the channel
    for (int i=0; i<int(chan->size()); i++) { // Pulse loop 
	 TotalArea[int(chan->at(i))] += area->at(i)/1000.;
	 if (ipulse->at(i) == 0 ) { // First pulse in this channel
	  if ((NeighborHit[int(chan->at(i))]) &&
	      (((area->at(i))/1000.>0.5*AThresholds[chan->at(i)][evtCategoryCode])) &&
	      (height->at(i)>0.5*VThresholds[chan->at(i)][evtCategoryCode]) &&
	      (duration->at(i)>DThresholds[chan->at(i)][evtCategoryCode])) { // Tagged by neighbor hit?
	    h_TaggedPulseArea[int(chan->at(i))]->Fill(float(area->at(i))/1000.,1.);
	    h_TaggedPulseHeight[int(chan->at(i))]->Fill(float(height->at(i)),1.);
	    h_TaggedDuration[int(chan->at(i))]->Fill(float(duration->at(i)),1.);
	  } // Tagged by neighbor hit?
	 } // First pulse in this channel?
    } // Pluse loop

    // Histograms of total pulse area for all channels
    for (int c=0; c<32; c++)
	 if (TotalArea[c]>1.E-5) h_TotalPulseArea[c]->Fill(TotalArea[c],1.);

    // Histograms of total pulse area for channels with a tagging neighbor hit
	for (int c=0; c<32; c++)
	 if (NeighborHit[c] && TotalArea[c]>0.5*AThresholds[c][evtCategoryCode])
 		h_TaggedTotalPulseArea[c]->Fill(TotalArea[c],1.);

    // Plot max pulse height for channels when they have hits in the two straight line interlayer partners.
    // The threshold for the other partners varies with evtCategoryCode as 10mV*evtCategoryCode.
    // There is no timing cuts on this comparison.
    bool partnersHit[32]; // True for each channel if there is a hit in its interlayer partners
    for (int c=0; c<32; c++)
      partnersHit[c] = false;
    float partnerThreshold = 10.*(1+evtCategoryCode); 
    if (maxSample[6]>partnerThreshold && maxSample[14]>partnerThreshold) partnersHit[0] = true;
    if (maxSample[7]>partnerThreshold && maxSample[15]>partnerThreshold) partnersHit[1] = true;
    if (maxSample[12]>partnerThreshold && maxSample[10]>partnerThreshold) partnersHit[2] = true;
    if (maxSample[13]>partnerThreshold && maxSample[11]>partnerThreshold) partnersHit[3] = true;
    if (maxSample[8]>partnerThreshold && maxSample[10]>partnerThreshold) partnersHit[4] = true;
    if (maxSample[9]>partnerThreshold && maxSample[11]>partnerThreshold) partnersHit[5] = true;
    if (maxSample[0]>partnerThreshold && maxSample[14]>partnerThreshold) partnersHit[6] = true;
    if (maxSample[1]>partnerThreshold && maxSample[15]>partnerThreshold) partnersHit[7] = true;
    if (maxSample[4]>partnerThreshold && maxSample[10]>partnerThreshold) partnersHit[8] = true;
    if (maxSample[5]>partnerThreshold && maxSample[11]>partnerThreshold) partnersHit[9] = true;
    if (maxSample[2]>partnerThreshold && maxSample[12]>partnerThreshold) partnersHit[10] = true;
    if (maxSample[3]>partnerThreshold && maxSample[13]>partnerThreshold) partnersHit[11] = true;
    if (maxSample[4]>partnerThreshold && maxSample[10]>partnerThreshold) partnersHit[12] = true;
    if (maxSample[3]>partnerThreshold && maxSample[11]>partnerThreshold) partnersHit[13] = true;
    if (maxSample[0]>partnerThreshold && maxSample[6]>partnerThreshold) partnersHit[14] = true;
    if (maxSample[1]>partnerThreshold && maxSample[7]>partnerThreshold) partnersHit[15] = true;

    // Plots of max sample for channels with 3-layer partner hits
    for (int c=0; c<32; c++)
     if (partnersHit[c] && maxSample[c]>7.) h_TripleLayerMax[c]->Fill(maxSample[c],1.);

    // Plot max pulse height for channels when they have hits in either of the nearest interlayer partners.
    // The threshold for the other partners varies with evtCategoryCode as 10mV*evtCategoryCode.
    // There is no timing cuts on this comparison.
    for (int c=0; c<32; c++)
      partnersHit[c] = false;
    partnerThreshold = 10.*(evtCategoryCode+1); 
    if (maxSample[6]>partnerThreshold || maxSample[14]>partnerThreshold) partnersHit[0] = true;
    if (maxSample[7]>partnerThreshold || maxSample[15]>partnerThreshold) partnersHit[1] = true;
    if (maxSample[12]>partnerThreshold || maxSample[10]>partnerThreshold) partnersHit[2] = true;
    if (maxSample[13]>partnerThreshold || maxSample[11]>partnerThreshold) partnersHit[3] = true;
    if (maxSample[8]>partnerThreshold || maxSample[10]>partnerThreshold) partnersHit[4] = true;
    if (maxSample[9]>partnerThreshold || maxSample[11]>partnerThreshold) partnersHit[5] = true;
    if (maxSample[0]>partnerThreshold || maxSample[14]>partnerThreshold) partnersHit[6] = true;
    if (maxSample[1]>partnerThreshold || maxSample[15]>partnerThreshold) partnersHit[7] = true;
    if (maxSample[4]>partnerThreshold || maxSample[10]>partnerThreshold) partnersHit[8] = true;
    if (maxSample[5]>partnerThreshold || maxSample[11]>partnerThreshold) partnersHit[9] = true;
    if (maxSample[2]>partnerThreshold || maxSample[12]>partnerThreshold) partnersHit[10] = true;
    if (maxSample[3]>partnerThreshold || maxSample[13]>partnerThreshold) partnersHit[11] = true;
    if (maxSample[4]>partnerThreshold || maxSample[10]>partnerThreshold) partnersHit[12] = true;
    if (maxSample[3]>partnerThreshold || maxSample[11]>partnerThreshold) partnersHit[13] = true;
    if (maxSample[0]>partnerThreshold || maxSample[6]>partnerThreshold) partnersHit[14] = true;
    if (maxSample[1]>partnerThreshold || maxSample[7]>partnerThreshold) partnersHit[15] = true;

    // Plots of max sample for channels with 2-layer partner hits
    for (int c=0; c<32; c++)
      if (partnersHit[c] && maxSample[c]>7.) h_DoubleLayerMax[c]->Fill(maxSample[c],1.);

 //    // Plot the time difference between hits in a channel and the time of hits in neighboring channels.
 //    // Fill time difference of vertical neighbors -- hard coded geometry.
 //    // Channel 1 and 3
 //    if (ChanHitIndex[1]>=0 && ChanHitIndex[3]>=0) {
	//  h_VertNeighborDeltaTime_1->Fill(float(time_module_calibrated->at(ChanHitIndex[1])-time_module_calibrated->at(ChanHitIndex[3])),1.);
	//  h_VertNeighborDeltaTime_3->Fill(float(time_module_calibrated->at(ChanHitIndex[3])-time_module_calibrated->at(ChanHitIndex[1])),1.);
 //    }
 //    // Channel 6 and 4
 //    if (ChanHitIndex[6]>=0 && ChanHitIndex[4]>=0) {
	//  h_VertNeighborDeltaTime_6->Fill(float(time_module_calibrated->at(ChanHitIndex[6])-time_module_calibrated->at(ChanHitIndex[4])),1.);
	//  h_VertNeighborDeltaTime_4->Fill(float(time_module_calibrated->at(ChanHitIndex[4])-time_module_calibrated->at(ChanHitIndex[6])),1.);
 //    }
 //      // Channel 7 and 5
 //      if (ChanHitIndex[7]>=0 && ChanHitIndex[5]>=0) {
	// h_VertNeighborDeltaTime_7->Fill(float(time_module_calibrated->at(ChanHitIndex[7])-time_module_calibrated->at(ChanHitIndex[5])),1.);
	// h_VertNeighborDeltaTime_5->Fill(float(time_module_calibrated->at(ChanHitIndex[5])-time_module_calibrated->at(ChanHitIndex[7])),1.);
 //      }
 //      // Channel 10 and 8
 //      if (ChanHitIndex[10]>=0 && ChanHitIndex[8]>=0) {
	// h_VertNeighborDeltaTime_10->Fill(float(time_module_calibrated->at(ChanHitIndex[10])-time_module_calibrated->at(ChanHitIndex[8])),1.);
	// h_VertNeighborDeltaTime_8->Fill(float(time_module_calibrated->at(ChanHitIndex[8])-time_module_calibrated->at(ChanHitIndex[10])),1.);
 //      }
 //      // Channel 11 and 9
 //      if (ChanHitIndex[11]>=0 && ChanHitIndex[9]>=0) {
	// h_VertNeighborDeltaTime_11->Fill(float(time_module_calibrated->at(ChanHitIndex[11])-time_module_calibrated->at(ChanHitIndex[9])),1.);
	// h_VertNeighborDeltaTime_9->Fill(float(time_module_calibrated->at(ChanHitIndex[9])-time_module_calibrated->at(ChanHitIndex[11])),1.);
 //      }
 //      // Fill difference of angled vertical neighbors -- hard coded geometry.
 //      // Channel 1 and 2
 //      if (ChanHitIndex[1]>=0 && ChanHitIndex[2]>=0) {
	// h_AngleNeighborDeltaTime_1->Fill(float(time_module_calibrated->at(ChanHitIndex[1])-time_module_calibrated->at(ChanHitIndex[2])),1.);
	// h_AngleNeighborDeltaTime_2->Fill(float(time_module_calibrated->at(ChanHitIndex[2])-time_module_calibrated->at(ChanHitIndex[1])),1.);
 //      }
 //      // Channel 3 and 15. Only record in the sidecar channel
 //      if (ChanHitIndex[3]>=0 && ChanHitIndex[15]>=0) {
	// h_AngleNeighborDeltaTime_15->Fill(float(time_module_calibrated->at(ChanHitIndex[15])-time_module_calibrated->at(ChanHitIndex[3])),1.);
 //      }
 //      // Channel 6 and 5
 //      if (ChanHitIndex[6]>=0 && ChanHitIndex[5]>=0) {
	// h_AngleNeighborDeltaTime_6->Fill(float(time_module_calibrated->at(ChanHitIndex[6])-time_module_calibrated->at(ChanHitIndex[5])),1.);
	// h_AngleNeighborDeltaTime_5->Fill(float(time_module_calibrated->at(ChanHitIndex[5])-time_module_calibrated->at(ChanHitIndex[6])),1.);
 //      }
 //      // Channel 7 and 4
 //      if (ChanHitIndex[7]>=0 && ChanHitIndex[4]>=0) {
	// h_AngleNeighborDeltaTime_7->Fill(float(time_module_calibrated->at(ChanHitIndex[7])-time_module_calibrated->at(ChanHitIndex[4])),1.);
	// h_AngleNeighborDeltaTime_4->Fill(float(time_module_calibrated->at(ChanHitIndex[4])-time_module_calibrated->at(ChanHitIndex[7])),1.);
 //      }
 //      // Channel 7 and 12. Only record in the sidecar channel
 //      if (ChanHitIndex[7]>=0 && ChanHitIndex[12]>=0) {
	// h_AngleNeighborDeltaTime_12->Fill(float(time_module_calibrated->at(ChanHitIndex[12])-time_module_calibrated->at(ChanHitIndex[7])),1.);
 //      }
 //      // Channel 10 and 9
 //      if (ChanHitIndex[10]>=0 && ChanHitIndex[9]>=0) {
	// h_AngleNeighborDeltaTime_10->Fill(float(time_module_calibrated->at(ChanHitIndex[10])-time_module_calibrated->at(ChanHitIndex[9])),1.);
	// h_AngleNeighborDeltaTime_9->Fill(float(time_module_calibrated->at(ChanHitIndex[9])-time_module_calibrated->at(ChanHitIndex[10])),1.);
 //      }
 //      // Channel 11 and 8
 //      if (ChanHitIndex[11]>=0 && ChanHitIndex[8]>=0) {
	// h_AngleNeighborDeltaTime_11->Fill(float(time_module_calibrated->at(ChanHitIndex[11])-time_module_calibrated->at(ChanHitIndex[8])),1.);
	// h_AngleNeighborDeltaTime_8->Fill(float(time_module_calibrated->at(ChanHitIndex[8])-time_module_calibrated->at(ChanHitIndex[11])),1.);
 //      }
 //      // Channel 11 and 14; Only record in the sidecar channel
 //      if (ChanHitIndex[11]>=0 && ChanHitIndex[14]>=0) {
	// h_AngleNeighborDeltaTime_14->Fill(float(time_module_calibrated->at(ChanHitIndex[14])-time_module_calibrated->at(ChanHitIndex[11])),1.);
 //      }

 //      // Fill difference of layer neighbors -- hard coded geometry.
 //      // Only fill in the lower layer channel to avoid ambiguity
 //      // Channel 2 and 4
 //      if (ChanHitIndex[2]>=0 && ChanHitIndex[4]>=0) {
	// h_LayerNeighborDeltaTime_2->Fill(float(time_module_calibrated->at(ChanHitIndex[2])-time_module_calibrated->at(ChanHitIndex[4])),1.);
 //      }
 //      // Channel 3 and 5
 //      if (ChanHitIndex[3]>=0 && ChanHitIndex[5]>=0) {
	// h_LayerNeighborDeltaTime_3->Fill(float(time_module_calibrated->at(ChanHitIndex[3])-time_module_calibrated->at(ChanHitIndex[5])),1.);
 //      }
 //      // Channel 1 and 7
 //      if (ChanHitIndex[1]>=0 && ChanHitIndex[7]>=0) {
	// h_LayerNeighborDeltaTime_1->Fill(float(time_module_calibrated->at(ChanHitIndex[1])-time_module_calibrated->at(ChanHitIndex[7])),1.);
 //      }
 //      // Channel 15 and 12
 //      if (ChanHitIndex[15]>=0 && ChanHitIndex[12]>=0) {
	// h_LayerNeighborDeltaTime_15->Fill(float(time_module_calibrated->at(ChanHitIndex[15])-time_module_calibrated->at(ChanHitIndex[12])),1.);
 //      }
 //      // Channel 6 and 8
 //      if (ChanHitIndex[6]>=0 && ChanHitIndex[8]>=0) {
	// h_LayerNeighborDeltaTime_6->Fill(float(time_module_calibrated->at(ChanHitIndex[6])-time_module_calibrated->at(ChanHitIndex[8])),1.);
 //      }
 //      // Channel 7 and 9
 //      if (ChanHitIndex[7]>=0 && ChanHitIndex[9]>=0) {
	// h_LayerNeighborDeltaTime_7->Fill(float(time_module_calibrated->at(ChanHitIndex[7])-time_module_calibrated->at(ChanHitIndex[9])),1.);
 //      }
 //      // Channel 4 and 10
 //      if (ChanHitIndex[4]>=0 && ChanHitIndex[10]>=0) {
	// h_LayerNeighborDeltaTime_4->Fill(float(time_module_calibrated->at(ChanHitIndex[4])-time_module_calibrated->at(ChanHitIndex[10])),1.);
 //      }
 //      // Channel 5 and 11
 //      if (ChanHitIndex[5]>=0 && ChanHitIndex[11]>=0) {
	// h_LayerNeighborDeltaTime_5->Fill(float(time_module_calibrated->at(ChanHitIndex[5])-time_module_calibrated->at(ChanHitIndex[11])),1.);
 //      }
 //      // Channel 12 and 14
 //      if (ChanHitIndex[12]>=0 && ChanHitIndex[14]>=0) {
	// h_LayerNeighborDeltaTime_12->Fill(float(time_module_calibrated->at(ChanHitIndex[12])-time_module_calibrated->at(ChanHitIndex[14])),1.);
 //      }

     // Plots of coincidence rates with other channels
     for (int c = 0; c<32; c++)
      for (int c2 = 0; c2<32; c2++) 
       if ( (c2 != c) && HitChan[c] && HitChan[c2] )
        h_NPairedHits[c]->Fill(float(c2),1.);
      
     // Plots of triplet hits
     //  1 = CH0-6-2
     //  2 = CH1-7-3
     //  3 = CH24-16-22
     //  4 = CH25-17-23
     //  5 = CH8-12-4
     //  6 = CH9-13-5
     //  7 = CH10-30-14
     //  8 = CH29-19-26
     //  9 = CH27-11-31
     // 10 = CH18-20-28
     // 11 = CH18-21-28
     if (HitChan[0] && HitChan[6] && HitChan[2]) h_NTripletHits->Fill(1.,1.);
     if (HitChan[1] && HitChan[7] && HitChan[3]) h_NTripletHits->Fill(2.,1.);
     if (HitChan[24] && HitChan[16] && HitChan[22]) h_NTripletHits->Fill(3.,1.);
     if (HitChan[25] && HitChan[17] && HitChan[23]) h_NTripletHits->Fill(4.,1.);
     if (HitChan[8] && HitChan[12] && HitChan[4]) h_NTripletHits->Fill(5.,1.);
     if (HitChan[9] && HitChan[13] && HitChan[5]) h_NTripletHits->Fill(6.,1.);
     if (HitChan[10] && HitChan[30] && HitChan[14]) h_NTripletHits->Fill(7.,1.);
     if (HitChan[29] && HitChan[19] && HitChan[26]) h_NTripletHits->Fill(8.,1.);
     if (HitChan[27] && HitChan[11] && HitChan[31]) h_NTripletHits->Fill(9.,1.);
     if (HitChan[18] && HitChan[20] && HitChan[28]) h_NTripletHits->Fill(10.,1.);
     if (HitChan[18] && HitChan[21] && HitChan[28]) h_NTripletHits->Fill(11.,1.);
     // Triplet hits with no other hits present
     if (nHits == 3) { // 3 hits?
      if (HitChan[0] && HitChan[6] && HitChan[2]) h_NTripletHits3->Fill(1.,1.);
      if (HitChan[1] && HitChan[7] && HitChan[3]) h_NTripletHits3->Fill(2.,1.);
      if (HitChan[24] && HitChan[16] && HitChan[22]) h_NTripletHits3->Fill(3.,1.);
      if (HitChan[25] && HitChan[17] && HitChan[23]) h_NTripletHits3->Fill(4.,1.);
      if (HitChan[8] && HitChan[12] && HitChan[4]) h_NTripletHits3->Fill(5.,1.);
      if (HitChan[9] && HitChan[13] && HitChan[5]) h_NTripletHits3->Fill(6.,1.);
      if (HitChan[10] && HitChan[30] && HitChan[14]) h_NTripletHits3->Fill(7.,1.);
      if (HitChan[29] && HitChan[19] && HitChan[26]) h_NTripletHits3->Fill(8.,1.);
      if (HitChan[27] && HitChan[11] && HitChan[31]) h_NTripletHits3->Fill(9.,1.);
      if (HitChan[18] && HitChan[20] && HitChan[28]) h_NTripletHits3->Fill(10.,1.);
      if (HitChan[18] && HitChan[21] && HitChan[28]) h_NTripletHits3->Fill(11.,1.);
     } // 3 hits?
     if (nHits == 4) { // 4 hits?
      if (HitChan[0] && HitChan[6] && HitChan[2]) h_NTripletHits4->Fill(1.,1.);
      if (HitChan[1] && HitChan[7] && HitChan[3]) h_NTripletHits4->Fill(2.,1.);
      if (HitChan[24] && HitChan[16] && HitChan[22]) h_NTripletHits4->Fill(3.,1.);
      if (HitChan[25] && HitChan[17] && HitChan[23]) h_NTripletHits4->Fill(4.,1.);
      if (HitChan[8] && HitChan[12] && HitChan[4]) h_NTripletHits4->Fill(5.,1.);
      if (HitChan[9] && HitChan[13] && HitChan[5]) h_NTripletHits4->Fill(6.,1.);
      if (HitChan[10] && HitChan[30] && HitChan[14]) h_NTripletHits4->Fill(7.,1.);
      if (HitChan[29] && HitChan[19] && HitChan[26]) h_NTripletHits4->Fill(8.,1.);
      if (HitChan[27] && HitChan[11] && HitChan[31]) h_NTripletHits4->Fill(9.,1.);
      if (HitChan[18] && HitChan[20] && HitChan[28]) h_NTripletHits4->Fill(10.,1.);
      if (HitChan[18] && HitChan[21] && HitChan[28]) h_NTripletHits4->Fill(11.,1.);
     } // 4 hits?

     // Fill counters of hit pairs within and between layers, which are used to determine rates and validate the trigger.
     if (HitChan[0]&&HitChan[8]) nCH0CH8++;
     if (HitChan[1]&&HitChan[9]) nCH1CH9++; 
     if (HitChan[0]&&HitChan[9]) nCH0CH9++;
     if (HitChan[1]&&HitChan[8]) nCH1CH8++;
     if (HitChan[6]&&HitChan[12]) nCH6CH12++; 
     if (HitChan[7]&&HitChan[13]) nCH7CH13++; 
     if (HitChan[6]&&HitChan[13]) nCH6CH13++; 
     if (HitChan[7]&&HitChan[12]) nCH7CH12++;
     if (HitChan[2]&&HitChan[4]) nCH2CH4++; 
     if (HitChan[3]&&HitChan[5]) nCH3CH5++; 
     if (HitChan[2]&&HitChan[5]) nCH2CH5++; 
     if (HitChan[3]&&HitChan[4]) nCH3CH4++; 
     if (HitChan[0]&&HitChan[6]) nCH0CH6++;
     if (HitChan[1]&&HitChan[7]) nCH1CH7++;
     if (HitChan[0]&&HitChan[7]) nCH0CH7++;
     if (HitChan[1]&&HitChan[6]) nCH1CH6++;
     if (HitChan[8]&&HitChan[12]) nCH8CH12++;
     if (HitChan[9]&&HitChan[13]) nCH9CH13++;
     if (HitChan[8]&&HitChan[13]) nCH8CH13++;
     if (HitChan[9]&&HitChan[12]) nCH9CH12++;
     if (HitChan[0]&&HitChan[12]) nCH0CH12++;
     if (HitChan[1]&&HitChan[13]) nCH1CH13++;
     if (HitChan[0]&&HitChan[13]) nCH0CH13++;
     if (HitChan[1]&&HitChan[12]) nCH1CH12++;
     if (HitChan[8]&&HitChan[6]) nCH8CH6++;
     if (HitChan[9]&&HitChan[7]) nCH9CH7++;
     if (HitChan[8]&&HitChan[7]) nCH8CH7++;
     if (HitChan[9]&&HitChan[6]) nCH9CH6++;
     if (HitChan[6]&&HitChan[2]) nCH6CH2++;
     if (HitChan[7]&&HitChan[3]) nCH7CH3++;
     if (HitChan[6]&&HitChan[3]) nCH6CH3++;
     if (HitChan[7]&&HitChan[2]) nCH7CH2++;
     if (HitChan[12]&&HitChan[4]) nCH12CH4++;
     if (HitChan[13]&&HitChan[5]) nCH13CH5++;
     if (HitChan[12]&&HitChan[5]) nCH12CH5++;
     if (HitChan[13]&&HitChan[4]) nCH13CH4++;
     if (HitChan[6]&&HitChan[4]) nCH6CH4++;
     if (HitChan[7]&&HitChan[5]) nCH7CH5++;
     if (HitChan[6]&&HitChan[5]) nCH6CH5++;
     if (HitChan[7]&&HitChan[4]) nCH7CH4++;
     if (HitChan[12]&&HitChan[2]) nCH12CH2++;
     if (HitChan[13]&&HitChan[3]) nCH13CH3++;
     if (HitChan[12]&&HitChan[3]) nCH12CH3++;
     if (HitChan[13]&&HitChan[2]) nCH13CH2++;
     if (HitChan[0]&&HitChan[8]&&HitChan[10]) nCH0CH8CH10++; // 3 hits per layer
     if (HitChan[1]&&HitChan[9]&&HitChan[10]) nCH1CH9CH10++;
     if (HitChan[0]&&HitChan[9]&&HitChan[10]) nCH0CH9CH10++;
     if (HitChan[1]&&HitChan[8]&&HitChan[10]) nCH1CH8CH10++;
     if (HitChan[6]&&HitChan[12]&&HitChan[30]) nCH6CH12CH30++;
     if (HitChan[7]&&HitChan[13]&&HitChan[30]) nCH7CH13CH30++;
     if (HitChan[2]&&HitChan[4]&&HitChan[14]) nCH2CH4CH14++;
     if (HitChan[3]&&HitChan[5]&&HitChan[14]) nCH3CH5CH14++;
     if (HitChan[2]&&HitChan[5]&&HitChan[14]) nCH2CH5CH14++;
     if (HitChan[3]&&HitChan[4]&&HitChan[14]) nCH3CH4CH14++;
     if (HitChan[0]&&HitChan[6]&&HitChan[2]) nCH0CH6CH2++; // Hits across 3 layers
     if (HitChan[1]&&HitChan[7]&&HitChan[3]) nCH1CH7CH3++;
     if (HitChan[0]&&HitChan[6]&&HitChan[3]) nCH0CH6CH3++;
     if (HitChan[1]&&HitChan[7]&&HitChan[2]) nCH1CH7CH2++;
     if (HitChan[0]&&HitChan[7]&&HitChan[2]) nCH0CH7CH2++;
     if (HitChan[1]&&HitChan[6]&&HitChan[3]) nCH1CH6CH3++;
     if (HitChan[8]&&HitChan[12]&&HitChan[4]) nCH8CH12CH4++;
     if (HitChan[9]&&HitChan[13]&&HitChan[5]) nCH9CH13CH5++;
     if (HitChan[8]&&HitChan[12]&&HitChan[5]) nCH8CH12CH5++;
     if (HitChan[9]&&HitChan[13]&&HitChan[4]) nCH9CH13CH4++;
     if (HitChan[8]&&HitChan[13]&&HitChan[4]) nCH8CH13CH4++;
     if (HitChan[9]&&HitChan[12]&&HitChan[5]) nCH9CH12CH5++;
     if (HitChan[0]&&HitChan[6]&&HitChan[4]) nCH0CH6CH4++;
     if (HitChan[0]&&HitChan[6]&&HitChan[5]) nCH0CH6CH5++;
     if (HitChan[1]&&HitChan[7]&&HitChan[5]) nCH1CH7CH5++;
     if (HitChan[1]&&HitChan[7]&&HitChan[4]) nCH1CH7CH4++;
     if (HitChan[11]&&HitChan[8]) nCH11CH8++; // Two hits in digitizer 1...
     if (HitChan[11]&&HitChan[9]) nCH11CH9++; 
     if (HitChan[13]&&HitChan[8]) nCH13CH8++; 
     if (HitChan[13]&&HitChan[9]) nCH13CH9++;
     if (HitChan[16]&&HitChan[30]) nCH16CH30++;
     if (HitChan[17]&&HitChan[30]) nCH17CH30++;
     if (HitChan[16]&&HitChan[19]) nCH16CH19++;
     if (HitChan[17]&&HitChan[19]) nCH17CH19++;
     if (HitChan[22]&&HitChan[15]) nCH22CH15++;
     if (HitChan[23]&&HitChan[15]) nCH23CH15++;
     if (HitChan[22]&&HitChan[10]) nCH22CH10++;
     if (HitChan[23]&&HitChan[10]) nCH23CH10++;
     if (HitChan[27]&&HitChan[24]&&HitChan[29]) nCH27CH24CH29++; // Three hits within a layer
     if (HitChan[27]&&HitChan[25]&&HitChan[29]) nCH27CH25CH29++;
     if (HitChan[30]&&HitChan[16]&&HitChan[19]) nCH30CH16CH19++;
     if (HitChan[30]&&HitChan[17]&&HitChan[19]) nCH30CH17CH19++;
     if (HitChan[31]&&HitChan[22]&&HitChan[26]) nCH31CH22CH26++;
     if (HitChan[31]&&HitChan[23]&&HitChan[26]) nCH31CH23CH26++;
     if (HitChan[24]&&HitChan[16]&&HitChan[22]) nCH24CH16CH22++; // Hits across 3 layers
     if (HitChan[24]&&HitChan[16]&&HitChan[23]) nCH24CH16CH23++;
     if (HitChan[24]&&HitChan[17]&&HitChan[23]) nCH24CH17CH23++;
     if (HitChan[24]&&HitChan[17]&&HitChan[22]) nCH24CH17CH22++;
     if (HitChan[25]&&HitChan[17]&&HitChan[23]) nCH25CH17CH23++;
     if (HitChan[25]&&HitChan[16]&&HitChan[23]) nCH25CH16CH23++;
     if (HitChan[25]&&HitChan[16]&&HitChan[22]) nCH25CH16CH22++;
     if (HitChan[18]&&HitChan[20]) nCH18CH20++; // Two hits between slabs
     if (HitChan[20]&&HitChan[28]) nCH20CH28++;
     if (HitChan[28]&&HitChan[21]) nCH28CH21++;
     if (HitChan[18]&&HitChan[28]) nCH18CH28++;
     if (HitChan[20]&&HitChan[21]) nCH20CH21++;
     if (HitChan[18]&&HitChan[21]) nCH18CH21++;
     if (HitChan[18]&&HitChan[20]&&HitChan[28]) nCH18CH20CH28++; // Hits across slab layers
     if (HitChan[20]&&HitChan[28]&&HitChan[21]) nCH20CH28CH21++;
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]) nCH18CH20CH28CH21++; 
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      HitChan[0]&&HitChan[6]&&HitChan[2]) nCH18CH20CH28CH21CH0CH6CH2++; 
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      HitChan[1]&&HitChan[7]&&HitChan[3]) nCH18CH20CH28CH21CH1CH7CH3++; 
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      HitChan[24]&&HitChan[16]&&HitChan[22]) nCH18CH20CH28CH21CH24CH16CH22++; 
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      HitChan[25]&&HitChan[17]&&HitChan[23]) nCH18CH20CH28CH21CH25CH17CH23++; 
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      HitChan[12]&&HitChan[4]) nCH18CH20CH28CH21CH12CH4++; 
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      HitChan[9]&&HitChan[13]&&HitChan[5]) nCH18CH20CH28CH21CH9CH13CH5++; 
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      (HitChan[0]||HitChan[1]||HitChan[24]||HitChan[25]||HitChan[8]||HitChan[9])&&
      (HitChan[6]||HitChan[7]||HitChan[16]||HitChan[17]||HitChan[12]||HitChan[13])&&
      (HitChan[2]||HitChan[3]||HitChan[22]||HitChan[23]||HitChan[4]||HitChan[5]) ) {
      nAllSlabsAllLayers++;
      if (evtCategoryCode==4) 
      	cout << "AllSlabsAllLayers event: makeDisplays.py "<<RunNum<<" 1 -s "<<'"'<<"file=="<<filenum<<"&&event=="<<event<<'"'<<" -r 0 640 --noBounds --tag AllSlabsAllLayers\n";
     }
     if (HitChan[18]&&HitChan[20]&&HitChan[28]&&HitChan[21]&&
      (HitChan[0]||HitChan[1]||HitChan[24]||HitChan[25]||HitChan[8]||HitChan[9]||
       HitChan[6]||HitChan[7]||HitChan[16]||HitChan[17]||HitChan[12]||HitChan[13]||
       HitChan[2]||HitChan[3]||HitChan[22]||HitChan[23]||HitChan[4]||HitChan[5]) )
      nAllSlabsAnyLayer++;
     if ((HitChan[18]||HitChan[20]||HitChan[28]||HitChan[21])&&
      (HitChan[0]||HitChan[1]||HitChan[24]||HitChan[25]||HitChan[8]||HitChan[9])&&
      (HitChan[6]||HitChan[7]||HitChan[16]||HitChan[17]||HitChan[12]||HitChan[13])&&
      (HitChan[2]||HitChan[3]||HitChan[22]||HitChan[23]||HitChan[4]||HitChan[5]) )
      nAllLayersAnySlab++;

     // Dump an ANT for mergine with the hodoscope data
     if (evtCategoryCode==0 && RangeCode == "A") { // TS1 ANT
       cout <<setprecision(4);
       cout << "\nTS1:";
       cout << run<<" "; // 1
       cout << filenum<<" "; // 2
       cout << event<<" "; // 3
       cout << ia<<" "; // 4
       cout <<beam<<" "; // 5
//       cout << long(10000*(event_time_fromTDC-1529000000.))<<" "; // 6; Event time in 100 microsecond steps since 152900000
       cout << long(100*(event_time_fromTDC))<<" "; // 6; Event time in 100 microsecond steps since 152900000
       cout << timeBetweenEvents<<" "; // 7; Time between events in microseconds
       cout << event_trigger_time_tag_b0 << " "; // 8; Time for digitizer 0
       cout << event_trigger_time_tag_b1 << " "; // 9; Time for digitizer 1
       cout << event_trigger_time_tag_b0-prev_event_trigger_time_tag_b0 << " "; // 10; Time between events for dig0
       cout << event_trigger_time_tag_b1-prev_event_trigger_time_tag_b1 << " "; // 11; Time between events for dig0
       prev_event_trigger_time_tag_b0 = event_trigger_time_tag_b0;
       prev_event_trigger_time_tag_b1 = event_trigger_time_tag_b1;
       cout << event_trigger_time_tag_b1-event_trigger_time_tag_b0<<" "; // 12; TTT mismatch between boards

       cout << groupTDC_b0->at(0) << " "; // 13; gTDC for digitizer 0
       cout << groupTDC_b1->at(0) << " "; // 14; gTDC for digitizer 1
       cout << groupTDC_b0->at(0)-prev_gTDC_b0 << " "; // 15; gTDC diff from prev event for dig 0
       cout << groupTDC_b1->at(0)-prev_gTDC_b1 << " "; // 16; gTDC diff from prev event for dig 0
       prev_gTDC_b0 = groupTDC_b0->at(0);
       prev_gTDC_b1 = groupTDC_b1->at(0);
       cout << groupTDC_b1->at(0)-groupTDC_b0->at(0) << " "; // 17; gTDC mismatch between boards

       for (int c=0; c<32; c++) {
         cout << " "<<maxSample[c];
         cout << " "<<TotalArea[c];
       }
       cout << "\n";
     } // TS1 ANT

    } // Good event?
    
  } //end over loop over all the events

  outFile->cd();
  outFile->Write();
  outFile->Close();

  if (RangeCode == "A") { // Report the rates for each combination of hits
   int nSeconds = int(0.5+runTime);
   int nX;
   nX = nCH0CH8; // Within a layer
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH8 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH9; 
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH9 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH9;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH9 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH8;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH8 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH6CH12; 
   cout << "CombinedHits: "<<EventCategory<<" CH6+CH12 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH7CH13; 
   cout << "CombinedHits: "<<EventCategory<<" CH7+CH13 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH6CH13; 
   cout << "CombinedHits: "<<EventCategory<<" CH6+CH13 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH7CH12;
   cout << "CombinedHits: "<<EventCategory<<" CH7+CH12 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH2CH4; 
   cout << "CombinedHits: "<<EventCategory<<" CH2+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH3CH5; 
   cout << "CombinedHits: "<<EventCategory<<" CH3+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH2CH5; 
   cout << "CombinedHits: "<<EventCategory<<" CH2+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH3CH4; 
   cout << "CombinedHits: "<<EventCategory<<" CH3+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH6;   // Between layers
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH6 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH7;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH7 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH7;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH7 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH6;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH6 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH8CH12;
   cout << "CombinedHits: "<<EventCategory<<" CH8+CH12 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH9CH13;
   cout << "CombinedHits: "<<EventCategory<<" CH9+CH13 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH8CH13;
   cout << "CombinedHits: "<<EventCategory<<" CH8+CH13 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH9CH12;
   cout << "CombinedHits: "<<EventCategory<<" CH9+CH12 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH12;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH12 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH13;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH13 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH13;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH13 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH12;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH12 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH8CH6;
   cout << "CombinedHits: "<<EventCategory<<" CH8+CH6 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH9CH7;
   cout << "CombinedHits: "<<EventCategory<<" CH9+CH7 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH8CH7;
   cout << "CombinedHits: "<<EventCategory<<" CH8+CH7 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH9CH6;
   cout << "CombinedHits: "<<EventCategory<<" CH9+CH6 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH6CH2;
   cout << "CombinedHits: "<<EventCategory<<" CH6+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH7CH3;
   cout << "CombinedHits: "<<EventCategory<<" CH7+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH6CH3;
   cout << "CombinedHits: "<<EventCategory<<" CH6+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH7CH2;
   cout << "CombinedHits: "<<EventCategory<<" CH7+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH12CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH12+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH13CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH13+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH12CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH12+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH13CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH13+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH6CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH6+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH7CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH7+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH6CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH6+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH7CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH7+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH12CH2;
   cout << "CombinedHits: "<<EventCategory<<" CH12+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH13CH3;
   cout << "CombinedHits: "<<EventCategory<<" CH13+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH12CH3;
   cout << "CombinedHits: "<<EventCategory<<" CH12+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH13CH2;
   cout << "CombinedHits: "<<EventCategory<<" CH13+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH8CH10; // 3 hits per layer
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH8+CH10 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH9CH10;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH9+CH10 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH9CH10;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH9+CH10 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH8CH10;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH8+CH10 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH6CH12CH30;
   cout << "CombinedHits: "<<EventCategory<<" CH6+CH12+CH30 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH7CH13CH30;
   cout << "CombinedHits: "<<EventCategory<<" CH7+CH13+CH30 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH2CH4CH14;
   cout << "CombinedHits: "<<EventCategory<<" CH2+CH4+CH14 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH3CH5CH14;
   cout << "CombinedHits: "<<EventCategory<<" CH3+CH5+CH14 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH2CH5CH14;
   cout << "CombinedHits: "<<EventCategory<<" CH2+CH5+CH14 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH3CH4CH14;
   cout << "CombinedHits: "<<EventCategory<<" CH3+CH4+CH14 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH6CH2; // Hits across 3 layers
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH6+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH7CH3;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH7+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH6CH3;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH6+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH7CH2;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH7+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH7CH2;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH7+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH6CH3;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH6+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH8CH12CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH8+CH12+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH9CH13CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH9+CH13+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH8CH12CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH8+CH12+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH9CH13CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH9+CH13+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH8CH13CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH8+CH13+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH9CH12CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH9+CH12+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH6CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH6+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH0CH6CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH0+CH6+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH7CH5;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH7+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH1CH7CH4;
   cout << "CombinedHits: "<<EventCategory<<" CH1+CH7+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH11CH8;
   cout << "CombinedHits: "<<EventCategory<<" CH11+CH8 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH11CH9; 
   cout << "CombinedHits: "<<EventCategory<<" CH11+CH9 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH13CH8; 
   cout << "CombinedHits: "<<EventCategory<<" CH13+CH8 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH13CH9;
   cout << "CombinedHits: "<<EventCategory<<" CH13+CH9 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH16CH30;
   cout << "CombinedHits: "<<EventCategory<<" CH16+CH30 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH17CH30;
   cout << "CombinedHits: "<<EventCategory<<" CH17+CH30 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH16CH19;
   cout << "CombinedHits: "<<EventCategory<<" CH16+CH19 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH17CH19;
   cout << "CombinedHits: "<<EventCategory<<" CH17+CH19 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH22CH15;
   cout << "CombinedHits: "<<EventCategory<<" CH22+CH15 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH23CH15;
   cout << "CombinedHits: "<<EventCategory<<" CH23+CH15 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH22CH10;
   cout << "CombinedHits: "<<EventCategory<<" CH22+CH10 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH23CH10;
   cout << "CombinedHits: "<<EventCategory<<" CH23+CH10 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH27CH24CH29; // Three hits within a layer
   cout << "CombinedHits: "<<EventCategory<<" CH27+CH24+CH29 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH27CH25CH29;
   cout << "CombinedHits: "<<EventCategory<<" CH27+CH25+CH29 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH30CH16CH19;
   cout << "CombinedHits: "<<EventCategory<<" CH30+CH16+CH19 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH30CH17CH19;
   cout << "CombinedHits: "<<EventCategory<<" CH30+CH17+CH19 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH31CH22CH26;
   cout << "CombinedHits: "<<EventCategory<<" CH31+CH22+CH26 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH31CH23CH26;
   cout << "CombinedHits: "<<EventCategory<<" CH31+CH23+CH26 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH24CH16CH22; // Hits across 3 layers
   cout << "CombinedHits: "<<EventCategory<<" CH24+CH16+CH22 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH24CH16CH23;
   cout << "CombinedHits: "<<EventCategory<<" CH24+CH16+CH23 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH24CH17CH23;
   cout << "CombinedHits: "<<EventCategory<<" CH24+CH17+CH23 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH24CH17CH22;
   cout << "CombinedHits: "<<EventCategory<<" CH24+CH17+CH22 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH25CH17CH23;
   cout << "CombinedHits: "<<EventCategory<<" CH25+CH17+CH23 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH25CH16CH23;
   cout << "CombinedHits: "<<EventCategory<<" CH25+CH16+CH23 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH25CH16CH22;
   cout << "CombinedHits: "<<EventCategory<<" CH25+CH16+CH22 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20; // Two hits between slabs
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH20CH28;
   cout << "CombinedHits: "<<EventCategory<<" CH20+CH28 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH28CH21;
   cout << "CombinedHits: "<<EventCategory<<" CH28+CH21 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH28;
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH28 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH20CH21;
   cout << "CombinedHits: "<<EventCategory<<" CH20+CH21 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH21;
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH21 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20CH28; // Hits across slab layers
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH20CH28CH21;
   cout << "CombinedHits: "<<EventCategory<<" CH20+CH28+CH21 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20CH28CH21; 
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28+CH21 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour. = "<<nX/totalLumi<<" +- "<<sqrt(nX)/totalLumi<<" per ipb.\n";
   nX = nCH18CH20CH28CH21CH0CH6CH2; 
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28+CH21+CH0+CH6+CH2 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20CH28CH21CH1CH7CH3; 
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28+CH21+CH1+CH7+CH3 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20CH28CH21CH24CH16CH22; 
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28+CH21+CH24+CH16+CH22 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20CH28CH21CH25CH17CH23; 
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28+CH21+CH25+CH17+CH23 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20CH28CH21CH12CH4; 
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28+CH21+CH12+CH4 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nCH18CH20CH28CH21CH9CH13CH5; 
   cout << "CombinedHits: "<<EventCategory<<" CH18+CH20+CH28+CH21+CH9+CH13+CH5 "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour.\n";
   nX = nAllSlabsAllLayers; 
   cout << "CombinedHits: "<<EventCategory<<" AllSlabsAllLayers "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour. = "<<nX/totalLumi<<" +- "<<sqrt(nX)/totalLumi<<" per ipb.\n";
   nX = nAllSlabsAnyLayer; 
   cout << "CombinedHits: "<<EventCategory<<" AllSlabsAnyLayer "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour. = "<<nX/totalLumi<<" +- "<<sqrt(nX)/totalLumi<<" per ipb.\n";
   nX = nAllLayersAnySlab; 
   cout << "CombinedHits: "<<EventCategory<<" AllLayersAnySlab "<<nX<<"/"<<nSeconds<<" = "<<nX*3600./nSeconds<<" +- "<<sqrt(nX)*3600./nSeconds<<" per hour. = "<<nX/totalLumi<<" +- "<<sqrt(nX)/totalLumi<<" per ipb.\n";
  }
}


