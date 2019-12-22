// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TSystemDirectory.h"


using namespace std;

const float tFire = 380;  // lowest limit of firing pulse window
const float dtMax = 100;    // max time difference btw pulses
const float dtMin = 15;     // max time difference btw calibrated times
const float tMin = tFire - dtMax;   //smallest time of 100ns window
const float tMax = tFire + dtMax;   //largest time of 100ns window
const float hMin = 100;
const float nPEminSlab = 150;
const float nPEminPanel = 15;
const float nPEminMuon = 750;
const float nPEminMuonR7725 = 350;
const float nPEmin = 1.;
const vector<float> nPEminVec = {nPEminMuon,nPEminMuon, nPEminMuon,nPEminMuon, nPEminMuon,nPEminMuonR7725, nPEminMuon,nPEminMuon, nPEminMuon,nPEminMuon, nPEminPanel,nPEminPanel, nPEminMuon,nPEminMuon, nPEminPanel, 0, nPEminMuon,nPEminMuon, nPEminSlab, nPEminPanel, nPEminSlab,nPEminSlab, nPEminMuonR7725,nPEminMuon, nPEminMuon,nPEminMuon, nPEminPanel, nPEminPanel, nPEminSlab, nPEminPanel, nPEminPanel, nPEminPanel};


// plotting options
const float minBin = -1;
const float maxBin = 3.3;
const int nBins = 10;

const float maxNPEratio = 10;



// trigger configurations for double-hit within same layer
// first layer
const vector<int> trig100248 = {10,0,24,8};
const vector<int> trig101259 = {10,1,25,9};
// second layer
const vector<int> trig1161612 = {11,6,16,12};
const vector<int> trig1171713 = {11,7,17,13};
// third layer
const vector<int> trig142224 = {14,2,22,4};
const vector<int> trig143235 = {14,3,23,5};

// layer 1
const vector<int> trigL1 = {0,1,24,25,8,9};
const vector<int> trigL2 = {6,7,16,17,12,13};
const vector<int> trigL3 = {2,3,22,23,4,5};

// panels
const vector<int> trigPanL1 = {27,29};
const vector<int> trigPanL2 = {19,30};
const vector<int> trigPanL3 = {26,31};


const vector<int> trigL1left = {0,24,8};
const vector<int> trigL1right = {1,25,9};
const vector<int> trigL2left = {6,16,12};
const vector<int> trigL2right = {7,17,13};
const vector<int> trigL3left = {2,22,4};
const vector<int> trigL3right = {3,23,5};


const vector<int>  trigL2L3 = {6,16,12,7,17,13,2,22,4,3,23,5};
const vector<int>  trigL1L3 = {0,24,8,1,25,9,2,22,4,3,23,5};
const vector<int>  trigL1L2 = {0,24,8,1,25,9,6,16,12,7,17,13};

const vector<int> trigBarsAndPanels = {0,1,24,25,8,9,6,7,16,17,12,13,2,3,22,23,4,5, 10,11,14,27,30,31,29,19,26};
const vector<int> trigBarsPanelsSlabs = {0,1,24,25,8,9,6,7,16,17,12,13,2,3,22,23,4,5, 10,11,14,27,30,31,29,19,26, 18,20,21,28};
const vector<int> trigBars = {0,1,24,25,8,9,6,7,16,17,12,13,2,3,22,23,4,5};
const vector<int> trigSlabs = {18,20,21,28};



