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
const float nPEminBar = 750;
const float nPEminBarR7725 = 350;
const float nPEminBarS = 0.1;
const float maxSPE = 5;


// trigger configurations for double-hit within same layer
// first layer
const vector<int> trig01 = {0,1};
const vector<int> trig2425 = {24,25};
const vector<int> trig89 = {8,9};

const vector<int> trig024 = {0,24};
const vector<int> trig124 = {1,24};
const vector<int> trig025 = {0,25};
const vector<int> trig125 = {1,25};

const vector<int> trig08 = {0,8};
const vector<int> trig18 = {1,8};
const vector<int> trig09 = {0,9};
const vector<int> trig19 = {1,9};

const vector<int> trig248 = {24,8};
const vector<int> trig258 = {25,8};
const vector<int> trig249 = {24,9};
const vector<int> trig259 = {25,9};


// second layer
const vector<int> trig616 = {6,16};
const vector<int> trig617 = {6,17};
const vector<int> trig716 = {7,16};
const vector<int> trig717 = {7,17};

const vector<int> trig1612 = {16,12};
const vector<int> trig1613 = {16,13};
const vector<int> trig1712 = {17,12};
const vector<int> trig1713 = {17,13};

const vector<int> trig612 = {6,12};
const vector<int> trig613 = {6,13};
const vector<int> trig712 = {7,12};
const vector<int> trig713 = {7,13};


// third layer
const vector<int> trig222 = {2,22};
const vector<int> trig223 = {2,23};
const vector<int> trig322 = {3,22};
const vector<int> trig323 = {3,23};

const vector<int> trig224 = {22,4};
const vector<int> trig225 = {22,5};
const vector<int> trig234 = {23,4};
const vector<int> trig235 = {23,5};

const vector<int> trig24 = {2,4};
const vector<int> trig25 = {2,5};
const vector<int> trig34 = {3,4};
const vector<int> trig35 = {3,5};

// trigger configurations for through going particles
const vector<int> trig062 = {0,2,6};
const vector<int> trig063 = {0,3,6};
const vector<int> trig162 = {1,2,6};
const vector<int> trig173 = {1,3,7};
const vector<int> trig172 = {1,2,7};
const vector<int> trig073 = {0,3,7};
const vector<int> trig241622 = {16,22,24};
const vector<int> trig241623 = {16,23,24};
const vector<int> trig251622 = {16,22,25};
const vector<int> trig251723 = {17,23,25};
const vector<int> trig251722 = {17,22,25};
const vector<int> trig241723 = {17,23,24};
const vector<int> trig8124 = {8,12,4};
const vector<int> trig8125 = {8,12,5};
const vector<int> trig9124 = {9,12,4};
const vector<int> trig9135 = {9,13,5};
const vector<int> trig9134 = {9,13,4};
const vector<int> trig8135 = {8,13,5};

 
bool sort_wrt_channel( vector<float> vec1, vector<float> vec2 );

bool sort_wrt_height( vector<float> vec1, vector<float> vec2 );

bool sort_wrt_time( vector<float> vec1, vector<float> vec2 );

void read_list_files(const char *dirname, TChain *chainptr);

vector<vector<float>> sameLayerSelectionNoHit( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<int> TrigChannels );

vector<vector<float>> sameLayerSelectionNoHit_v2( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<int> TrigChannels );

vector<vector<vector<float>>> sameLayerSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels );

vector<vector<vector<float>>> pulseBarSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels );

vector<vector<vector<float>>> pulseSlabSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels );
