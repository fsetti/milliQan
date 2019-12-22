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
const float nPEmin = 150;
// trigger configurations for through going particles
const vector<int> trig062 = {0,2,6};
const vector<int> trig173 = {1,3,7};
const vector<int> trig241622 = {16,22,24};
const vector<int> trig251723 = {17,23,25};
const vector<int> trig8124 = {8,12,4};
const vector<int> trig9135 = {9,13,5};

 
bool sort_wrt_channel( vector<float> vec1, vector<float> vec2 );

bool sort_wrt_height( vector<float> vec1, vector<float> vec2 );

bool sort_wrt_time( vector<float> vec1, vector<float> vec2 );

void read_list_files(const char *dirname, TChain *chainptr);

vector<vector<float>> pulsePreSelection( vector<vector<float>> Pulses, vector<unsigned int> TrigChannels );

vector<vector<float>> pulseSelection( vector<vector<float>> Pulses, unsigned int nTriggerChannels );
