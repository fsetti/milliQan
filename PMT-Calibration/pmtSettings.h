// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include <iostream>
#include <vector>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"
#include "TLegend.h"
#include "TSystemDirectory.h"

using namespace std;


// stuff for binning between 1 and 1e3 nPEs on a log scale
const unsigned int minBin = 0;
const unsigned int maxBin = 3;
const unsigned int nBins = 1000;
const float increment = float(maxBin-minBin)/nBins;

const TString pmtTemplate = "template_1p6GHz";

const float spe_fit_mean_ET = 25.82;
const float spe_fit_mean_R7725 = 151.49;
const float spe_fit_mean_R878 = 69.15;

const float threshold_ET = 3.;
const float threshold_R7725 = 3.;
const float threshold_R878 = 3.;

void timeDelay(TString path1, TString pmtType, float threshold, TString outputFileName );
