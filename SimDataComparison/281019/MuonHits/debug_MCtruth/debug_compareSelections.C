// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/MCntupleHandler.h"


void selection( TChain *chain ) 
{
  double weight = 1;
  double weight_qcd = 0.014301851 * totLumi;
  double weight_qcd_nonbc = 0.0033202620 * totLumi;
  double weight_DY = 0.00011057182 * totLumi;
  double weight_W = 0.00026728950 * totLumi;
  TFile *currentFile;
  TString currentFileName;

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


//  ofstream outfile;
//  outfile.open("MCcomparison.txt");
   ifstream infile;
   infile.open("outEvtsBenNoFran.txt");
   vector<int> fileIDs, orig_events, channels;
   int fID, oEvts, chanN;
   while ( infile >> fID >> oEvts >> chanN ){
	fileIDs.push_back(fID);
	orig_events.push_back(oEvts);
	channels.push_back(chanN);
   }



  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);


    bool passID = false;
    int currentChannel = -1;
    for ( unsigned int i=0; i<fileIDs.size(); i++){
	if ( fileID == fileIDs.at(i) && orig_evt == orig_events.at(i) ){
		passID = true;
		currentChannel = channels.at(i);
	}
    }
    if ( !passID ) continue;

	for (unsigned int i=0; i<chan->size(); i++){
		if ( chan->at(i) == currentChannel ){
			cout << "reco nPE: " << nPE->at(i) << ", truth MC: " << chan_trueNPE[chan->at(i)] << "   fileID: " << fileID << " event " << orig_evt << ", chan: " << chan->at(i) << endl;
		}
	}

    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > maxSlabdT ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

   
/////////////////////////////////////////////////////////////////////////////
//	Muon Selection - hits in straight line 

    vector<vector<float>> nonNeighPulses, neighPulses;
    nonNeighPulses = NonNeighBarSelection ( chan,height,ptime,nPE,time_module_calibrated, trigBars);
    neighPulses = NeighBarSelection ( chan,height,ptime,nPE,time_module_calibrated, trigBars);

	for ( unsigned int i=0; i<nonNeighPulses.size(); i++){
	}	
	
	//	repeat for neighbouring channels
	for ( unsigned int i=0; i<neighPulses.size(); i++){
//		if ( (int)neighPulses.at(i).at(0) == currentChannel ) cout << "reco nPE: " << neighPulses.at(i).at(3) << ", MCtruth: " << chan_trueNPE[(int)neighPulses.at(i).at(0)] << ", channel: " << (int)neighPulses.at(i).at(0) << endl ;
		
	}

//	end of main loop
  }


}
