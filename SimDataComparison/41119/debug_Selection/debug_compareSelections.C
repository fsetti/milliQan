// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain ) 
{

   ifstream infile;
   infile.open("outEventsNewMinusOld.txt");
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
	if ( filenum == fileIDs.at(i) && event == orig_events.at(i) ){
		passID = true;
		currentChannel = channels.at(i);
	}
    }
    if ( !passID ) continue;


    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > maxSlabdT ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

   
/////////////////////////////////////////////////////////////////////////////
//	Muon Selection - hits in straight line 

    vector<vector<vector<float>>> nonNeighPulses, neighPulses;
    neighPulses = SameLayerSelection_debug ( chan,height,ptime,nPE,time_module_calibrated, {trig01,trig024,trig125,trig2425,trig248,trig259,trig89,trig67,trig616,trig717,trig1617,trig1612,trig1713,trig1213,trig23,trig222,trig323,trig2223,trig224,trig235,trig45});


cout << "---------------------------------- \n current ch: " << currentChannel << endl;
	for ( unsigned int i=0; i<neighPulses.size(); i++){
		cout << "MUON   - reco nPE: " << neighPulses.at(i).at(0).at(3) << ", channel: " << (int)neighPulses.at(i).at(0).at(0) << endl ;
		cout << "reco nPE: " << neighPulses.at(i).at(1).at(3) << ", channel: " << (int)neighPulses.at(i).at(1).at(0) << endl ;
		cout << "\n " << endl;
	}
cout << "----------------------------------- \n " << endl;

//	end of main loop
  }


}
