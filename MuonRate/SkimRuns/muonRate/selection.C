// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "muonAnalysis.cc"
#include "milliHists.h"


void selection( TChain *chain ) 
{

  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();

  int nSlabs=0;
  int nBars=0;

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);

//	skip non-synchronised events && require beam On
    bool synchBoards = false;			// only consider events with synched boards
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true; 
    if (!synchBoards ) continue;
    if ( !beam ) continue;



// require one hit in each slab, with dT btw first,last slab less than 12ns
    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > 12 ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

    nSlabs++;


/////////////////////////////////////////////////////////////////////////////
////      Layer selection
////      Bar selection I, require through going particles in upper, lower channels 
    vector<vector<vector<float>>> throughPulses; // give trigger configurations for straight line paths trig062 means Ch0,6,2, etc.
    throughPulses = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig241622,trig251723,trig8124,trig9135} );

    if ( throughPulses.size() != 0 ) nBars++;
//	end of main loop
  }


cout << " number of muons hitting the slabs: " << nSlabs << "\n number of muons hitting the slabs and bars in straight line: " << nBars << endl;

}
