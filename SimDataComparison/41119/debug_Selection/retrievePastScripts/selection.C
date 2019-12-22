// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename, bool Data = true ) 
{
  ofstream outfile;
  outfile.open("oldSelection.txt");


  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);

    bool synchBoards = false;                   // only consider events with synched boards
    if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true;
    if ( (!synchBoards  && Data) || ( !beam && Data ) ) continue;



    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > maxSlabdT ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

   
/////////////////////////////////////////////////////////////////////////////
//	Muon Selection - hits in straight line 

    vector<vector<vector<float>>> neighPulses, nonNeighPulses;
    neighPulses = SameLayerSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig01,trig024,trig125,trig2425,trig248,trig259,trig89,trig67,trig616,trig717,trig1617,trig1612,trig1713,trig1213,trig23,trig222,trig323,trig2223,trig224,trig235,trig45});


    for ( unsigned int i=0; i<neighPulses.size(); i++){
        outfile << filenum << " " << event << " " << (int)neighPulses.at(i).at(1).at(0) << endl;
    }


//	end of main loop
  }

  outfile.close();

}
