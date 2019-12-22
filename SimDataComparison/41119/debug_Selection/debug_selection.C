// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain ) 
{

  TFile *currentFile;
  TString currentFileName;

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


  ofstream outfile;
  outfile.open("newSelection.txt");
//   ifstream infile;
//   infile.open("outEvts.txt");
//   vector<int> fileIDs, orig_events, channels;
//   int fID, oEvts, chanN;
//   while ( infile >> fID >> oEvts >> chanN ){
//	fileIDs.push_back(fID);
//	orig_events.push_back(oEvts);
//	channels.push_back(chanN);
//   }

  TH1::SetDefaultSumw2();
  TH1D *h_nPE, *h_nNnPE, *h_NnPE, *h_dT;
  double nPEbins[23] = { 0., 10., 20., 30., 50., 70., 90., 120., 150., 180., 210., 240., 270., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750. };
  h_nNnPE = new TH1D("SimData_nNnPE","nPE",  22, nPEbins  );
  h_NnPE = new TH1D("SimData_NnPE","nPE",22, nPEbins  );
  h_dT = new TH1D("dT_slabs","#Delta_{t} slabs", 26,-13,13);


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
    if ( !synchBoards  ||  !beam  ) continue;


    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > maxSlabdT ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

   
/////////////////////////////////////////////////////////////////////////////
//	Muon Selection - hits in straight line 

    vector<vector<float>> nonNeighPulses, neighPulses;
    neighPulses = NeighBarSelection ( chan,height,ptime,nPE,time_module_calibrated, trigBars);

/////////////////////////////////////////////////////////////////////////////
//	Organise the plotting  
    for ( unsigned int i=0; i<neighPulses.size(); i++){
    	outfile << filenum << " " << event << " " << (int)neighPulses.at(i).at(0) << endl;
    }

//	end of main loop
  }

  outfile.close();

}
