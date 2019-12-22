// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"
//#include "/home/users/fsetti/milliQan/Headers/MCntupleHandler.h"


void selection( TChain *chain, TString output_filename, bool Data = true ) 
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


  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

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
    if ( (!synchBoards  && Data) || ( !beam && Data ) ) continue;



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

/////////////////////////////////////////////////////////////////////////////
//	Organise the plotting  

    if ( Data ) {
	for ( unsigned int i=0; i<nonNeighPulses.size(); i++){
		h_nNnPE->Fill ( nonNeighPulses.at(i).at(3) );
	}
	//	repeat for neighbouring channels
	for ( unsigned int i=0; i<neighPulses.size(); i++){
		h_NnPE->Fill ( neighPulses.at(i).at(3) );
    	}
    }


    if ( !Data ) {
        currentFile = chain->GetFile();
        currentFileName = currentFile->GetName();
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcd.root" ) weight = weight_qcd;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_qcd_nonbc.root" ) weight = weight_qcd_nonbc;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_dy.root" ) weight = weight_DY;
        if ( currentFileName == "/nfs-7/userdata/fsetti/sim_run3_w.root" ) weight = weight_W;

	for ( unsigned int i=0; i<nonNeighPulses.size(); i++){
		h_nNnPE->Fill ( nonNeighPulses.at(i).at(3) , weight );
//		h_nNnPE->Fill ( chan_trueNPE[(int)nonNeighPulses.at(i).at(0)], weight );
	}	
	
	//	repeat for neighbouring channels
	for ( unsigned int i=0; i<neighPulses.size(); i++){
		h_NnPE->Fill ( neighPulses.at(i).at(3) , weight );
//		h_NnPE->Fill ( chan_trueNPE[(int)neighPulses.at(i).at(0)], weight );
	}
    }

//	end of main loop
  }
  outFile->cd();

  h_NnPE->Write();
  h_nNnPE->Write();

  outFile->Write();
  outFile->Close();

}
