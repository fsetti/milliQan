// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.

#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"


void selection( TChain *chain, TString output_filename ) 
{


  vector<int> runNum;
  vector<float> runLumi; // in pb^-1
  FILE *fillFile;
  fillFile = fopen("runLumi.txt","r");
  if (!fillFile) {
    cerr << "Could not open fill file.\n";
  }
  else { // Opened file successfully?  
    int ncol = 0;
    while (ncol >= 0) {
      int intVal; float fVal;
      ncol=fscanf(fillFile,"%i", &intVal);
      runNum.push_back(intVal);
      ncol=fscanf(fillFile,"%f", &fVal);
      runLumi.push_back(fVal);
//	cout <<"run: " << intVal << ", and total lumi: " << fVal << endl;
      char c;
      c = '0'; while (c != '\n' && ncol>0) ncol=fscanf(fillFile,"%c",&c);
      }
    fclose(fillFile);
  } // Opened file successfully?



	TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
	outFile->cd();

	const unsigned int nBins = runNum.size();	
	TH1D *h_NnPE, *h_Nunscaled, *h_nNnPE, *h_nNunscaled;
	h_NnPE = new TH1D("neigh_nPEs","1<nPE<10 / pb^{-1}",nBins+1,0,nBins+1);
	h_Nunscaled = new TH1D("Nunscaled","1<nPE<10 unscaled",nBins+1,0,nBins+1);
	h_nNnPE = new TH1D("nonNeigh_nPEs","1<nPE<10 / pb^{-1}",nBins+1,0,nBins+1);
	h_nNunscaled = new TH1D("nNunscaled","1<nPE<10 unscaled",nBins+1,0,nBins+1);

	InitializeChain(chain);
	//Number of events to loop over
	Int_t nentries = (Int_t)chain->GetEntries();
	
	int curRun, prevRun;
	float lumi;
	curRun = 0;
	prevRun = 0;
	lumi = -1;
	int idx = 0;
	vector<unsigned int> nNnPEs(nBins,0);
	vector<unsigned int> NnPEs(nBins,0);
	vector<float> weights(nBins,0);
///////////////////////////////////////////////////////////////////////////
// Main Event Loop 
	for(int ia = 0; ia<nentries; ia++){
	
		chain->GetEntry(ia);
		curRun = run;
		if ( curRun != prevRun ){
			for ( unsigned int i=0; i<runNum.size(); i++){
				if ( run == runNum.at(i) ) {
					lumi = runLumi.at(i);
					idx = i;
					if ( lumi != 0 ) weights.at(idx) = 1./lumi;
				}
			}
		}
		if ( lumi == 0 ) continue;

		bool synchBoards = false;                   // only consider events with synched boards
		if ((event_trigger_time_tag_b1==event_trigger_time_tag_b0) || (groupTDC_b1->at(0) == groupTDC_b0->at(0))) synchBoards = true;
		if ( !synchBoards  ||  !beam  ) continue;
		
		vector<vector<float>> slabPulses;
		slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
		if ( slabPulses.size() == 0 ) continue;
		std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
		float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
		if ( fabs( dt_cal ) > 12 ) continue;        // require dt_cal between -12 and +12 ns from dt_cal distribution


		vector<vector<vector<float>>> nonNeighPulses, neighPulses;
		nonNeighPulses = SameLayerSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {trig18,trig19,trig08,trig09,trig612,trig613,trig712,trig713,trig24,trig25,trig34,trig35}, false );
		neighPulses = SameLayerSelection_v2 ( chan,height,ptime,nPE,time_module_calibrated, {trig01,trig024,trig125,trig2425,trig248,trig259,trig89,trig67,trig616,trig717,trig1617,trig1612,trig1713,trig1213,trig23,trig222,trig323,trig2223,trig224,trig235,trig45} );


		for ( unsigned int i=0; i<neighPulses.size(); i++){
			if ( lumi != 0 && neighPulses.at(i).at(1).at(3) > 1 && neighPulses.at(i).at(1).at(3) < 10 ) NnPEs.at(idx)+=1;
		}
		for ( unsigned int i=0; i<nonNeighPulses.size(); i++){
			if ( lumi != 0 && nonNeighPulses.at(i).at(1).at(3) > 1 && nonNeighPulses.at(i).at(1).at(3) < 10 ) nNnPEs.at(idx)+=1;
		}
	
		prevRun = curRun;
		
//	end of main loop
	}
				
 for (unsigned int i=0; i<nBins; i++){
	if ( NnPEs.at(i) == 0 ) continue;
	TString label;
	if ( i%2 == 0 ) label = std::to_string(runNum.at(i));
	if ( i%2 != 0 ) label = std::to_string(runNum.at(i))+"    ";
	h_NnPE->GetXaxis()->SetBinLabel(i+1,label );	
	h_Nunscaled->GetXaxis()->SetBinLabel(i+1,label );	
	h_NnPE->SetBinContent( i+1, NnPEs.at(i)*weights.at(i) );
	h_Nunscaled->SetBinContent( i+1, NnPEs.at(i));
 }
 for (unsigned int i=0; i<nBins; i++){
	if ( nNnPEs.at(i) == 0 ) continue;
	TString label;
	if ( i%2 == 0 ) label = std::to_string(runNum.at(i));
	if ( i%2 != 0 ) label = std::to_string(runNum.at(i))+"    ";
	h_nNnPE->GetXaxis()->SetBinLabel(i+1,label );	
	h_nNunscaled->GetXaxis()->SetBinLabel(i+1,label );	
	h_nNnPE->SetBinContent( i+1, nNnPEs.at(i)*weights.at(i) );
	h_nNunscaled->SetBinContent( i+1, nNnPEs.at(i));
 }

 h_NnPE->LabelsOption("v","X");	
 h_Nunscaled->LabelsOption("v","X");	
 h_nNnPE->LabelsOption("v","X");	
 h_nNunscaled->LabelsOption("v","X");	

 outFile->cd();

 h_NnPE->Write();
 h_Nunscaled->Write();
 h_nNnPE->Write();

 h_nNunscaled->Write();
 outFile->Write();
 outFile->Close();

}
