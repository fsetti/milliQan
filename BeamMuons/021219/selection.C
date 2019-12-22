// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"
#include "/home/users/fsetti/milliQan/Headers/milliHists.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"


void selection( TChain *chain, TString output_filename ) 
{

  TString ChStr[32];
  ChStr[0] = "0"; ChStr[1] = "1"; ChStr[2] = "2"; ChStr[3] = "3"; ChStr[4] = "4"; ChStr[5] = "5"; ChStr[6] = "6"; ChStr[7] = "7"; ChStr[8] = "8"; ChStr[9] = "9";
  ChStr[10] = "10"; ChStr[11] = "11"; ChStr[12] = "12"; ChStr[13] = "13"; ChStr[14] = "14"; ChStr[15] = "15"; ChStr[16] = "16"; ChStr[17] = "17"; ChStr[18] = "18"; ChStr[19] = "19"; ChStr[10] = "10";
  ChStr[20] = "20"; ChStr[21] = "21"; ChStr[22] = "22"; ChStr[23] = "23"; ChStr[24] = "24"; ChStr[25] = "25"; ChStr[26] = "26"; ChStr[27] = "27"; ChStr[28] = "28"; ChStr[29] = "29";
  ChStr[30] = "30"; ChStr[31] = "31";


  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");
  outFile->cd();

  TH1D *h_nPE, *h_dT878, *h_dTET, *h_dT7725, *h_dT;
  float increment = (maxBin-minBin)/nBins;
  Double_t edges[nBins+1];
  float x=minBin;
  for ( unsigned int i=0; i<nBins+1; i++){
	edges[i]= pow(10,x);
	x+=increment;
  }
  TH1D *h_RMS;
  h_RMS = new TH1D ("h_RMS", "RMS on |#Delta_{T}|", nBins, edges );
  h_dT = new TH1D( "dT"," #Delta_{t} after pulses" , 40,0, 200 );
  h_dT878 = new TH1D( "dT878"," #Delta_{t} 878 after pulses" , 40, 0 , 200 );
  h_dTET = new TH1D( "dTET"," #Delta_{t} ET after pulses" , 40, 0 , 200 );
  h_dT7725 = new TH1D( "dT7725"," #Delta_{t} 7725 after pulses" , 40, 0 , 200 );

  TH2D *h_nAfterPulses, *h_nAP878, *h_nAPET, *h_nAP7725;
  h_nAfterPulses = new TH2D ("h_nAfterPulses","# afterpulses",nBins,edges, 10,0,10);
  h_nAP878  = new TH2D ("h_nAP878","# afterpulses R878",nBins,edges, 10,0,10);
  h_nAPET   = new TH2D ("h_nAPET","# afterpulses ET ",nBins,edges, 10,0,10);
  h_nAP7725 = new TH2D ("h_nAP7725","# afterpulses R7725",nBins,edges, 10,0,10);

  TH2D *h_nPEdT, *h_nPEdTsameLayer, *h_nPEdT2;
  h_nPEdT = new TH2D("nPEdT","nPE vs dT", 80, -100, 100, nBins, edges );
  h_nPEdT2 = new TH2D("nPEdT2","max nPE vs |dT|", 80, -100, 100, nBins, edges );
  h_nPEdTsameLayer = new TH2D("nPEdT_sameLayer","nPE vs dT", 80, -100, 100, nBins, edges );
  TProfile *hprof, *hprofX;
  hprof  = new TProfile("hprof","nPE vs dT",nBins, edges, "s" );
  hprofX  = new TProfile("hprofX","nPE vs dT",nBins, edges );

  InitializeChain(chain);

  //Number of events to loop over
  Int_t nentries = (Int_t)chain->GetEntries();


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//	Main Event Loop 
  for(int ia = 0; ia<nentries; ia++){

    chain->GetEntry(ia);

    vector<vector<float>> slabPulses;
    slabPulses = SlabSelection( chan,height,ptime,nPE,time_module_calibrated );
    if ( slabPulses.size() == 0 ) continue;
    std::sort(slabPulses.begin(), slabPulses.end(), sort_wrt_channel);
    float dt_cal = slabPulses.at(2)[4] - slabPulses.at(0)[4];
    if ( fabs( dt_cal ) > 12 ) continue;	// require dt_cal between -12 and +12 ns from dt_cal distribution

   
/////////////////////////////////////////////////////////////////////////////
//	Muon Selection - hits in straight line 
    vector<vector<vector<float>>> muonPulses;
    muonPulses = pulseBarSelection ( chan,height,ptime,nPE,time_module_calibrated, {trig062,trig173,trig241622,trig251723,trig8124,trig9135} );
    if ( muonPulses.size() == 0 ) continue;

    vector<vector<float>> neighPulsesL1, neighPulsesL2, neighPulsesL3;
    neighPulsesL1 = NeighBarSelection_diag ( chan,height,ptime,nPE,time_module_calibrated, trigL1 );
    neighPulsesL2 = NeighBarSelection_diag ( chan,height,ptime,nPE,time_module_calibrated, trigL2 );
    neighPulsesL3 = NeighBarSelection_diag ( chan,height,ptime,nPE,time_module_calibrated, trigL3 );




// store dT of after pulses
    for ( unsigned int i=0; i<neighPulsesL1.size(); i++){
	for (unsigned int j=0; j<chan->size(); j++){
		if ( (int)neighPulsesL1.at(i).at(0) == chan->at(j) && neighPulsesL1.at(i).at(3) != nPE->at(j) && nPE->at(j) < nPEminVec.at(chan->at(j)) && nPE->at(j) > nPEmin  && ptime->at(j) > tMin )
		{
			h_dT->Fill ( - neighPulsesL1.at(i).at(4) + time_module_calibrated->at(j) );
			if ( chan->at(j) == 5 || chan->at(j) == 22 ) h_dT7725->Fill ( - neighPulsesL1.at(i).at(4) + time_module_calibrated->at(j) );
			if ( chan->at(j) == 9 || chan->at(j) == 17 || chan->at(j) == 24 || chan->at(j) == 25 ) h_dTET->Fill ( - neighPulsesL1.at(i).at(4) + time_module_calibrated->at(j) );
			else h_dT878->Fill ( - neighPulsesL1.at(i).at(4) + time_module_calibrated->at(j) );
		}
	}
    }

    for ( unsigned int i=0; i<neighPulsesL2.size(); i++){
	for (unsigned int j=0; j<chan->size(); j++){
		if ( (int)neighPulsesL2.at(i).at(0) == chan->at(j) && neighPulsesL2.at(i).at(3) != nPE->at(j) && nPE->at(j) < nPEminVec.at(chan->at(j)) && nPE->at(j) > nPEmin && ptime->at(j) > tMin )
		{
			h_dT->Fill ( - neighPulsesL2.at(i).at(4) + time_module_calibrated->at(j) );
			if ( chan->at(j) == 5 || chan->at(j) == 22 ) h_dT7725->Fill ( - neighPulsesL2.at(i).at(4) + time_module_calibrated->at(j) );
			if ( chan->at(j) == 9 || chan->at(j) == 17 || chan->at(j) == 24 || chan->at(j) == 25 ) h_dTET->Fill ( - neighPulsesL2.at(i).at(4) + time_module_calibrated->at(j) );
			else h_dT878->Fill ( - neighPulsesL2.at(i).at(4) + time_module_calibrated->at(j) );
		}
	}
    }

    for ( unsigned int i=0; i<neighPulsesL3.size(); i++){
	for (unsigned int j=0; j<chan->size(); j++){
		if ( (int)neighPulsesL3.at(i).at(0) == chan->at(j) && neighPulsesL3.at(i).at(3) != nPE->at(j) && nPE->at(j) < nPEminVec.at(chan->at(j)) && nPE->at(j) > nPEmin && ptime->at(j) > tMin )
		{
			h_dT->Fill ( - neighPulsesL3.at(i).at(4) + time_module_calibrated->at(j) );
			if ( chan->at(j) == 5 || chan->at(j) == 22 ) h_dT7725->Fill ( - neighPulsesL3.at(i).at(4) + time_module_calibrated->at(j) );
			if ( chan->at(j) == 9 || chan->at(j) == 17 || chan->at(j) == 24 || chan->at(j) == 25 ) h_dTET->Fill ( - neighPulsesL3.at(i).at(4) + time_module_calibrated->at(j) );
			else h_dT878->Fill ( - neighPulsesL3.at(i).at(4) + time_module_calibrated->at(j) );
		}
	}
    }



//	different layer activity
    for (unsigned int i=0; i<neighPulsesL1.size(); i++){
	for (unsigned int j=0; j<neighPulsesL2.size(); j++){
	  if (std::max(neighPulsesL1.at(i).at(3), neighPulsesL2.at(j).at(3)) / std::min(neighPulsesL1.at(i).at(3), neighPulsesL2.at(j).at(3)) > maxNPEratio ){
		continue;  }
	  h_nPEdT->Fill( neighPulsesL2.at(j).at(4) - neighPulsesL1.at(i).at(4),  std::max(neighPulsesL1.at(i).at(3), neighPulsesL2.at(j).at(3)) );
	  h_nPEdT2->Fill( fabs(neighPulsesL2.at(j).at(4) - neighPulsesL1.at(i).at(4)),  std::max(neighPulsesL1.at(i).at(3), neighPulsesL2.at(j).at(3)) );
	  hprof->Fill(   std::max(neighPulsesL1.at(i).at(3), neighPulsesL2.at(j).at(3)), neighPulsesL2.at(j).at(4) - neighPulsesL1.at(i).at(4));
	  hprofX->Fill(   std::max(neighPulsesL1.at(i).at(3), neighPulsesL2.at(j).at(3)), fabs(neighPulsesL2.at(j).at(4) - neighPulsesL1.at(i).at(4)) );
	}
	for (unsigned int j=0; j<neighPulsesL3.size(); j++){
	  if (std::max(neighPulsesL1.at(i).at(3), neighPulsesL3.at(j).at(3)) / std::min(neighPulsesL1.at(i).at(3), neighPulsesL3.at(j).at(3)) > maxNPEratio ){
		continue;  }
	  h_nPEdT->Fill( neighPulsesL3.at(j).at(4) - neighPulsesL1.at(i).at(4), std::max(neighPulsesL1.at(i).at(3), neighPulsesL3.at(j).at(3)) );
	  h_nPEdT2->Fill( fabs(neighPulsesL3.at(j).at(4) - neighPulsesL1.at(i).at(4)), std::max(neighPulsesL1.at(i).at(3), neighPulsesL3.at(j).at(3)) );
	  hprof->Fill( std::max(neighPulsesL1.at(i).at(3), neighPulsesL3.at(j).at(3)) , neighPulsesL3.at(j).at(4) - neighPulsesL1.at(i).at(4) );
	  hprofX->Fill( std::max(neighPulsesL1.at(i).at(3), neighPulsesL3.at(j).at(3)) , fabs(neighPulsesL3.at(j).at(4) - neighPulsesL1.at(i).at(4)) );
	}
    }
    for (unsigned int i=0; i<neighPulsesL2.size(); i++){
	for (unsigned int j=0; j<neighPulsesL3.size(); j++){
	  if (std::max(neighPulsesL2.at(i).at(3), neighPulsesL3.at(j).at(3)) / std::min(neighPulsesL2.at(i).at(3), neighPulsesL3.at(j).at(3)) > maxNPEratio ){
		continue;  }
	  h_nPEdT->Fill( neighPulsesL3.at(j).at(4) - neighPulsesL2.at(i).at(4), std::max(neighPulsesL2.at(i).at(3), neighPulsesL3.at(j).at(3)) );
	  h_nPEdT2->Fill( fabs(neighPulsesL3.at(j).at(4) - neighPulsesL2.at(i).at(4)), std::max(neighPulsesL2.at(i).at(3), neighPulsesL3.at(j).at(3)) );
	  hprof->Fill(  std::max(neighPulsesL2.at(i).at(3), neighPulsesL3.at(j).at(3)), neighPulsesL3.at(j).at(4) - neighPulsesL2.at(i).at(4) );
	  hprofX->Fill(  std::max(neighPulsesL2.at(i).at(3), neighPulsesL3.at(j).at(3)), fabs(neighPulsesL3.at(j).at(4) - neighPulsesL2.at(i).at(4)) );
	}
    }




 
//	same layer activity & after pulses 
    for (unsigned int i=0; i<neighPulsesL1.size(); i++){	

// sort by PMT type
    if ( (int)neighPulsesL1.at(i).at(0) == 5 || (int)neighPulsesL1.at(i).at(0) == 22 ) h_nAP7725->Fill(neighPulsesL1.at(i).at(3), neighPulsesL1.at(i).at(5));
    if ( (int)neighPulsesL1.at(i).at(0) == 9 || (int)neighPulsesL1.at(i).at(0) == 17 || (int)neighPulsesL1.at(i).at(0) == 24 || (int)neighPulsesL1.at(i).at(0) == 25 ) h_nAPET->Fill(neighPulsesL1.at(i).at(3), neighPulsesL1.at(i).at(5));
    if ( (int)neighPulsesL1.at(i).at(0) != 5 && (int)neighPulsesL1.at(i).at(0) != 22 && (int)neighPulsesL1.at(i).at(0) != 9 && (int)neighPulsesL1.at(i).at(0) != 17 && (int)neighPulsesL1.at(i).at(0) != 24 && (int)neighPulsesL1.at(i).at(0) != 25 ) h_nAP878->Fill(neighPulsesL1.at(i).at(3), neighPulsesL1.at(i).at(5));

    h_nAfterPulses->Fill ( neighPulsesL1.at(i).at(3), neighPulsesL1.at(i).at(5) );
    	for (unsigned int j=i+1; j<neighPulsesL1.size(); j++){
	  if (std::max(neighPulsesL1.at(i).at(3), neighPulsesL1.at(j).at(3)) / std::min(neighPulsesL1.at(i).at(3), neighPulsesL1.at(j).at(3)) < maxNPEratio ){
	    h_nPEdTsameLayer->Fill( neighPulsesL1.at(j).at(4) - neighPulsesL1.at(i).at(4),  std::max(neighPulsesL1.at(i).at(3), neighPulsesL1.at(j).at(3)) );
	  }
	}
    }
    for (unsigned int i=0; i<neighPulsesL2.size(); i++){

// sort by PMT type
    if ( (int)neighPulsesL2.at(i).at(0) == 5 || (int)neighPulsesL2.at(i).at(0) == 22 ) h_nAP7725->Fill(neighPulsesL2.at(i).at(3), neighPulsesL2.at(i).at(5));
    if ( (int)neighPulsesL2.at(i).at(0) == 9 || (int)neighPulsesL2.at(i).at(0) == 17 || (int)neighPulsesL2.at(i).at(0) == 24 || (int)neighPulsesL2.at(i).at(0) == 25 ) h_nAPET->Fill(neighPulsesL2.at(i).at(3), neighPulsesL2.at(i).at(5));
    if ( (int)neighPulsesL2.at(i).at(0) != 5 && (int)neighPulsesL2.at(i).at(0) != 22 && (int)neighPulsesL2.at(i).at(0) != 9 && (int)neighPulsesL2.at(i).at(0) != 17 && (int)neighPulsesL2.at(i).at(0) != 24 && (int)neighPulsesL2.at(i).at(0) != 25 ) h_nAP878->Fill(neighPulsesL2.at(i).at(3), neighPulsesL2.at(i).at(5));


    h_nAfterPulses->Fill ( neighPulsesL2.at(i).at(3), neighPulsesL2.at(i).at(5) );
    	for (unsigned int j=i+2; j<neighPulsesL2.size(); j++){
	  if (std::max(neighPulsesL2.at(i).at(3), neighPulsesL2.at(j).at(3)) / std::min(neighPulsesL2.at(i).at(3), neighPulsesL2.at(j).at(3)) < maxNPEratio ){
	    h_nPEdTsameLayer->Fill( neighPulsesL2.at(j).at(4) - neighPulsesL2.at(i).at(4),  std::max(neighPulsesL2.at(i).at(3), neighPulsesL2.at(j).at(3)) );
	  }
	}
    }
    for (unsigned int i=0; i<neighPulsesL3.size(); i++){

// sort by PMT type
    if ( (int)neighPulsesL3.at(i).at(0) == 5 || (int)neighPulsesL3.at(i).at(0) == 22 ) h_nAP7725->Fill(neighPulsesL3.at(i).at(3), neighPulsesL3.at(i).at(5));
    if ( (int)neighPulsesL3.at(i).at(0) == 9 || (int)neighPulsesL3.at(i).at(0) == 17 || (int)neighPulsesL3.at(i).at(0) == 24 || (int)neighPulsesL3.at(i).at(0) == 25 ) h_nAPET->Fill(neighPulsesL3.at(i).at(3), neighPulsesL3.at(i).at(5));
    if ( (int)neighPulsesL3.at(i).at(0) != 5 && (int)neighPulsesL3.at(i).at(0) != 22 && (int)neighPulsesL3.at(i).at(0) != 9 && (int)neighPulsesL3.at(i).at(0) != 17 && (int)neighPulsesL3.at(i).at(0) != 24 && (int)neighPulsesL3.at(i).at(0) != 25 ) h_nAP878->Fill(neighPulsesL3.at(i).at(3), neighPulsesL3.at(i).at(5));


    h_nAfterPulses->Fill ( neighPulsesL3.at(i).at(3), neighPulsesL3.at(i).at(5) );
    	for (unsigned int j=i+3; j<neighPulsesL3.size(); j++){
	  if (std::max(neighPulsesL3.at(i).at(3), neighPulsesL3.at(j).at(3)) / std::min(neighPulsesL3.at(i).at(3), neighPulsesL3.at(j).at(3)) < maxNPEratio ){
	    h_nPEdTsameLayer->Fill( neighPulsesL3.at(j).at(4) - neighPulsesL3.at(i).at(4),  std::max(neighPulsesL3.at(i).at(3), neighPulsesL3.at(j).at(3)) );
	  }
	}
    }
//	end of main loop
  }


  for (unsigned int i=0; i<nBins; i++){
	h_RMS->SetBinContent( i+1, hprof->GetBinError(i+1) );
  }




  outFile->cd();
  outFile->Write();
  outFile->Close();

}



void Projections (TString path1, TString h1Name, TString outPath ){


        TFile* file1 = TFile::Open(path1);
        TH2D *h1_Copy = (TH2D*)file1->Get(h1Name);
        TH2D *h1 = (TH2D*)h1_Copy->Clone(h1Name);
        h1->SetDirectory(0);
        file1->Close();

	TFile *outFile = new TFile(outPath.Data(),"RECREATE");
	outFile->cd();

	TH1D *h_Proj[nBins], *h_Sigmas;
	h_Sigmas = new TH1D( "sigmas","sigmas", nBins, 0 ,nBins );
	for (unsigned int i=3; i<nBins-1; i++){
		TString name = "projX_"+std::to_string(i+1);
		h_Proj[i] = h1->ProjectionX(name, i+1,i+1 );
		TFitResultPtr r = h_Proj[i]->Fit("gaus","S","",-40,40);
		Int_t fitStatus = r;
		if ( fitStatus == 0 ){
			Double_t sigma = r->Parameter(2);
			h_Sigmas->SetBinContent( i+1 , sigma );
		}
	}

	outFile->Write();
	outFile->Close();
}




void Quadrature (TString path1, TString path2, TString h1Name, TString outPath ){

        TFile* file1 = TFile::Open(path1);
        TH2D *h1_Copy = (TH2D*)file1->Get(h1Name);
        TH2D *h1 = (TH2D*)h1_Copy->Clone(h1Name);
        h1->SetDirectory(0);
        file1->Close();

        TFile* file2 = TFile::Open(path2);
        TH2D *h2_Copy = (TH2D*)file2->Get(h1Name);
        TH2D *h2 = (TH2D*)h2_Copy->Clone(h1Name);
        h2->SetDirectory(0);
        file2->Close();

	TFile *outFile = new TFile(outPath.Data(),"RECREATE");
	outFile->cd();

	TH1D *h_Quad;
	h_Quad = new TH1D( "quad","#sqrt{x_{Data}^{2}-x_{MC}^{2}}", nBins, 0 ,nBins );
	for (unsigned int i=0; i<nBins; i++){
		double xData = h1->GetBinContent(i+1);
		double xSim  = h2->GetBinContent(i+1);
		h_Quad->SetBinContent(i+1, TMath::Sqrt( xData*xData - xSim*xSim ) );
	}

	outFile->Write();
	outFile->Close();
}
