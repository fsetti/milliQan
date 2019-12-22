{
  gROOT->ProcessLine(".L /homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/SkimRuns/muonSelectionSkim.C+");
  gROOT->ProcessLine(".L /homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/plotHistsMuon.C+");

  bool Synch = true;
  TChain *ch = new TChain("t");
  ch->Add("/net/cms26/cms26r0/milliqan/haddedTrees/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  TString output_file_name = "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_muonStudy.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  muonSelectionSkim( ch , output_file_name, "Skim", {18,20,21,28} , Synch );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;


  TString nTimeDiffLayer1, nTimeDiffLayer2, nTimeDiffLayer3, nTimeDiffSecondPulsesL, nTimeDiffSecondPulses;
  nTimeDiffLayer = "RunSkim_TimeDiffSecondPulses";
  nTimeDiffLayerL = "RunSkim_TimeDiffSecondPulsesL";
  nTimeDiffLayer1 = "RunSkim_TimeDiffLayer1";
  nTimeDiffLayer2 = "RunSkim_TimeDiffLayer2";
  nTimeDiffLayer3 = "RunSkim_TimeDiffLayer3";
  compareHists( output_file_name, nTimeDiffLayer, nTimeDiffLayerL );
  stackHists( output_file_name, nTimeDiffLayer1, nTimeDiffLayer2, nTimeDiffLayer3 );
}
