{
  gROOT->ProcessLine(".L /homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/SkimRuns/SimComparison/compareSimData_horizontal.C+");

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/net/cms26/cms26r0/milliqan/haddedTrees/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_name = "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_SimComparison.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  compareSimData( ch , output_file_name, "Skim", {18,20,21,28} );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

}
