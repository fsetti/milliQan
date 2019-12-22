{
  gROOT->ProcessLine(".L /homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/SkimRuns/muonSelectionSkim.C+");
  gROOT->ProcessLine(".L /homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/plotHistsMuon.C+");

  bool Synch = true;
  bool MuonBarHit = true;
  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/net/cms26/cms26r0/milliqan/haddedTrees/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  if ( MuonBarHit )  output_file_name = "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_muonStudy.root";
  if ( !MuonBarHit ) output_file_name = "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_muonStudyNoBar.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  muonSelectionSkim( ch , output_file_name, "Skim", {18,20,21,28} , Synch, MuonBarHit );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

//  compareHists_v2 (  output_file_name, "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_muonStudyNoBar.root" , "RunSkim_tDiffSlab" );
//  compareHists_v2 (  output_file_name, "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_muonStudyNoBar.root" , "RunSkim_MuHitsSL" );
//  compareHists_v2 (  output_file_name, "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_muonStudyNoBar.root" , "RunSkim_nPE", true );
//  compareProfiles(   output_file_name, "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Skim_muonStudyNoBar.root" , "profile" );
//  compareHists (  output_file_name,  "RunSkim_tDiffSlab", "RunSkim_tDiffAll" );
//  plotHist ( output_file_name, "RunSkim_tDiffSlab" );
//  plotHist ( output_file_name, "RunSkim_MuHitsSL" );
}
