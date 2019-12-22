{
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/SimDataComparison/beamMuons/selection.C+");
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C+");

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_name = "/home/users/fsetti/milliQan/MuonRate/MuonFiles/SimDataComparison.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( ch , output_file_name );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

//  plotHist( output_file_name, "SimData_nPEdT", true, true );
}
