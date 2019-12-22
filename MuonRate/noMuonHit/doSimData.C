{
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/noMuonHit/selection.C+");

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_name = "/home/users/fsetti/milliQan/OutputFiles/Skim_noMuonHit.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( ch , output_file_name );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

}
