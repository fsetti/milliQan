{
//  gROOT->ProcessLine(".L debug_compareSelections.C+");
  gROOT->ProcessLine(".L debug_selection.C+");

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( ch  );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

}
