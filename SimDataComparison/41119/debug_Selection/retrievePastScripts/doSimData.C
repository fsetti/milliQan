{
  gROOT->ProcessLine(".L selection.C+");


  TString output_file_nameData;
  TChain *chD = new TChain("t");
  chD->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_nameData = "/home/users/fsetti/milliQan/OutputFiles/41119/Data_Bar_retrieveOld.root";

  cout << "-------------------------" << endl;
  cout << "Skim of all Runs" << endl;
  selection( chD , output_file_nameData );
  cout << "-------------------------" << endl;
  cout << "\n \n \n " << endl;

}
