{
  gROOT->ProcessLine(".L selection.C+");

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/sim_run3_*.root");
  output_file_name = "/home/users/fsetti/milliQan/OutputFiles/Sim_Bar.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( ch , output_file_name );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

}
