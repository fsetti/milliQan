{
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C+");

  output_file_nameNew =    "/home/users/fsetti/milliQan/OutputFiles/281019/Data_Bar.root";
  output_file_name    =    "/home/users/fsetti/milliQan/OutputFiles/281019/Data_Bar_retrieveOld.root";
//  output_file_nameNew =    "/home/users/fsetti/milliQan/OutputFiles/281019/Data_Bar_retrieveOld.root";
//  output_file_name    =    "/home/users/bemarsh/analysis/milliqan/geantDemoSim/slim_ntupler/datacomps/franny_output.root";

  TString NoutPath = "/home/users/fsetti/public_html/milliQan/SimDataComparison/281019/MuonHit/neigh_";

//  compareHists_v2 ( output_file_nameNew, output_file_name, "SimData_NnPE", true );
  compareHists_v2 ( output_file_nameNew, output_file_name, "SimData_nNnPE",  true );
}
