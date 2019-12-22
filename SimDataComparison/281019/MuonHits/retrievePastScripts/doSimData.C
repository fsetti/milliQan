{
  gROOT->ProcessLine(".L selection.C+");
//  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C++");


  TString output_file_nameData;
  TChain *chD = new TChain("t");
  chD->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_nameData = "/home/users/fsetti/milliQan/OutputFiles/281019/Data_Bar_retrieveOld.root";

  cout << "-------------------------" << endl;
  cout << "Skim of all Runs" << endl;
  selection( chD , output_file_nameData );
  cout << "-------------------------" << endl;
  cout << "\n \n \n " << endl;

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/sim_run3_*.root");
  output_file_name = "/home/users/fsetti/milliQan/OutputFiles/281019/Sim_Bar_retrieveOld.root";

  cout << "-------------------------" << endl;
  cout << "Skim of all Runs" << endl;
//  selection( ch , output_file_name, false );    
  cout << "-------------------------" << endl;
  cout << "\n \n \n " << endl;

  TString outPath =  "/home/users/fsetti/public_html/milliQan/SimDataComparison/281019/MuonHit/Old_nonNeigh_";
  TString NoutPath = "/home/users/fsetti/public_html/milliQan/SimDataComparison/281019/MuonHit/Old_neigh_";

//  ratioHists ( output_file_nameData, output_file_name, "SimData_nPE"  , outPath );
//  ratioHists ( output_file_nameData, output_file_name, "SimData_NnPE" , NoutPath);
}
