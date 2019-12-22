{
  gROOT->ProcessLine(".L selection.C+");
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C++");

  TString output_file_nameData;
  TChain *chD = new TChain("t");
  chD->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_nameData = "/home/users/fsetti/milliQan/OutputFiles/Data_Bar.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( chD , output_file_nameData );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/sim_run3_*.root");
  output_file_name = "/home/users/fsetti/milliQan/OutputFiles/Sim_Bar.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
//  selection( ch , output_file_name, true );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

  TString outPath = "/home/users/fsetti/public_html/milliQan/SimDataComparison/beamMuons/nNBars/";
  TString NoutPath = "/home/users/fsetti/public_html/milliQan/SimDataComparison/beamMuons/NBars/";



//  ratioHists ( output_file_nameData, output_file_name, "SimData_nPE"  , outPath );
//  ratioHists ( output_file_nameData, output_file_name, "SimData_nPEET"  , outPath );
//  ratioHists ( output_file_nameData, output_file_name, "SimData_nPER878"  , outPath );
//  ratioHists ( output_file_nameData, output_file_name, "SimData_nPER7725"  , outPath );
//
//  ratioHists ( output_file_nameData, output_file_name, "SimData_NnPE" , NoutPath);
//  ratioHists ( output_file_nameData, output_file_name, "SimData_NnPEET" , NoutPath);
//  ratioHists ( output_file_nameData, output_file_name, "SimData_NnPER878" , NoutPath);
//  ratioHists ( output_file_nameData, output_file_name, "SimData_NnPER7725" , NoutPath);
//  plotHist ( output_file_name, "SimData_nPEdTall" , true, true );
//  plotHist ( output_file_name, "SimData_NnPEdTall" , true, true );
//  plotHist ( output_file_nameData, "SimData_nPEdTall" , true, true );
//  plotHist ( output_file_nameData, "SimData_NnPEdTall" , true, true );

}
