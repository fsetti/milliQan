{
  gROOT->ProcessLine(".L selection.C+");
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C++");

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/sim_run3_*.root");
  output_file_name = "/home/users/fsetti/milliQan/OutputFiles/Sim_noBar.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( ch , output_file_name, true );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

  TString output_file_nameData;
  TChain *chD = new TChain("t");
  chD->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_nameData = "/home/users/fsetti/milliQan/OutputFiles/Data_noBar.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( chD , output_file_nameData );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

  TString outputPath = "/home/users/fsetti/public_html/milliQan/MuonRate/SimDataComparison/beamMuons/SlabNoBar/";

  ratioHists ( output_file_nameData, output_file_name, "SimData_nPE"     , outputPath);
  ratioHists ( output_file_nameData, output_file_name, "SimData_nPEET"   , outputPath);
  ratioHists ( output_file_nameData, output_file_name, "SimData_nPER878" , outputPath);
  ratioHists ( output_file_nameData, output_file_name, "SimData_nPER7725", outputPath);
//  plotHist ( output_file_name,"SimData_nPEdTall", true, true );
//  plotHist ( output_file_name,"SimData_nPEdTET", true, true );
//  plotHist ( output_file_name,"SimData_nPEdTR878", true, true );
//  plotHist ( output_file_name,"SimData_nPEdTR7725", true, true );
//  plotHist ( output_file_nameData, "SimData_nPEdTall", true, true );
//  plotHist ( output_file_nameData, "SimData_nPEdTET", true, true );
//  plotHist ( output_file_nameData, "SimData_nPEdTR878", true, true );
//  plotHist ( output_file_nameData, "SimData_nPEdTR7725", true, true );
}
