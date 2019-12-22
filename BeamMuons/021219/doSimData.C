{
  gROOT->ProcessLine(".L selection.C+");
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C+");

  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/sim_run3_*.root");
  output_file_name = "/home/users/fsetti/milliQan/OutputFiles/021219/Sim_ratio10.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( ch , output_file_name );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

  TString output_file_nameData;
  TChain *chD = new TChain("t");
  chD->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_threeLayerHit_191204.root");
  output_file_nameData = "/home/users/fsetti/milliQan/OutputFiles/021219/Data_ratio10.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( chD , output_file_nameData );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

//  compareHists_v2 ( output_file_nameData, output_file_name, "h_RMS" );

  TString outProjData, outProjMC;
  outProjData = "/home/users/fsetti/milliQan/OutputFiles/021219/Data_ratio10_Projections.root";
  outProjMC   = "/home/users/fsetti/milliQan/OutputFiles/021219/Sim_ratio10_Projections.root" ;

//  Projections( output_file_nameData , "nPEdT", outProjData );
//  Projections( output_file_name     , "nPEdT", outProjMC   );
//  Quadrature ( outProjData,outProjMC, "sigmas", "/home/users/fsetti/milliQan/OutputFiles/021219/QuadratureHist.root" );

}
