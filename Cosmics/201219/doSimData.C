{
  gROOT->ProcessLine(".L selection.C+");

  TString output_file_nameData;
  TChain *chD = new TChain("t");
  chD->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_cosmicSkim_191210.root");
  output_file_nameData = "/home/users/fsetti/milliQan/OutputFiles/201219/Cosmics_ratio10.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
  selection( chD , output_file_nameData );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

//  compareHists_v2 ( output_file_nameData, output_file_name, "h_RMS" );

  TString outProjData, outProjMC;
  outProjData = "/home/users/fsetti/milliQan/OutputFiles/201219/Cosmics_ratio10_Projections.root";

//  Projections( output_file_nameData , "nPEdT", outProjData );
//  Quadrature ( outProjData,outProjMC, "sigmas", "/home/users/fsetti/milliQan/OutputFiles/021219/QuadratureHist.root" );

}
