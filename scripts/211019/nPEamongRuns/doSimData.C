{
  gROOT->ProcessLine(".L selection.C+");
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C+");

	
  TChain *ch = new TChain("t");
  TString fileName = "/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root";
  ch->Add(fileName);
  TString output_file_name = "/home/users/fsetti/milliQan/OutputFiles/211019/nPErunComparison.root";
  
  cout << "-------------------------" << endl; 
  selection( ch , output_file_name );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

  plotHist( output_file_name, "neigh_nPEs",     "/home/users/fsetti/public_html/milliQan/lowNPErunComparison/",false, true );
  plotHist( output_file_name, "Nunscaled", "/home/users/fsetti/public_html/milliQan/lowNPErunComparison/" );
  plotHist( output_file_name, "nonNeigh_nPEs",     "/home/users/fsetti/public_html/milliQan/lowNPErunComparison/",false, true );
  plotHist( output_file_name, "nNunscaled", "/home/users/fsetti/public_html/milliQan/lowNPErunComparison/" );
}
