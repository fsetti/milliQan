{
  gROOT->ProcessLine(".L selection.C+");
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/MuonRate/plotHistsMuon.C++");

// Data
  TString output_file_name;
  TChain *ch = new TChain("t");
  ch->Add("/nfs-7/userdata/fsetti/allPhysicsAndTripleChannelSinceTS1_slabSkim_181017.root");
  output_file_name = "/home/users/fsetti/milliQan/OutputFiles/Data_slabTime.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of all Runs" << endl; 
//  selection( ch , output_file_name );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

// MC
  TString output_file_nameMC;
  TChain *chMC = new TChain("t");
  chMC->Add("/nfs-7/userdata/fsetti/sim_run3_*.root");
  output_file_nameMC = "/home/users/fsetti/milliQan/OutputFiles/MC_slabTime.root";

  cout << "-------------------------" << endl; 
  cout << "Skim of MC" << endl; 
//  selection( chMC , output_file_nameMC, false );	
  cout << "-------------------------" << endl; 
  cout << "\n \n \n " << endl;

  plotHist( output_file_name, "dT_slabs0_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs1_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs2_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabsNoBar_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 

  plotHist( output_file_name, "dT_slabs0_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs1_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs2_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabsNoBar_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 

  plotHist( output_file_name, "dT_slabs0_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs1_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs2_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabsNoBar_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 

  plotHist( output_file_name, "dT_slabs0_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs1_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabs2_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 
  plotHist( output_file_name, "dT_slabsNoBar_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/Data_" ); 


  plotHist( output_file_nameMC, "dT_slabs0_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs1_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs2_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabsNoBar_Ch18", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 

  plotHist( output_file_nameMC, "dT_slabs0_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs1_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs2_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabsNoBar_Ch21", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 

  plotHist( output_file_nameMC, "dT_slabs0_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs1_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs2_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabsNoBar_Ch28", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 

  plotHist( output_file_nameMC, "dT_slabs0_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs1_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabs2_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
  plotHist( output_file_nameMC, "dT_slabsNoBar_Ch20", "/home/users/fsetti/public_html/milliQan/MuonRate/14Oct2019/MC_" ); 
}
