{
  gROOT->ProcessLine(".L selection.C+");
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/Headers/MCplotting.C+");

  std::ifstream file("MCfiles.txt");
  std::string line;
  vector<TString> plotInfo;
  while (std::getline(file, line) ){
         int iter=0;
         std::stringstream ss(line);
	TString output_file;
	TString input_file;
         while( ss.good() ){
                string subs;
                ss >> subs;
                if ( subs == "#" || subs == "//" ) break;
                if ( subs == ""  ) continue;
                if (iter==0) {  input_file = subs ;}
                if (iter!=0 ){ output_file = subs ;}
                iter++;
        }
	if ( iter < 1 ) continue;
	
	TChain *ch = new TChain("t");
	TString fileName = "/hadoop/cms/store/user/bemarsh/milliqan/milliq_mcgen/mixed_trees/mcp_v7_v1_save2m_skim0p25m_simmcp_v1_v1/"+input_file+"/*.root";
	ch->Add(fileName);
	TString output_file_name = "/home/users/fsetti/milliQan/OutputFiles/14Oct2019/MC_"+output_file+".root";
	
	cout << "-------------------------" << endl; 
	cout << "Analysing MC for mCP:  " << output_file << endl; 
	selection( ch , output_file_name );	
	cout << "-------------------------" << endl; 
	cout << "\n \n \n " << endl;

	plotInfo.push_back( output_file_name );
   }

   compareHists( plotInfo.at(0),   plotInfo.at(1),  plotInfo.at(2),  plotInfo.at(3),  plotInfo.at(4), "nPE_adjBars" , false, true);
   compareHists( plotInfo.at(0),   plotInfo.at(1),  plotInfo.at(2),  plotInfo.at(3),  plotInfo.at(4), "nPE_non_adjBars", false, true );
   compareHists( plotInfo.at(0),   plotInfo.at(1),  plotInfo.at(2),  plotInfo.at(3),  plotInfo.at(4), "nPE_Panels" );
//   compareHists_v2( plotInfo.at(0),   plotInfo.at(1),  plotInfo.at(2),  plotInfo.at(3),  plotInfo.at(4), "trueMC_nPE_bars" );
//   compareHists( plotInfo.at(0),   plotInfo.at(1),  plotInfo.at(2),  plotInfo.at(3),  plotInfo.at(4), "dT_slabs_Bar" );
//   compareHists( plotInfo.at(0),   plotInfo.at(1),  plotInfo.at(2),  plotInfo.at(3),  plotInfo.at(4), "dT_slab18" );
//   compareHists( plotInfo.at(0),   plotInfo.at(1),  plotInfo.at(2),  plotInfo.at(3),  plotInfo.at(4), "dT_slab21" );
}
