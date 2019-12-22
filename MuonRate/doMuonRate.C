{
  gROOT->ProcessLine(".L /homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/muonSelection.C++");
  gROOT->ProcessLine(".L /homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/plotHistsMuon.C+");


  std::ifstream file("/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/listRunsMuonRate.txt");
  std::string line;
  while (std::getline(file, line) ){
         int iter=0;
         vector<TString> runVec;
         vector<float> trigList;
         vector<TString> chList;
         std::stringstream ss(line);
         while( ss.good() ){
             	string subs;
              	ss >> subs; 
	      	if ( subs == "#" || subs == "//" ) break;
	      	if ( subs == "" || subs == " " ) continue;
              	if (iter<1){ runVec.push_back(subs) ; }
              	if (iter>=1 ){ 
			trigList.push_back(std::stof(subs));
			chList.push_back(subs);
		}
	      	iter++;
        }
	if ( runVec.size() == 0 ) continue;

	TString runName = runVec[0];
	cout <<"Processing runs: " << runName << ", with trigger set on channels ";
        for (unsigned int i=0; i<trigList.size(); i++){ cout << trigList.at(i) << " " ;}
	cout << endl;

	bool Synch = true;
  	TChain *ch = new TChain("t");
	TString path = "/net/cms26/cms26r0/milliqan/milliqanOffline/trees/Run"+runName;
	read_list_files( path, ch );
	TString output_file_name = "/homes/fsetti/CMSSW_8_1_0/src/milliQan/MuonRate/MuonFiles/Run" + runName + "_muonStudy.root";

	cout << "\n \n \n RUN " << runName << endl;
	cout << "-------------------------" << endl; 
  	muonSelection( ch , output_file_name, runName, trigList, Synch);	
	cout << "-------------------------" << endl; 
	cout << "\n \n \n " << endl;
 
	TString hName1 = "Run"+runName+"_CH"+chList[0]+"_Height";
	TString hName2 = "Run"+runName+"_CH"+chList[1]+"_Height";
	TString hName3 = "Run"+runName+"_CH"+chList[2]+"_Height";
	TString hName4 = "Run"+runName+"_CH"+chList[3]+"_Height";
//	comparePlots( output_file_name, hName1, hName2, hName3, hName4 );

	TString nName1 = "Run"+runName+"_CH"+chList[0]+"_nPE";
	TString nName2 = "Run"+runName+"_CH"+chList[1]+"_nPE";
	TString nName3 = "Run"+runName+"_CH"+chList[2]+"_nPE";
	TString nName4 = "Run"+runName+"_CH"+chList[3]+"_nPE";
//	comparePlots( output_file_name, nName1, nName2, nName3, nName4 );

	TString tName1 = "Run"+runName+"_CH"+chList[0]+"_Time";
	TString tName2 = "Run"+runName+"_CH"+chList[1]+"_Time";
	TString tName3 = "Run"+runName+"_CH"+chList[2]+"_Time";
	TString tName4 = "Run"+runName+"_CH"+chList[3]+"_Time";
//	comparePlots( output_file_name, tName1, tName2, tName3, tName4 );

	TString ctName1 = "Run"+runName+"_CH"+chList[0]+"_CalTime";
	TString ctName2 = "Run"+runName+"_CH"+chList[1]+"_CalTime";
	TString ctName3 = "Run"+runName+"_CH"+chList[2]+"_CalTime";
	TString ctName4 = "Run"+runName+"_CH"+chList[3]+"_CalTime";
//	comparePlots( output_file_name, ctName1, ctName2, ctName3, ctName4 );
  }
}
