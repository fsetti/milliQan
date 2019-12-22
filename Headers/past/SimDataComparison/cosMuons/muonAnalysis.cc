// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/SimDataComparison/cosMuons/muonAnalysis.h"

 
bool sort_wrt_channel( vector<float> vec1, vector<float> vec2 ){
        return vec1[0] < vec2[0];
}

bool sort_wrt_height( vector<float> vec1, vector<float> vec2 ){
        return vec1[1] > vec2[1];
}

bool sort_wrt_time( vector<float> vec1, vector<float> vec2 ){
        return vec1[2] < vec2[2];
}

bool sort_wrt_nPE( vector<float> vec1, vector<float> vec2 ){
        return vec1[3] > vec2[3];
}

void read_list_files(const char *dirname, TChain *chainptr )
{  
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files) {
       TSystemFile *file;
       TString fname;
       TIter next(files);
       while ((file=(TSystemFile*)next())) {
          fname = file->GetName();
          if ( fname.Length()>2 ) {
               TString path;
               path = dirname + TString("/") + fname ;
               chainptr->Add(path);
          }
        }
    }
}


vector<vector<vector<float>>> CosmicSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
	vector<vector<vector<float>>> outputPulses;
	for (unsigned int i=0; i<TrigChannels.size(); i++){
	        vector<vector<float>> preSelPulses;
	        for (unsigned int j=0; j<pulseChan->size(); j++){
	                for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
	                        if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
	                         && pulsenPE->at(j)  >  nPEminVec.at(pulseChan->at(j))
	                         && pulseTime->at(j)  <  tMax
	                         && pulseTime->at(j) >  tMin
	                                ){
	                                preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
	                        }
	                }
	        }
	
	        std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
	        for (unsigned int j=0; j<preSelPulses.size(); j++){
	                for (unsigned int k=j+1; k<preSelPulses.size(); k++){
	                        if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
	//                                preSelPulses.at(j).at(3)+= preSelPulses.at(k).at(3);
	                                preSelPulses.erase( preSelPulses.begin() + k );
	                                k -=1;
	                        }
	                }
	        }
		if ( preSelPulses.size() == TrigChannels.at(i).size() ) outputPulses.push_back( preSelPulses );
	}
        return outputPulses;
}


vector<vector<vector<float>>> SameLayerActivity ( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<vector<float>>> cosmicPulses ){
	vector<vector<vector<float>>> outputPulses;
	for (unsigned int i=0; i<cosmicPulses.size(); i++){
		vector<int> TrigChannels;
		std::sort( cosmicPulses.at(i).begin(), cosmicPulses.at(i).end(), sort_wrt_channel );
		if ( cosmicPulses.at(i).at(0).at(0) == 0 ) TrigChannels = trigL1right;
		if ( cosmicPulses.at(i).at(0).at(0) == 1 ) TrigChannels = trigL1left;
		if ( cosmicPulses.at(i).at(0).at(0) == 6 ) TrigChannels = trigL2right;
		if ( cosmicPulses.at(i).at(0).at(0) == 7 ) TrigChannels = trigL2left;
		if ( cosmicPulses.at(i).at(0).at(0) == 2 ) TrigChannels = trigL3right;
		if ( cosmicPulses.at(i).at(0).at(0) == 3 ) TrigChannels = trigL3left;
	        vector<vector<float>> preSelPulses;
	        for (unsigned int j=0; j<pulseChan->size(); j++){
	                for (unsigned int k=0; k<TrigChannels.size(); k++){
	                        if( pulseChan->at(j)   == TrigChannels.at(k)
	                         && pulsenPE->at(j)  >  nPEmin
	                         && pulseTime->at(j)  <  tMax
	                         && pulseTime->at(j) >  tMin
	                                ){
	                                preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
	                        }
	                }
	        }
		if ( preSelPulses.size() == 0 ) continue;

	        std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
	        for (unsigned int j=0; j<preSelPulses.size(); j++){
	                for (unsigned int k=j+1; k<preSelPulses.size(); k++){
	                        if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
	                                preSelPulses.at(j).at(3)+= preSelPulses.at(k).at(3);
	                                preSelPulses.erase( preSelPulses.begin() + k );
	                                k -=1;
	                        }
	                }
	        }
		if ( preSelPulses.size() != 0 ) outputPulses.push_back( preSelPulses );
	}
	return outputPulses;
}

vector<vector<float>> SlabSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal ){
	vector<vector<float>> nul;
	vector<int> TrigChannels = {18,20,21,28 }; //	set trig channels for slab selection
	vector<vector<float>> preSelPulses;
	for (unsigned int j=0; j<pulseChan->size(); j++){
        	for (unsigned int k=0; k<TrigChannels.size(); k++){
	                if( pulseChan->at(j)   == TrigChannels.at(k)
	                 && pulseTime->at(j)  <  tMax
	                 && pulseTime->at(j) >  tMin
	                 && pulsenPE->at(j)   >  nPEminSlab
	                        ){
	                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
	                        break;
	                }
       		}
     	}

//	check one good hit per channel, keep earliest pulses
	std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
	for (unsigned int j=0; j<preSelPulses.size(); j++){
		for (unsigned int k=j+1; k<preSelPulses.size(); k++){
			if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
                                preSelPulses.at(j).at(3)+= preSelPulses.at(k).at(3);
				preSelPulses.erase( preSelPulses.begin() + k );
				k -=1;
			}
		}
	}
	if ( preSelPulses.size() == TrigChannels.size() ) return preSelPulses; 
	else return nul; 
}



