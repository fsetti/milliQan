// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.h"

 
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

bool sort_wrt_calTime( vector<float> vec1, vector<float> vec2 ){
        return vec1[4] < vec2[4];
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



vector<vector<vector<float>>> mCPulses( vector<int> *pulseChan, int pulseTrueNPE[32], float pulseMuDist[32], vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
        vector<vector<vector<float>>> outputPulses;
        for (unsigned int i=0; i<TrigChannels.size(); i++){
                vector<vector<float>> preSelPulses;
                for (unsigned int j=0; j<pulseChan->size(); j++){
                        for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
                                if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
                                        ){
                                        preSelPulses.push_back({(float)pulseChan->at(j), (float)pulseTrueNPE[pulseChan->at(j)], pulseMuDist[pulseChan->at(j)], pulsenPE->at(j), pulseTimeCal->at(j)});
                                        break;
                                }
                        }
                }

//	check one good hit per channel, keep earliest pulses
		std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_calTime );
		for (unsigned int j=0; j<preSelPulses.size(); j++){
			for (unsigned int k=j+1; k<preSelPulses.size(); k++){
				if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
	                                preSelPulses.at(j).at(3)+= preSelPulses.at(k).at(3);
					preSelPulses.erase( preSelPulses.begin() + k );
					k -=1;
				}
			}
		}

//	check at least one mCP hit per trig configuration
		bool mCPhit = false;
		for (unsigned int j=0; j<preSelPulses.size(); j++){
			if ( preSelPulses.at(j).at(2) > 0 ){ 
				mCPhit = true; 
				preSelPulses.erase( preSelPulses.begin() + j );
				j -=1;
			}
		}
		if ( !mCPhit ) continue;
		else outputPulses.push_back( preSelPulses );
        }
        return outputPulses;
}




