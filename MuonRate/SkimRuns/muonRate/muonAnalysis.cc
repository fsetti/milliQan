// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "muonAnalysis.h"

 
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





bool MuonInBars (  vector<int> *pulseChan, vector<float> *pulseTime, vector<float> *pulsenPE,  vector<int> TrigChannels ){
	bool muonInBars = false;
	for (unsigned int j=0; j<pulseChan->size(); j++){
		for (unsigned int i=0; i<TrigChannels.size(); i++){
			if( pulseChan->at(j)   == TrigChannels.at(i)
		    	&& pulseTime->at(j)  <  tMax
     	            	&& pulseTime->at(j) >  tMin
     	            	){
				if ( pulseChan->at(j) != 5
     	                	  && pulseChan->at(j) != 22
     	                	  && pulsenPE->at(j)  >  nPEminMuon ){
     	   				muonInBars = true;
     	   				break; 
     	                   	}
     	               		if ( ( pulseChan->at(j) == 5 || pulseChan->at(j) == 22 )
     	                	    && pulsenPE->at(j)  >  nPEminMuonR7725  ){
     	   				muonInBars = true;
     	   				break; 
     	                	}
     	   		}
		}
     	}
	return muonInBars;
}




vector<vector<vector<float>>> SameLayerSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels, bool Neighbouring = true  ){
        vector<vector<vector<float>>> outputPulses;
        for (unsigned int i=0; i<TrigChannels.size(); i++){
                vector<vector<float>> preSelPulses;
		
                for (unsigned int j=0; j<pulseChan->size(); j++){
                        for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
                                if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
                                 && pulsenPE->at(j)  >  nPEmin
                                 && pulseTime->at(j) <  tMax
                                 && pulseTime->at(j) >  tMin
                                        ){
                                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
                                        break;
                                }
                        }
                }

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
		if ( preSelPulses.size() < TrigChannels.at(i).size() ) continue;

                std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_nPE );
                if ( preSelPulses.at(0).at(0) != 5
                  && preSelPulses.at(0).at(0) != 22
                  && preSelPulses.at(0).at(3) < nPEminMuon ){
                        continue;
                }
                if ( ( preSelPulses.at(0).at(0) == 5 || preSelPulses.at(0).at(0) == 22 )
                  && preSelPulses.at(0).at(3) < nPEminMuonR7725 ){
                        continue;
                }

                if ( preSelPulses.at(1).at(0) != 5
                  && preSelPulses.at(1).at(0) != 22
                  && preSelPulses.at(1).at(3) > nPEminMuon ){
                        continue;
                }
                if ( ( preSelPulses.at(1).at(0) == 5 || preSelPulses.at(1).at(0) == 22 )
                  && preSelPulses.at(1).at(3) > nPEminMuonR7725 ){
                        continue;
                }
//	make sure no muons in bar close to secondary pulse
		if ( !Neighbouring ){
			bool nearbyMuon = false ; 
			if ( preSelPulses.at(1).at(0) == 0  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {1,24,25} );
			if ( preSelPulses.at(1).at(0) == 1  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {0,24,25} );
			if ( preSelPulses.at(1).at(0) == 8  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {9,24,25} );
			if ( preSelPulses.at(1).at(0) == 9  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {8,24,25} );
			if ( preSelPulses.at(1).at(0) == 6  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {7,16,17} );
			if ( preSelPulses.at(1).at(0) == 7  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {6,16,17} );
			if ( preSelPulses.at(1).at(0) == 12 )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE,{13,16,17} );
			if ( preSelPulses.at(1).at(0) == 13 )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE,{12,16,17} );
			if ( preSelPulses.at(1).at(0) == 2  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {3,22,23} );
			if ( preSelPulses.at(1).at(0) == 3  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {2,22,23} );
			if ( preSelPulses.at(1).at(0) == 4  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {5,22,23} );
			if ( preSelPulses.at(1).at(0) == 5  )  nearbyMuon =  MuonInBars ( pulseChan, pulseTime, pulsenPE, {4,22,23} );
	
			if ( nearbyMuon ) continue;
		}

                if ( preSelPulses.size() == TrigChannels.at(i).size()  ) outputPulses.push_back( preSelPulses );  
        }
        return outputPulses;
}





vector<vector<float>> selectPulses( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<int> TrigChannels ){
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

        std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_calTime );
//        for (unsigned int j=0; j<preSelPulses.size(); j++){
//                for (unsigned int k=j+1; k<preSelPulses.size(); k++){
//                        if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
//                                preSelPulses.at(j).at(3)+= preSelPulses.at(k).at(3);
//                                preSelPulses.erase( preSelPulses.begin() + k );
//                                k -=1;
//                        }
//                }
//        }
        return preSelPulses;
}


//same as above but w/o nPE threshold
vector<vector<float>> selectPulses_v2( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<int> TrigChannels ){
        vector<vector<float>> preSelPulses;
        for (unsigned int j=0; j<pulseChan->size(); j++){
                for (unsigned int k=0; k<TrigChannels.size(); k++){
                        if( pulseChan->at(j)   == TrigChannels.at(k)
                         && pulseTime->at(j)  <  tMax
                         && pulseTime->at(j) >  tMin
                                ){
                                preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
                        }
                }
        }

        std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_calTime );
//        for (unsigned int j=0; j<preSelPulses.size(); j++){
//                for (unsigned int k=j+1; k<preSelPulses.size(); k++){
//                        if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
//                                preSelPulses.at(j).at(3)+= preSelPulses.at(k).at(3);
//                                preSelPulses.erase( preSelPulses.begin() + k );
//                                k -=1;
//                        }
//                }
//        }
        return preSelPulses;
}






vector<vector<vector<float>>> pulseBarSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
        vector<vector<vector<float>>> outputPulses;
        for (unsigned int i=0; i<TrigChannels.size(); i++){
                vector<vector<float>> preSelPulses;
                for (unsigned int j=0; j<pulseChan->size(); j++){
                        for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
                                if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
                                 && pulseTime->at(j)  <  tMax
                                 && pulseTime->at(j)  >  tMin
                                        ){
                                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
                                        break;
                                }
                        }
                }

// only retain earliest pulse
                std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_calTime );
                for (unsigned int j=0; j<preSelPulses.size(); j++){
                        for (unsigned int k=j+1; k<preSelPulses.size(); k++){
                                if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
                                        preSelPulses.erase( preSelPulses.begin() + k );
                                        k -=1;
                                }
                        }
                }
		if ( preSelPulses.size() != TrigChannels.at(i).size() ) continue;


                for (unsigned int j=0; j<preSelPulses.size(); j++){
                        if ( preSelPulses.at(j).at(0) != 5
                          && preSelPulses.at(j).at(0) != 22
                          && preSelPulses.at(j).at(3)   <  nPEminMuon ){
                                preSelPulses.erase( preSelPulses.begin() + j );
                                j -=1;
                        }
                }
                for (unsigned int j=0; j<preSelPulses.size(); j++){
                        if ( ( preSelPulses.at(j).at(0) == 5 ||  preSelPulses.at(j).at(0) == 22 )
                          && preSelPulses.at(j).at(3)   <  nPEminMuonR7725 ){
                                preSelPulses.erase( preSelPulses.begin() + j );
                                j -=1;
                        }
                }

                if ( preSelPulses.size() == TrigChannels.at(i).size() ) outputPulses.push_back( preSelPulses );
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
	if ( preSelPulses.size() == TrigChannels.size() ) return preSelPulses; 
	else return nul; 
}



