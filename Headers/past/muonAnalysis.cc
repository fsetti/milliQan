// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/past/muonAnalysis.h"

 
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



vector<vector<vector<float>>> sameLayerSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
	vector<vector<vector<float>>> outputPulses;
	for (unsigned int i=0; i<TrigChannels.size(); i++){
		vector<int> trigChStatus(TrigChannels.at(i).size(),0);  // set no pulses in tag channel at the beginning of every event
		vector<vector<float>> preSelPulses;
		for (unsigned int j=0; j<pulseChan->size(); j++){
	        	for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
		                if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
	                 	 && pulsenPE->at(j)  >  nPEminBarS
		                 && pulseTime->at(j)  <  tMax
		                 && pulseTime->at(j) >  tMin
		                        ){
		                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
		                        trigChStatus.at(k) = 1;
		                        break;
		                }
	       		}
	     	}
	//	check 1 hit in all channels of TrigChannels
		unsigned int nHits = std::accumulate( trigChStatus.begin(), trigChStatus.end(), 0 );
		if ( nHits != TrigChannels.at(i).size() ){ 
			continue;
		}

//	only retain largest pulse in each channel, require largest pulse_nPE > nPEmin
		std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
		for (unsigned int j=0; j<preSelPulses.size(); j++){
			for (unsigned int k=j+1; k<preSelPulses.size(); k++){
				if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
					preSelPulses.erase( preSelPulses.begin() + k );
					k -=1;
				}
			}
		}

		std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_nPE );
//		require minimum nPE on muon pulse
		if ( preSelPulses.at(0).at(0) != 5
		  && preSelPulses.at(0).at(0) != 22
		  && preSelPulses.at(0).at(3) < nPEminBar ){
			continue;
		}
		if ( ( preSelPulses.at(0).at(0) == 5 || preSelPulses.at(0).at(0) == 22 )
		  && preSelPulses.at(0).at(3) < nPEminBarR7725 ){
			continue;
		}

//		require maximum nPE on secondary pulse
		if ( preSelPulses.at(1).at(0) != 5
		  && preSelPulses.at(1).at(0) != 22
		  && preSelPulses.at(1).at(3) > nPEminBar ){
			continue;
		}
		if ( ( preSelPulses.at(1).at(0) == 5 || preSelPulses.at(1).at(0) == 22 )
		  && preSelPulses.at(1).at(3) > nPEminBarR7725 ){
			continue;
		}

		if ( preSelPulses.size() == TrigChannels.at(i).size()  ){ outputPulses.push_back( preSelPulses );  }
	}
	return outputPulses;
}




vector<vector<float>> sameLayerSelectionNoHit( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<int> TrigChannels ){
	vector<vector<float>> preSelPulses;
	for (unsigned int j=0; j<pulseChan->size(); j++){
        	for (unsigned int k=0; k<TrigChannels.size(); k++){
	                if( pulseChan->at(j)   == TrigChannels.at(k)
	                 && pulsenPE->at(j)  >  nPEminBarS
	                 && pulseTime->at(j)  <  tMax
	                 && pulseTime->at(j) >  tMin
	                        ){
		                preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
		                break;
			}
       		}
	}

//	only retain earliest pulse in each channel
	std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
	for (unsigned int j=0; j<preSelPulses.size(); j++){
		for (unsigned int k=j+1; k<preSelPulses.size(); k++){
			if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
				preSelPulses.erase( preSelPulses.begin() + k );
				k -=1;
			}
		}
	}

	for (unsigned int i=0; i<preSelPulses.size(); i++){ 
		if ( preSelPulses.at(i).at(0) != 5
		  && preSelPulses.at(i).at(0) != 22
		  && preSelPulses.at(i).at(3) > nPEminBar ){
			preSelPulses.erase( preSelPulses.begin() + i );
			i -=1;
		}
	}
	for (unsigned int i=0; i<preSelPulses.size(); i++){ 
		if ( ( preSelPulses.at(i).at(0) == 5 || preSelPulses.at(i).at(0) == 22 )
		  && preSelPulses.at(i).at(3) > nPEminBarR7725 ){
			preSelPulses.erase( preSelPulses.begin() + i );
			i -=1;
		}
	}

	return preSelPulses;
}



vector<vector<vector<float>>> sameLayerSelectionNoHit_v2( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
	vector<vector<vector<float>>> outputPulses;
	for (unsigned int i=0; i<TrigChannels.size(); i++){
		vector<int> trigChStatus(TrigChannels.at(i).size(),0);  // set no pulses in tag channel at the beginning of every event
		vector<vector<float>> preSelPulses;
		for (unsigned int j=0; j<pulseChan->size(); j++){
	        	for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
		                if( pulseChan->at(j)  == TrigChannels.at(i).at(k)
	                 	 && pulsenPE->at(j)   >  nPEminBarS
//	                 	 && pulsenPE->at(j)   <  maxSPE
	                 	 && pulseTime->at(j)  <  tMax
	                 	 && pulseTime->at(j)  >  tMin
		                        ){
		                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
		                        trigChStatus.at(k) = 1;
		                        break;
		                }
	       		}
	     	}
	//	check 1 hit in all channels of TrigChannels
		unsigned int nHits = std::accumulate( trigChStatus.begin(), trigChStatus.end(), 0 );
		if ( nHits != TrigChannels.at(i).size() ){ 
			continue;
		}

//	only retain largest pulse in each channel, require largest pulse_nPE > nPEmin
		std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
		for (unsigned int j=0; j<preSelPulses.size(); j++){
			for (unsigned int k=j+1; k<preSelPulses.size(); k++){
				if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
					preSelPulses.erase( preSelPulses.begin() + k );
					k -=1;
				}
			}
		}

//		make sure first pulse in channel is not above nPE max threshold
		bool acceptPulse = true;
		for (unsigned int j=0; j<preSelPulses.size(); j++){
                 	 if ( preSelPulses.at(j).at(3)   >  maxSPE ){ 
			 	acceptPulse = false;
				continue;
			 }
		}
		if ( !acceptPulse ) continue;

		std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_channel );
		if ( preSelPulses.size() != 0 ) outputPulses.push_back( preSelPulses ); 
	}
	return outputPulses;
}




vector<vector<vector<float>>> pulseBarSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
	vector<vector<vector<float>>> outputPulses;
	for (unsigned int i=0; i<TrigChannels.size(); i++){
		vector<int> trigChStatus(TrigChannels.at(i).size(),0);  // set no pulses in tag channel at the beginning of every event
		vector<vector<float>> preSelPulses;
		for (unsigned int j=0; j<pulseChan->size(); j++){
	        	for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
		                if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
		                 && pulseTime->at(j)  <  tMax
		                 && pulseTime->at(j)  >  tMin
		                        ){
		                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
		                        trigChStatus[k] = 1;
		                        break;
		                }
	       		}
	     	}
	//	check 1 hit in all channels of TrigChannels
		unsigned int nHits = std::accumulate( trigChStatus.begin(), trigChStatus.end(), 0 );
		if ( nHits != TrigChannels.at(i).size() ){ 
			continue;
		}

//		check one good hit per channel, keep earliest pulses
		std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
		for (unsigned int j=0; j<preSelPulses.size(); j++){
			for (unsigned int k=j+1; k<preSelPulses.size(); k++){
				if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
					preSelPulses.erase( preSelPulses.begin() + k );
					k -=1;
				}
			}
		}
		for (unsigned int j=0; j<preSelPulses.size(); j++){
	         	if ( preSelPulses.at(j).at(0) != 5
			  && preSelPulses.at(j).at(0) != 22  
			  && preSelPulses.at(j).at(3)   <  nPEminBar ){
				preSelPulses.erase( preSelPulses.begin() + j );
				j -=1;
			}
		}
		for (unsigned int j=0; j<preSelPulses.size(); j++){
	         	if ( ( preSelPulses.at(j).at(0) == 5 ||  preSelPulses.at(j).at(0) == 22 )
			  && preSelPulses.at(j).at(3)   <  nPEminBarR7725 ){
				preSelPulses.erase( preSelPulses.begin() + j );
				j -=1;
			}
		}

		if ( preSelPulses.size() == TrigChannels.at(i).size() ) outputPulses.push_back( preSelPulses );
	}
	return outputPulses;
}



vector<vector<float>> pulseBarSelection_v2( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<int> TrigChannels ){

	vector<vector<float>> preSelPulses;
	for (unsigned int j=0; j<pulseChan->size(); j++){
        	for (unsigned int k=0; k<TrigChannels.size(); k++){
	                if( pulseChan->at(j)   == TrigChannels.at(k)
	                 && pulseTime->at(j)  <  tMax
	                 && pulseTime->at(j)  >  tMin
	                        ){
		                preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
				break;
	                }
       		}
     	}

	std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
	for (unsigned int j=0; j<preSelPulses.size(); j++){
		for (unsigned int k=j+1; k<preSelPulses.size(); k++){
			if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
				preSelPulses.erase( preSelPulses.begin() + k );
				k -=1;
			}
		}
	}

	for (unsigned int j=0; j<preSelPulses.size(); j++){
         	if ( preSelPulses.at(j).at(0) != 5
		  && preSelPulses.at(j).at(0) != 22  
		  && preSelPulses.at(j).at(3)   <  nPEminBar ){
			preSelPulses.erase( preSelPulses.begin() + j );
			j -=1;
		}
	}
	for (unsigned int j=0; j<preSelPulses.size(); j++){
         	if ( ( preSelPulses.at(j).at(0) == 5 ||  preSelPulses.at(j).at(0) == 22 )
		  && preSelPulses.at(j).at(3)   <  nPEminBarR7725 ){
			preSelPulses.erase( preSelPulses.begin() + j );
			j -=1;
		}
	}

	return preSelPulses;
}



vector<vector<vector<float>>> pulseSlabSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
	vector<vector<vector<float>>> outputPulses;
	for (unsigned int i=0; i<TrigChannels.size(); i++){
		vector<int> trigChStatus(TrigChannels.at(i).size(),0);  // set no pulses in tag channel at the beginning of every event
		vector<vector<float>> preSelPulses;
		for (unsigned int j=0; j<pulseChan->size(); j++){
	        	for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
		                if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
		                 && pulseTime->at(j)  <  tMax
		                 && pulseTime->at(j) >  tMin
		                 && pulsenPE->at(j)   >  nPEminSlab
		                        ){
		                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
		                        trigChStatus[k] = 1;
		                        break;
		                }
	       		}
	     	}
	//	check 1 hit in all channels of TrigChannels
		unsigned int nHits = std::accumulate( trigChStatus.begin(), trigChStatus.end(), 0 );
		if ( nHits != TrigChannels.at(i).size() ){ 
			continue;
		}

//		check one good hit per channel, keep earliest pulses
		std::sort( preSelPulses.begin(), preSelPulses.end(), sort_wrt_time );
		for (unsigned int j=0; j<preSelPulses.size(); j++){
			for (unsigned int k=j+1; k<preSelPulses.size(); k++){
				if ( preSelPulses.at(j)[0] == preSelPulses.at(k)[0] ){
					preSelPulses.erase( preSelPulses.begin() + k );
					k -=1;
				}
			}
		}
		outputPulses.push_back( preSelPulses );
	}
	return outputPulses;
}



