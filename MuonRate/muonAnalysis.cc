// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/homes/fsetti/CMSSW_8_1_0/src/milliQan/Headers/muonAnalysis.h"

 
bool sort_wrt_channel( vector<float> vec1, vector<float> vec2 ){
        return vec1[0] < vec2[0];
}

bool sort_wrt_height( vector<float> vec1, vector<float> vec2 ){
        return vec1[1] > vec2[1];
}

bool sort_wrt_time( vector<float> vec1, vector<float> vec2 ){
        return vec1[2] < vec2[2];
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


vector<vector<vector<float>>> pulsePreSelection( vector<int> *pulseChan, vector<float> *pulseHeight, vector<float> *pulseTime, vector<float> *pulsenPE, vector<float> *pulseTimeCal, vector<vector<int>> TrigChannels ){
	vector<vector<vector<float>>> outputPulses;
	for (unsigned int i=0; i<TrigChannels.size(); i++){
		vector<int> trigChStatus(TrigChannels.at(i).size(),0);  // set no pulses in tag channel at the beginning of every event
		vector<vector<float>> preSelPulses;
		for (unsigned int j=0; j<pulseChan->size(); j++){
	        	for (unsigned int k=0; k<TrigChannels.at(i).size(); k++){
		                if( pulseChan->at(j)   == TrigChannels.at(i).at(k)
		                 && pulseTime->at(j)  <  tMax
		                 && pulseTime->at(j) >  tMin
		                 && pulsenPE->at(j)   >  nPEmin
		                        ){
		                        preSelPulses.push_back({(float)pulseChan->at(j), pulseHeight->at(j), pulseTime->at(j), pulsenPE->at(j), pulseTimeCal->at(j)});
		                        trigChStatus[k] = 1;
		                        break;
		                }
	       		}
	     	}
	//	check 1 hit in all channels of TrigChannels
		unsigned int nHits = std::accumulate( trigChStatus.begin(), trigChStatus.end(), 0 );
		if ( nHits == TrigChannels.at(i).size() ) outputPulses.push_back( preSelPulses );
		if ( nHits != TrigChannels.at(i).size() ) outputPulses.push_back( {{0.}} );
	}
	return outputPulses;
}




vector<vector<float>> pulseSelection( vector<vector<float>> Pulses, unsigned int nTriggerChannels ){ 
    vector<vector<float>> outputPulses;
    std::sort( Pulses.begin(), Pulses.end(), sort_wrt_time );
    for (unsigned int i=0; i<Pulses.size(); i++){
//	require 3 hits within 100ns in 3 different channels	
	for (unsigned int j=i+1; j<Pulses.size(); j++){
		if ( Pulses.at(i)[0] == Pulses.at(j)[0]
//		  || fabs(Pulses.at(i)[2]-Pulses.at(j)[2])>dtMax   
//		  || fabs(Pulses.at(i)[4]-Pulses.at(j)[4])>dtMin   
		){ 
			continue;
		}
		if ( nTriggerChannels == 3 ){
			for (unsigned int k=j+1; k<Pulses.size(); k++){
				if ( Pulses.at(i)[0] == Pulses.at(k)[0] || Pulses.at(j)[0] == Pulses.at(k)[0]
//				  || fabs(Pulses.at(i)[2]-Pulses.at(k)[2])>dtMax || fabs(Pulses.at(j)[2]-Pulses.at(k)[2])>dtMax 
//				  || fabs(Pulses.at(i)[4]-Pulses.at(k)[4])>dtMin || fabs(Pulses.at(j)[4]-Pulses.at(k)[4])>dtMin 
				){ 
					continue;
				}
				if ( outputPulses.size() == 0 ){
					outputPulses.push_back( Pulses.at(i) );
					outputPulses.push_back( Pulses.at(j) );
					outputPulses.push_back( Pulses.at(k) );
					break;
				}
			}
			if ( outputPulses.size() != 0 ) break;
		}
//	require 4 hits within 100ns in 4 different channels	
		if ( nTriggerChannels == 4 ){
			for (unsigned int k=j+1; k<Pulses.size(); k++){
				if ( Pulses.at(i)[0] == Pulses.at(k)[0] || Pulses.at(j)[0] == Pulses.at(k)[0]
//				  || fabs(Pulses.at(i)[2]-Pulses.at(k)[2])>dtMax || fabs(Pulses.at(j)[2]-Pulses.at(k)[2])>dtMax 
//				  || fabs(Pulses.at(i)[4]-Pulses.at(k)[4])>dtMin || fabs(Pulses.at(j)[4]-Pulses.at(k)[4])>dtMin 
				){ 
					continue;
				}
				for (unsigned int l=k+1; l<Pulses.size(); l++){
					if ( Pulses.at(i)[0] == Pulses.at(l)[0] || Pulses.at(j)[0] == Pulses.at(l)[0] || Pulses.at(k)[0] == Pulses.at(l)[0] 
//				  	  || fabs(Pulses.at(i)[2]-Pulses.at(l)[2])>dtMax || fabs(Pulses.at(j)[2]-Pulses.at(l)[2])>dtMax || fabs(Pulses.at(k)[2]-Pulses.at(l)[2])>dtMax 
//				  	  || fabs(Pulses.at(i)[4]-Pulses.at(l)[4])>dtMin || fabs(Pulses.at(j)[4]-Pulses.at(l)[4])>dtMin || fabs(Pulses.at(k)[4]-Pulses.at(l)[4])>dtMin 
					){ 
						continue;
					}
					if ( outputPulses.size() == 0 ){
						outputPulses.push_back( Pulses.at(i) );
						outputPulses.push_back( Pulses.at(j) );
						outputPulses.push_back( Pulses.at(k) );
						outputPulses.push_back( Pulses.at(l) );
						break;
					}
				}
				if ( outputPulses.size() != 0 ) break;
			}
			if ( outputPulses.size() != 0 ) break;
		}
	}	
	if ( outputPulses.size() != 0 ) break;
    }
    return outputPulses;
}
