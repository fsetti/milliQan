//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  6 22:01:07 2017 by ROOT version 5.34/36
// from TTree t/t
// found on file: UX5MilliQan_Run25.root
//////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

unsigned int nChannels = 32;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Int_t           filenum;
   Int_t           fillnum;
   Bool_t          beam;
   Float_t         fillAvgLumi;
   Float_t         fillTotalLumi;
   Int_t           nRollOvers;
   Bool_t          present_b0;
   Bool_t          present_b1;
   Long64_t        event_time_b0;
   Long64_t        event_time_b1;
   Double_t        event_time_fromTDC;
   Int_t           t_since_fill_start;
   string          *event_t_string;
   Long64_t        event_trigger_time_tag_b0;
   Long64_t        event_trigger_time_tag_b1;
   vector<int>     *chan;
   vector<float>   *triggerCandidates;
   vector<float>   *triggerCandidatesEnd;
   vector<int>     *triggerCandidatesChannel;
   vector<int>     *layer;
   vector<int>     *row;
   vector<int>     *column;
   vector<int>     *type;
   vector<float>   *height;
   vector<float>   *ptime;
   vector<float>   *time_module_calibrated;
   vector<float>   *delay;
   vector<float>   *area;
   vector<float>   *nPE;
   vector<int>     *ipulse;
   vector<int>     *npulses;
   vector<float>   *duration;
   vector<bool>    *quiet;
   vector<float>   *presample_mean;
   vector<float>   *presample_RMS;
   vector<float>   *sideband_mean;
   vector<float>   *sideband_RMS;
   vector<Long64_t> *groupTDC_b0;
   vector<Long64_t> *groupTDC_b1;
   vector<float>   *bx;
   vector<float>   *by;
   vector<float>   *bz;
   vector<float>   *max_sample;
   vector<float>   *min_sample;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_filenum;   //!
   TBranch        *b_fillnum;   //!
   TBranch        *b_nRollOvers;   //!
   TBranch        *b_beam;   //!
   TBranch        *b_fillAvgLumi;   //!
   TBranch        *b_fillTotalLumi;   //!
   TBranch        *b_present_b0;   //!
   TBranch        *b_present_b1;   //!
   TBranch        *b_event_time_b0;   //!
   TBranch        *b_event_time_b1;   //!
   TBranch        *b_event_time_fromTDC;   //!
   TBranch        *b_t_since_fill_start;   //!
   TBranch        *b_event_t_string;   //!
   TBranch        *b_event_trigger_time_tag_b0;   //!
   TBranch        *b_event_trigger_time_tag_b1;   //!
   TBranch        *b_chan;   //!
   TBranch        *b_triggerCandidates;   //!
   TBranch        *b_triggerCandidatesEnd;   //!
   TBranch        *b_triggerCandidatesChannel;   //!
   TBranch        *b_layer;   //!
   TBranch        *b_row;   //!
   TBranch        *b_column;   //!
   TBranch        *b_type;   //!
   TBranch        *b_height;   //!
   TBranch        *b_ptime;   //!
   TBranch        *b_time_module_calibrated;   //!
   TBranch        *b_delay;   //!
   TBranch        *b_area;   //!
   TBranch        *b_nPE;   //!
   TBranch        *b_ipulse;   //!
   TBranch        *b_npulses;   //!
   TBranch        *b_duration;   //!
   TBranch        *b_quiet;   //!
   TBranch        *b_presample_mean;   //!
   TBranch        *b_presample_RMS;   //!
   TBranch        *b_sideband_mean;   //!
   TBranch        *b_sideband_RMS;   //!
   TBranch        *b_groupTDC_b0;   //!
   TBranch        *b_groupTDC_b1;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_by;   //!
   TBranch        *b_bz;   //!
   TBranch        *b_max;   //!
   TBranch        *b_min;   //!

void InitializeChain(TChain *fChain)
{

   // Set object pointer
   event_t_string = 0;
   chan = 0;
   height = 0;
   ptime = 0;
   area = 0;
   ipulse = 0;
   npulses = 0;
   duration = 0;
   sideband_mean = 0;
   sideband_RMS = 0;
   bx = 0;
   by = 0;
   bz = 0;

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("file", &filenum, &b_filenum);
   fChain->SetBranchAddress("fill", &fillnum, &b_fillnum);
   fChain->SetBranchAddress("nRollOvers", &nRollOvers, &b_nRollOvers);
   fChain->SetBranchAddress("beam", &beam, &b_beam);
   fChain->SetBranchAddress("fillAvgLumi", &fillAvgLumi, &b_fillAvgLumi);
   fChain->SetBranchAddress("fillTotalLumi", &fillTotalLumi, &b_fillTotalLumi);
   fChain->SetBranchAddress("present_b0", &present_b0, &b_present_b0);
   fChain->SetBranchAddress("present_b1", &present_b1, &b_present_b1);
   fChain->SetBranchAddress("event_time_b0", &event_time_b0, &b_event_time_b0);
   fChain->SetBranchAddress("event_time_b1", &event_time_b1, &b_event_time_b1);
   fChain->SetBranchAddress("event_time_fromTDC", &event_time_fromTDC, &b_event_time_fromTDC);
   fChain->SetBranchAddress("t_since_fill_start", &t_since_fill_start, &b_t_since_fill_start);
   fChain->SetBranchAddress("event_t_string", &event_t_string, &b_event_t_string);
   fChain->SetBranchAddress("event_trigger_time_tag_b0", &event_trigger_time_tag_b0, &b_event_trigger_time_tag_b0);
   fChain->SetBranchAddress("event_trigger_time_tag_b1", &event_trigger_time_tag_b1, &b_event_trigger_time_tag_b1);
   fChain->SetBranchAddress("chan", &chan, &b_chan);
   fChain->SetBranchAddress("triggerCandidates", &triggerCandidates, &b_triggerCandidates);
   fChain->SetBranchAddress("triggerCandidatesEnd", &triggerCandidatesEnd, &b_triggerCandidatesEnd);
   fChain->SetBranchAddress("triggerCandidatesChannel", &triggerCandidatesChannel, &b_triggerCandidatesChannel);
   fChain->SetBranchAddress("layer", &layer, &b_layer);
   fChain->SetBranchAddress("row", &row, &b_row);
   fChain->SetBranchAddress("column", &column, &b_column);
   fChain->SetBranchAddress("type", &type, &b_type);
   fChain->SetBranchAddress("height", &height, &b_height);
   fChain->SetBranchAddress("time", &ptime, &b_ptime);
   fChain->SetBranchAddress("time_module_calibrated", &time_module_calibrated, &b_time_module_calibrated);
   fChain->SetBranchAddress("delay", &delay, &b_delay);
   fChain->SetBranchAddress("area", &area, &b_area);
   fChain->SetBranchAddress("nPE", &nPE, &b_nPE);
   fChain->SetBranchAddress("ipulse", &ipulse, &b_ipulse);
   fChain->SetBranchAddress("npulses", &npulses, &b_npulses);
   fChain->SetBranchAddress("duration", &duration, &b_duration);
   fChain->SetBranchAddress("quiet", &quiet, &b_quiet);
   fChain->SetBranchAddress("presample_mean", &presample_mean, &b_presample_mean);
   fChain->SetBranchAddress("presample_RMS", &presample_RMS, &b_presample_RMS);
   fChain->SetBranchAddress("sideband_mean", &sideband_mean, &b_sideband_mean);
   fChain->SetBranchAddress("sideband_RMS", &sideband_RMS, &b_sideband_RMS);
   fChain->SetBranchAddress("groupTDC_b0", &groupTDC_b0, &b_groupTDC_b0);
   fChain->SetBranchAddress("groupTDC_b1", &groupTDC_b1, &b_groupTDC_b1);
   fChain->SetBranchAddress("max", &max_sample, &b_max);
   fChain->SetBranchAddress("min", &min_sample, &b_min);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("by", &by, &b_by);
   fChain->SetBranchAddress("bz", &bz, &b_bz);
}

