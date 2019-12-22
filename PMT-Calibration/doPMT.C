{
  gROOT->ProcessLine(".L /home/users/fsetti/milliQan/PMT-Calibration/pmtFunctions.cc+");

  TString output_file_name = "/home/users/fsetti/milliQan/MuonRate/MuonFiles/pmtCalibration";
  TString pathET = "/home/users/fsetti/milliQan/PMT-Calibration/pmt-calibration/photon_template_generator/assets/templates_et.root";
  TString pathR878 = "/home/users/fsetti/milliQan/PMT-Calibration/pmt-calibration/photon_template_generator/assets/templates_r878.root";
  TString pathR7725 = "/home/users/fsetti/milliQan/PMT-Calibration/pmt-calibration/photon_template_generator/assets/templates_r7725.root";

  timeDelay( pathET, "ET", 10. , output_file_name );
  timeDelay( pathR878, "R878", 10. , output_file_name );
  timeDelay( pathR7725, "R7725", 10. , output_file_name );

  plotAllPmts( output_file_name+"R878.root", "R878", output_file_name+"R7725.root", "R7725", output_file_name+"ET.root", "ET" ); 
}
