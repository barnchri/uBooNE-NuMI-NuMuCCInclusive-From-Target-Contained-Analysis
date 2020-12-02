#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"

void Making_Cut_By_Cut_Plots() {
  
  TFile* file = new TFile("all_plots_contained_only.root", "RECREATE");

  // Data
  double data_normalized_events;
  double data_muon_candidate_kinetic_energy;
  double data_muon_candidate_kinetic_energy_range;
  double data_muon_candidate_kinetic_energy_mcs;
  double data_muon_candidate_momentum;
  double data_muon_candidate_length;
  int    data_muon_candidate_is_contained;
  double data_reco_vtx_x;
  double data_reco_vtx_y;
  double data_reco_vtx_z;
  double data_other_end_location_x;
  double data_other_end_location_y;
  double data_other_end_location_z;
  double data_track_resultant;
  double data_track_x_component;
  double data_track_y_component;
  double data_track_z_component;
  double data_directional_cosine;
  double data_flash_z_difference;
  double data_flash_y_difference;
  double data_flash_y;
  double data_flash_z;
  double data_flash_PEs;
  double data_Muon_Candidate_NuScore_;
  double data_Muon_Candidate_Muon_Track_Score_;
  double data_Muon_Candidate_Proton_PID_;
  int    data_fails_fiducial_volume_req;
  int    data_fails_residuals_req;
  int    data_fails_perc_used_hits_in_cluster;

  // Signal 
  double signal_spline_fix_mcweight;
  double signal_rootino_fix_mcweight;
  double signal_central_value_mcweight;
  double signal_ppfx_value_mcweight;
  double signal_normalized_events;
  double signal_muon_candidate_kinetic_energy;
  double signal_muon_candidate_kinetic_energy_range;
  double signal_muon_candidate_kinetic_energy_mcs;
  double signal_muon_candidate_momentum;
  double signal_muon_candidate_length;
  int    signal_muon_candidate_is_contained;
  double signal_parent_fvz;
  int    signal_NC_channel;
  int    signal_num_muminus_tracks;
  int    signal_num_muplus_tracks;
  int    signal_event_is_in_AV_in_truth;
  double signal_reco_vtx_x;
  double signal_reco_vtx_y;
  double signal_reco_vtx_z;
  double signal_other_end_location_x;
  double signal_other_end_location_y;
  double signal_other_end_location_z;
  double signal_track_resultant;
  double signal_track_x_component;
  double signal_track_y_component;
  double signal_track_z_component;
  double signal_directional_cosine;
  double signal_flash_z_difference;
  double signal_flash_y_difference;
  double signal_flash_y;
  double signal_flash_z;
  double signal_flash_PEs;
  double signal_Muon_Candidate_NuScore_;
  double signal_Muon_Candidate_Muon_Track_Score_;
  double signal_Muon_Candidate_Proton_PID_;
  int    signal_fails_fiducial_volume_req;
  int    signal_fails_residuals_req;
  int    signal_fails_perc_used_hits_in_cluster;

  // EXT Background
  double ext_background_normalized_events;
  double ext_background_muon_candidate_kinetic_energy;
  double ext_background_muon_candidate_kinetic_energy_range;
  double ext_background_muon_candidate_kinetic_energy_mcs;
  double ext_background_muon_candidate_momentum;
  double ext_background_muon_candidate_length;
  int    ext_background_muon_candidate_is_contained;
  double ext_background_reco_vtx_x;
  double ext_background_reco_vtx_y;
  double ext_background_reco_vtx_z;
  double ext_background_other_end_location_x;
  double ext_background_other_end_location_y;
  double ext_background_other_end_location_z;
  double ext_background_track_resultant;
  double ext_background_track_x_component;
  double ext_background_track_y_component;
  double ext_background_track_z_component;
  double ext_background_directional_cosine;
  double ext_background_flash_z_difference;
  double ext_background_flash_y_difference;
  double ext_background_flash_y;
  double ext_background_flash_z;
  double ext_background_flash_PEs;
  double ext_background_Muon_Candidate_NuScore_;
  double ext_background_Muon_Candidate_Muon_Track_Score_;
  double ext_background_Muon_Candidate_Proton_PID_;
  int    ext_background_fails_fiducial_volume_req;
  int    ext_background_fails_residuals_req;
  int    ext_background_fails_perc_used_hits_in_cluster;

  // Neutrino Background
  double nu_background_spline_fix_mcweight;
  double nu_background_rootino_fix_mcweight;
  double nu_background_central_value_mcweight;
  double nu_background_ppfx_value_mcweight;
  double background_numi_ext_normalized_events;
  double background_numi_nu_normalized_events;
  double nu_background_muon_candidate_kinetic_energy;
  double nu_background_muon_candidate_kinetic_energy_range;
  double nu_background_muon_candidate_kinetic_energy_mcs;
  double nu_background_muon_candidate_momentum;
  double nu_background_muon_candidate_length;
  int    nu_background_muon_candidate_is_contained;
  double nu_background_parent_fvz;
  int    nu_background_NC_channel;
  int    nu_background_num_muminus_tracks;
  int    nu_background_num_muplus_tracks;
  int    nu_background_event_is_in_AV_in_truth;
  double nu_background_reco_vtx_x;
  double nu_background_reco_vtx_y;
  double nu_background_reco_vtx_z;
  double nu_background_other_end_location_x;
  double nu_background_other_end_location_y;
  double nu_background_other_end_location_z;
  double nu_background_track_resultant;
  double nu_background_track_x_component;
  double nu_background_track_y_component;
  double nu_background_track_z_component;
  double nu_background_directional_cosine;
  double nu_background_flash_z_difference;
  double nu_background_flash_y_difference;
  double nu_background_flash_y;
  double nu_background_flash_z;
  double nu_background_flash_PEs;
  double nu_background_Muon_Candidate_NuScore_;
  double nu_background_Muon_Candidate_Muon_Track_Score_;
  double nu_background_Muon_Candidate_Proton_PID_;
  int    nu_background_fails_fiducial_volume_req;
  int    nu_background_fails_residuals_req;
  int    nu_background_fails_perc_used_hits_in_cluster;
  
  // Dirt
  double dirt_spline_fix_mcweight;
  double dirt_rootino_fix_mcweight;
  double dirt_central_value_mcweight;
  double dirt_ppfx_value_mcweight;
  double dirt_normalized_events;
  double dirt_muon_candidate_kinetic_energy;
  double dirt_muon_candidate_kinetic_energy_range;
  double dirt_muon_candidate_kinetic_energy_mcs;
  double dirt_muon_candidate_momentum;
  double dirt_muon_candidate_length;
  int    dirt_muon_candidate_is_contained;
  double dirt_reco_vtx_x;
  double dirt_reco_vtx_y;
  double dirt_reco_vtx_z;
  double dirt_other_end_location_x;
  double dirt_other_end_location_y;
  double dirt_other_end_location_z;
  double dirt_track_resultant;
  double dirt_track_x_component;
  double dirt_track_y_component;
  double dirt_track_z_component;
  double dirt_directional_cosine;
  double dirt_flash_y_difference;
  double dirt_flash_z_difference;
  double dirt_flash_y;
  double dirt_flash_z;
  double dirt_flash_PEs;
  double dirt_Muon_Candidate_NuScore_;
  double dirt_Muon_Candidate_Muon_Track_Score_;
  double dirt_Muon_Candidate_Proton_PID_;
  int    dirt_fails_fiducial_volume_req;
  int    dirt_fails_residuals_req;
  int    dirt_fails_perc_used_hits_in_cluster;

  double NuMI_target_to_MicroBooNE[3];

  NuMI_target_to_MicroBooNE[0] = 0.462926; 
  NuMI_target_to_MicroBooNE[1] = 0.049604;
  NuMI_target_to_MicroBooNE[2] = 0.885008;
  

  // Declare the histograms.
  TH1F*    muon_candidate_length_signal                                      = new TH1F("muon_candidate_length_signal", "Signal Muon Candidate Length", 20, 0., 800.);
  TH1F*    muon_candidate_length_dirt                                        = new TH1F("muon_candidate_length_dirt", "Dirt Muon Candidate Length", 20, 0., 800.);
  TH1F*    muon_candidate_length_data                                        = new TH1F("muon_candidate_length_data", "Data Muon Candidate Length", 20, 0., 800.);
  TH1F*    muon_candidate_length_ext_background                              = new TH1F("muon_candidate_length_ext_background", "NuMI EXT Muon Candidate Length", 20, 0., 800.);
  TH1F*    muon_candidate_length_nu_background                               = new TH1F("muon_candidate_length_nu_background", "NuMI Neutrino Background Muon Candidate Length", 20, 0., 800.);
  THStack* muon_candidate_length_stack                                       = new THStack("muon_candidate_length_stack", "");
  TH1F*    muon_candidate_length_ratio                                       = new TH1F("muon_candidate_length_ratio", "Muon Candidate Length Ratio", 20, 0., 800.);
  
  TH1F*    muon_candidate_kinetic_energy_signal                              = new TH1F("muon_candidate_kinetic_energy_signal", "Signal Muon Candidate Kinetic_Energy", 100, 0., 2500.);
  TH1F*    muon_candidate_kinetic_energy_dirt                                = new TH1F("muon_candidate_kinetic_energy_dirt", "Dirt Muon Candidate Kinetic_Energy", 100, 0., 2500.);
  TH1F*    muon_candidate_kinetic_energy_data                                = new TH1F("muon_candidate_kinetic_energy_data", "Data Muon Candidate Kinetic_Energy", 100, 0., 2500.);
  TH1F*    muon_candidate_kinetic_energy_ext_background                      = new TH1F("muon_candidate_kinetic_energy_ext_background", "NuMI EXT Muon Candidate Kinetic_Energy", 100, 0., 2500.);
  TH1F*    muon_candidate_kinetic_energy_nu_background                       = new TH1F("muon_candidate_kinetic_energy_nu_background", "NuMI Neutrino Background Muon Candidate Kinetic_Energy", 100, 0., 2500.);
  THStack* muon_candidate_kinetic_energy_stack                               = new THStack("muon_candidate_kinetic_energy_stack", "");
  TH1F*    muon_candidate_kinetic_energy_ratio                               = new TH1F("muon_candidate_kinetic_energy_ratio", "Muon Candidate Kinetic Energy", 100, 0., 2500.);

  TH1F*    flash_y_difference_signal                                         = new TH1F("flash_y_difference_signal", "Signal Flash Y Difference", 25, -150., 150.);
  TH1F*    flash_y_difference_dirt                                           = new TH1F("flash_y_difference_dirt", "Dirt Flash Y Difference", 25, -150., 150.);
  TH1F*    flash_y_difference_data                                           = new TH1F("flash_y_difference_data", "Data Flash Y Difference", 25, -150., 150.);
  TH1F*    flash_y_difference_ext_background                                 = new TH1F("flash_y_difference_ext_background", "NuMI EXT Flash Y Difference", 25, -150., 150.);
  TH1F*    flash_y_difference_nu_background                                  = new TH1F("flash_y_difference_nu_background", "NuMI Neutrino Background Flash Y Difference", 25, -150., 150.);
  THStack* flash_y_difference_stack                                          = new THStack("flash_y_difference_stack", "");
  TH1F*    flash_y_difference_ratio                                          = new TH1F("flash_y_difference_ratio", "Flash Y Difference Ratio", 25, -150., 150.);
  
  TH1F*    flash_z_difference_signal                                         = new TH1F("flash_z_difference_signal", "Signal Flash Z Difference", 20, -200.0, 200.0);
  TH1F*    flash_z_difference_dirt                                           = new TH1F("flash_z_difference_dirt", "Dirt Flash Z Difference", 20, -200.0, 200.0);
  TH1F*    flash_z_difference_data                                           = new TH1F("flash_z_difference_data", "Data Flash Z Difference", 20, -200.0, 200.0);
  TH1F*    flash_z_difference_ext_background                                 = new TH1F("flash_z_difference_ext_background", "NuMI EXT Flash Z Difference", 20, -200.0, 200.0);
  TH1F*    flash_z_difference_nu_background                                  = new TH1F("flash_z_difference_nu_background", "NuMI Neutrino Background Flash Z Difference", 20, -200.0, 200.0);
  THStack* flash_z_difference_stack                                          = new THStack("flash_z_difference_stack", "");
  TH1F*    flash_z_difference_ratio                                          = new TH1F("flash_z_difference_ratio", "Flash Z Difference Ratio", 20, -200.0, 200.0);
  
  TH1F*    reco_vtx_x_signal                                                 = new TH1F("reco_vtx_x_signal", "Signal Vtx X", 25, 20., 236.35);
  TH1F*    reco_vtx_x_dirt                                                   = new TH1F("reco_vtx_x_dirt", "Dirt Vtx X", 25, 20., 236.35);
  TH1F*    reco_vtx_x_nu_background                                          = new TH1F("reco_vtx_x_nu_background", "Nu Background Vtx X", 25, 20., 236.35);
  TH1F*    reco_vtx_x_ext_background                                         = new TH1F("reco_vtx_x_ext_background", "EXT Background Vtx X", 25, 20., 236.35);
  TH1F*    reco_vtx_x_data                                                   = new TH1F("reco_vtx_x_data", "Data Vtx X", 25, 20., 236.35);
  THStack* reco_vtx_x_stack                                                  = new THStack("reco_vtx_x_stack", "");
  TH1F*    reco_vtx_x_ratio                                                  = new TH1F("reco_vtx_x_ratio", "Vertex X Ratio", 25, 20., 236.35);
  
  TH1F*    reco_vtx_y_signal                                                 = new TH1F("reco_vtx_y_signal", "Signal Vertex Y", 25, -96.5, 96.5);
  TH1F*    reco_vtx_y_dirt                                                   = new TH1F("reco_vtx_y_dirt", "Dirt Vertex Y", 25, -96.5, 96.5);
  TH1F*    reco_vtx_y_nu_background                                          = new TH1F("reco_vtx_y_nu_background", "Nu Background Vertex Y", 25, -96.5, 96.5);
  TH1F*    reco_vtx_y_ext_background                                         = new TH1F("reco_vtx_y_ext_background", "EXT Background Vertex Y", 25, -96.5, 96.5);
  TH1F*    reco_vtx_y_data                                                   = new TH1F("reco_vtx_y_data", "Data Vertex Y", 25, -96.5, 96.5);
  THStack* reco_vtx_y_stack                                                  = new THStack("reco_vtx_y_stack", "");
  TH1F*    reco_vtx_y_ratio                                                  = new TH1F("reco_vtx_y_ratio", "Vertex Y Ratio", 25, -96.5, 96.5);
  
  TH1F*    reco_vtx_z_signal                                                 = new TH1F("reco_vtx_z_signal", "Signal Vertex Z", 50, 20.0, 1016.8);
  TH1F*    reco_vtx_z_dirt                                                   = new TH1F("reco_vtx_z_dirt", "Dirt Vertex Z", 50, 20.0, 1016.8);
  TH1F*    reco_vtx_z_nu_background                                          = new TH1F("reco_vtx_z_nu_background", "Nu Background Vertex Z", 50, 20.0, 1016.8);
  TH1F*    reco_vtx_z_ext_background                                         = new TH1F("reco_vtx_z_ext_background", "EXT Background Vertex Z", 50, 20.0, 1016.8);
  TH1F*    reco_vtx_z_data                                                   = new TH1F("reco_vtx_z_data", "Data Vertex Z", 50, 20.0, 1016.8);
  THStack* reco_vtx_z_stack                                                  = new THStack("reco_vtx_z_stack", "");
  TH1F*    reco_vtx_z_ratio                                                  = new TH1F("reco_vtx_z_ratio", "Vertex Z Ratio", 50, 20.0, 1016.8);
  
  TH1F*    Muon_Candidate_Proton_PID__signal                                 = new TH1F("Muon_Candidate_Proton_PID__signal", "Signal Muon_Candidate_Proton_PID_", 40, 0., 500.);
  TH1F*    Muon_Candidate_Proton_PID__dirt                                   = new TH1F("Muon_Candidate_Proton_PID__dirt", "Dirt Muon_Candidate_Proton_PID_", 40, 0., 500.);
  TH1F*    Muon_Candidate_Proton_PID__nu_background                          = new TH1F("Muon_Candidate_Proton_PID__nu_background", "Nu Background Muon_Candidate_Proton_PID_", 40, 0., 500.);
  TH1F*    Muon_Candidate_Proton_PID__ext_background                         = new TH1F("Muon_Candidate_Proton_PID__ext_background", "EXT Background Muon_Candidate_Proton_PID_", 40, 0., 500.);
  TH1F*    Muon_Candidate_Proton_PID__data                                   = new TH1F("Muon_Candidate_Proton_PID__data", "Data Muon_Candidate_Proton_PID_", 40, 0., 500.);
  THStack* Muon_Candidate_Proton_PID__stack                                  = new THStack("Muon_Candidate_Proton_PID__stack", "");
  TH1F*    Muon_Candidate_Proton_PID__ratio                                  = new TH1F("Muon_Candidate_Proton_PID__ratio", "Muon Candidate Proton PID", 40, 0., 500.);
  
  TH1F*    Muon_Candidate_Muon_Track_Score__signal                           = new TH1F("Muon_Candidate_Muon_Track_Score__signal", "Signal Muon_Candidate_Muon_Track_Score_", 20, 0., 1.);
  TH1F*    Muon_Candidate_Muon_Track_Score__dirt                             = new TH1F("Muon_Candidate_Muon_Track_Score__dirt", "Dirt Muon_Candidate_Muon_Track_Score_", 20, 0., 1.);
  TH1F*    Muon_Candidate_Muon_Track_Score__nu_background                    = new TH1F("Muon_Candidate_Muon_Track_Score__nu_background", "Nu Background Muon_Candidate_Muon_Track_Score_", 20, 0., 1.);
  TH1F*    Muon_Candidate_Muon_Track_Score__ext_background                   = new TH1F("Muon_Candidate_Muon_Track_Score__ext_background", "EXT Background Muon_Candidate_Muon_Track_Score_", 20, 0., 1.);
  TH1F*    Muon_Candidate_Muon_Track_Score__data                             = new TH1F("Muon_Candidate_Muon_Track_Score__data", "Data Muon_Candidate_Muon_Track_Score_", 20, 0., 1.);
  THStack* Muon_Candidate_Muon_Track_Score__stack                            = new THStack("Muon_Candidate_Muon_Track_Score__stack", "");
  TH1F*    Muon_Candidate_Muon_Track_Score__ratio                            = new TH1F("Muon_Candidate_Muon_Track_Score__ratio", "Muon Candidate Muon Track Score ratio", 20, 0., 1.);
  
  TH1F*    Muon_Candidate_NuScore__signal                                    = new TH1F("Muon_Candidate_NuScore__signal", "Signal Muon_Candidate_NuScore_", 20, 0., 1.);
  TH1F*    Muon_Candidate_NuScore__dirt                                      = new TH1F("Muon_Candidate_NuScore__dirt", "Dirt Muon_Candidate_NuScore_", 20, 0., 1.);
  TH1F*    Muon_Candidate_NuScore__nu_background                             = new TH1F("Muon_Candidate_NuScore__nu_background", "Nu Background Muon_Candidate_NuScore_", 20, 0., 1.);
  TH1F*    Muon_Candidate_NuScore__ext_background                            = new TH1F("Muon_Candidate_NuScore__ext_background", "EXT Background Muon_Candidate_NuScore_", 20, 0., 1.);
  TH1F*    Muon_Candidate_NuScore__data                                      = new TH1F("Muon_Candidate_NuScore__data", "Data Muon_Candidate_NuScore_", 20, 0., 1.);
  THStack* Muon_Candidate_NuScore__stack                                     = new THStack("Muon_Candidate_NuScore__stack", "");
  TH1F*    Muon_Candidate_NuScore__ratio                                     = new TH1F("Muon_Candidate_NuScore__ratio", "Muon Candidate NuScore ratio", 20, 0., 1.);

  TH1F*    muon_candidate_directional_cosine_signal                          = new TH1F("muon_candidate_directional_cosine_signal", "Signal Muon Candidate Directional Cosine", 20, -1., 1. );
  TH1F*    muon_candidate_directional_cosine_dirt                            = new TH1F("muon_candidate_directional_cosine_dirt", "Dirt Muon Candidate Directional Cosine", 20, -1., 1. );
  TH1F*    muon_candidate_directional_cosine_data                            = new TH1F("muon_candidate_directional_cosine_data", "Data Muon Candidate Directional Cosine", 20, -1., 1. );
  TH1F*    muon_candidate_directional_cosine_ext_background                  = new TH1F("muon_candidate_directional_cosine_ext_background", "NuMI EXT Muon Candidate Directional Cosine", 20, -1., 1. );
  TH1F*    muon_candidate_directional_cosine_nu_background                   = new TH1F("muon_candidate_directional_cosine_nu_background", "NuMI Neutrino Background Muon Candidate Directional Cosine", 20, -1., 1. ); 
  THStack* muon_candidate_directional_cosine_stack                           = new THStack("muon_candidate_directional_cosine_stack", "");
  TH1F*    muon_candidate_directional_cosine_ratio                           = new TH1F("muon_candidate_directional_cosine_ratio", "Muon Candidate Directional Cosine Ratio", 20, -1., 1. );
  
  TChain* data  = new TChain("UBXSec/_passing_events_tree");
  data->Add("/home/barnchri/Making_Cut_By_Cut_Plots_In_Analysis/NuMI_Run1_OnBeam_Output_For_Going_Cut_By_Cut.root");
  data->SetBranchAddress("muon_candidate_kinetic_energy_range", &data_muon_candidate_kinetic_energy_range);
  data->SetBranchAddress("muon_candidate_kinetic_energy_mcs", &data_muon_candidate_kinetic_energy_mcs);
  data->SetBranchAddress("muon_candidate_length", &data_muon_candidate_length);
  data->SetBranchAddress("muon_candidate_is_contained", &data_muon_candidate_is_contained);
  data->SetBranchAddress("vertex_location_x", &data_reco_vtx_x);
  data->SetBranchAddress("vertex_location_y", &data_reco_vtx_y);
  data->SetBranchAddress("vertex_location_z", &data_reco_vtx_z);
  data->SetBranchAddress("other_end_location_x", &data_other_end_location_x);
  data->SetBranchAddress("other_end_location_y", &data_other_end_location_y);
  data->SetBranchAddress("other_end_location_z", &data_other_end_location_z);
  data->SetBranchAddress("flash_z", &data_flash_z);
  data->SetBranchAddress("flash_y", &data_flash_y);
  data->SetBranchAddress("flash_PEs", &data_flash_PEs);
  data->SetBranchAddress("Muon_Candidate_NuScore_", &data_Muon_Candidate_NuScore_);
  data->SetBranchAddress("Muon_Candidate_Muon_Track_Score_", &data_Muon_Candidate_Muon_Track_Score_);
  data->SetBranchAddress("Muon_Candidate_Proton_PID_", &data_Muon_Candidate_Proton_PID_);
  data->SetBranchAddress("fails_fiducial_volume_req", &data_fails_fiducial_volume_req);
  data->SetBranchAddress("fails_residuals_req", &data_fails_residuals_req);
  data->SetBranchAddress("fails_perc_used_hits_in_cluster", &data_fails_perc_used_hits_in_cluster);
  
  TChain* signal = new TChain("UBXSec/_passing_events_tree");
  signal->Add("/home/barnchri/Making_Cut_By_Cut_Plots_In_Analysis/NuMI_Run1_Overlay_Output_For_Going_Cut_By_Cut.root");
  signal->SetBranchAddress("spline_fix_mcweight", &signal_spline_fix_mcweight);
  signal->SetBranchAddress("rootino_fix_mcweight", &signal_rootino_fix_mcweight);
  signal->SetBranchAddress("central_value_mcweight", &signal_central_value_mcweight);
  signal->SetBranchAddress("ppfx_value_mcweight", &signal_ppfx_value_mcweight );
  signal->SetBranchAddress("muon_candidate_kinetic_energy_range", &signal_muon_candidate_kinetic_energy_range);
  signal->SetBranchAddress("muon_candidate_kinetic_energy_mcs", &signal_muon_candidate_kinetic_energy_mcs);
  signal->SetBranchAddress("muon_candidate_length", &signal_muon_candidate_length);
  signal->SetBranchAddress("muon_candidate_is_contained", &signal_muon_candidate_is_contained);
  signal->SetBranchAddress("parent_fvz", &signal_parent_fvz);
  signal->SetBranchAddress("NC_channel", &signal_NC_channel);
  signal->SetBranchAddress("num_muminus_tracks", &signal_num_muminus_tracks);
  signal->SetBranchAddress("num_muplus_tracks", &signal_num_muplus_tracks);
  signal->SetBranchAddress("event_is_in_AV_in_truth", &signal_event_is_in_AV_in_truth);
  signal->SetBranchAddress("vertex_location_x", &signal_reco_vtx_x);
  signal->SetBranchAddress("vertex_location_y", &signal_reco_vtx_y);
  signal->SetBranchAddress("vertex_location_z", &signal_reco_vtx_z);
  signal->SetBranchAddress("other_end_location_x", &signal_other_end_location_x);
  signal->SetBranchAddress("other_end_location_y", &signal_other_end_location_y);
  signal->SetBranchAddress("other_end_location_z", &signal_other_end_location_z);
  signal->SetBranchAddress("flash_z", &signal_flash_z);
  signal->SetBranchAddress("flash_y", &signal_flash_y);
  signal->SetBranchAddress("flash_PEs", &signal_flash_PEs);
  signal->SetBranchAddress("Muon_Candidate_NuScore_", &signal_Muon_Candidate_NuScore_);
  signal->SetBranchAddress("Muon_Candidate_Muon_Track_Score_", &signal_Muon_Candidate_Muon_Track_Score_);
  signal->SetBranchAddress("Muon_Candidate_Proton_PID_", &signal_Muon_Candidate_Proton_PID_);
  signal->SetBranchAddress("fails_fiducial_volume_req", &signal_fails_fiducial_volume_req);
  signal->SetBranchAddress("fails_residuals_req", &signal_fails_residuals_req);
  signal->SetBranchAddress("fails_perc_used_hits_in_cluster", &signal_fails_perc_used_hits_in_cluster);

  TChain* ext_background  = new TChain("UBXSec/_passing_events_tree");
  ext_background->Add("/home/barnchri/Making_Cut_By_Cut_Plots_In_Analysis/NuMI_Run1_OffBeam_Output_For_Going_Cut_By_Cut.root");
  ext_background->SetBranchAddress("muon_candidate_kinetic_energy_range", &ext_background_muon_candidate_kinetic_energy_range);
  ext_background->SetBranchAddress("muon_candidate_kinetic_energy_mcs", &ext_background_muon_candidate_kinetic_energy_mcs);
  ext_background->SetBranchAddress("muon_candidate_length", &ext_background_muon_candidate_length);
  ext_background->SetBranchAddress("muon_candidate_is_contained", &ext_background_muon_candidate_is_contained);
  ext_background->SetBranchAddress("flash_z", &ext_background_flash_z);
  ext_background->SetBranchAddress("flash_y", &ext_background_flash_y);
  ext_background->SetBranchAddress("flash_PEs", &ext_background_flash_PEs);
  ext_background->SetBranchAddress("vertex_location_x", &ext_background_reco_vtx_x);
  ext_background->SetBranchAddress("vertex_location_y", &ext_background_reco_vtx_y);
  ext_background->SetBranchAddress("vertex_location_z", &ext_background_reco_vtx_z);
  ext_background->SetBranchAddress("other_end_location_x", &ext_background_other_end_location_x);
  ext_background->SetBranchAddress("other_end_location_y", &ext_background_other_end_location_y);
  ext_background->SetBranchAddress("other_end_location_z", &ext_background_other_end_location_z);
  ext_background->SetBranchAddress("Muon_Candidate_NuScore_", &ext_background_Muon_Candidate_NuScore_);
  ext_background->SetBranchAddress("Muon_Candidate_Muon_Track_Score_", &ext_background_Muon_Candidate_Muon_Track_Score_);
  ext_background->SetBranchAddress("Muon_Candidate_Proton_PID_", &ext_background_Muon_Candidate_Proton_PID_);
  ext_background->SetBranchAddress("fails_fiducial_volume_req", &ext_background_fails_fiducial_volume_req);
  ext_background->SetBranchAddress("fails_residuals_req", &ext_background_fails_residuals_req);
  ext_background->SetBranchAddress("fails_perc_used_hits_in_cluster", &ext_background_fails_perc_used_hits_in_cluster);
  
  TChain* nu_background = new TChain("UBXSec/_passing_events_tree");
  nu_background->Add("/home/barnchri/Making_Cut_By_Cut_Plots_In_Analysis/NuMI_Run1_Overlay_Output_For_Going_Cut_By_Cut.root");
  nu_background->SetBranchAddress("spline_fix_mcweight", &nu_background_spline_fix_mcweight);
  nu_background->SetBranchAddress("rootino_fix_mcweight", &nu_background_rootino_fix_mcweight);
  nu_background->SetBranchAddress("central_value_mcweight", &nu_background_central_value_mcweight);
  nu_background->SetBranchAddress("ppfx_value_mcweight", &nu_background_ppfx_value_mcweight );
  nu_background->SetBranchAddress("muon_candidate_kinetic_energy_range", &nu_background_muon_candidate_kinetic_energy_range);
  nu_background->SetBranchAddress("muon_candidate_kinetic_energy_mcs", &nu_background_muon_candidate_kinetic_energy_mcs);
  nu_background->SetBranchAddress("muon_candidate_length", &nu_background_muon_candidate_length);
  nu_background->SetBranchAddress("muon_candidate_is_contained", &nu_background_muon_candidate_is_contained);
  nu_background->SetBranchAddress("parent_fvz", &nu_background_parent_fvz);
  nu_background->SetBranchAddress("NC_channel", &nu_background_NC_channel);
  nu_background->SetBranchAddress("num_muminus_tracks", &nu_background_num_muminus_tracks);
  nu_background->SetBranchAddress("num_muplus_tracks", &nu_background_num_muplus_tracks);
  nu_background->SetBranchAddress("event_is_in_AV_in_truth", &nu_background_event_is_in_AV_in_truth);
  nu_background->SetBranchAddress("vertex_location_x", &nu_background_reco_vtx_x);
  nu_background->SetBranchAddress("vertex_location_y", &nu_background_reco_vtx_y);
  nu_background->SetBranchAddress("vertex_location_z", &nu_background_reco_vtx_z);
  nu_background->SetBranchAddress("other_end_location_x", &nu_background_other_end_location_x);
  nu_background->SetBranchAddress("other_end_location_y", &nu_background_other_end_location_y);
  nu_background->SetBranchAddress("other_end_location_z", &nu_background_other_end_location_z);
  nu_background->SetBranchAddress("flash_z", &nu_background_flash_z);
  nu_background->SetBranchAddress("flash_y", &nu_background_flash_y);
  nu_background->SetBranchAddress("flash_PEs", &nu_background_flash_PEs);
  nu_background->SetBranchAddress("Muon_Candidate_NuScore_", &nu_background_Muon_Candidate_NuScore_);
  nu_background->SetBranchAddress("Muon_Candidate_Muon_Track_Score_", &nu_background_Muon_Candidate_Muon_Track_Score_);
  nu_background->SetBranchAddress("Muon_Candidate_Proton_PID_", &nu_background_Muon_Candidate_Proton_PID_);
  nu_background->SetBranchAddress("fails_fiducial_volume_req", &nu_background_fails_fiducial_volume_req);
  nu_background->SetBranchAddress("fails_residuals_req", &nu_background_fails_residuals_req);
  nu_background->SetBranchAddress("fails_perc_used_hits_in_cluster", &nu_background_fails_perc_used_hits_in_cluster);
  
  // Dirt
  TChain* dirt = new TChain("UBXSec/_passing_events_tree");
  dirt->Add("/home/barnchri/Making_Cut_By_Cut_Plots_In_Analysis/NuMI_Run1_Dirt_Output_For_Going_Cut_By_Cut.root");
  dirt->SetBranchAddress("spline_fix_mcweight", &dirt_spline_fix_mcweight);
  dirt->SetBranchAddress("rootino_fix_mcweight", &dirt_rootino_fix_mcweight);
  dirt->SetBranchAddress("central_value_mcweight", &dirt_central_value_mcweight);
  dirt->SetBranchAddress("ppfx_value_mcweight", &dirt_ppfx_value_mcweight );
  dirt->SetBranchAddress("muon_candidate_kinetic_energy_range", &dirt_muon_candidate_kinetic_energy_range);
  dirt->SetBranchAddress("muon_candidate_kinetic_energy_mcs", &dirt_muon_candidate_kinetic_energy_mcs);
  dirt->SetBranchAddress("muon_candidate_length", &dirt_muon_candidate_length);
  dirt->SetBranchAddress("muon_candidate_is_contained", &dirt_muon_candidate_is_contained);
  dirt->SetBranchAddress("vertex_location_x", &dirt_reco_vtx_x);
  dirt->SetBranchAddress("vertex_location_y", &dirt_reco_vtx_y);
  dirt->SetBranchAddress("vertex_location_z", &dirt_reco_vtx_z);
  dirt->SetBranchAddress("other_end_location_x", &dirt_other_end_location_x);
  dirt->SetBranchAddress("other_end_location_y", &dirt_other_end_location_y);
  dirt->SetBranchAddress("other_end_location_z", &dirt_other_end_location_z);
  dirt->SetBranchAddress("flash_z", &dirt_flash_z);
  dirt->SetBranchAddress("flash_y", &dirt_flash_y);
  dirt->SetBranchAddress("flash_PEs", &dirt_flash_PEs);
  dirt->SetBranchAddress("Muon_Candidate_NuScore_", &dirt_Muon_Candidate_NuScore_);
  dirt->SetBranchAddress("Muon_Candidate_Muon_Track_Score_", &dirt_Muon_Candidate_Muon_Track_Score_);
  dirt->SetBranchAddress("Muon_Candidate_Proton_PID_", &dirt_Muon_Candidate_Proton_PID_);
  dirt->SetBranchAddress("fails_fiducial_volume_req", &dirt_fails_fiducial_volume_req);
  dirt->SetBranchAddress("fails_residuals_req", &dirt_fails_residuals_req);
  dirt->SetBranchAddress("fails_perc_used_hits_in_cluster", &dirt_fails_perc_used_hits_in_cluster);
  
  int num_data_entries                        = data->GetEntries();
  double num_data_events_above_score_threshold = 0.;

  std::cout << "Making it to the top of the event loops." << std::endl;
  
  for ( int data_iter = 0; data_iter < num_data_entries; data_iter++ ) {

    data->GetEntry( data_iter );

    if ( data_muon_candidate_is_contained == 1 )
      data_muon_candidate_kinetic_energy = data_muon_candidate_kinetic_energy_range; 

    else
      continue;

    if ( data_Muon_Candidate_Muon_Track_Score_ < 0.8 || data_Muon_Candidate_Proton_PID_ < 78. || data_Muon_Candidate_NuScore_ < 0.1 )
      continue;

    if ( data_flash_PEs < 50. )
      continue;

    if ( data_fails_fiducial_volume_req == 1 )
      continue;

    if ( data_fails_residuals_req == 1 )
      continue;

    if ( data_fails_perc_used_hits_in_cluster == 1 )
      continue;
    
    muon_candidate_kinetic_energy_data->Fill( data_muon_candidate_kinetic_energy );
    muon_candidate_length_data->Fill( data_muon_candidate_length );
    
    data_flash_y_difference = data_reco_vtx_y - data_flash_y;
    data_flash_z_difference = data_reco_vtx_z - data_flash_z;
    
    flash_y_difference_data->Fill( data_flash_y_difference );
    flash_z_difference_data->Fill( data_flash_z_difference );
    reco_vtx_x_data->Fill( data_reco_vtx_x );
    reco_vtx_y_data->Fill( data_reco_vtx_y );
    reco_vtx_z_data->Fill( data_reco_vtx_z );
    Muon_Candidate_Proton_PID__data->Fill( data_Muon_Candidate_Proton_PID_ );
    Muon_Candidate_Muon_Track_Score__data->Fill( data_Muon_Candidate_Muon_Track_Score_ );
    Muon_Candidate_NuScore__data->Fill( data_Muon_Candidate_NuScore_ );

    // Calculate the directional cosine.
    data_track_resultant = TMath::Sqrt( ( data_other_end_location_x - data_reco_vtx_x ) * ( data_other_end_location_x - data_reco_vtx_x ) + ( data_other_end_location_y - data_reco_vtx_y ) * ( data_other_end_location_y - data_reco_vtx_y ) + ( data_other_end_location_z - data_reco_vtx_z ) * ( data_other_end_location_z - data_reco_vtx_z ) );

    data_track_x_component = ( data_other_end_location_x - data_reco_vtx_x ) / data_track_resultant;
    data_track_y_component = ( data_other_end_location_y - data_reco_vtx_y ) / data_track_resultant;
    data_track_z_component = ( data_other_end_location_z - data_reco_vtx_z ) / data_track_resultant;

    data_directional_cosine = ( NuMI_target_to_MicroBooNE[0] * data_track_x_component + NuMI_target_to_MicroBooNE[1] * data_track_y_component + NuMI_target_to_MicroBooNE[2] * data_track_z_component );

    muon_candidate_directional_cosine_data->Fill( data_directional_cosine );
    
    num_data_events_above_score_threshold += 1.0;

  }

  std::cout << "The number of data events = " << num_data_events_above_score_threshold << "." << std::endl;
  
  int num_signal_entries = signal->GetEntries();

  double   num_unweighted_signal_events_above_score_threshold = 0;
  double   num_weighted_signal_events_above_score_threshold   = 0;

  std::cout << "Total number signal events = " << num_signal_entries << "." << std::endl;
  std::cout << "Outside the signal loop." << std::endl;
  
  for ( int signal_iter = 0; signal_iter < num_signal_entries; signal_iter++ ) {

    std::cout << "Looping over signal entry #" << signal_iter << "." << std::endl;
    
    signal->GetEntry( signal_iter );

    if ( signal_NC_channel == 1 || ( signal_num_muminus_tracks == 0 && signal_num_muplus_tracks == 0 ) || signal_event_is_in_AV_in_truth == 0 || signal_parent_fvz > 3000. )
      continue;
    
    if ( signal_muon_candidate_is_contained == 1 )
      signal_muon_candidate_kinetic_energy = signal_muon_candidate_kinetic_energy_range;
      
    else
      continue;

    if ( std::isnan( signal_rootino_fix_mcweight ) || std::isinf( signal_rootino_fix_mcweight ) )
      signal_rootino_fix_mcweight = 1.;

    if ( std::isnan( signal_central_value_mcweight ) || std::isinf( signal_central_value_mcweight ) )
      signal_central_value_mcweight = 1.;

    if ( std::isnan( signal_spline_fix_mcweight ) || std::isinf( signal_spline_fix_mcweight ) )
      signal_spline_fix_mcweight = 1.;

    if ( std::isnan( signal_ppfx_value_mcweight ) || std::isinf( signal_ppfx_value_mcweight ) )
      signal_ppfx_value_mcweight = 1.;

    if ( signal_central_value_mcweight > 2000. )
      signal_central_value_mcweight = 1.;

    if ( signal_Muon_Candidate_Muon_Track_Score_ < 0.8 || signal_Muon_Candidate_Proton_PID_ < 78. || signal_Muon_Candidate_NuScore_ < 0.1 )
      continue;

    if ( signal_flash_PEs < 50. )
      continue;
    
    if ( signal_fails_fiducial_volume_req == 1 )
      continue;

    if ( signal_fails_residuals_req == 1 )
      continue;

    if ( signal_fails_perc_used_hits_in_cluster == 1 )
      continue;

    muon_candidate_length_signal->Fill( signal_muon_candidate_length, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight);
    muon_candidate_kinetic_energy_signal->Fill( signal_muon_candidate_kinetic_energy, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );

    signal_flash_y_difference = signal_reco_vtx_y - signal_flash_y;
    signal_flash_z_difference = signal_reco_vtx_z - signal_flash_z;
    
    flash_y_difference_signal->Fill( signal_flash_y_difference, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    flash_z_difference_signal->Fill( signal_flash_z_difference, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    reco_vtx_x_signal->Fill( signal_reco_vtx_x, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    reco_vtx_y_signal->Fill( signal_reco_vtx_y, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    reco_vtx_z_signal->Fill( signal_reco_vtx_z, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    Muon_Candidate_Proton_PID__signal->Fill( signal_Muon_Candidate_Proton_PID_, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    Muon_Candidate_Muon_Track_Score__signal->Fill( signal_Muon_Candidate_Muon_Track_Score_, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    Muon_Candidate_NuScore__signal->Fill( signal_Muon_Candidate_NuScore_, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );

    // Calculate the directional cosine.
    signal_track_resultant = TMath::Sqrt( ( signal_other_end_location_x - signal_reco_vtx_x ) * ( signal_other_end_location_x - signal_reco_vtx_x ) + ( signal_other_end_location_y - signal_reco_vtx_y ) * ( signal_other_end_location_y - signal_reco_vtx_y ) + ( signal_other_end_location_z - signal_reco_vtx_z ) * ( signal_other_end_location_z - signal_reco_vtx_z ) );

    signal_track_x_component = ( signal_other_end_location_x - signal_reco_vtx_x ) / signal_track_resultant;
    signal_track_y_component = ( signal_other_end_location_y - signal_reco_vtx_y ) / signal_track_resultant;
    signal_track_z_component = ( signal_other_end_location_z - signal_reco_vtx_z ) / signal_track_resultant;

    signal_directional_cosine = ( NuMI_target_to_MicroBooNE[0] * signal_track_x_component + NuMI_target_to_MicroBooNE[1] * signal_track_y_component + NuMI_target_to_MicroBooNE[2] * signal_track_z_component );

    muon_candidate_directional_cosine_signal->Fill( signal_directional_cosine, signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );
    
    num_unweighted_signal_events_above_score_threshold++;
    num_weighted_signal_events_above_score_threshold += ( signal_ppfx_value_mcweight * signal_spline_fix_mcweight * signal_rootino_fix_mcweight * signal_central_value_mcweight );

  }

  int num_ext_entries                                     = ext_background->GetEntries();
  double num_unweighted_ext_events_above_score_threshold  = 0.;
  
  for ( int ext_background_iter = 0; ext_background_iter < num_ext_entries; ext_background_iter++ ) {

    ext_background->GetEntry( ext_background_iter );

    if ( ext_background_muon_candidate_is_contained == 1 )
      ext_background_muon_candidate_kinetic_energy = ext_background_muon_candidate_kinetic_energy_range; 

    else
      continue;

    if ( ext_background_Muon_Candidate_Muon_Track_Score_ < 0.8 || ext_background_Muon_Candidate_Proton_PID_ < 78. || ext_background_Muon_Candidate_NuScore_ < 0.1 )
      continue;

    if ( ext_background_flash_PEs < 50. )
      continue;

    if ( ext_background_fails_fiducial_volume_req == 1 )
      continue;

    if ( ext_background_fails_residuals_req == 1 )
      continue;

    if ( ext_background_fails_perc_used_hits_in_cluster == 1 )
      continue;
    
    ext_background_flash_y_difference = ext_background_reco_vtx_y - ext_background_flash_y;
    ext_background_flash_z_difference = ext_background_reco_vtx_z - ext_background_flash_z;

    flash_y_difference_ext_background->Fill( ext_background_flash_y_difference );
    flash_z_difference_ext_background->Fill( ext_background_flash_z_difference );
    reco_vtx_x_ext_background->Fill( ext_background_reco_vtx_x );
    reco_vtx_y_ext_background->Fill( ext_background_reco_vtx_y );
    reco_vtx_z_ext_background->Fill( ext_background_reco_vtx_z );
    Muon_Candidate_Proton_PID__ext_background->Fill( ext_background_Muon_Candidate_Proton_PID_ );
    Muon_Candidate_Muon_Track_Score__ext_background->Fill( ext_background_Muon_Candidate_Muon_Track_Score_ );
    Muon_Candidate_NuScore__ext_background->Fill( ext_background_Muon_Candidate_NuScore_ );

    muon_candidate_length_ext_background->Fill( ext_background_muon_candidate_length );
    muon_candidate_kinetic_energy_ext_background->Fill( ext_background_muon_candidate_kinetic_energy );

    // Calculate the directional cosine.
    ext_background_track_resultant = TMath::Sqrt( ( ext_background_other_end_location_x - ext_background_reco_vtx_x ) * ( ext_background_other_end_location_x - ext_background_reco_vtx_x ) + ( ext_background_other_end_location_y - ext_background_reco_vtx_y ) * ( ext_background_other_end_location_y - ext_background_reco_vtx_y ) + ( ext_background_other_end_location_z - ext_background_reco_vtx_z ) * ( ext_background_other_end_location_z - ext_background_reco_vtx_z ) );

    ext_background_track_x_component = ( ext_background_other_end_location_x - ext_background_reco_vtx_x ) / ext_background_track_resultant;
    ext_background_track_y_component = ( ext_background_other_end_location_y - ext_background_reco_vtx_y ) / ext_background_track_resultant;
    ext_background_track_z_component = ( ext_background_other_end_location_z - ext_background_reco_vtx_z ) / ext_background_track_resultant;

    ext_background_directional_cosine = ( NuMI_target_to_MicroBooNE[0] * ext_background_track_x_component + NuMI_target_to_MicroBooNE[1] * ext_background_track_y_component + NuMI_target_to_MicroBooNE[2] * ext_background_track_z_component );

    muon_candidate_directional_cosine_ext_background->Fill( ext_background_directional_cosine );
        
    num_unweighted_ext_events_above_score_threshold++;

  }


  std::cout << "The number of unweighted KDAR signal events with a score = " << num_unweighted_signal_events_above_score_threshold << "." << std::endl;

  int num_nu_background_entries                                                  = nu_background->GetEntries();
  double num_unweighted_nu_events_above_score_threshold                          = 0.;
  double num_weighted_nu_events_above_score_threshold                            = 0.;

  std::cout << "Total number nu background events = " << num_nu_background_entries <<"." << std::endl;
  std::cout << "Outside the nu background loop." << std::endl;
  
  for ( int nu_background_iter = 0; nu_background_iter < num_nu_background_entries; nu_background_iter++ ) {

    std::cout << "nu_background event #" << nu_background_iter << "." << std::endl;
    
    nu_background->GetEntry( nu_background_iter );
    
    if ( nu_background_muon_candidate_is_contained == 1 )
      nu_background_muon_candidate_kinetic_energy = nu_background_muon_candidate_kinetic_energy_range;

    else
      continue;

    if ( nu_background_Muon_Candidate_Muon_Track_Score_ < 0.8 || nu_background_Muon_Candidate_Proton_PID_ < 78. || nu_background_Muon_Candidate_NuScore_ < 0.1 )
      continue;

    if ( nu_background_flash_PEs < 50. )
      continue;

    if ( nu_background_fails_fiducial_volume_req == 1 )
      continue;
    
    if ( nu_background_fails_residuals_req == 1 )
      continue;

    if ( nu_background_fails_perc_used_hits_in_cluster == 1 )
      continue;
    
    if ( nu_background_NC_channel == 0 && ( nu_background_num_muminus_tracks > 0 || nu_background_num_muplus_tracks > 0 ) && nu_background_event_is_in_AV_in_truth == 1 && nu_background_parent_fvz < 3000. )
      continue;
      
    if ( std::isnan( nu_background_rootino_fix_mcweight ) || std::isinf( nu_background_rootino_fix_mcweight ) )
      nu_background_rootino_fix_mcweight = 1.;

    if ( std::isnan( nu_background_central_value_mcweight ) || std::isinf( nu_background_central_value_mcweight ) )
      nu_background_central_value_mcweight = 1.;

    if ( std::isnan( nu_background_spline_fix_mcweight ) || std::isinf( nu_background_spline_fix_mcweight ) )
      nu_background_spline_fix_mcweight = 1.;

    if ( std::isnan( nu_background_ppfx_value_mcweight ) || std::isinf( nu_background_ppfx_value_mcweight ) )
      nu_background_ppfx_value_mcweight = 1.;

    muon_candidate_length_nu_background->Fill( nu_background_muon_candidate_length, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight);
    muon_candidate_kinetic_energy_nu_background->Fill( nu_background_muon_candidate_kinetic_energy, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    
    nu_background_flash_y_difference = nu_background_reco_vtx_y - nu_background_flash_y;
    nu_background_flash_z_difference = nu_background_reco_vtx_z - nu_background_flash_z;
    
    flash_y_difference_nu_background->Fill( nu_background_flash_y_difference, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    flash_z_difference_nu_background->Fill( nu_background_flash_z_difference, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    reco_vtx_x_nu_background->Fill( nu_background_reco_vtx_x, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    reco_vtx_y_nu_background->Fill( nu_background_reco_vtx_y, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    reco_vtx_z_nu_background->Fill( nu_background_reco_vtx_z, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    Muon_Candidate_Proton_PID__nu_background->Fill( nu_background_Muon_Candidate_Proton_PID_, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    Muon_Candidate_Muon_Track_Score__nu_background->Fill( nu_background_Muon_Candidate_Muon_Track_Score_, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );
    Muon_Candidate_NuScore__nu_background->Fill( nu_background_Muon_Candidate_NuScore_, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );

    // Calculate the directional cosine.
    nu_background_track_resultant = TMath::Sqrt( ( nu_background_other_end_location_x - nu_background_reco_vtx_x ) * ( nu_background_other_end_location_x - nu_background_reco_vtx_x ) + ( nu_background_other_end_location_y - nu_background_reco_vtx_y ) * ( nu_background_other_end_location_y - nu_background_reco_vtx_y ) + ( nu_background_other_end_location_z - nu_background_reco_vtx_z ) * ( nu_background_other_end_location_z - nu_background_reco_vtx_z ) );

    nu_background_track_x_component = ( nu_background_other_end_location_x - nu_background_reco_vtx_x ) / nu_background_track_resultant;
    nu_background_track_y_component = ( nu_background_other_end_location_y - nu_background_reco_vtx_y ) / nu_background_track_resultant;
    nu_background_track_z_component = ( nu_background_other_end_location_z - nu_background_reco_vtx_z ) / nu_background_track_resultant;

    nu_background_directional_cosine = ( NuMI_target_to_MicroBooNE[0] * nu_background_track_x_component + NuMI_target_to_MicroBooNE[1] * nu_background_track_y_component + NuMI_target_to_MicroBooNE[2] * nu_background_track_z_component );

    muon_candidate_directional_cosine_nu_background->Fill( nu_background_directional_cosine, nu_background_ppfx_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight * nu_background_central_value_mcweight );  
    
    num_weighted_nu_events_above_score_threshold += ( nu_background_ppfx_value_mcweight * nu_background_central_value_mcweight * nu_background_spline_fix_mcweight * nu_background_rootino_fix_mcweight );
    
  }

  // Loop through each of the systematic universes and fill the histograms.
  // Dirt
  int num_dirt_entries                                  = dirt->GetEntries();
  double num_weighted_dirt_events_above_score_threshold = 0.;

  for ( int dirt_iter = 0; dirt_iter < num_dirt_entries; dirt_iter++ ) {

    dirt->GetEntry( dirt_iter );
    
    if ( dirt_muon_candidate_is_contained == 1 )
      dirt_muon_candidate_kinetic_energy = dirt_muon_candidate_kinetic_energy_range; 
    
    else
      continue;

    if ( std::isnan( dirt_rootino_fix_mcweight ) || std::isinf( dirt_rootino_fix_mcweight ) )
      dirt_rootino_fix_mcweight = 1.;
    
    if ( std::isnan( dirt_central_value_mcweight ) || std::isinf( dirt_central_value_mcweight ) )
      dirt_central_value_mcweight = 1.;

    if ( std::isnan( dirt_spline_fix_mcweight ) || std::isinf( dirt_spline_fix_mcweight ) )
      dirt_spline_fix_mcweight = 1.;

    if ( std::isnan( dirt_ppfx_value_mcweight ) || std::isinf( dirt_ppfx_value_mcweight ) )
      dirt_ppfx_value_mcweight = 1.;

    if ( dirt_Muon_Candidate_Muon_Track_Score_ < 0.8 || dirt_Muon_Candidate_Proton_PID_ < 78. || dirt_Muon_Candidate_NuScore_ < 0.1 )
      continue;
    
    if ( dirt_flash_PEs < 50. )
      continue;

    if ( dirt_fails_fiducial_volume_req == 1 )
      continue;

    if ( dirt_fails_residuals_req == 1 )
      continue;

    if ( dirt_fails_perc_used_hits_in_cluster == 1 )
      continue;

    muon_candidate_length_dirt->Fill( dirt_muon_candidate_length, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight);
    muon_candidate_kinetic_energy_dirt->Fill( dirt_muon_candidate_kinetic_energy, dirt_ppfx_value_mcweight * dirt_central_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight );

    dirt_flash_y_difference = dirt_reco_vtx_y - dirt_flash_y;
    dirt_flash_z_difference = dirt_reco_vtx_z - dirt_flash_z;

    flash_y_difference_dirt->Fill( dirt_flash_y_difference, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    flash_z_difference_dirt->Fill( dirt_flash_z_difference, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    reco_vtx_x_dirt->Fill( dirt_reco_vtx_x, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    reco_vtx_y_dirt->Fill( dirt_reco_vtx_y, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    reco_vtx_z_dirt->Fill( dirt_reco_vtx_z, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    Muon_Candidate_Proton_PID__dirt->Fill( dirt_Muon_Candidate_Proton_PID_, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    Muon_Candidate_Muon_Track_Score__dirt->Fill( dirt_Muon_Candidate_Muon_Track_Score_, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    Muon_Candidate_NuScore__dirt->Fill( dirt_Muon_Candidate_NuScore_, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );

    // Calculate the directional cosine.
    dirt_track_resultant = TMath::Sqrt( ( dirt_other_end_location_x - dirt_reco_vtx_x ) * ( dirt_other_end_location_x - dirt_reco_vtx_x ) + ( dirt_other_end_location_y - dirt_reco_vtx_y ) * ( dirt_other_end_location_y - dirt_reco_vtx_y ) + ( dirt_other_end_location_z - dirt_reco_vtx_z ) * ( dirt_other_end_location_z - dirt_reco_vtx_z ) );
    
    dirt_track_x_component = ( dirt_other_end_location_x - dirt_reco_vtx_x ) / dirt_track_resultant;
    dirt_track_y_component = ( dirt_other_end_location_y - dirt_reco_vtx_y ) / dirt_track_resultant;
    dirt_track_z_component = ( dirt_other_end_location_z - dirt_reco_vtx_z ) / dirt_track_resultant;

    dirt_directional_cosine = ( NuMI_target_to_MicroBooNE[0] * dirt_track_x_component + NuMI_target_to_MicroBooNE[1] * dirt_track_y_component + NuMI_target_to_MicroBooNE[2] * dirt_track_z_component );
    muon_candidate_directional_cosine_dirt->Fill( dirt_directional_cosine, dirt_ppfx_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight * dirt_central_value_mcweight );
    
    
    num_weighted_dirt_events_above_score_threshold += ( dirt_ppfx_value_mcweight * dirt_central_value_mcweight * dirt_spline_fix_mcweight * dirt_rootino_fix_mcweight );

  }
  
  std::cout << "The number of unweighted EXT events with a score = " << num_unweighted_ext_events_above_score_threshold << "." << std::endl;
  std::cout << "The number of unweighted NuMI background events = " << num_unweighted_nu_events_above_score_threshold << "." << std::endl;

  double signal_normalization_factor         = 0.0378;
  double numi_ext_normalization_factor       = 0.678;
  double numi_nu_normalization_factor        = 0.0378;
  double dirt_normalization_factor           = 0.0518;
  
  double numi_ext_extra_normalization_factor = 0.98;
  double dirt_extra_normalization_factor     = 0.35;
  
  signal_normalized_events              = ( num_weighted_signal_events_above_score_threshold * signal_normalization_factor );
  background_numi_ext_normalized_events = ( num_unweighted_ext_events_above_score_threshold * numi_ext_normalization_factor * numi_ext_extra_normalization_factor ); 
  background_numi_nu_normalized_events  = ( num_weighted_nu_events_above_score_threshold * numi_nu_normalization_factor );
  dirt_normalized_events                = ( num_weighted_dirt_events_above_score_threshold * dirt_normalization_factor * dirt_extra_normalization_factor );
  data_normalized_events                = num_data_events_above_score_threshold;
  
  std::cout << "The number of normalized signal events = " << signal_normalized_events << "." << std::endl;
  std::cout << "The number of normalized NuMI EXT events = " << background_numi_ext_normalized_events << "." << std::endl;
  std::cout << "The number of normalized background nu events = " << background_numi_nu_normalized_events << "." << std::endl;
  std::cout << "The number of normalized dirt events = " << dirt_normalized_events << "." << std::endl;
  std::cout << "The total number of normalized background events (NuMI EXT + background nu) = " << ( background_numi_ext_normalized_events + background_numi_nu_normalized_events ) << "." << std::endl;
  std::cout << "The number of data events = " << data_normalized_events << "." << std::endl;

  std::cout << "The significance ( Signal / Sqrt( Background + Signal ) ) = " << ( signal_normalized_events / TMath::Sqrt( ( signal_normalized_events + background_numi_ext_normalized_events + background_numi_nu_normalized_events ) ) ) << "." << std::endl;  

  muon_candidate_kinetic_energy_signal->Scale( signal_normalized_events / muon_candidate_kinetic_energy_signal->Integral() );
  muon_candidate_kinetic_energy_ext_background->Scale(  background_numi_ext_normalized_events / muon_candidate_kinetic_energy_ext_background->Integral() );
  muon_candidate_kinetic_energy_nu_background->Scale( background_numi_nu_normalized_events / muon_candidate_kinetic_energy_nu_background->Integral() );
  muon_candidate_kinetic_energy_dirt->Scale( dirt_normalized_events / muon_candidate_kinetic_energy_dirt->Integral() );
  muon_candidate_kinetic_energy_data->Scale( data_normalized_events / muon_candidate_kinetic_energy_data->Integral() );
  
  // Declare new histograms that will still contain all of the information but can be manipulated properly.
  TH1D* muon_candidate_kinetic_energy_ext_background_copy                      = new TH1D("muon_candidate_kinetic_energy_ext_background_copy", "NuMI EXT Muon Candidate Kinetic_Energy", 100, 0., 2500.);
  TH1D* muon_candidate_kinetic_energy_nu_background_copy                       = new TH1D("muon_candidate_kinetic_energy_nu_background_copy", "NuMI Nu Background Muon Candidate Kinetic_Energy", 100, 0., 2500.);
  TH1D* muon_candidate_kinetic_energy_dirt_copy                                = new TH1D("muon_candidate_kinetic_energy_dirt_copy", "NuMI Dirt Muon Candidate Kinetic Energy", 100, 0., 2500.);
  TH1D* muon_candidate_kinetic_energy_signal_copy                              = new TH1D("muon_candidate_kinetic_energy_signal_copy", "Signal Muon Candidate Kinetic_Energy", 100, 0., 2500.);

  double maximum_bin_stacked_hist = -1.;
  double maximum_bin_data         = -1.;
  
  for ( size_t bin_num = 1; bin_num < 101; bin_num++ ) {

    muon_candidate_kinetic_energy_ext_background_copy->SetBinContent( bin_num, muon_candidate_kinetic_energy_ext_background->GetBinContent( bin_num ) );
    muon_candidate_kinetic_energy_nu_background_copy->SetBinContent( bin_num, muon_candidate_kinetic_energy_nu_background->GetBinContent( bin_num ) );
    muon_candidate_kinetic_energy_dirt_copy->SetBinContent( bin_num, muon_candidate_kinetic_energy_dirt->GetBinContent( bin_num ) );
    muon_candidate_kinetic_energy_signal_copy->SetBinContent( bin_num, muon_candidate_kinetic_energy_signal->GetBinContent( bin_num ) );

    double sum = ( muon_candidate_kinetic_energy_ext_background->GetBinContent( bin_num ) + muon_candidate_kinetic_energy_nu_background->GetBinContent( bin_num ) + muon_candidate_kinetic_energy_dirt->GetBinContent( bin_num ) + muon_candidate_kinetic_energy_signal->GetBinContent( bin_num ) );

    if ( sum > 0. )
      muon_candidate_kinetic_energy_ratio->SetBinContent( bin_num, muon_candidate_kinetic_energy_data->GetBinContent( bin_num ) / sum );
    
    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( muon_candidate_kinetic_energy_data->GetBinContent( bin_num ) > maximum_bin_data )
      maximum_bin_data = muon_candidate_kinetic_energy_data->GetBinContent( bin_num );
    
  }

  // Print out the integral of the KDAR signal plot.
  std::cout << "The integral of the KDAR signal plot copy = " << muon_candidate_kinetic_energy_signal_copy->Integral() << "." << std::endl;

  // Fill the "stacked" plot.
  muon_candidate_kinetic_energy_nu_background_copy->SetFillColor(kOrange);
  muon_candidate_kinetic_energy_ext_background_copy->SetFillColor(kGreen);
  muon_candidate_kinetic_energy_dirt_copy->SetFillColor(kBlue);
  muon_candidate_kinetic_energy_signal_copy->SetFillColor(kRed);
  muon_candidate_kinetic_energy_stack->Add(muon_candidate_kinetic_energy_nu_background_copy);
  muon_candidate_kinetic_energy_stack->Add(muon_candidate_kinetic_energy_ext_background_copy);
  muon_candidate_kinetic_energy_stack->Add(muon_candidate_kinetic_energy_dirt_copy);
  muon_candidate_kinetic_energy_stack->Add(muon_candidate_kinetic_energy_signal_copy);

  TCanvas* muon_candidate_kinetic_energy_canvas = new TCanvas("muon_candidate_kinetic_energy_canvas", "muon_candidate_kinetic_energy_canvas", 600, 600);
  muon_candidate_kinetic_energy_canvas->cd();
  muon_candidate_kinetic_energy_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 20.0 );
  gStyle->SetOptStat(0);
  muon_candidate_kinetic_energy_stack->Draw();
  muon_candidate_kinetic_energy_stack->SetTitle("Muon Candidate Kinetic Energy");
  muon_candidate_kinetic_energy_stack->GetXaxis()->SetTitle("Muon Candidate Kinetic Energy [MeV]");
  muon_candidate_kinetic_energy_stack->GetYaxis()->SetTitle("Events Per Bin");
  muon_candidate_kinetic_energy_stack->GetXaxis()->SetTitleOffset(1.10);
  muon_candidate_kinetic_energy_stack->GetYaxis()->SetTitleOffset(1.40);
  muon_candidate_kinetic_energy_stack->GetXaxis()->SetTickLength(0.015);
  muon_candidate_kinetic_energy_stack->GetYaxis()->SetTickLength(0.015);
  muon_candidate_kinetic_energy_stack->GetXaxis()->SetNdivisions( 8 );
  muon_candidate_kinetic_energy_stack->GetYaxis()->SetNdivisions( 10 );
  muon_candidate_kinetic_energy_stack->SetTitle("Reconstructed Muon Candidate MCS Kinetic Energy: Contained Tracks");
  gPad->Modified();
  muon_candidate_kinetic_energy_data->SetLineColor(kBlack);
  muon_candidate_kinetic_energy_data->Draw("Sames");
  
  muon_candidate_kinetic_energy_canvas->cd();
  auto* muon_candidate_kinetic_energy_legend = new TLegend(0.48, 0.48, 0.88, 0.88);
  muon_candidate_kinetic_energy_legend->AddEntry(muon_candidate_kinetic_energy_nu_background_copy,"NuMI Monte Carlo Neutrino Background","f");
  muon_candidate_kinetic_energy_legend->AddEntry(muon_candidate_kinetic_energy_ext_background_copy,"NuMI EXT","f");
  muon_candidate_kinetic_energy_legend->AddEntry(muon_candidate_kinetic_energy_dirt_copy, "NuMI Dirt", "f");
  muon_candidate_kinetic_energy_legend->AddEntry(muon_candidate_kinetic_energy_signal_copy,"NuMI Monte Carlo Signal","f");
  muon_candidate_kinetic_energy_legend->AddEntry(muon_candidate_kinetic_energy_data, "NuMI Data W/ Full Statistal Uncertainty", "pel");
  muon_candidate_kinetic_energy_legend->SetLineColor(kWhite);
  muon_candidate_kinetic_energy_legend->Draw("Sames");
  muon_candidate_kinetic_energy_canvas->Write();

  // Muon Candidate Length
  muon_candidate_length_signal->Scale( signal_normalized_events / muon_candidate_length_signal->Integral() );
  muon_candidate_length_ext_background->Scale(  background_numi_ext_normalized_events / muon_candidate_length_ext_background->Integral() );
  muon_candidate_length_dirt->Scale( dirt_normalized_events / muon_candidate_length_dirt->Integral() );
  muon_candidate_length_nu_background->Scale( background_numi_nu_normalized_events / muon_candidate_length_nu_background->Integral() );
  muon_candidate_length_data->Scale( data_normalized_events / muon_candidate_length_data->Integral() );

  TH1F*    muon_candidate_length_signal_copy                        = new TH1F("muon_candidate_length_signal_copy", "Signal Muon Candidate Length", 20, 0., 800.);
  TH1F*    muon_candidate_length_nu_background_copy                 = new TH1F("muon_candidate_length_nu_background_copy", "Nu Background Muon Candidate Length", 20, 0., 800.);
  TH1F*    muon_candidate_length_ext_background_copy                = new TH1F("muon_candidate_length_ext_background_copy", "EXT Background Muon Candidate Length", 20, 0., 800.);
  TH1F*    muon_candidate_length_dirt_copy                          = new TH1F("muon_candidate_length_dirt_copy", "Dirt Muon Candidate Length", 20, 0., 800.);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;

  for ( size_t copy_iter = 1; copy_iter < 21; copy_iter++ ) {

    muon_candidate_length_signal_copy->SetBinContent( copy_iter, muon_candidate_length_signal->GetBinContent( copy_iter ) );
    muon_candidate_length_nu_background_copy->SetBinContent( copy_iter, muon_candidate_length_nu_background->GetBinContent( copy_iter ) );
    muon_candidate_length_ext_background_copy->SetBinContent( copy_iter, muon_candidate_length_ext_background->GetBinContent( copy_iter ) );
    muon_candidate_length_dirt_copy->SetBinContent( copy_iter, muon_candidate_length_dirt->GetBinContent( copy_iter ) );

    muon_candidate_length_ratio->SetBinContent( copy_iter, 0. );

    double sum = ( muon_candidate_length_signal->GetBinContent( copy_iter ) + muon_candidate_length_nu_background->GetBinContent( copy_iter ) + muon_candidate_length_ext_background->GetBinContent( copy_iter ) + muon_candidate_length_dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      muon_candidate_length_ratio->SetBinContent( copy_iter, muon_candidate_length_data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( muon_candidate_length_data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = muon_candidate_length_data->GetBinContent( copy_iter );

  }

  // Fill the "stacked" plot.
  muon_candidate_length_nu_background_copy->SetFillColor(kOrange);
  muon_candidate_length_ext_background_copy->SetFillColor(kGreen);
  muon_candidate_length_signal_copy->SetFillColor(kRed);
  muon_candidate_length_dirt_copy->SetFillColor(kBlue);
  muon_candidate_length_stack->Add(muon_candidate_length_nu_background_copy);
  muon_candidate_length_stack->Add(muon_candidate_length_ext_background_copy);
  muon_candidate_length_stack->Add(muon_candidate_length_signal_copy);
  muon_candidate_length_stack->Add(muon_candidate_length_dirt_copy);

  TCanvas* muon_candidate_length_canvas = new TCanvas("muon_candidate_length_canvas", "muon_candidate_length_canvas", 600, 600);
  muon_candidate_length_canvas->cd();
  muon_candidate_length_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 50.0 );
  gStyle->SetOptStat(0);
  muon_candidate_length_stack->Draw();
  muon_candidate_length_stack->GetXaxis()->SetTitle("Muon Candidate Length [cm]");
  muon_candidate_length_stack->GetYaxis()->SetTitle("Events Per Bin");
  muon_candidate_length_stack->GetXaxis()->SetTitleOffset(1.10);
  muon_candidate_length_stack->GetYaxis()->SetTitleOffset(1.48);
  muon_candidate_length_stack->GetXaxis()->SetTickLength(0.015);
  muon_candidate_length_stack->GetYaxis()->SetTickLength(0.015);
  muon_candidate_length_stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  muon_candidate_length_data->SetLineColor(kBlack);
  muon_candidate_length_data->Draw("SamesE1");
  muon_candidate_length_canvas->Write();

  // Flash Y Difference
  flash_y_difference_signal->Scale( signal_normalized_events / flash_y_difference_signal->Integral() );
  flash_y_difference_ext_background->Scale(  background_numi_ext_normalized_events / flash_y_difference_ext_background->Integral() );
  flash_y_difference_dirt->Scale( dirt_normalized_events / flash_y_difference_dirt->Integral() );
  flash_y_difference_nu_background->Scale( background_numi_nu_normalized_events / flash_y_difference_nu_background->Integral() );
  flash_y_difference_data->Scale( data_normalized_events / flash_y_difference_data->Integral() );

  TH1F*    flash_y_difference_signal_copy                        = new TH1F("flash_y_difference_signal_copy", "Signal Flash Y Difference", 25, -150.0, 150.0);
  TH1F*    flash_y_difference_nu_background_copy                 = new TH1F("flash_y_difference_nu_background_copy", "Nu Background Flash Y Difference", 25, -150.0, 150.0);
  TH1F*    flash_y_difference_ext_background_copy                = new TH1F("flash_y_difference_ext_background_copy", "EXT Background Flash Y Difference", 25, -150.0, 150.0);
  TH1F*    flash_y_difference_dirt_copy                          = new TH1F("flash_y_difference_dirt_copy", "Dirt Flash Y Difference", 25, -150.0, 150.0);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
  
  for ( size_t copy_iter = 1; copy_iter < 26; copy_iter++ ) {

    flash_y_difference_signal_copy->SetBinContent( copy_iter, flash_y_difference_signal->GetBinContent( copy_iter ) );
    flash_y_difference_nu_background_copy->SetBinContent( copy_iter, flash_y_difference_nu_background->GetBinContent( copy_iter ) );
    flash_y_difference_ext_background_copy->SetBinContent( copy_iter, flash_y_difference_ext_background->GetBinContent( copy_iter ) );
    flash_y_difference_dirt_copy->SetBinContent( copy_iter, flash_y_difference_dirt->GetBinContent( copy_iter ) );

    flash_y_difference_ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( flash_y_difference_signal->GetBinContent( copy_iter ) + flash_y_difference_nu_background->GetBinContent( copy_iter ) + flash_y_difference_ext_background->GetBinContent( copy_iter ) + flash_y_difference_dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      flash_y_difference_ratio->SetBinContent( copy_iter, flash_y_difference_data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( flash_y_difference_data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = flash_y_difference_data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  flash_y_difference_nu_background_copy->SetFillColor(kOrange);
  flash_y_difference_ext_background_copy->SetFillColor(kGreen);
  flash_y_difference_signal_copy->SetFillColor(kRed);
  flash_y_difference_dirt_copy->SetFillColor(kBlue);
  flash_y_difference_stack->Add(flash_y_difference_nu_background_copy);
  flash_y_difference_stack->Add(flash_y_difference_ext_background_copy);
  flash_y_difference_stack->Add(flash_y_difference_signal_copy);
  flash_y_difference_stack->Add(flash_y_difference_dirt_copy);
  
  TCanvas* flash_y_difference_canvas = new TCanvas("flash_y_difference_canvas", "flash_y_difference_canvas", 600, 600);
  flash_y_difference_canvas->cd();
  flash_y_difference_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 50.0 );
  gStyle->SetOptStat(0);
  flash_y_difference_stack->Draw();
  flash_y_difference_stack->GetXaxis()->SetTitle("Flash Y Difference [cm]");
  flash_y_difference_stack->GetYaxis()->SetTitle("Events Per Bin");
  flash_y_difference_stack->GetXaxis()->SetTitleOffset(1.10);
  flash_y_difference_stack->GetYaxis()->SetTitleOffset(1.48);
  flash_y_difference_stack->GetXaxis()->SetTickLength(0.015);
  flash_y_difference_stack->GetYaxis()->SetTickLength(0.015);
  flash_y_difference_stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  flash_y_difference_data->SetLineColor(kBlack);
  flash_y_difference_data->Draw("SamesE1");
  flash_y_difference_canvas->Write();

  // Flash Z Difference
  flash_z_difference_signal->Scale( signal_normalized_events / flash_z_difference_signal->Integral() );
  flash_z_difference_ext_background->Scale(  background_numi_ext_normalized_events / flash_z_difference_ext_background->Integral() );
  flash_z_difference_nu_background->Scale( background_numi_nu_normalized_events / flash_z_difference_nu_background->Integral() );
  flash_z_difference_dirt->Scale( dirt_normalized_events / flash_z_difference_dirt->Integral() );
  flash_z_difference_data->Scale( data_normalized_events / flash_z_difference_data->Integral() );

  TH1F*    flash_z_difference_signal_copy                        = new TH1F("flash_z_difference_signal_copy", "Signal Flash Z Difference", 20, -200.0, 200.0);
  TH1F*    flash_z_difference_nu_background_copy                 = new TH1F("flash_z_difference_nu_background_copy", "Nu Background Flash Z Difference", 20, -200.0, 200.0);
  TH1F*    flash_z_difference_dirt_copy                          = new TH1F("flash_z_difference_dirt_copy", "Dirt Flash Z Difference", 20, -200.0, 200.0);
  TH1F*    flash_z_difference_ext_background_copy                = new TH1F("flash_z_difference_ext_background_copy", "EXT Background Flash Z Difference", 20, -200.0, 200.0);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
  
  for ( size_t copy_iter = 1; copy_iter < 21; copy_iter++ ) {

    flash_z_difference_signal_copy->SetBinContent( copy_iter, flash_z_difference_signal->GetBinContent( copy_iter ) );
    flash_z_difference_nu_background_copy->SetBinContent( copy_iter, flash_z_difference_nu_background->GetBinContent( copy_iter ) );
    flash_z_difference_ext_background_copy->SetBinContent( copy_iter, flash_z_difference_ext_background->GetBinContent( copy_iter ) );
    flash_z_difference_dirt_copy->SetBinContent( copy_iter, flash_z_difference_dirt->GetBinContent( copy_iter ) );

    flash_z_difference_ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( flash_z_difference_signal->GetBinContent( copy_iter ) + flash_z_difference_nu_background->GetBinContent( copy_iter ) + flash_z_difference_ext_background->GetBinContent( copy_iter ) + flash_z_difference_dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      flash_z_difference_ratio->SetBinContent( copy_iter, flash_z_difference_data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( flash_z_difference_data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = flash_z_difference_data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  flash_z_difference_nu_background_copy->SetFillColor(kOrange);
  flash_z_difference_ext_background_copy->SetFillColor(kGreen);
  flash_z_difference_signal_copy->SetFillColor(kRed);
  flash_z_difference_dirt_copy->SetFillColor(kBlue);
  flash_z_difference_stack->Add(flash_z_difference_nu_background_copy);
  flash_z_difference_stack->Add(flash_z_difference_ext_background_copy);
  flash_z_difference_stack->Add(flash_z_difference_signal_copy);
  flash_z_difference_stack->Add(flash_z_difference_dirt_copy);
  
  TCanvas* flash_z_difference_canvas = new TCanvas("flash_z_difference_canvas", "flash_z_difference_canvas", 600, 600);
  flash_z_difference_canvas->cd();
  flash_z_difference_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 50.0 );
  gStyle->SetOptStat(0);
  flash_z_difference_stack->Draw();
  flash_z_difference_stack->GetXaxis()->SetTitle("Flash Z Difference [cm]");
  flash_z_difference_stack->GetYaxis()->SetTitle("Events Per Bin");
  flash_z_difference_stack->GetXaxis()->SetTitleOffset(1.10);
  flash_z_difference_stack->GetYaxis()->SetTitleOffset(1.48);
  flash_z_difference_stack->GetXaxis()->SetTickLength(0.015);
  flash_z_difference_stack->GetYaxis()->SetTickLength(0.015);
  flash_z_difference_stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  flash_z_difference_data->SetLineColor(kBlack);
  flash_z_difference_data->Draw("SamesE1");
  flash_z_difference_canvas->Write();
  
  // Reco Vtx X
  reco_vtx_x_signal->Scale( signal_normalized_events / reco_vtx_x_signal->Integral() );
  reco_vtx_x_ext_background->Scale(  background_numi_ext_normalized_events / reco_vtx_x_ext_background->Integral() );
  reco_vtx_x_nu_background->Scale( background_numi_nu_normalized_events / reco_vtx_x_nu_background->Integral() );
  reco_vtx_x_dirt->Scale( dirt_normalized_events / reco_vtx_x_dirt->Integral() );
  reco_vtx_x_data->Scale( data_normalized_events / reco_vtx_x_data->Integral() );

  TH1F*    reco_vtx_x_signal_copy                        = new TH1F("reco_vtx_x_signal_copy", "Signal Reco Vtx X", 25, 20., 236.35);
  TH1F*    reco_vtx_x_nu_background_copy                 = new TH1F("reco_vtx_x_nu_background_copy", "Nu Background Reco Vtx X", 25, 20., 236.35);
  TH1F*    reco_vtx_x_ext_background_copy                = new TH1F("reco_vtx_x_ext_background_copy", "EXT Background Reco Vtx X", 25, 20., 236.35);
  TH1F*    reco_vtx_x_dirt_copy                          = new TH1F("reco_vtx_x_dirt_copy", "Dirt Reco Vtx X", 25, 20., 236.35);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
  
  for ( size_t copy_iter = 1; copy_iter < 26; copy_iter++ ) {

    reco_vtx_x_signal_copy->SetBinContent( copy_iter, reco_vtx_x_signal->GetBinContent( copy_iter ) );
    reco_vtx_x_nu_background_copy->SetBinContent( copy_iter, reco_vtx_x_nu_background->GetBinContent( copy_iter ) );
    reco_vtx_x_ext_background_copy->SetBinContent( copy_iter, reco_vtx_x_ext_background->GetBinContent( copy_iter ) );
    reco_vtx_x_dirt_copy->SetBinContent( copy_iter, reco_vtx_x_dirt->GetBinContent( copy_iter ) );

    reco_vtx_x_ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( reco_vtx_x_signal->GetBinContent( copy_iter ) + reco_vtx_x_nu_background->GetBinContent( copy_iter ) + reco_vtx_x_ext_background->GetBinContent( copy_iter ) + reco_vtx_x_dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      reco_vtx_x_ratio->SetBinContent( copy_iter, reco_vtx_x_data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( reco_vtx_x_data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = reco_vtx_x_data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  reco_vtx_x_nu_background_copy->SetFillColor(kOrange);
  reco_vtx_x_ext_background_copy->SetFillColor(kGreen);
  reco_vtx_x_signal_copy->SetFillColor(kRed);
  reco_vtx_x_dirt_copy->SetFillColor(kBlue);
  reco_vtx_x_stack->Add(reco_vtx_x_nu_background_copy);
  reco_vtx_x_stack->Add(reco_vtx_x_ext_background_copy);
  reco_vtx_x_stack->Add(reco_vtx_x_signal_copy);
  reco_vtx_x_stack->Add(reco_vtx_x_dirt_copy);

  TCanvas* reco_vtx_x_canvas = new TCanvas("reco_vtx_x_canvas", "reco_vtx_x_canvas", 600, 600);
  reco_vtx_x_canvas->cd();
  reco_vtx_x_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 20.0 );
  gStyle->SetOptStat(0);
  reco_vtx_x_stack->Draw();
  reco_vtx_x_stack->SetTitle("Reconstructed Vertex #it{x}-Coordinate");
  reco_vtx_x_stack->GetXaxis()->SetTitle("Reco Vtx X [cm]");
  reco_vtx_x_stack->GetYaxis()->SetTitle("Events Per Bin");
  reco_vtx_x_stack->GetXaxis()->SetTitleOffset(1.10);
  reco_vtx_x_stack->GetYaxis()->SetTitleOffset(1.40);
  reco_vtx_x_stack->GetXaxis()->SetTickLength(0.015);
  reco_vtx_x_stack->GetYaxis()->SetTickLength(0.015);
  reco_vtx_x_stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  reco_vtx_x_data->SetLineColor(kBlack);
  reco_vtx_x_data->Draw("SamesE1");
  reco_vtx_x_canvas->Write();

  // Reco Vtx Y
  reco_vtx_y_signal->Scale( signal_normalized_events / reco_vtx_y_signal->Integral() );
  reco_vtx_y_ext_background->Scale(  background_numi_ext_normalized_events / reco_vtx_y_ext_background->Integral() );
  reco_vtx_y_nu_background->Scale( background_numi_nu_normalized_events / reco_vtx_y_nu_background->Integral() );
  reco_vtx_y_dirt->Scale( dirt_normalized_events / reco_vtx_y_dirt->Integral() );
  reco_vtx_y_data->Scale( data_normalized_events / reco_vtx_y_data->Integral() );

  TH1F*    reco_vtx_y_signal_copy                        = new TH1F("reco_vtx_y_signal_copy", "Signal Reco Vtx Y", 25, -96.5, 96.5);
  TH1F*    reco_vtx_y_nu_background_copy                 = new TH1F("reco_vtx_y_nu_background_copy", "Nu Background Reco Vtx Y", 25, -96.5, 96.5);
  TH1F*    reco_vtx_y_ext_background_copy                = new TH1F("reco_vtx_y_ext_background_copy", "EXT Background Reco Vtx Y", 25, -96.5, 96.5);
  TH1F*    reco_vtx_y_dirt_copy                          = new TH1F("reco_vtx_y_dirt_copy", "Dirt Reco Vtx Y", 25, -96.5, 96.5);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
 
  for ( size_t copy_iter = 1; copy_iter < 26; copy_iter++ ) {

    reco_vtx_y_signal_copy->SetBinContent( copy_iter, reco_vtx_y_signal->GetBinContent( copy_iter ) );
    reco_vtx_y_nu_background_copy->SetBinContent( copy_iter, reco_vtx_y_nu_background->GetBinContent( copy_iter ) );
    reco_vtx_y_ext_background_copy->SetBinContent( copy_iter, reco_vtx_y_ext_background->GetBinContent( copy_iter ) );
    reco_vtx_y_dirt_copy->SetBinContent( copy_iter, reco_vtx_y_dirt->GetBinContent( copy_iter ) );

    reco_vtx_y_ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( reco_vtx_y_signal->GetBinContent( copy_iter ) + reco_vtx_y_nu_background->GetBinContent( copy_iter ) + reco_vtx_y_ext_background->GetBinContent( copy_iter ) + reco_vtx_y_dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      reco_vtx_y_ratio->SetBinContent( copy_iter, reco_vtx_y_data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( reco_vtx_y_data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = reco_vtx_y_data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  reco_vtx_y_nu_background_copy->SetFillColor(kOrange);
  reco_vtx_y_ext_background_copy->SetFillColor(kGreen);
  reco_vtx_y_signal_copy->SetFillColor(kRed);
  reco_vtx_y_dirt_copy->SetFillColor(kBlue);
  reco_vtx_y_stack->Add(reco_vtx_y_nu_background_copy);
  reco_vtx_y_stack->Add(reco_vtx_y_ext_background_copy);
  reco_vtx_y_stack->Add(reco_vtx_y_signal_copy);
  reco_vtx_y_stack->Add(reco_vtx_y_dirt_copy);

  TCanvas* reco_vtx_y_canvas = new TCanvas("reco_vtx_y_canvas", "reco_vtx_y_canvas", 600, 600);
  reco_vtx_y_canvas->cd();
  reco_vtx_y_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 20.0 );
  gStyle->SetOptStat(0);
  reco_vtx_y_stack->Draw();
  reco_vtx_y_stack->SetTitle("Reconstructed Vertex #it{y}-Coordinate");
  reco_vtx_y_stack->GetXaxis()->SetTitle("Reco Vtx Y [cm]");
  reco_vtx_y_stack->GetYaxis()->SetTitle("Events Per Bin");
  reco_vtx_y_stack->GetXaxis()->SetTitleOffset(1.10);
  reco_vtx_y_stack->GetYaxis()->SetTitleOffset(1.40);
  reco_vtx_y_stack->GetXaxis()->SetTickLength(0.015);
  reco_vtx_y_stack->GetYaxis()->SetTickLength(0.015);
  reco_vtx_y_stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  reco_vtx_y_data->SetLineColor(kBlack);
  reco_vtx_y_data->Draw("SamesE1");
  reco_vtx_y_canvas->Write();

  // Reco Vtx Z
  reco_vtx_z_signal->Scale( signal_normalized_events / reco_vtx_z_signal->Integral() );
  reco_vtx_z_ext_background->Scale(  background_numi_ext_normalized_events / reco_vtx_z_ext_background->Integral() );
  reco_vtx_z_nu_background->Scale( background_numi_nu_normalized_events / reco_vtx_z_nu_background->Integral() );
  reco_vtx_z_dirt->Scale( dirt_normalized_events / reco_vtx_z_dirt->Integral() );
  reco_vtx_z_data->Scale( data_normalized_events / reco_vtx_z_data->Integral() );

  TH1F*    reco_vtx_z_signal_copy                        = new TH1F("reco_vtx_z_signal_copy", "Signal Reco Vtx Z", 50, 20., 1016.8);
  TH1F*    reco_vtx_z_nu_background_copy                 = new TH1F("reco_vtx_z_nu_background_copy", "Nu Background Reco Vtx Z", 50, 20., 1016.8);
  TH1F*    reco_vtx_z_ext_background_copy                = new TH1F("reco_vtx_z_ext_background_copy", "EXT Background Reco Vtx Z", 50, 20., 1016.8);
  TH1F*    reco_vtx_z_dirt_copy                          = new TH1F("reco_vtx_z_dirt_copy", "Dirt Reco Vtx z", 50, 20., 1016.8);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
  
  for ( size_t copy_iter = 1; copy_iter < 51; copy_iter++ ) {

    reco_vtx_z_signal_copy->SetBinContent( copy_iter, reco_vtx_z_signal->GetBinContent( copy_iter ) );
    reco_vtx_z_nu_background_copy->SetBinContent( copy_iter, reco_vtx_z_nu_background->GetBinContent( copy_iter ) );
    reco_vtx_z_ext_background_copy->SetBinContent( copy_iter, reco_vtx_z_ext_background->GetBinContent( copy_iter ) );
    reco_vtx_z_dirt_copy->SetBinContent( copy_iter, reco_vtx_z_dirt->GetBinContent( copy_iter ) );

    reco_vtx_z_ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( reco_vtx_z_signal->GetBinContent( copy_iter ) + reco_vtx_z_nu_background->GetBinContent( copy_iter ) + reco_vtx_z_ext_background->GetBinContent( copy_iter ) + reco_vtx_z_dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      reco_vtx_z_ratio->SetBinContent( copy_iter, reco_vtx_z_data->GetBinContent( copy_iter ) / sum );
    
    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( reco_vtx_z_data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = reco_vtx_z_data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  reco_vtx_z_nu_background_copy->SetFillColor(kOrange);
  reco_vtx_z_ext_background_copy->SetFillColor(kGreen);
  reco_vtx_z_signal_copy->SetFillColor(kRed);
  reco_vtx_z_dirt_copy->SetFillColor(kBlue);
  reco_vtx_z_stack->Add(reco_vtx_z_nu_background_copy);
  reco_vtx_z_stack->Add(reco_vtx_z_ext_background_copy);
  reco_vtx_z_stack->Add(reco_vtx_z_signal_copy);
  reco_vtx_z_stack->Add(reco_vtx_z_dirt_copy);

  TCanvas* reco_vtx_z_canvas = new TCanvas("reco_vtx_z_canvas", "reco_vtx_z_canvas", 600, 600);
  reco_vtx_z_canvas->cd();
  reco_vtx_z_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 20.0 );
  gStyle->SetOptStat(0);
  reco_vtx_z_stack->Draw();
  reco_vtx_z_stack->SetTitle("Reconstructed Vertex #it{z}-Coordinate");
  reco_vtx_z_stack->GetXaxis()->SetTitle("Reco Vtx Z [cm]");
  reco_vtx_z_stack->GetYaxis()->SetTitle("Events Per Bin");
  reco_vtx_z_stack->GetXaxis()->SetTitleOffset(1.10);
  reco_vtx_z_stack->GetYaxis()->SetTitleOffset(1.40);
  reco_vtx_z_stack->GetXaxis()->SetTickLength(0.015);
  reco_vtx_z_stack->GetYaxis()->SetTickLength(0.015);
  reco_vtx_z_stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  reco_vtx_z_data->SetLineColor(kBlack);
  reco_vtx_z_data->Draw("SamesE1");
  reco_vtx_z_canvas->Write();

  // Muon_Candidate_Proton_PID_
  Muon_Candidate_Proton_PID__signal->Scale( signal_normalized_events / Muon_Candidate_Proton_PID__signal->Integral() );
  Muon_Candidate_Proton_PID__ext_background->Scale(  background_numi_ext_normalized_events / Muon_Candidate_Proton_PID__ext_background->Integral() );
  Muon_Candidate_Proton_PID__nu_background->Scale( background_numi_nu_normalized_events / Muon_Candidate_Proton_PID__nu_background->Integral() );
  Muon_Candidate_Proton_PID__data->Scale( data_normalized_events / Muon_Candidate_Proton_PID__data->Integral() );
  Muon_Candidate_Proton_PID__dirt->Scale( dirt_normalized_events / Muon_Candidate_Proton_PID__dirt->Integral() );
  
  TH1F*    Muon_Candidate_Proton_PID__signal_copy                        = new TH1F("Muon_Candidate_Proton_PID__signal_copy", "Signal Muon Candidate Proton PID", 40, 0., 500.);
  TH1F*    Muon_Candidate_Proton_PID__nu_background_copy                 = new TH1F("Muon_Candidate_Proton_PID__nu_background_copy", "Nu Background Muon Candidate Proton PID", 40, 0., 500.);
  TH1F*    Muon_Candidate_Proton_PID__ext_background_copy                = new TH1F("Muon_Candidate_Proton_PID__ext_background_copy", "EXT Background Muon Candidate Proton PID", 40, 0., 500.);
  TH1F*    Muon_Candidate_Proton_PID__dirt_copy                          = new TH1F("Muon_Candidate_Proton_PID__dirt_copy", "Dirt Muon Candidate Proton PID", 40, 0., 500.);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
  
  for ( size_t copy_iter = 1; copy_iter < 41; copy_iter++ ) {

    Muon_Candidate_Proton_PID__signal_copy->SetBinContent( copy_iter, Muon_Candidate_Proton_PID__signal->GetBinContent( copy_iter ) );
    Muon_Candidate_Proton_PID__nu_background_copy->SetBinContent( copy_iter, Muon_Candidate_Proton_PID__nu_background->GetBinContent( copy_iter ) );
    Muon_Candidate_Proton_PID__ext_background_copy->SetBinContent( copy_iter, Muon_Candidate_Proton_PID__ext_background->GetBinContent( copy_iter ) );
    Muon_Candidate_Proton_PID__dirt_copy->SetBinContent( copy_iter, Muon_Candidate_Proton_PID__dirt->GetBinContent( copy_iter ) );

    Muon_Candidate_Proton_PID__ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( Muon_Candidate_Proton_PID__signal->GetBinContent( copy_iter ) + Muon_Candidate_Proton_PID__nu_background->GetBinContent( copy_iter ) + Muon_Candidate_Proton_PID__ext_background->GetBinContent( copy_iter ) + Muon_Candidate_Proton_PID__dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      Muon_Candidate_Proton_PID__ratio->SetBinContent( copy_iter, Muon_Candidate_Proton_PID__data->GetBinContent( copy_iter ) / sum );
    
    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( Muon_Candidate_Proton_PID__data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = Muon_Candidate_Proton_PID__data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  Muon_Candidate_Proton_PID__nu_background_copy->SetFillColor(kOrange);
  Muon_Candidate_Proton_PID__ext_background_copy->SetFillColor(kGreen);
  Muon_Candidate_Proton_PID__signal_copy->SetFillColor(kRed);
  Muon_Candidate_Proton_PID__dirt_copy->SetFillColor(kBlue);
  Muon_Candidate_Proton_PID__stack->Add(Muon_Candidate_Proton_PID__nu_background_copy);
  Muon_Candidate_Proton_PID__stack->Add(Muon_Candidate_Proton_PID__ext_background_copy);
  Muon_Candidate_Proton_PID__stack->Add(Muon_Candidate_Proton_PID__signal_copy);
  Muon_Candidate_Proton_PID__stack->Add(Muon_Candidate_Proton_PID__dirt_copy);
  
  TCanvas* Muon_Candidate_Proton_PID__canvas = new TCanvas("Muon_Candidate_Proton_PID__canvas", "Muon_Candidate_Proton_PID__canvas", 600, 600);
  Muon_Candidate_Proton_PID__canvas->cd();
  Muon_Candidate_Proton_PID__stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 50.0 );
  gStyle->SetOptStat(0);
  Muon_Candidate_Proton_PID__stack->Draw();
  Muon_Candidate_Proton_PID__stack->SetTitle("Muon Candidate Proton PID Value");
  Muon_Candidate_Proton_PID__stack->GetXaxis()->SetTitle("Muon Candidate Proton PID [Unitless]");
  Muon_Candidate_Proton_PID__stack->GetYaxis()->SetTitle("Events Per Bin");
  Muon_Candidate_Proton_PID__stack->GetXaxis()->SetTitleOffset(1.10);
  Muon_Candidate_Proton_PID__stack->GetYaxis()->SetTitleOffset(1.48);
  Muon_Candidate_Proton_PID__stack->GetXaxis()->SetTickLength(0.015);
  Muon_Candidate_Proton_PID__stack->GetYaxis()->SetTickLength(0.015);
  Muon_Candidate_Proton_PID__stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  Muon_Candidate_Proton_PID__data->SetLineColor(kBlack);
  Muon_Candidate_Proton_PID__data->Draw("SamesE1");
  Muon_Candidate_Proton_PID__canvas->Write();

  // Muon_Candidate_Muon_Track_Score_
  Muon_Candidate_Muon_Track_Score__signal->Scale( signal_normalized_events / Muon_Candidate_Muon_Track_Score__signal->Integral() );
  Muon_Candidate_Muon_Track_Score__ext_background->Scale(  background_numi_ext_normalized_events / Muon_Candidate_Muon_Track_Score__ext_background->Integral() );
  Muon_Candidate_Muon_Track_Score__nu_background->Scale( background_numi_nu_normalized_events / Muon_Candidate_Muon_Track_Score__nu_background->Integral() );
  Muon_Candidate_Muon_Track_Score__data->Scale( data_normalized_events / Muon_Candidate_Muon_Track_Score__data->Integral() );
  Muon_Candidate_Muon_Track_Score__dirt->Scale( dirt_normalized_events / Muon_Candidate_Muon_Track_Score__dirt->Integral() );
  
  TH1F*    Muon_Candidate_Muon_Track_Score__signal_copy                        = new TH1F("Muon_Candidate_Muon_Track_Score__signal_copy", "Signal Muon Candidate Track Score", 20, 0., 1.);
  TH1F*    Muon_Candidate_Muon_Track_Score__nu_background_copy                 = new TH1F("Muon_Candidate_Muon_Track_Score__nu_background_copy", "Nu Background Muon Candidate Track Score", 20, 0., 1.);
  TH1F*    Muon_Candidate_Muon_Track_Score__ext_background_copy                = new TH1F("Muon_Candidate_Muon_Track_Score__ext_background_copy", "EXT Background Muon Candidate Track Score", 20, 0., 1.);
  TH1F*    Muon_Candidate_Muon_Track_Score__dirt_copy                          = new TH1F("Muon_Candidate_Muon_Track_Score__dirt_copy", "Dirt Muon Candidate Track Score", 20, 0., 1.);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
  
  for ( size_t copy_iter = 1; copy_iter < 21; copy_iter++ ) {

    Muon_Candidate_Muon_Track_Score__signal_copy->SetBinContent( copy_iter, Muon_Candidate_Muon_Track_Score__signal->GetBinContent( copy_iter ) );
    Muon_Candidate_Muon_Track_Score__nu_background_copy->SetBinContent( copy_iter, Muon_Candidate_Muon_Track_Score__nu_background->GetBinContent( copy_iter ) );
    Muon_Candidate_Muon_Track_Score__ext_background_copy->SetBinContent( copy_iter, Muon_Candidate_Muon_Track_Score__ext_background->GetBinContent( copy_iter ) );
    Muon_Candidate_Muon_Track_Score__dirt_copy->SetBinContent( copy_iter, Muon_Candidate_Muon_Track_Score__dirt->GetBinContent( copy_iter ) );

    Muon_Candidate_Muon_Track_Score__ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( Muon_Candidate_Muon_Track_Score__signal->GetBinContent( copy_iter ) + Muon_Candidate_Muon_Track_Score__nu_background->GetBinContent( copy_iter ) + Muon_Candidate_Muon_Track_Score__ext_background->GetBinContent( copy_iter ) + Muon_Candidate_Muon_Track_Score__dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      Muon_Candidate_Muon_Track_Score__ratio->SetBinContent( copy_iter, Muon_Candidate_Muon_Track_Score__data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( Muon_Candidate_Muon_Track_Score__data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = Muon_Candidate_Muon_Track_Score__data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  Muon_Candidate_Muon_Track_Score__nu_background_copy->SetFillColor(kOrange);
  Muon_Candidate_Muon_Track_Score__ext_background_copy->SetFillColor(kGreen);
  Muon_Candidate_Muon_Track_Score__signal_copy->SetFillColor(kRed);
  Muon_Candidate_Muon_Track_Score__dirt_copy->SetFillColor(kBlue);
  Muon_Candidate_Muon_Track_Score__stack->Add(Muon_Candidate_Muon_Track_Score__nu_background_copy);
  Muon_Candidate_Muon_Track_Score__stack->Add(Muon_Candidate_Muon_Track_Score__ext_background_copy);
  Muon_Candidate_Muon_Track_Score__stack->Add(Muon_Candidate_Muon_Track_Score__signal_copy);
  Muon_Candidate_Muon_Track_Score__stack->Add(Muon_Candidate_Muon_Track_Score__dirt_copy);
  
  TCanvas* Muon_Candidate_Muon_Track_Score__canvas = new TCanvas("Muon_Candidate_Muon_Track_Score__canvas", "Muon_Candidate_Muon_Track_Score__canvas", 600, 600);
  Muon_Candidate_Muon_Track_Score__canvas->cd();
  Muon_Candidate_Muon_Track_Score__stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 50.0 );
  gStyle->SetOptStat(0);
  Muon_Candidate_Muon_Track_Score__stack->Draw();
  Muon_Candidate_Muon_Track_Score__stack->SetTitle("Muon Candidate Muon Track Score");
  Muon_Candidate_Muon_Track_Score__stack->GetXaxis()->SetTitle("Muon Candidate Track Score [Unitless]");
  Muon_Candidate_Muon_Track_Score__stack->GetYaxis()->SetTitle("Events Per Bin");
  Muon_Candidate_Muon_Track_Score__stack->GetXaxis()->SetTitleOffset(1.10);
  Muon_Candidate_Muon_Track_Score__stack->GetYaxis()->SetTitleOffset(1.48);
  Muon_Candidate_Muon_Track_Score__stack->GetXaxis()->SetTickLength(0.015);
  Muon_Candidate_Muon_Track_Score__stack->GetYaxis()->SetTickLength(0.015);
  Muon_Candidate_Muon_Track_Score__stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  Muon_Candidate_Muon_Track_Score__data->SetLineColor(kBlack);
  Muon_Candidate_Muon_Track_Score__data->Draw("SamesE1");
  Muon_Candidate_Muon_Track_Score__canvas->Write();

  // Muon_Candidate_NuScore_
  Muon_Candidate_NuScore__signal->Scale( signal_normalized_events / Muon_Candidate_NuScore__signal->Integral() );
  Muon_Candidate_NuScore__ext_background->Scale(  background_numi_ext_normalized_events / Muon_Candidate_NuScore__ext_background->Integral() );
  Muon_Candidate_NuScore__nu_background->Scale( background_numi_nu_normalized_events / Muon_Candidate_NuScore__nu_background->Integral() );
  Muon_Candidate_NuScore__data->Scale( data_normalized_events / Muon_Candidate_NuScore__data->Integral() );
  Muon_Candidate_NuScore__dirt->Scale( dirt_normalized_events / Muon_Candidate_NuScore__dirt->Integral() );

  TH1F*    Muon_Candidate_NuScore__signal_copy                        = new TH1F("Muon_Candidate_NuScore__signal_copy", "Signal Muon Candidate NuScore", 20, 0., 1.);
  TH1F*    Muon_Candidate_NuScore__nu_background_copy                 = new TH1F("Muon_Candidate_NuScore__nu_background_copy", "Nu Background Muon Candidate NuScore", 20, 0., 1.);
  TH1F*    Muon_Candidate_NuScore__ext_background_copy                = new TH1F("Muon_Candidate_NuScore__ext_background_copy", "EXT Background Muon Candidate NuScore", 20, 0., 1.);
  TH1F*    Muon_Candidate_NuScore__dirt_copy                          = new TH1F("Muon_Candidate_NuScore__dirt_copy", "Dirt Muon Candidate NuScore", 20, 0., 1.);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;
  
  for ( size_t copy_iter = 1; copy_iter < 21; copy_iter++ ) {

    Muon_Candidate_NuScore__signal_copy->SetBinContent( copy_iter, Muon_Candidate_NuScore__signal->GetBinContent( copy_iter ) );
    Muon_Candidate_NuScore__nu_background_copy->SetBinContent( copy_iter, Muon_Candidate_NuScore__nu_background->GetBinContent( copy_iter ) );
    Muon_Candidate_NuScore__ext_background_copy->SetBinContent( copy_iter, Muon_Candidate_NuScore__ext_background->GetBinContent( copy_iter ) );
    Muon_Candidate_NuScore__dirt_copy->SetBinContent( copy_iter, Muon_Candidate_NuScore__dirt->GetBinContent( copy_iter ) );

    Muon_Candidate_NuScore__ratio->SetBinContent( copy_iter, 0. );
    
    double sum = ( Muon_Candidate_NuScore__signal->GetBinContent( copy_iter ) + Muon_Candidate_NuScore__nu_background->GetBinContent( copy_iter ) + Muon_Candidate_NuScore__ext_background->GetBinContent( copy_iter ) + Muon_Candidate_NuScore__dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      Muon_Candidate_NuScore__ratio->SetBinContent( copy_iter, Muon_Candidate_NuScore__data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( Muon_Candidate_NuScore__data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = Muon_Candidate_NuScore__data->GetBinContent( copy_iter );
    
  }

  // Fill the "stacked" plot.
  Muon_Candidate_NuScore__nu_background_copy->SetFillColor(kOrange);
  Muon_Candidate_NuScore__ext_background_copy->SetFillColor(kGreen);
  Muon_Candidate_NuScore__signal_copy->SetFillColor(kRed);
  Muon_Candidate_NuScore__dirt_copy->SetFillColor(kBlue);
  Muon_Candidate_NuScore__stack->Add(Muon_Candidate_NuScore__nu_background_copy);
  Muon_Candidate_NuScore__stack->Add(Muon_Candidate_NuScore__ext_background_copy);
  Muon_Candidate_NuScore__stack->Add(Muon_Candidate_NuScore__signal_copy);
  Muon_Candidate_NuScore__stack->Add(Muon_Candidate_NuScore__dirt_copy);

  TCanvas* Muon_Candidate_NuScore__canvas = new TCanvas("Muon_Candidate_NuScore__canvas", "Muon_Candidate_NuScore__canvas", 600, 600);
  Muon_Candidate_NuScore__canvas->cd();
  Muon_Candidate_NuScore__stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 50.0 );
  gStyle->SetOptStat(0);
  Muon_Candidate_NuScore__stack->Draw();
  Muon_Candidate_NuScore__stack->SetTitle("Muon Candidate NuScore");
  Muon_Candidate_NuScore__stack->GetXaxis()->SetTitle("Muon Candidate NuScore [Unitless]");
  Muon_Candidate_NuScore__stack->GetYaxis()->SetTitle("Events Per Bin");
  Muon_Candidate_NuScore__stack->GetXaxis()->SetTitleOffset(1.10);
  Muon_Candidate_NuScore__stack->GetYaxis()->SetTitleOffset(1.48);
  Muon_Candidate_NuScore__stack->GetXaxis()->SetTickLength(0.015);
  Muon_Candidate_NuScore__stack->GetYaxis()->SetTickLength(0.015);
  Muon_Candidate_NuScore__stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  Muon_Candidate_NuScore__data->SetLineColor(kBlack);
  Muon_Candidate_NuScore__data->Draw("SamesE1");
  Muon_Candidate_NuScore__canvas->Write();

  // Muon Directional Cosine.
  muon_candidate_directional_cosine_signal->Scale( signal_normalized_events / muon_candidate_directional_cosine_signal->Integral() );
  muon_candidate_directional_cosine_ext_background->Scale(  background_numi_ext_normalized_events / muon_candidate_directional_cosine_ext_background->Integral() );
  muon_candidate_directional_cosine_dirt->Scale( dirt_normalized_events / muon_candidate_directional_cosine_dirt->Integral() );
  muon_candidate_directional_cosine_nu_background->Scale( background_numi_nu_normalized_events / muon_candidate_directional_cosine_nu_background->Integral() );
  muon_candidate_directional_cosine_data->Scale( data_normalized_events / muon_candidate_directional_cosine_data->Integral() );

  TH1F*    muon_candidate_directional_cosine_signal_copy                        = new TH1F("muon_candidate_directional_cosine_signal_copy", "Signal Muon Candidate Directional Cosine", 20, -1., 1.);
  TH1F*    muon_candidate_directional_cosine_nu_background_copy                 = new TH1F("muon_candidate_directional_cosine_nu_background_copy", "Nu Background Muon Candidate Directional Cosine", 20, -1., 1.);
  TH1F*    muon_candidate_directional_cosine_ext_background_copy                = new TH1F("muon_candidate_directional_cosine_ext_background_copy", "EXT Background Muon Candidate Directional Cosine", 20, -1., 1.);
  TH1F*    muon_candidate_directional_cosine_dirt_copy                          = new TH1F("muon_candidate_directional_cosine_dirt_copy", "Dirt Muon Candidate Directional Cosine", 20, -1., 1.);

  maximum_bin_stacked_hist = -1.;
  maximum_bin_data         = -1.;

  for ( size_t copy_iter = 1; copy_iter < 21; copy_iter++ ) {

    muon_candidate_directional_cosine_signal_copy->SetBinContent( copy_iter, muon_candidate_directional_cosine_signal->GetBinContent( copy_iter ) );
    muon_candidate_directional_cosine_nu_background_copy->SetBinContent( copy_iter, muon_candidate_directional_cosine_nu_background->GetBinContent( copy_iter ) );
    muon_candidate_directional_cosine_ext_background_copy->SetBinContent( copy_iter, muon_candidate_directional_cosine_ext_background->GetBinContent( copy_iter ) );
    muon_candidate_directional_cosine_dirt_copy->SetBinContent( copy_iter, muon_candidate_directional_cosine_dirt->GetBinContent( copy_iter ) );

    muon_candidate_directional_cosine_ratio->SetBinContent( copy_iter, 0. );

    double sum = ( muon_candidate_directional_cosine_signal->GetBinContent( copy_iter ) + muon_candidate_directional_cosine_nu_background->GetBinContent( copy_iter ) + muon_candidate_directional_cosine_ext_background->GetBinContent( copy_iter ) + muon_candidate_directional_cosine_dirt->GetBinContent( copy_iter ) );

    if ( sum > 0. )
      muon_candidate_directional_cosine_ratio->SetBinContent( copy_iter, muon_candidate_directional_cosine_data->GetBinContent( copy_iter ) / sum );

    if ( sum > maximum_bin_stacked_hist )
      maximum_bin_stacked_hist = sum;

    if ( muon_candidate_directional_cosine_data->GetBinContent( copy_iter ) > maximum_bin_data )
      maximum_bin_data = muon_candidate_directional_cosine_data->GetBinContent( copy_iter );

  }

  // Fill the "stacked" plot.
  muon_candidate_directional_cosine_nu_background_copy->SetFillColor(kOrange);
  muon_candidate_directional_cosine_ext_background_copy->SetFillColor(kGreen);
  muon_candidate_directional_cosine_signal_copy->SetFillColor(kRed);
  muon_candidate_directional_cosine_dirt_copy->SetFillColor(kBlue);
  muon_candidate_directional_cosine_stack->Add(muon_candidate_directional_cosine_nu_background_copy);
  muon_candidate_directional_cosine_stack->Add(muon_candidate_directional_cosine_ext_background_copy);
  muon_candidate_directional_cosine_stack->Add(muon_candidate_directional_cosine_signal_copy);
  muon_candidate_directional_cosine_stack->Add(muon_candidate_directional_cosine_dirt_copy);

  TCanvas* muon_candidate_directional_cosine_canvas = new TCanvas("muon_candidate_directional_cosine_canvas", "muon_candidate_directional_cosine_canvas", 600, 600);
  muon_candidate_directional_cosine_canvas->cd();
  muon_candidate_directional_cosine_stack->SetMaximum( std::max( maximum_bin_stacked_hist, maximum_bin_data ) + 50.0 );
  gStyle->SetOptStat(0);
  muon_candidate_directional_cosine_stack->Draw();
  muon_candidate_directional_cosine_stack->GetXaxis()->SetTitle("Muon Candidate Directional Cosine [Unitless]");
  muon_candidate_directional_cosine_stack->GetYaxis()->SetTitle("Events Per Bin");
  muon_candidate_directional_cosine_stack->GetXaxis()->SetTitleOffset(1.10);
  muon_candidate_directional_cosine_stack->GetYaxis()->SetTitleOffset(1.48);
  muon_candidate_directional_cosine_stack->GetXaxis()->SetTickLength(0.015);
  muon_candidate_directional_cosine_stack->GetYaxis()->SetTickLength(0.015);
  muon_candidate_directional_cosine_stack->GetYaxis()->SetNdivisions( 10 );
  gPad->Modified();
  muon_candidate_directional_cosine_data->SetLineColor(kBlack);
  muon_candidate_directional_cosine_data->Draw("SamesE1");
  muon_candidate_directional_cosine_canvas->Write();
  
  file->Write();
  file->Close();

}
