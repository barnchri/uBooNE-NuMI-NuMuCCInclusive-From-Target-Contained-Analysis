////////////////////////////////////////////////////////////////////////
// Class:       UBXSec
// Plugin Type: producer
// File:        UBXSec_module.cc
//
// Generated at Monday August 3 12:00:00 2020 by Christopher Barnes using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class UBXSec
 *
 * \ingroup UBXSec
 *
 * \brief Art producer module
 *
 *
 * \author Christopher Barnes <cbarnes4@fnal.gov>
 *
 * \version producer
 *
 * \date 2020/08/03
 *
 * Contact: cbarnes4@fnal.gov
 *
 * Created on: Monday August 3rd 12:00:00 by Christopher Barnes
 *
 */

// Art include
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "ubobj/UBXSec/FlashMatch.h"
#include "ubobj/UBXSec/MCGhost.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "ubobj/UBXSec/TPCObject.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "ubobj/UBXSec/UBXSecEvent.h"
#include "ubobj/UBXSec/SelectionResult.h"

#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "ubevt/Utilities/SignalShapingServiceMicroBooNE.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


// LArSoft include
#include "ubreco/UBFlashFinder/PECalib.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/MCBase/MCDataHolder.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// Algorithms include
#include "ubana/UBXSec/Algorithms/UBXSecHelper.h"
#include "ubana/UBXSec/Algorithms/VertexCheck.h"
#include "ubana/UBXSec/Algorithms/FindDeadRegions.h"
#include "ubana/UBXSec/Algorithms/MuonCandidateFinder.h"
#include "ubana/UBXSec/Algorithms/FiducialVolume.h"
#include "ubana/UBXSec/Algorithms/NuMuCCEventSelection.h"
#include "ubana/UBXSec/Algorithms/TrackQuality.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

// Include statements for BackTracker.
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "ubana/myClasses/BackTrackerTruthMatch.h"
#include "ubana/myClasses/TruncMean.h"

// Include the modules for the space charge correction.
#include "ubevt/SpaceCharge/SpaceChargeMicroBooNE.h"
#include "ubevt/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"

// Include the file that we need for the momentum calculation.                                                                                                                                          
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

// event weight include file.                                                                                                                                                                            
#include "larsim/EventWeight/Base/MCEventWeight.h"

// Helper files.
#include "helpers/PandoraInterfaceHelper.h"
#include "helpers/TrackHelper.h"

// PID files.
#include "PID.h"

// Root include
#include "TString.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <fstream>
#include <string>
#include <iterator>

// Geometry package.
#include "larcore/Geometry/Geometry.h"

// Additional geometry package.
#include "larcorealg/Geometry/geo_vectors_utils.h"

namespace ubxsec {
  struct Hit3D_t {
    double x;
    double y;
    double z;
    double q;
  };
}

using namespace recob::tracking;

class UBXSec;


class UBXSec : public art::EDProducer {
public:
  explicit UBXSec(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBXSec(UBXSec const &) = delete;
  UBXSec(UBXSec &&) = delete;
  UBXSec & operator = (UBXSec const &) = delete;
  UBXSec & operator = (UBXSec &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
  void endSubRun(art::SubRun &sr) override;

private:

  /// Prints MC particles from GENIE on the screen
  void PrintMC(std::vector<art::Ptr<simb::MCTruth>> mclist);

  /// Calculates flash position
  void GetFlashLocation(std::vector<double>, double&, double&, double&, double&);
  float GetTrackShowerScore(art::Ptr<recob::PFParticle> pfParticle, const lar_pandora::PFParticlesToMetadata pfParticleToMetadata);
  FindDeadRegions deadRegionsFinder;
void GetTaggedPFP(art::Event const & e, std::string cosmictag_producer, double score_cut, lar_pandora::PFParticleVector & pfpTaggedOut,std::vector<int> & tagid_v);
int PFPInCommon(lar_pandora::PFParticleVector first, lar_pandora::PFParticleVector second);
  //ubxsec::McPfpMatch mcpfpMatcher;
  ::ubana::FiducialVolume _fiducial_volume;
  ::ubana::MuonCandidateFinder _muon_finder;
  ::ubana::NuMuCCEventSelection _event_selection;
  ::pmtana::PECalib _pecalib;
  ::trkf::TrackMomentumCalculator  _trk_mom_calculator{0.1};

  // Database to understand particle pdg
  const TDatabasePDG* _database_pdg = TDatabasePDG::Instance();

  // Services
  ::detinfo::DetectorProperties const* _detector_properties;
  ::detinfo::DetectorClocks const* _detector_clocks;
  spacecharge::SpaceCharge const* _SCE;

  // To be set via fcl parameters
  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_tag_producer;
  std::string _neutrino_flash_match_producer;
  std::string _cosmic_flash_match_producer;
  std::string _opflash_producer_beam;

 std::string _cosmic_flash_tag_producer;
  std::string _cosmic_geo_tag_producer;
  std::string _cosmic_acpt_tag_producer;
  std::string _cosmic_stopmu_tag_producer;

  std::string _tpcobject_producer;
  std::string _potsum_producer;
  std::string _particle_id_producer;
  std::string _mc_ghost_producer;
  std::string _geocosmictag_producer;
  std::string _candidateconsistency_producer;
  std::string _mcsfitresult_mu_producer;
  std::string _calorimetry_producer;
  std::string _eventweight_producer;
  std::string _genie_eventweight_pm1_producer;
  std::string _genie_eventweight_multisim_producer;
  std::string _flux_eventweight_multisim_producer;
  bool _debug = false;                   ///< Debug mode
  bool _debug_cr = false;                   ///< Debug mode

  int _minimumHitRequirement;           ///< Minimum number of hits in at least a plane for a track
  double _minimumDistDeadReg;           ///< Minimum distance the track end points can have to a dead region
  bool _use_genie_info;                 ///< Turn this off if looking at cosmic only files
  double _beam_spill_start;             ///< Start time of the beam spill (us)
  double _beam_spill_end;               ///< Start time of the beam spill (us)
  double _total_pe_cut;                 ///< PE cut to be applied to beam flash
  double _geo_cosmic_score_cut;         ///< Cut on the score of the pandoraNu geo cosmic tagger
  double _tolerance_track_multiplicity; ///< Tolerance to consider a track coming from the nu reco vertex
  double _min_track_len;                ///< Min track length for momentum calculation
  bool _make_ophit_csv;                 ///< If true makea a csv file with ophit info
  bool _make_pida_csv;                  ///< If true makea a csv file with pida/tracklength info

  bool _do_opdet_swap;                  ///< If true swaps reconstructed OpDets according to _opdet_swap_map
  std::vector<int> _opdet_swap_map;     ///< The OpDet swap map for reco flashes
  double _cosmic_flash_tag_score_cut; ///< Score cut used in the analysis to consider the PFP as cosmic (applied to flash tagger)
  double _cosmic_geo_tag_score_cut;   ///< Score cut used in the analysis to consider the PFP as cosmic (applied to geo tagger)
  double _cosmic_acpt_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to acpt tagger)
  double _cosmic_stopmu_tag_score_cut;  ///< Score cut used in the analysis to consider the PFP as cosmic (applied to stopmu tagger)
  // Constants
  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;

  // To be filled within module
  bool _is_data, _is_mc;
  //double _candidate_flash_time;
  double _drift_velocity;

  // Outputs trees
  TTree* _tree1;
  UBXSecEvent *ubxsec_event = new UBXSecEvent();

  TH2F * _deadRegion2P;
  TH2F * _deadRegion3P;

  // PIDA related variables
  TH1D * _h_pida_proton,     * _h_pida_muon,     * _h_pida_pion,     * _h_pida_kaon;
  TH2D * _h_pida_len_proton, * _h_pida_len_muon, * _h_pida_len_pion, * _h_pida_len_kaon;

  // Momentum Related Variables
  TH2D* _h_mom_true_mcs; ///< 2D histogram of true muon momentum VS reconstructed (using MCS)
  TH2D* _h_mom_true_mcs_contained; ///< 2D histogram of true muon momentum VS reconstructed (using MCS) (contained tracks)
  TH2D* _h_mom_true_mcs_uncontained; ///< 2D histogram of true muon momentum VS reconstructed (using MCS) (uncontained tracks)
  TH2D* _h_mom_true_range_contained; ///< 2D histogram of true muon momentum VS reconstructed (using Length) (contained tracks)
  TH2D* _h_mom_range_mcs_contained;  ///< 2D histogram of reconstructed (using MCS) muon momentum VS reconstructed (using Length) (contained tracks)
  TH1D* _h_mcs_cosmic_track_direction; ///< Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward)
  TH2D* _h_mcs_cosmic_track_direction_deltall; /// Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward) VS delta LL from MCS fit
  TH2D* _h_mcs_cosmic_track_direction_ratioll; /// Track direction from cosmic origin TPCObjects as given by mcs (0: downward, 1: upward) VS ratio LL from MCS fit
  TTree *_mom_tree_contained, *_mom_tree_uncontained;
  double _mom_true_contained;
  double _mom_mcs_contained;
  double _mom_range_contained;
  double _mom_true_uncontained;
  double _mom_mcs_uncontained;
  TTree *_mcs_cosmic_track_direction_tree;
  double _mcs_cosmic_track_direction;
  double _mcs_cosmic_track_downll;
  double _mcs_cosmic_track_upll;
  TTree *_mom_cosmic_tree;
  double _mom_cosmic_true;
  double _mom_cosmic_mcs;
  double _mom_cosmic_mcs_downforced;
  double _mom_cosmic_range;
  bool _mom_cosmic_down;

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot;

  TTree* _passing_events_tree;
  TTree* Truth_Tree;
  int    _run;
  int    _subrun;
  int    _event;
  int    num_of_neutrinos;
  double neutrino_energy;
  int    truth_vertex_is_in_AV;
  int    truth_vertex_is_in_FV;
  double truth_theta_angle_for_weighting;
  double parent_fvx;
  double parent_fvy;
  double parent_fvz;
  int    NC_channel;
  int    neutrino_PDG;
  int    neutrino_interaction_mode;
  double nu_vtx_x_truth;
  double nu_vtx_y_truth;
  double nu_vtx_z_truth;
  double neutrino_px;
  double neutrino_py;
  double neutrino_pz;
  double spline_fix_mcweight;
  double central_value_mcweight;
  double rootino_fix_mcweight;
  double ppfx_value_mcweight;
  double flux_other_universe_mcweights[600];
  double other_universe_mcweights[100];
  double axialff_mcweights[2];
  double rpaccqe_mcweights[2];
  double xsecshape_mcweights[2];
  int    num_muminus_tracks;
  int    num_muplus_tracks;
  int    num_piplus_tracks;
  int    num_piminus_tracks;
  int    num_pi0_tracks;
  int    num_proton_tracks;
  int    num_electron_showers;
  int    num_positron_showers;
  int    num_photon_showers;
  int    num_mctrack_points;
  double truth_muon_starting_x_coordinate;
  double truth_muon_starting_y_coordinate;
  double truth_muon_starting_z_coordinate;
  double truth_muon_ending_x_coordinate;
  double truth_muon_ending_y_coordinate;
  double truth_muon_ending_z_coordinate;
  double truth_muon_length;
  int    truth_muon_track_is_contained;
  double truth_muon_px;
  double truth_muon_py;
  double truth_muon_pz;
  double truth_muon_kinetic_energy;
  double vertex_location_x;
  double vertex_location_y;
  double vertex_location_z;
  double other_end_location_x;
  double other_end_location_y;
  double other_end_location_z;
  double pandora_vertex_x;
  double pandora_vertex_y;
  double pandora_vertex_z;
  double truth_reco_vtx_distance;
  double distance_between_truth_vertex_and_pandora_vertex;
  double distance_between_truth_vertex_and_closest_track_terminal_point;
  int    at_least_one_terminal_point_in_FV;
  int    _nflashes_in_beamgate;
  int    _nflashes_in_beamspill;
  int    _nflashes_in_beamspill_window_passing_filter_PE_cut;
  double flash_time;
  double flash_PEs;
  double flash_z;
  double flash_y;
  double flash_z_width;
  double flash_y_width;
  int    nslices_in_event;
  int    number_of_tracks_in_TPCObject;
  int    num_of_tracks_originating_from_vertex;
  int    num_of_tracks_originating_from_end;
  double average_track_length_originating_from_vertex_other_than_muon;
  double average_track_length_originating_from_end;
  double muon_candidate_length;
  double longest_TPCObject_track_length;
  int    muon_candidate_is_contained;
  int    entire_event_is_contained;
  double muon_candidate_kinetic_energy_range;
  double muon_candidate_kinetic_energy_mcs;
  double sum_of_TPCObject_track_lengths;
  double total_length_of_tracks_originating_from_vertex;
  double Muon_Candidate_Muon_Track_Score_;
  double Muon_Candidate_Proton_PID_;
  double Muon_Candidate_NuScore_;
  double Muon_Candidate_PDGCode_;
  double ubxsec_muon_phi;
  double ubxsec_muon_cos_theta;
  int    fails_fiducial_volume_req;
  double muon_candidate_residuals;
  int    fails_residuals_req;
  double muon_candidate_percent_used_hits_in_cluster;
  int    fails_perc_used_hits_in_cluster;

  std::ofstream _csvfile, _csvfile2;

  // Declare the PID variables.
  double this_chimuon_;
  double this_chiproton_;
  double this_chipion_;
  double Proton_PID_;
  double NuScore_;
  double Muon_Track_Score_;
  double muon_proton_ratio;
  double muon_pion_ratio;

  // Declare the variables for the track orientation algorithm.
  bool u_plane_has_more_than_20_points;
  bool v_plane_has_more_than_20_points;
  bool y_plane_has_more_than_20_points;
  bool u_plane_has_track_forward;
  bool v_plane_has_track_forward;
  bool y_plane_has_track_forward;

  // Variables used to measure track quantities.
  double muon_track_first_point_x;
  double muon_track_first_point_y;
  double muon_track_first_point_z;
  double muon_track_last_point_x;
  double muon_track_last_point_y;
  double muon_track_last_point_z;

  // Count the total number of events, the number of events passing, and the number of events failing.                                                                                                 
  int    total_num_events;
  int    total_num_events_failing;
  int    total_num_events_passing;
  double total_num_events_weighted;
  double total_num_events_failing_weighted;
  double total_num_events_passing_weighted;

  double total_num_events_failed_beam_disc_flashes_weighted;
  double total_num_events_failed_beam_spill_flash_weighted;
  double total_num_events_failed_has_slices_weighted;
  double total_num_events_failed_fiducial_volume_weighted;
  double total_num_events_failed_ntrack_weighted;
  double total_num_events_failed_residuals_std_up_weighted;
  double total_num_events_failed_perc_used_hits_in_cluster_weighted;

  // The number of events failing each of the different cuts.
  int    total_num_events_failed_beam_disc_flashes;
  int    total_num_events_failed_beam_spill_flash;
  int    total_num_events_failed_has_slices;
  int    total_num_events_failed_track_length;
  int    total_num_events_failed_has_slice_tagged_as_neutrino;
  int    total_num_events_failed_fiducial_volume;
  int    total_num_events_failed_ntrack;
  int    total_num_events_failed_residuals_std_up;
  int    total_num_events_failed_perc_used_hits_in_cluster;
  
  int  event_counter;

};


UBXSec::UBXSec(fhicl::ParameterSet const & p) 
{

  ::art::ServiceHandle<geo::Geometry> geo;

  _pfp_producer                        = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel                      = p.get<std::string>("HitProducer");
  _geantModuleLabel                    = p.get<std::string>("GeantModule");
  _spacepointLabel                     = p.get<std::string>("SpacePointProducer");
  _neutrino_flash_match_producer       = p.get<std::string>("NeutrinoFlashMatchProducer");
  _cosmic_flash_match_producer         = p.get<std::string>("CosmicFlashMatchProducer");
  _opflash_producer_beam               = p.get<std::string>("OpFlashBeamProducer");
  _tpcobject_producer                  = p.get<std::string>("TPCObjectProducer");
  _potsum_producer                     = p.get<std::string>("POTSummaryProducer");
  _particle_id_producer                = p.get<std::string>("ParticleIDProducer");
  _mc_ghost_producer                   = p.get<std::string>("MCGhostProducer");
  _geocosmictag_producer               = p.get<std::string>("GeoCosmicTaggerProducer");
  _candidateconsistency_producer       = p.get<std::string>("CandidateConsistencyProducer");
  _mcsfitresult_mu_producer            = p.get<std::string>("MCSFitResultMuProducer");
  _calorimetry_producer                = p.get<std::string>("CalorimetryProducer");
  _eventweight_producer                = p.get<std::string>("EventWeightProducer");
  _genie_eventweight_pm1_producer      = p.get<std::string>("GenieEventWeightPMOneProducer");
  _genie_eventweight_multisim_producer = p.get<std::string>("GenieEventWeightMultisimProducer");
  _flux_eventweight_multisim_producer  = p.get<std::string>("FluxEventWeightMultisimProducer");
  _cosmic_flash_tag_producer           = p.get<std::string>("CosmicFlashTagProducer");
  _cosmic_geo_tag_producer             = p.get<std::string>("CosmicGeoTagProducer");
  _cosmic_acpt_tag_producer            = p.get<std::string>("CosmicACPTTagProducer");
  _cosmic_stopmu_tag_producer          = p.get<std::string>("CosmicStopMuTagProducer");
  _mc_ghost_producer                   = p.get<std::string>("MCGhostProducer");
  _cosmic_flash_tag_score_cut          = p.get<double>("CosmicFlashTagScoreCut",0.99);
  _cosmic_geo_tag_score_cut            = p.get<double>("CosmicGeoTagScoreCut",0.6);
  _cosmic_acpt_tag_score_cut           = p.get<double>("CosmicACPTTagScoreCut",0.99);
  _cosmic_stopmu_tag_score_cut         = p.get<double>("CosmicStopMuTagScoreCut",0.99);

  _use_genie_info                      = p.get<bool>("UseGENIEInfo", false);
  _minimumHitRequirement               = p.get<int>("MinimumHitRequirement", 3);
  _minimumDistDeadReg                  = p.get<double>("MinimumDistanceToDeadRegion", 5.);

  _beam_spill_start                    = p.get<double>("BeamSpillStart", 3.2);
  _beam_spill_end                      = p.get<double>("BeamSpillEnd",   4.8);
  _total_pe_cut                        = p.get<double>("TotalPECut",     50);

  _do_opdet_swap                       = p.get<bool>("DoOpDetSwap", false);
  _opdet_swap_map                      = p.get<std::vector<int> >("OpDetSwapMap");

  _geo_cosmic_score_cut                = p.get<double>("GeoCosmicScoreCut", 0.6);
  _tolerance_track_multiplicity        = p.get<double>("ToleranceTrackMultiplicity", 5.);

  _make_ophit_csv                      = p.get<bool>("MakeOpHitCSV", false);
  _make_pida_csv                       = p.get<bool>("MakePIDACSV", false);

  _pecalib.Configure(p.get<fhicl::ParameterSet>("PECalib"));

  _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
                             geo->DetHalfHeight(),
                             2.*geo->DetHalfWidth(),
                             geo->DetLength());

  _fiducial_volume.PrintConfig();

  _muon_finder.Configure(p.get<fhicl::ParameterSet>("MuonCandidateFinderSettings"));

  _muon_finder.PrintConfig();

  _event_selection.Configure(p.get<fhicl::ParameterSet>("NuMuCCSelectionSettings"));

  _event_selection.PrintConfig();

  //_mcs_fitter.Configure( p.get<fhicl::ParameterSet>("MCSFitter"));

  _detector_properties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _detector_clocks = lar::providerFrom<detinfo::DetectorClocksService>();
  _SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  std::cout << "E Field: "            << _detector_properties->Efield() << std::endl;
  std::cout << "Temperature: "        << _detector_properties->Temperature() << std::endl;
  std::cout << "Drift Velocity: "     << _detector_properties->DriftVelocity(_detector_properties->Efield(), _detector_properties->Temperature())<< std::endl;
  std::cout << "Sampling Rate: "      << _detector_properties->SamplingRate() << std::endl;
  std::cout << "Beam Spill Start = "  << _beam_spill_start << " us." << std::endl;
  std::cout << "Beam Spill Ends = "   << _beam_spill_end << " us." << std::endl;

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");

  int bufsize    = 16000;
  int splitlevel = 99;
  _tree1->Branch("ubxsec_event_split", &ubxsec_event, bufsize, splitlevel);

  _deadRegion2P                         = fs->make<TH2F>("deadRegion2P","deadRegion2P", 10350,0.0,1035.0,2300,-115.0,115.0);
  _deadRegion3P                         = fs->make<TH2F>("deadRegion3P","deadRegion3P", 10350,0.0,1035.0,2300,-115.0,115.0);

  _h_pida_muon                          = fs->make<TH1D>("h_pida_muon", "Muon tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_proton                        = fs->make<TH1D>("h_pida_proton", "Proton tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_pion                          = fs->make<TH1D>("h_pida_pion", "Pion tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);
  _h_pida_kaon                          = fs->make<TH1D>("h_pida_kaon", "Kaon tracks;PIDa [MeV/cm^{1.42}];", 50, 0, 20);

  _h_pida_len_muon                      = fs->make<TH2D>("h_pida_len_muon", "Muon tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_proton                    = fs->make<TH2D>("h_pida_len_proton", "Proton tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_pion                      = fs->make<TH2D>("h_pida_len_pion", "Pion tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);
  _h_pida_len_kaon                      = fs->make<TH2D>("h_pida_len_kaon", "Kaon tracks;PIDa [MeV/cm^{1.42}];Track length [cm];", 50, 0, 20, 100, 0, 700);

  _h_mom_true_mcs                       = fs->make<TH2D>("h_mom_true_mcs", ";True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_mcs_contained             = fs->make<TH2D>("h_mom_true_mcs_contained", "Contained;True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_mcs_uncontained           = fs->make<TH2D>("h_mom_true_mcs_uncontained", "Uncontained;True Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_true_range_contained           = fs->make<TH2D>("h_mom_true_range_contained", "Contained;True Muon Momentum [GeV];Reconstructed (via Length) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);
  _h_mom_range_mcs_contained            = fs->make<TH2D>("h_mom_range_mcs_contained", "Contained;Reconstructed (via Length) Muon Momentum [GeV];Reconstructed (via MCS) Muon Momentum [GeV];", 80, 0, 2, 80, 0, 2);

  _h_mcs_cosmic_track_direction         = fs->make<TH1D>("h_mcs_cosmic_track_direction", "0: down, 1: up;;", 2, 0, 2);
  _h_mcs_cosmic_track_direction_deltall = fs->make<TH2D>("h_mcs_cosmic_track_direction_deltall", ";0: down, 1: up;(FWD - BWD) LL;", 2, 0, 2, 500, -1e-6, 1e-6);
  _h_mcs_cosmic_track_direction_ratioll = fs->make<TH2D>("h_mcs_cosmic_track_direction_ratioll", ";0: down, 1: up;(FWD / BWD) LL;", 2, 0, 2, 500, 0, 2);

  _mom_tree_contained              = fs->make<TTree>("mom_tree_contained","");
  _mom_tree_contained->Branch("run",                                                             &_run,                                               "run/I");
  _mom_tree_contained->Branch("subrun",                                                          &_subrun,                                            "subrun/I");
  _mom_tree_contained->Branch("event",                                                           &_event,                                             "event/I");
  _mom_tree_contained->Branch("mom_true_contained",                                              &_mom_true_contained,                                "mom_true_contained/D");
  _mom_tree_contained->Branch("mom_mcs_contained",                                               &_mom_mcs_contained,                                 "mom_mcs_contained/D");
  _mom_tree_contained->Branch("mom_range_contained",                                             &_mom_range_contained,                               "mom_range_contained/D");

  _mom_tree_uncontained            = fs->make<TTree>("mom_tree_uncontained","");
  _mom_tree_uncontained->Branch("run",                                                           &_run,                                                "run/I");
  _mom_tree_uncontained->Branch("subrun",                                                        &_subrun,                                             "subrun/I");
  _mom_tree_uncontained->Branch("event",                                                         &_event,                                              "event/I");
  _mom_tree_uncontained->Branch("mom_true_uncontained",                                          &_mom_true_uncontained,                               "mom_true_uncontained/D");
  _mom_tree_uncontained->Branch("mom_mcs_uncontained",                                           &_mom_mcs_uncontained,                                "mom_mcs_uncontained/D");

  _mcs_cosmic_track_direction_tree = fs->make<TTree>("mcs_cosmic_track_direction_tree","");
  _mcs_cosmic_track_direction_tree->Branch("run",                                                &_run,                                                "run/I");
  _mcs_cosmic_track_direction_tree->Branch("subrun",                                             &_subrun,                                             "subrun/I");
  _mcs_cosmic_track_direction_tree->Branch("event",                                              &_event,                                              "event/I");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_direction",                         &_mcs_cosmic_track_direction,                         "mcs_cosmic_track_direction/D");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_downll",                            &_mcs_cosmic_track_downll,                            "mcs_cosmic_track_downll/D");
  _mcs_cosmic_track_direction_tree->Branch("mcs_cosmic_track_upll",                              &_mcs_cosmic_track_upll,                              "mcs_cosmic_track_upll/D");

  _mom_cosmic_tree                 = fs->make<TTree>("mom_cosmic_tree","");
  _mom_cosmic_tree->Branch("run",                                                                &_run,                                                "run/I");
  _mom_cosmic_tree->Branch("subrun",                                                             &_subrun,                                             "subrun/I");
  _mom_cosmic_tree->Branch("event",                                                              &_event,                                              "event/I");
  _mom_cosmic_tree->Branch("mom_cosmic_true",                                                    &_mom_cosmic_true,                                    "mom_cosmic_true/D");
  _mom_cosmic_tree->Branch("mom_cosmic_mcs",                                                     &_mom_cosmic_mcs,                                     "mom_cosmic_mcs/D");
  _mom_cosmic_tree->Branch("mom_cosmic_mcs_downforced",                                          &_mom_cosmic_mcs_downforced,                          "mom_cosmic_mcs_downforced/D");
  _mom_cosmic_tree->Branch("mom_cosmic_range",                                                   &_mom_cosmic_range,                                   "mom_cosmic_range/D");
  _mom_cosmic_tree->Branch("mom_cosmic_down",                                                    &_mom_cosmic_down,                                    "mom_cosmic_down/O");

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run",                                                                        &_sr_run,                                             "run/I");
  _sr_tree->Branch("subrun",                                                                     &_sr_subrun,                                          "subrun/I");
  _sr_tree->Branch("begintime",                                                                  &_sr_begintime,                                       "begintime/D");
  _sr_tree->Branch("endtime",                                                                    &_sr_endtime,                                         "endtime/D");
  _sr_tree->Branch("pot",                                                                        &_sr_pot,                                             "pot/D");

  _passing_events_tree = fs->make<TTree>("_passing_events_tree", "A tree that contains muon quantities for the events that pass the NuMuCCInclusive Filter");
  _passing_events_tree->Branch("run",                                                            &_run,                                                "run/I");
  _passing_events_tree->Branch("subrun",                                                         &_subrun,                                             "subrun/I");
  _passing_events_tree->Branch("event",                                                          &_event,                                              "event/I");
  _passing_events_tree->Branch("num_of_neutrinos",                                               &num_of_neutrinos,                                    "num_of_neutrinos/I");
  _passing_events_tree->Branch("neutrino_energy",                                                &neutrino_energy,                                     "neutrino_energy/D");
  _passing_events_tree->Branch("truth_vertex_is_in_AV",                                          &truth_vertex_is_in_AV,                               "truth_vertex_is_in_AV/I");
  _passing_events_tree->Branch("truth_vertex_is_in_FV",                                          &truth_vertex_is_in_FV,                               "truth_vertex_is_in_FV/I");
  _passing_events_tree->Branch("truth_theta_angle_for_weighting",                                &truth_theta_angle_for_weighting,                     "truth_theta_angle_for_weighting/D");
  _passing_events_tree->Branch("parent_fvx",                                                     &parent_fvx,                                          "parent_fvx/D");
  _passing_events_tree->Branch("parent_fvy",                                                     &parent_fvy,                                          "parent_fvy/D");
  _passing_events_tree->Branch("parent_fvz",                                                     &parent_fvz,                                          "parent_fvz/D");
  _passing_events_tree->Branch("NC_channel",                                                     &NC_channel,                                          "NC_channel/I");
  _passing_events_tree->Branch("neutrino_PDG",                                                   &neutrino_PDG,                                        "neutrino_PDG/I");
  _passing_events_tree->Branch("neutrino_interaction_mode",                                      &neutrino_interaction_mode,                           "neutrino_interaction_mode/I");
  _passing_events_tree->Branch("nu_vtx_x_truth",                                                 &nu_vtx_x_truth,                                      "nu_vtx_x_truth/D");
  _passing_events_tree->Branch("nu_vtx_y_truth",                                                 &nu_vtx_y_truth,                                      "nu_vtx_y_truth/D");
  _passing_events_tree->Branch("nu_vtx_z_truth",                                                 &nu_vtx_z_truth,                                      "nu_vtx_z_truth/D");
  _passing_events_tree->Branch("neutrino_px",                                                    &neutrino_px,                                         "neutrino_px/D");
  _passing_events_tree->Branch("neutrino_py",                                                    &neutrino_py,                                         "neutrino_py/D");
  _passing_events_tree->Branch("neutrino_pz",                                                    &neutrino_pz,                                         "neutrino_pz/D");
  _passing_events_tree->Branch("spline_fix_mcweight",                                            &spline_fix_mcweight,                                 "spline_fix_mcweight/D");
  _passing_events_tree->Branch("central_value_mcweight",                                         &central_value_mcweight,                              "central_value_mcweight/D");
  _passing_events_tree->Branch("rootino_fix_mcweight",                                           &rootino_fix_mcweight,                                "rootino_fix_mcweight/D");
  _passing_events_tree->Branch("ppfx_value_mcweight",                                            &ppfx_value_mcweight,                                 "ppfx_value_mcweight/D");
  _passing_events_tree->Branch("flux_other_universe_mcweights",                                  &flux_other_universe_mcweights,                       "flux_other_universe_mcweights[600]/D");
  _passing_events_tree->Branch("other_universe_mcweights",                                       &other_universe_mcweights,                            "other_universe_mcweights[100]/D");
  _passing_events_tree->Branch("axialff_mcweights",                                              &axialff_mcweights,                                   "axialff_mcweights[2]/D");
  _passing_events_tree->Branch("rpaccqe_mcweights",                                              &rpaccqe_mcweights,                                   "rpaccqe_mcweights[2]/D");
  _passing_events_tree->Branch("xsecshape_mcweights",                                            &xsecshape_mcweights,                                 "xsecshape_mcweights[2]/D");
  _passing_events_tree->Branch("num_muminus_tracks",                                             &num_muminus_tracks,                                  "num_muminus_tracks/I");
  _passing_events_tree->Branch("num_muplus_tracks",                                              &num_muplus_tracks,                                   "num_muplus_tracks/I");
  _passing_events_tree->Branch("num_piplus_tracks",                                              &num_piplus_tracks,                                   "num_piplus_tracks/I");
  _passing_events_tree->Branch("num_piminus_tracks",                                             &num_piminus_tracks,                                  "num_piminus_tracks/I");
  _passing_events_tree->Branch("num_pi0_tracks",                                                 &num_pi0_tracks,                                      "num_pi0_tracks/I");
  _passing_events_tree->Branch("num_proton_tracks",                                              &num_proton_tracks,                                   "num_proton_tracks/I");
  _passing_events_tree->Branch("num_electron_showers",                                           &num_electron_showers,                                "num_electron_showers/I");
  _passing_events_tree->Branch("num_positron_showers",                                           &num_positron_showers,                                "num_positron_showers/I");
  _passing_events_tree->Branch("num_photon_showers",                                             &num_photon_showers,                                  "num_photon_showers/I");
  _passing_events_tree->Branch("num_mctrack_points",                                             &num_mctrack_points,                                  "num_mctrack_points/I");
  _passing_events_tree->Branch("truth_muon_starting_x_coordinate",                               &truth_muon_starting_x_coordinate,                    "truth_muon_starting_x_coordinate/D");
  _passing_events_tree->Branch("truth_muon_starting_y_coordinate",                               &truth_muon_starting_y_coordinate,                    "truth_muon_starting_y_coordinate/D");
  _passing_events_tree->Branch("truth_muon_starting_z_coordinate",                               &truth_muon_starting_z_coordinate,                    "truth_muon_starting_z_coordinate/D");
  _passing_events_tree->Branch("truth_muon_ending_x_coordinate",                                 &truth_muon_ending_x_coordinate,                      "truth_muon_ending_x_coordinate/D");
  _passing_events_tree->Branch("truth_muon_ending_y_coordinate",                                 &truth_muon_ending_y_coordinate,                      "truth_muon_ending_y_coordinate/D");
  _passing_events_tree->Branch("truth_muon_ending_z_coordinate",                                 &truth_muon_ending_z_coordinate,                      "truth_muon_ending_z_coordinate/D");
  _passing_events_tree->Branch("truth_muon_length",                                              &truth_muon_length,                                   "truth_muon_length/D");
  _passing_events_tree->Branch("truth_muon_track_is_contained",                                  &truth_muon_track_is_contained,                       "truth_muon_track_is_contained/I");
  _passing_events_tree->Branch("truth_muon_px",                                                  &truth_muon_px,                                       "truth_muon_px/D");
  _passing_events_tree->Branch("truth_muon_py",                                                  &truth_muon_py,                                       "truth_muon_py/D");
  _passing_events_tree->Branch("truth_muon_pz",                                                  &truth_muon_pz,                                       "truth_muon_pz/D");
  _passing_events_tree->Branch("truth_muon_kinetic_energy",                                      &truth_muon_kinetic_energy,                           "truth_muon_kinetic_energy/D");
  _passing_events_tree->Branch("vertex_location_x",                                              &vertex_location_x,                                   "vertex_location_x/D");
  _passing_events_tree->Branch("vertex_location_y",                                              &vertex_location_y,                                   "vertex_location_y/D");
  _passing_events_tree->Branch("vertex_location_z",                                              &vertex_location_z,                                   "vertex_location_z/D");
  _passing_events_tree->Branch("other_end_location_x",                                           &other_end_location_x,                                "other_end_location_x/D");
  _passing_events_tree->Branch("other_end_location_y",                                           &other_end_location_y,                                "other_end_location_y/D");
  _passing_events_tree->Branch("other_end_location_z",                                           &other_end_location_z,                                "other_end_location_z/D");
  _passing_events_tree->Branch("pandora_vertex_x",                                               &pandora_vertex_x,                                    "pandora_vertex_x/D");
  _passing_events_tree->Branch("pandora_vertex_y",                                               &pandora_vertex_y,                                    "pandora_vertex_y/D");
  _passing_events_tree->Branch("pandora_vertex_z",                                               &pandora_vertex_z,                                    "pandora_vertex_z/D");
  _passing_events_tree->Branch("truth_reco_vtx_distance",                                        &truth_reco_vtx_distance,                             "truth_reco_vtx_distance/D");
  _passing_events_tree->Branch("distance_between_truth_vertex_and_pandora_vertex",               &distance_between_truth_vertex_and_pandora_vertex,    "distance_between_truth_vertex_and_pandora_vertex/D\
");
  _passing_events_tree->Branch("distance_between_truth_vertex_and_closest_track_terminal_point", &distance_between_truth_vertex_and_closest_track_terminal_point,                                            "distance_between_truth_vertex_and_closest_track_terminal_point/D");
  _passing_events_tree->Branch("at_least_one_terminal_point_in_FV",                              &at_least_one_terminal_point_in_FV,                   "at_least_one_terminal_point_in_FV/I");
  _passing_events_tree->Branch("num_beam_flashes",                                               &_nflashes_in_beamgate,                               "num_beam_flashes/I"); 
  _passing_events_tree->Branch("_nflashes_in_beamspill",                                         &_nflashes_in_beamspill,                              "_nflashes_in_beamspill/I");
  _passing_events_tree->Branch("_nflashes_in_beamspill_window_passing_filter_PE_cut",            &_nflashes_in_beamspill_window_passing_filter_PE_cut,                                                       "_nflashes_in_beamspill_window_passing_filter_PE_cut/I");
  _passing_events_tree->Branch("flash_time",                                                     &flash_time,                                           "flash_time/D");
  _passing_events_tree->Branch("flash_PEs",                                                      &flash_PEs,                                            "flash_PEs/D");
  _passing_events_tree->Branch("flash_z",                                                        &flash_z,                                              "flash_z/D");
  _passing_events_tree->Branch("flash_y",                                                        &flash_y,                                              "flash_y/D");
  _passing_events_tree->Branch("flash_z_width",                                                  &flash_z_width,                                        "flash_z_width/D");
  _passing_events_tree->Branch("flash_y_width",                                                  &flash_y_width,                                        "flash_y_width/D");
  _passing_events_tree->Branch("nslices_in_event",                                               &nslices_in_event,                                     "nslices_in_event/I");
  _passing_events_tree->Branch("number_of_tracks_in_TPCObject",                                  &number_of_tracks_in_TPCObject,                        "number_of_tracks_in_TPCObject/I");
  _passing_events_tree->Branch("num_of_tracks_originating_from_vertex",                          &num_of_tracks_originating_from_vertex,                "num_of_tracks_originating_from_vertex/I");
  _passing_events_tree->Branch("num_of_tracks_originating_from_end",                             &num_of_tracks_originating_from_end,                   "num_of_tracks_originating_from_end/I");
  _passing_events_tree->Branch("average_track_length_originating_from_vertex_other_than_muon",   &average_track_length_originating_from_vertex_other_than_muon,                                              "average_track_length_originating_from_vertex_other_than_muon/D");
  _passing_events_tree->Branch("average_track_length_originating_from_end",                      &average_track_length_originating_from_end,                                                                 "average_track_length_originating_from_end/D");
  _passing_events_tree->Branch("muon_candidate_length",                                          &muon_candidate_length,                                "muon_candidate_length/D");
  _passing_events_tree->Branch("longest_TPCObject_track_length",                                 &longest_TPCObject_track_length,                       "longest_TPCObject_track_length/D");
  _passing_events_tree->Branch("muon_candidate_is_contained",                                    &muon_candidate_is_contained,                          "muon_candidate_is_contained/I");
  _passing_events_tree->Branch("entire_event_is_contained",                                      &entire_event_is_contained,                            "entire_event_is_contained/I");
  _passing_events_tree->Branch("muon_candidate_kinetic_energy_range",                            &muon_candidate_kinetic_energy_range,                  "muon_candidate_kinetic_energy_range/D");
  _passing_events_tree->Branch("muon_candidate_kinetic_energy_mcs",                              &muon_candidate_kinetic_energy_mcs,                    "muon_candidate_kinetic_energy_mcs/D");
  _passing_events_tree->Branch("sum_of_TPCObject_track_lengths",                                 &sum_of_TPCObject_track_lengths,                       "sum_of_TPCObject_track_lengths/D");
  _passing_events_tree->Branch("total_length_of_tracks_originating_from_vertex",                 &total_length_of_tracks_originating_from_vertex,       "total_length_of_tracks_originating_from_vertex/D");
  _passing_events_tree->Branch("Muon_Candidate_Muon_Track_Score_",                               &Muon_Candidate_Muon_Track_Score_,                     "Muon_Candidate_Muon_Track_Score_/D");
  _passing_events_tree->Branch("Muon_Candidate_Proton_PID_",                                     &Muon_Candidate_Proton_PID_,                           "Muon_Candidate_Proton_PID_/D");
  _passing_events_tree->Branch("Muon_Candidate_NuScore_",                                        &Muon_Candidate_NuScore_,                              "Muon_Candidate_NuScore_/D");
  _passing_events_tree->Branch("Muon_Candidate_PDGCode_",                                        &Muon_Candidate_PDGCode_,                              "Muon_Candidate_PDGCode_/D");
  _passing_events_tree->Branch("ubxsec_muon_phi",                                                &ubxsec_muon_phi,                                     "ubxsec_muon_phi/D");
  _passing_events_tree->Branch("ubxsec_muon_cos_theta",                                          &ubxsec_muon_cos_theta,                               "ubxsec_muon_cos_theta/D");
  _passing_events_tree->Branch("fails_fiducial_volume_req",                                      &fails_fiducial_volume_req,                            "fails_fiducial_volume_req/I");
  _passing_events_tree->Branch("muon_candidate_residuals",                                       &muon_candidate_residuals,                             "muon_candidate_residuals/D");
  _passing_events_tree->Branch("fails_residuals_req",                                            &fails_residuals_req,                                  "fails_residuals_req/I");
  _passing_events_tree->Branch("muon_candidate_percent_used_hits_in_cluster",                    &muon_candidate_percent_used_hits_in_cluster,          "muon_candidate_percent_used_hits_in_cluster/D");
  _passing_events_tree->Branch("fails_perc_used_hits_in_cluster",                                &fails_perc_used_hits_in_cluster,                      "fails_perc_used_hits_in_cluster/I");

  Truth_Tree = fs->make<TTree>("Truth_Tree", "A tree for the study of the truth distributions of signal events before any cuts");
  Truth_Tree->Branch("run",                                                            &_run,                                                "run/I");
  Truth_Tree->Branch("subrun",                                                         &_subrun,                                             "subrun/I");
  Truth_Tree->Branch("event",                                                          &_event,                                              "event/I");
  Truth_Tree->Branch("spline_fix_mcweight",                                            &spline_fix_mcweight,                                 "spline_fix_mcweight/D");
  Truth_Tree->Branch("central_value_mcweight",                                         &central_value_mcweight,                              "central_value_mcweight/D");
  Truth_Tree->Branch("rootino_fix_mcweight",                                           &rootino_fix_mcweight,                                "rootino_fix_mcweight/D");
  Truth_Tree->Branch("ppfx_value_mcweight",                                            &ppfx_value_mcweight,                                 "ppfx_value_mcweight/D");
  Truth_Tree->Branch("flux_other_universe_mcweights",                                  &flux_other_universe_mcweights,                       "flux_other_universe_mcweights[600]/D");
  Truth_Tree->Branch("other_universe_mcweights",                                       &other_universe_mcweights,                            "other_universe_mcweights[100]/D");
  Truth_Tree->Branch("axialff_mcweights",                                              &axialff_mcweights,                                   "axialff_mcweights[2]/D");
  Truth_Tree->Branch("rpaccqe_mcweights",                                              &rpaccqe_mcweights,                                   "rpaccqe_mcweights[2]/D");
  Truth_Tree->Branch("xsecshape_mcweights",                                            &xsecshape_mcweights,                                 "xsecshape_mcweights[2]/D");
  Truth_Tree->Branch("parent_fvx",                                                     &parent_fvx,                                          "parent_fvx/D");
  Truth_Tree->Branch("parent_fvy",                                                     &parent_fvy,                                          "parent_fvy/D");
  Truth_Tree->Branch("parent_fvz",                                                     &parent_fvz,                                          "parent_fvz/D");
  Truth_Tree->Branch("neutrino_PDG",                                                   &neutrino_PDG,                                        "neutrino_PDG/I");
  Truth_Tree->Branch("NC_channel",                                                     &NC_channel,                                          "NC_channel/I");
  Truth_Tree->Branch("neutrino_interaction_mode",                                      &neutrino_interaction_mode,                           "neutrino_interaction_mode/I");
  Truth_Tree->Branch("nu_vtx_x_truth",                                                 &nu_vtx_x_truth,                                      "nu_vtx_x_truth/D");
  Truth_Tree->Branch("nu_vtx_y_truth",                                                 &nu_vtx_y_truth,                                      "nu_vtx_y_truth/D");
  Truth_Tree->Branch("nu_vtx_z_truth",                                                 &nu_vtx_z_truth,                                      "nu_vtx_z_truth/D");
  Truth_Tree->Branch("neutrino_energy",                                                &neutrino_energy,                                     "neutrino_energy/D");
  Truth_Tree->Branch("neutrino_px",                                                    &neutrino_px,                                         "neutrino_px/D");
  Truth_Tree->Branch("neutrino_py",                                                    &neutrino_py,                                         "neutrino_py/D");
  Truth_Tree->Branch("neutrino_pz",                                                    &neutrino_pz,                                         "neutrino_pz/D");
  Truth_Tree->Branch("num_muminus_tracks",                                             &num_muminus_tracks,                                  "num_muminus_tracks/I");
  Truth_Tree->Branch("num_muplus_tracks",                                              &num_muplus_tracks,                                   "num_muplus_tracks/I");
  Truth_Tree->Branch("num_piplus_tracks",                                              &num_piplus_tracks,                                   "num_piplus_tracks/I");
  Truth_Tree->Branch("num_piminus_tracks",                                             &num_piminus_tracks,                                  "num_piminus_tracks/I");
  Truth_Tree->Branch("num_pi0_tracks",                                                 &num_pi0_tracks,                                      "num_pi0_tracks/I");
  Truth_Tree->Branch("num_proton_tracks",                                              &num_proton_tracks,                                   "num_proton_tracks/I");
  Truth_Tree->Branch("num_electron_showers",                                           &num_electron_showers,                                "num_electron_showers/I");
  Truth_Tree->Branch("num_positron_showers",                                           &num_positron_showers,                                "num_positron_showers/I");
  Truth_Tree->Branch("num_photon_showers",                                             &num_photon_showers,                                  "num_photon_showers/I");
  Truth_Tree->Branch("num_mctrack_points",                                             &num_mctrack_points,                                  "num_mctrack_points/I");
  Truth_Tree->Branch("truth_muon_starting_x_coordinate",                               &truth_muon_starting_x_coordinate,                    "truth_muon_starting_x_coordinate/D");
  Truth_Tree->Branch("truth_muon_starting_y_coordinate",                               &truth_muon_starting_y_coordinate,                    "truth_muon_starting_y_coordinate/D");
  Truth_Tree->Branch("truth_muon_starting_z_coordinate",                               &truth_muon_starting_z_coordinate,                    "truth_muon_starting_z_coordinate/D");
  Truth_Tree->Branch("truth_muon_ending_x_coordinate",                                 &truth_muon_ending_x_coordinate,                      "truth_muon_ending_x_coordinate/D");
  Truth_Tree->Branch("truth_muon_ending_y_coordinate",                                 &truth_muon_ending_y_coordinate,                      "truth_muon_ending_y_coordinate/D");
  Truth_Tree->Branch("truth_muon_ending_z_coordinate",                                 &truth_muon_ending_z_coordinate,                      "truth_muon_ending_z_coordinate/D");
  Truth_Tree->Branch("truth_muon_length",                                              &truth_muon_length,                                   "truth_muon_length/D");
  Truth_Tree->Branch("truth_muon_track_is_contained",                                  &truth_muon_track_is_contained,                       "truth_muon_track_is_contained/I");
  Truth_Tree->Branch("truth_muon_px",                                                  &truth_muon_px,                                       "truth_muon_px/D");
  Truth_Tree->Branch("truth_muon_py",                                                  &truth_muon_py,                                       "truth_muon_py/D");
  Truth_Tree->Branch("truth_muon_pz",                                                  &truth_muon_pz,                                       "truth_muon_pz/D");
  Truth_Tree->Branch("truth_muon_kinetic_energy",                                      &truth_muon_kinetic_energy,                           "truth_muon_kinetic_energy/D");

  if(_make_pida_csv) _csvfile.open ("pida_trklen.csv", std::ofstream::out | std::ofstream::trunc);
  if(_make_pida_csv) _csvfile << "pida,trklen,y" << std::endl;

  if(_make_ophit_csv) _csvfile2.open("ophit.csv", std::ofstream::out | std::ofstream::trunc);
  if(_make_ophit_csv) _csvfile2 << "ophit,opdet,time,pe" << std::endl;

  // Set the variables for the total numbers of events to 0.                                                                                                                                             
  total_num_events                                                        = 0;
  total_num_events_failing                                                = 0;
  total_num_events_passing                                                = 0;

  total_num_events_weighted                                               = 0.;
  total_num_events_passing_weighted                                       = 0.;
  total_num_events_failing_weighted                                       = 0.;

  total_num_events_failed_beam_disc_flashes_weighted                      = 0.;
  total_num_events_failed_beam_spill_flash_weighted                       = 0.;
  total_num_events_failed_has_slices_weighted                             = 0.;
  total_num_events_failed_fiducial_volume_weighted                        = 0.;
  total_num_events_failed_ntrack_weighted                                 = 0.;
  total_num_events_failed_residuals_std_up_weighted                       = 0.;
  total_num_events_failed_perc_used_hits_in_cluster_weighted              = 0.;
  
  // Set all of these variables equal to 0.                                                                                                                                                             
  total_num_events_failed_beam_disc_flashes                               = 0;
  total_num_events_failed_beam_spill_flash                                = 0;
  total_num_events_failed_has_slices                                      = 0;
  total_num_events_failed_track_length                                    = 0;
  total_num_events_failed_has_slice_tagged_as_neutrino                    = 0;
  total_num_events_failed_fiducial_volume                                 = 0;
  total_num_events_failed_ntrack                                          = 0;
  total_num_events_failed_residuals_std_up                                = 0;
  total_num_events_failed_perc_used_hits_in_cluster                       = 0;

  produces<std::vector<ubana::SelectionResult>>();
  produces<art::Assns<ubana::SelectionResult, ubana::TPCObject>>();

  // For the neutrino id filter
  produces<art::Assns<recob::Vertex, recob::Track>>();
  produces< art::Assns<recob::Vertex, recob::PFParticle>>();

  event_counter  = 0;

}



void UBXSec::produce(art::Event & e) {

  std::cout << "Currently looping over event #" << event_counter << "." << std::endl;

  // Reset the event counter;
  event_counter++;

  // Reset the channel info.
  NC_channel                                                 = 0;

  // Declare the tool for the track momentum calculator.
  trkf::TrackMomentumCalculator p_calculator_from_length;

  // Include the object for the Space Charge effect service.
  spacecharge::SpaceChargeMicroBooNE const* sce =
    reinterpret_cast<spacecharge::SpaceChargeMicroBooNE const*>(lar::providerFrom<spacecharge::SpaceChargeService>());

  // Declare the object for using the offsets that are determined by data.
  auto const* SCE_data = lar::providerFrom<spacecharge::SpaceChargeService>();

  // Load the wire objects that I will need to use to find charge in the vicinity of the vertex.                                                                                                       
  art::Handle<std::vector<recob::Hit> > hit_h;
  e.getByLabel("gaushit",hit_h);

  // make sure the hit objects look good                                                                                                                                                                  
  if(!hit_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Hit Objects!"<<std::endl;
    throw std::exception();
  }

  // Load the wire objects that I will need to use to find charge in the vicinity of the vertex.
  art::Handle<std::vector<recob::Wire> > wire_h;
  e.getByLabel("butcher", wire_h);

  // make sure wire objects look good                                                                                                                                                                     
  if(!wire_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Wire Objects!"<<std::endl;
    throw std::exception();
  }

  // Load the tracks from pandora.
  art::Handle<std::vector<recob::Track> > pandora_track_h;
  e.getByLabel("pandora",pandora_track_h);

  // make sure pandora tracks look good                                                                                                                                                                    
  if(!pandora_track_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate pandora track objects!"<<std::endl;
    throw std::exception();
  }

  std::vector< art::Ptr<recob::Track> > pandora_track_pointers_h;
  art::fill_ptr_vector(pandora_track_pointers_h, pandora_track_h);

  // Declare the vector for the MCS objects with the SCE correction.
  art::Handle<std::vector<recob::MCSFitResult> > pandora_MCSResult_h;
  e.getByLabel("pandoraMCSMu", pandora_MCSResult_h);

  // make sure pandora tracks look good                                                                                                                                                                    
  if(!pandora_MCSResult_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate pandora MCSFitResult objects!"<<std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<recob::MCSFitResult>> pandora_MCSResult_pointers_h;
  art::fill_ptr_vector(pandora_MCSResult_pointers_h, pandora_MCSResult_h);

  // Declare the handle full of track hits.                                                                                                                                                              
  art::FindManyP<recob::Hit> hits_per_track(pandora_track_h, e, "pandora");

  // Load the hit information for these tracks.
  // Declare the association between tracks and hits as well.                                                                                                                                             
  art::FindMany<recob::Hit> trk_hit_assn_v(pandora_track_h, e, "pandora");

  // Load the calorimetry information for these tracks.
  art::FindMany<anab::Calorimetry> trk_calo_assn_v( pandora_track_h, e, "pandoracaliSCE" );

  // Use this to find the calorimetry objects that belong to the muon track in the event but without the space charge correction.                                                                          
  art::FindMany<anab::Calorimetry> trk_calo_assn_v_no_SCE_corrections( pandora_track_h, e, "pandoracali" );

  // Add in the information for removing the bad hits from the calorimetry information.
  auto const & track_list_ptr = e.getValidHandle<std::vector <recob::Track> >("pandora");
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(track_list_ptr, e, "pandora");

  // Include the variables needed to fill track PID.
  lar_pandora::LArPandoraHelper    larpandora;
  lar_pandora::PFParticlesToTracks particlesToTracks;
  lar_pandora::TrackVector         pftracks;
  PandoraInterfaceHelper           pandoraInterfaceHelper;
  TrackHelper                      trackHelper;
  std::vector<double> TrackScore_;
  TrackScore_.clear();
  std::vector<double> TrackPID_chiproton_;
  TrackPID_chiproton_.clear();

  std::vector<double> TrackStart_x_sce_;
  std::vector<double> TrackStart_y_sce_;
  std::vector<double> TrackStart_z_sce_;
  std::vector<double> TrackEnd_x_sce_;
  std::vector<double> TrackEnd_y_sce_;
  std::vector<double> TrackEnd_z_sce_;

  TrackStart_x_sce_.clear();
  TrackStart_y_sce_.clear();
  TrackStart_z_sce_.clear();
  TrackEnd_x_sce_.clear();
  TrackEnd_y_sce_.clear();
  TrackEnd_z_sce_.clear();

  // Collect the tracks using the larpandora method.
  larpandora.CollectTracks(e, "pandora", pftracks, particlesToTracks);
  
  const art::FindManyP<anab::ParticleID> trackPIDAssn(pandora_track_h, e, "pandoracalipidSCE");
  if (!trackPIDAssn.isValid()){
    std::cout << "[NumuCCana::getVertex] Event failed: PID is invalid" << std::endl;
  }

  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleVector pfdaughters;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  
  larpandora.CollectPFParticleMetadata(e, "pandora", pfparticles, particlesToMetadata);

  larpandora.BuildPFParticleMap(pfparticles, particleMap);

  if ( pfparticles.size() > 0 ) {

    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);

    }

  if ( pfneutrinos.size() > 0 ) {

    art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();

    pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);

  }

  // Repeat this same procedure with the other systematic weights.
  art::Handle<std::vector<evwgh::MCEventWeight> > genie_eventweight_h;
  e.getByLabel( "eventweight::EventWeightMar18" , genie_eventweight_h );

  if(!genie_eventweight_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate 'eventweight' EventWeight!"<<std::endl;
    throw std::exception();
  }

  std::cout << "The size of genie_eventweight_h = " << genie_eventweight_h->size() << "." << std::endl;
  
  for ( size_t i = 0; i < genie_eventweight_h->size(); i++ ) {

    auto const& mc_weight = genie_eventweight_h->at(i);

    std::map<std::string, std::vector<double> > weight_map = mc_weight.fWeight;

    std::cout << "Above printing out all of the values in the eventweight map."<< std::endl;

    for ( std::map<std::string, std::vector<double>>::iterator it = weight_map.begin(); it != weight_map.end(); it++ ) 
      {
	std::cout << it->first  // string (key)
	  //		  << ':'
	  //  << it.second   // string's value 
		  << std::endl;

      }

    std::cout << "Below printing out all of the values in the eventweight map." << std::endl;

    const std::vector<double>& spline_fix_weights = weight_map.at("splines_general_Spline");

    spline_fix_mcweight = spline_fix_weights.front();

    if ( std::isnan( spline_fix_mcweight ) || std::isinf( spline_fix_mcweight ) )
      spline_fix_mcweight = 1.;

    const std::vector<double>& central_value_weights = weight_map.at("TunedCentralValue_UBGenie");

    central_value_mcweight = central_value_weights.front();

    if ( std::isnan( central_value_mcweight ) || std::isinf( central_value_mcweight ) )
      central_value_mcweight = 1.;

    const std::vector<double>& rootino_fix_weights = weight_map.at("RootinoFix_UBGenie");

    rootino_fix_mcweight = rootino_fix_weights.front();

    if ( std::isnan( rootino_fix_mcweight ) || std::isinf( rootino_fix_mcweight ) )
      rootino_fix_mcweight = 1.;

    std::cout << "About to implement the ppfx weight." << std::endl;

    const std::vector<double>& ppfx_value_weights = weight_map.at("ppfx_cv_UBPPFXCV");
    
    ppfx_value_mcweight  = ppfx_value_weights.front();

    if ( std::isnan( ppfx_value_mcweight ) || std::isinf( ppfx_value_mcweight ) )
      ppfx_value_mcweight = 1.;

    std::cout << "After implementing the ppfx weight." << std::endl;

    std::cout << "About to implement the other flux universes weight." << std::endl;

    const std::vector<double>& other_flux_universes_weight = weight_map.at("ppfx_ms_UBPPFX");

    std::cout << "The size of 'other_flux_universes_weight' = " << other_flux_universes_weight.size() << "." << std::endl;

    for ( size_t i = 0; i < 100; i++ ) {
      
      flux_other_universe_mcweights[i] = other_flux_universes_weight.at( i );

      std::cout << "The value in flux other universe #" << i << " = " << flux_other_universe_mcweights[i] << "." << std::endl;

    }

    const std::vector<double>& universe_weights = weight_map.at("All_UBGenie");

    for ( size_t i = 0; i < 100; i++ ) {

      other_universe_mcweights[i] = universe_weights.at( i );

    }
    
    const std::vector<double>& axialff_weights = weight_map.at("AxFFCCQEshape_UBGenie");

    for( size_t i = 0;i < axialff_weights.size(); i++ ) {

      axialff_mcweights[i] = axialff_weights.at( i );

    }   

    const std::vector<double>& rpaccqe_weights = weight_map.at("RPA_CCQE_UBGenie");

    std::cout << "Size of 'RPA_CCQE' = " << rpaccqe_weights.size() << "." << std::endl;

    for( size_t i = 0;i < rpaccqe_weights.size(); i++ ) {

      rpaccqe_mcweights[i] = rpaccqe_weights.at( i );

    }   

    const std::vector<double>& xsecshape_weights = weight_map.at("XSecShape_CCMEC_UBGenie");

    for( size_t i = 0; i < xsecshape_weights.size(); i++ ) {

      xsecshape_mcweights[i] = xsecshape_weights.at( i );

    }   

  }

  // Include the dk2nu weights at this point in the loop.  
  total_num_events_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
  total_num_events++;

  // load neutrino information.                                                                                                                                                                        
  if (_debug) { std::cout << "loading neutrino from producer generator." << std::endl; }
  art::Handle<std::vector<simb::MCTruth> > neutrino_h;
  e.getByLabel("generator", neutrino_h);

  // make sure MCTruth info looks good                                                                                                                                                                   
  if(!neutrino_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Neutrino!"<<std::endl;
    throw std::exception();
  }

  // Declare a vector of pointers for the neutrino object.                   
  auto neutrino = neutrino_h->at(0).GetNeutrino();
  auto nu = neutrino.Nu();                      

  // Find the direction of the truth track to be used in the weighting. 
  double truth_px_for_theta_for_weighting = nu.Px();
  double truth_py_for_theta_for_weighting = nu.Py();
  double truth_pz_for_theta_for_weighting = nu.Pz();

  // Variables
  TRotation RotDet2Beam;             // Rotations
  TVector3  detxyz, BeamCoords;      // Translations
  std::vector<double> rotmatrix;     // Inputs

  // input detector coordinates to translate
  detxyz = {truth_px_for_theta_for_weighting, truth_py_for_theta_for_weighting, truth_pz_for_theta_for_weighting};     

  // From beam to detector rotation matrix
  rotmatrix = {
    0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
    4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
    -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

  // Return the TRotation
  TVector3 newX, newY, newZ;
  newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
  newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
  newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

  RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam

  // Rotate to beam coords
  BeamCoords = RotDet2Beam * detxyz;

  TVector3 beam_dir = {0 , 0 , 1};
  truth_theta_angle_for_weighting = BeamCoords.Angle(beam_dir) * 180 / 3.1415926;

  std::vector<art::Ptr<simb::MCTruth> >NeutrinoVec;
  art::fill_ptr_vector(NeutrinoVec, neutrino_h);

  // load the truth MCFlux information.                                                                                                                                                                
  if (_debug) { std::cout << "loading MCFlux object from producer generator." << std::endl; }
  art::Handle<std::vector<simb::MCFlux> > parent_h;
  e.getByLabel("generator", parent_h);

  // make sure MCFlux info looks good                                                                                                                                                                  
  if(!parent_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Parent!" << std::endl;
    throw std::exception();
  }

  // Declare a vector of pointers with the 'MCFlux' object.                                                                                                                                            
  std::vector<art::Ptr<simb::MCFlux> > ParentVec;
  art::fill_ptr_vector(ParentVec, parent_h);

  for (auto& parent : ParentVec){

    parent_fvx   = parent->fvx;
    parent_fvy   = parent->fvy;
    parent_fvz   = parent->fvz;

  }

  num_of_neutrinos = 0;
  
  // loop through neutrinos themselves.                                                                                                                                                                 
  for (auto& neutrino : NeutrinoVec ) {
      
    // Unpack the neutrino object to find an MCParticle.                                                                                                                                                
    const simb::MCNeutrino& truth_neutrino = neutrino->GetNeutrino();
    const simb::MCParticle& truth_particle = truth_neutrino.Nu();

    // Unpack the coordinates for the vertex as well.                                                                                                                                                  
    nu_vtx_x_truth                         = truth_particle.Vx(0);
    nu_vtx_y_truth                         = truth_particle.Vy(0);
    nu_vtx_z_truth                         = truth_particle.Vz(0);
    
    neutrino_PDG                           = truth_particle.PdgCode(); 
    neutrino_energy                        = ( truth_particle.E(0) * 1000. );
    NC_channel                             = truth_neutrino.CCNC();
    neutrino_interaction_mode              = truth_neutrino.Mode();

    truth_vertex_is_in_AV = 0;
      
    if ( nu_vtx_x_truth > 0. && nu_vtx_x_truth < 256.35 && nu_vtx_y_truth > -116.5 && nu_vtx_y_truth < 116.5 && nu_vtx_z_truth > 0. && nu_vtx_z_truth < 1036.8 )
      truth_vertex_is_in_AV = 1;

    truth_vertex_is_in_FV = 0;

    if ( nu_vtx_x_truth > 20. && nu_vtx_x_truth < 236.35 && nu_vtx_y_truth > -96.5 && nu_vtx_y_truth < 96.5 && nu_vtx_z_truth > 20. && nu_vtx_z_truth < 1016.8 )
      truth_vertex_is_in_FV = 1;

    num_of_neutrinos++;
      
  } // End of the loop over the neutrino candidates in the event.                         

  // Include info about the MCParticles.
  // Load in the MCParticle information.                                                                                                                                                                  
  art::Handle<std::vector<simb::MCParticle> > mcparticle_h;
  e.getByLabel("largeant", mcparticle_h);

  // make sure MCParticle info looks good                                                                                                                                                                
  if(!mcparticle_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MCParticle!"<<std::endl;
    throw std::exception();
  }

  for ( size_t mcparticle_iter = 0; mcparticle_iter < mcparticle_h->size(); mcparticle_iter++ ) {

    if ( mcparticle_h->at( mcparticle_iter ).PdgCode() != 14 && mcparticle_h->at( mcparticle_iter ).PdgCode() != -14 )
      continue;

    if ( mcparticle_h->at( mcparticle_iter ).PdgCode() == 14 || mcparticle_h->at( mcparticle_iter ).PdgCode() == -14 ) 
      std::cout << "MCParticle PDG = " << mcparticle_h->at( mcparticle_iter ).PdgCode() << "." << std::endl;

    neutrino_px                            = mcparticle_h->at( mcparticle_iter ).Px(0);
    neutrino_py                            = mcparticle_h->at( mcparticle_iter ).Py(0);
    neutrino_pz                            = mcparticle_h->at( mcparticle_iter ).Pz(0);

    double normalized_vector               = TMath::Sqrt( neutrino_px * neutrino_px + neutrino_py * neutrino_py + neutrino_pz * neutrino_pz );

    neutrino_px                            = neutrino_px / normalized_vector;
    neutrino_py                            = neutrino_py / normalized_vector;
    neutrino_pz                            = neutrino_pz / normalized_vector;

  }

  // Load in the MCTrack information.
  art::Handle<std::vector<sim::MCTrack> > mctrack_h;
  e.getByLabel("mcreco", mctrack_h);

  // make sure MCTrack info looks good                                                                                                                                                                 
  if(!mctrack_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MCTrack!"<<std::endl;
    throw std::exception();
  }

  num_muminus_tracks        = 0;
  num_muplus_tracks         = 0;
  num_piplus_tracks         = 0;
  num_piminus_tracks        = 0;
  num_pi0_tracks            = 0;
  num_proton_tracks         = 0;

  truth_muon_kinetic_energy = 0.;

  std::cout << "Truth neutrino energy = " << neutrino_energy * 1000. << " MeV." << std::endl;

  for ( size_t mctrack_iter = 0; mctrack_iter < mctrack_h->size(); mctrack_iter++ ) {

    // At the top of the loop, find out which type of track this is.
    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 13 ) {

      num_muminus_tracks++;

      num_mctrack_points = mctrack_h->at( mctrack_iter ).size();

      std::cout << "Number of muminus mctrack points = " << num_mctrack_points << "." << std::endl;

      if ( num_mctrack_points > 0 ) {

	if ( fabs( nu_vtx_x_truth - mctrack_h->at( mctrack_iter ).at( 0 ).X() ) < 0.001 && fabs( nu_vtx_y_truth - mctrack_h->at( mctrack_iter ).at( 0 ).Y() ) < 0.001 && fabs( nu_vtx_z_truth - mctrack_h->at( mctrack_iter ).at( 0 ).Z() ) < 0.001 ) {

	  truth_muon_kinetic_energy        = ( mctrack_h->at(mctrack_iter ).at( 0 ).E() - 105.7 );

	  truth_muon_px                    = mctrack_h->at( mctrack_iter ).at( 0 ).Px();
	  truth_muon_py                    = mctrack_h->at( mctrack_iter ).at( 0 ).Py();
	  truth_muon_pz                    = mctrack_h->at( mctrack_iter ).at( 0 ).Pz();
	
	  truth_muon_starting_x_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).X();
	  truth_muon_starting_y_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).Y();
	  truth_muon_starting_z_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).Z();
	  truth_muon_ending_x_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).X();
	  truth_muon_ending_y_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).Y();
	  truth_muon_ending_z_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).Z();
	  
	  truth_muon_length                = TMath::Sqrt( ( truth_muon_starting_x_coordinate - truth_muon_ending_x_coordinate ) * ( truth_muon_starting_x_coordinate - truth_muon_ending_x_coordinate ) + ( truth_muon_starting_y_coordinate - truth_muon_ending_y_coordinate ) * ( truth_muon_starting_y_coordinate - truth_muon_ending_y_coordinate ) + ( truth_muon_starting_z_coordinate - truth_muon_ending_z_coordinate ) * ( truth_muon_starting_z_coordinate - truth_muon_ending_z_coordinate ) );
	
	  truth_muon_track_is_contained    = 1;

	  if ( truth_muon_starting_x_coordinate < 3.0 || truth_muon_starting_x_coordinate > 253.35 || truth_muon_starting_y_coordinate < -113.5 || truth_muon_starting_y_coordinate > 113.5 || truth_muon_starting_z_coordinate < 3.0 || truth_muon_starting_z_coordinate > 1033.8 || truth_muon_ending_x_coordinate < 3.0 || truth_muon_ending_x_coordinate > 253.35 || truth_muon_ending_y_coordinate < -113.5 ||truth_muon_ending_y_coordinate > 113.5 || truth_muon_ending_z_coordinate < 3.0 || truth_muon_ending_z_coordinate > 1033.8 )
	    truth_muon_track_is_contained = 0;

	  double p_resultant = TMath::Sqrt( truth_muon_px * truth_muon_px + truth_muon_py * truth_muon_py + truth_muon_pz * truth_muon_pz );

	  std::cout << "Truth lepton PDG Code = 13." << std::endl;
	  std::cout << "Truth muon energy = " << mctrack_h->at(mctrack_iter ).at( 0 ).E() << " MeV." << std::endl;
	  std::cout << "Truth muon kinetic energy = " << truth_muon_kinetic_energy << " MeV." << std::endl;
	  std::cout << "Truth momentum resultant = " << p_resultant << " MeV/c." << std::endl;
	  std::cout << "Truth energy from momentum = " << TMath::Sqrt( p_resultant * p_resultant + 105.7 * 105.7 ) << " MeV." << std::endl;

	  std::cout << "Number of mctrack points = " << num_mctrack_points << "." << std::endl;
	  std::cout << "Truth muminus starting coordinates: x = " << truth_muon_starting_x_coordinate << " cm y = " << truth_muon_starting_y_coordinate << " cm z = " << truth_muon_starting_z_coordinate << " cm." << std::endl;
	  std::cout << "Truth muminus ending coordinates: x = " << truth_muon_ending_x_coordinate << " cm y = " << truth_muon_ending_y_coordinate << " cm z = " << truth_muon_ending_z_coordinate << " cm." << std::endl;
	  std::cout << "Truth muminus length = " << truth_muon_length << " cm." << std::endl;
	  std::cout << "Truth muon track is contained = " << truth_muon_track_is_contained << "." << std::endl;

	}

      }

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == -13 ) {

      num_muplus_tracks++;
      
      num_mctrack_points = mctrack_h->at( mctrack_iter ).size();

      if ( num_mctrack_points > 0 ) {

	// Use this conditional for the possibility of a pion decay to a muplus.
	if ( ( mctrack_h->at(mctrack_iter ).at( 0 ).E() - 105.7 ) > 5.0 && fabs( nu_vtx_x_truth - mctrack_h->at( mctrack_iter ).at( 0 ).X() ) < 0.001  && fabs( nu_vtx_y_truth - mctrack_h->at( mctrack_iter ).at( 0 ).Y() ) < 0.001 && fabs( nu_vtx_z_truth - mctrack_h->at( mctrack_iter ).at( 0 ).Z() ) < 0.001 ) {

	  truth_muon_kinetic_energy = ( mctrack_h->at( mctrack_iter ).at( 0 ).E() - 105.7 );

	  truth_muon_px             = mctrack_h->at( mctrack_iter ).at( 0 ).Px();
	  truth_muon_py             = mctrack_h->at( mctrack_iter ).at( 0 ).Py();
	  truth_muon_pz             = mctrack_h->at( mctrack_iter ).at( 0 ).Pz();

	  double p_resultant        = TMath::Sqrt( truth_muon_px * truth_muon_px + truth_muon_py * truth_muon_py + truth_muon_pz * truth_muon_pz );
	
	  std::cout << "Truth lepton PDG Code = -13." << std::endl;
	  std::cout << "Truth muon energy = " << mctrack_h->at(mctrack_iter ).at( 0 ).E() << " MeV." << std::endl;
	  std::cout << "Truth muon kinetic energy = " << truth_muon_kinetic_energy << " MeV." << std::endl;
	  std::cout << "Truth momentum resultant = " << p_resultant << " MeV/c." << std::endl;
	  std::cout << "Truth energy from momentum = " <<TMath::Sqrt( p_resultant * p_resultant + 105.7 * 105.7 ) << " MeV." << std::endl;
	  
	  truth_muon_starting_x_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).X();
	  truth_muon_starting_y_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).Y();
	  truth_muon_starting_z_coordinate = mctrack_h->at( mctrack_iter ).at( 0 ).Z();
	  truth_muon_ending_x_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).X();
	  truth_muon_ending_y_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).Y();
	  truth_muon_ending_z_coordinate   = mctrack_h->at( mctrack_iter ).at( mctrack_h->at( mctrack_iter ).size() - 1 ).Z();

	  truth_muon_length                = TMath::Sqrt( ( truth_muon_starting_x_coordinate - truth_muon_ending_x_coordinate ) * ( truth_muon_starting_x_coordinate - truth_muon_ending_x_coordinate ) + ( truth_muon_starting_y_coordinate - truth_muon_ending_y_coordinate ) * ( truth_muon_starting_y_coordinate - truth_muon_ending_y_coordinate ) + ( truth_muon_starting_z_coordinate - truth_muon_ending_z_coordinate ) * ( truth_muon_starting_z_coordinate - truth_muon_ending_z_coordinate ) );								    

	  truth_muon_track_is_contained    = 1;

	  if ( truth_muon_starting_x_coordinate < 3.0 || truth_muon_starting_x_coordinate > 253.35 || truth_muon_starting_y_coordinate < -113.5 || truth_muon_starting_y_coordinate > 113.5 || truth_muon_starting_z_coordinate < 3.0 || truth_muon_starting_z_coordinate > 1033.8 || truth_muon_ending_x_coordinate < 3.0 || truth_muon_ending_x_coordinate > 253.35 || truth_muon_ending_y_coordinate < -113.5 ||truth_muon_ending_y_coordinate > 113.5 || truth_muon_ending_z_coordinate < 3.0 || truth_muon_ending_z_coordinate > 1033.8 )
	    truth_muon_track_is_contained = 0;

	  std::cout << "Number of mctrack points = " << num_mctrack_points << "."<< std::endl;
	  std::cout << "Truth muplus starting coordinates: x = " << truth_muon_starting_x_coordinate << " cm y = " << truth_muon_starting_y_coordinate << " cm z = " << truth_muon_starting_z_coordinate << " cm." << std::endl;
	  std::cout << "Truth muplus ending coordinates: x = " << truth_muon_ending_x_coordinate << " cm y = " << truth_muon_ending_y_coordinate << " cm z = " << truth_muon_ending_z_coordinate << " cm." << std::endl;
	  std::cout << "Truth muplus length = " << truth_muon_length << " cm." << std::endl;
	  std::cout << "Truth muon track is contained = " << truth_muon_track_is_contained << "." << std::endl;

	}
	  
      }

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 111 ) {

      num_pi0_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 211 ) {

      num_piplus_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == -211 ) {

      num_piminus_tracks++;

    }

    if ( mctrack_h->at( mctrack_iter ).PdgCode() == 2212 ) {

      num_proton_tracks++;

    }

  }

  num_electron_showers = 0;
  num_positron_showers = 0;
  num_photon_showers   = 0;

  // Load in the MCShower information.
  art::Handle<std::vector<sim::MCShower> > mcshower_h;
  e.getByLabel("mcreco", mcshower_h);

  // make sure MCShower info looks good                                                                                                                                                                 
  if(!mcshower_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate MCShower!"<<std::endl;
    throw std::exception();
  }

  for ( size_t mcshower_iter = 0; mcshower_iter < mcshower_h->size(); mcshower_iter++ ) {

    if ( mcshower_h->at( mcshower_iter ).PdgCode() == 11 ) {

      num_electron_showers++;

    }

    if ( mcshower_h->at( mcshower_iter ).PdgCode() == -11 ) {

      num_positron_showers++;

    }

    if ( mcshower_h->at( mcshower_iter ).PdgCode() == 22  ) {

      num_photon_showers++;

    }

  }

  if(_debug) std::cout << "********** UBXSec starts" << std::endl;
  if(_debug) std::cout << "[UBXSec] Run: "           << e.id().run()    <<
                          ", subRun: "               << e.id().subRun() <<
                          ", event: "                << e.id().event()  << std::endl;

  if (_do_opdet_swap && e.isRealData()) {
    std::cout << "[UBXSec] WARNING!!! Swapping OpDets. I hope you know what you are doing." << std::endl;
  }

  // Instantiate the output
  std::unique_ptr< std::vector<ubana::SelectionResult>>                   selectionResultVector           (new std::vector<ubana::SelectionResult>);
  std::unique_ptr< art::Assns<ubana::SelectionResult, ubana::TPCObject>>  assnOutSelectionResultTPCObject (new art::Assns<ubana::SelectionResult, ubana::TPCObject>);

  std::unique_ptr<art::Assns<recob::Vertex, recob::Track>>      vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
  std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);

  // Initialize the UBXSecEvent
  ubxsec_event->Init();

  _run    = ubxsec_event->run    = e.id().run();
  _subrun = ubxsec_event->subrun = e.id().subRun();
  _event  = ubxsec_event->event  = e.id().event();

  _is_data = e.isRealData();
  _is_mc   = !_is_data;

  if ( NC_channel == 0 && ( num_muplus_tracks > 0 || num_muminus_tracks > 0 ) && truth_vertex_is_in_AV == 1 )
    Truth_Tree->Fill();

  //::art::ServiceHandle<cheat::BackTracker> bt;
  ::art::ServiceHandle<geo::Geometry> geo;

  // Prepare the dead region finder
  //std::cout << "[UBXSec] Recreate channel status map" << std::endl;
  const lariov::ChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  for (unsigned int ch = 0; ch < 8256; ch++) {
    deadRegionsFinder.SetChannelStatus(ch, chanFilt.Status(ch));
  }
  //std::cout << "[UBXSec] Now force reload BWires" << std::endl;
  deadRegionsFinder.CreateBWires();

  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  _drift_velocity = _detp -> DriftVelocity(efield,temp);
  if (_debug) std::cout << "[UBXSec] Using drift velocity = " << _drift_velocity << " cm/us, with E = " << efield << ", and T = " << temp << std::endl;

 lar_pandora::PFParticleVector pfParticleList;              //vector of PFParticles
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, pfParticleList);

  // Collect tracks
  lar_pandora::TrackVector            allPfParticleTracks;
  lar_pandora::PFParticlesToTracks    pfParticleToTrackMap;
  lar_pandora::TracksToHits           trackToHitsMap;
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, pfParticleToTrackMap);
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, allPfParticleTracks, trackToHitsMap);

  // Collect showers
  lar_pandora::ShowerVector           _shower_v;
  lar_pandora::PFParticlesToShowers   _pfp_to_shower_map;
  lar_pandora::PFParticleVector _pfp_v;
  lar_pandora::LArPandoraHelper::CollectShowers(e, _pfp_producer, _shower_v, _pfp_to_shower_map);
  lar_pandora::PFParticlesToMetadata pfParticlesToMetadata;
  lar_pandora::LArPandoraHelper::CollectPFParticleMetadata(e, _pfp_producer, _pfp_v, pfParticlesToMetadata);
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kUseDaughters, true);

  // Get TPCObjects from the Event
  art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
  e.getByLabel(_tpcobject_producer, tpcobj_h);
  if (!tpcobj_h.isValid()) {
    std::cout << "[UBXSec] Cannote locate ubana::TPCObject." << std::endl;
  }
  art::FindManyP<ubana::FlashMatch> tpcobjToFlashMatchAssns(tpcobj_h, e, _neutrino_flash_match_producer);
  art::FindManyP<recob::Track>      tpcobjToTrackAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Shower>     tpcobjToShowerAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::PFParticle> tpcobjToPFPAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<recob::Vertex>     tpcobjToVertexAssns(tpcobj_h, e, _tpcobject_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToCosmicTagAssns(tpcobj_h, e, _geocosmictag_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToConsistency(tpcobj_h, e, _candidateconsistency_producer);
  art::FindManyP<anab::CosmicTag>   tpcobjToStopMu(tpcobj_h, e, _cosmic_stopmu_tag_producer);

  // Get Tracks
  art::Handle<std::vector<recob::Track>> track_h;
  e.getByLabel(_pfp_producer,track_h);
  if (!track_h.isValid() || track_h->empty()) {
    std::cout << "[UBXSec] Track handle is not valid or empty." << std::endl;
    //throw std::exception();
  }

  std::vector<art::Ptr<recob::Track>> track_p_v;

  art::fill_ptr_vector(track_p_v, track_h);

  art::FindManyP<recob::PFParticle> pfp_from_track(track_h, e, _pfp_producer);

  art::FindManyP<anab::Calorimetry> calos_from_track(track_h, e, _calorimetry_producer);

  // Get PFP
  art::Handle<std::vector<recob::PFParticle> > pfp_h;
  e.getByLabel(_pfp_producer,pfp_h);
  art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfp_h,e, _pfp_producer);

  if(!pfp_h.isValid()){
    std::cout << "[UBXSec] PFP product " << _pfp_producer << " not found..." << std::endl;
  }
  if(pfp_h->empty()) {
    std::cout << "[UBXSec] PFP "         << _pfp_producer << " is empty."    << std::endl;
  }
  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfp_h);
  art::FindManyP<recob::Track> tracks_from_pfp(pfp_h, e, _pfp_producer);

  ubxsec_event->n_pfp = ubxsec_event->n_pfp_primary = 0;
  for (size_t i = 0; i < pfp_h->size(); i++) {
    ubxsec_event->n_pfp++;
    if ((*pfp_h)[i].IsPrimary()&&abs((*pfp_h)[i].PdgCode()==14))
      ubxsec_event->n_pfp_primary++;

  }
 // Get Ghosts
  art::Handle<std::vector<ubana::MCGhost> > ghost_h;
  e.getByLabel(_mc_ghost_producer,ghost_h);
  if(!ghost_h.isValid()){
    std::cout << "[UBXSec] MCGhost product " << _mc_ghost_producer << " not found..." << std::endl;
  }
  art::FindManyP<ubana::MCGhost>   mcghost_from_pfp   (pfp_h,   e, _mc_ghost_producer);
  art::FindManyP<simb::MCParticle> mcpar_from_mcghost (ghost_h, e, _mc_ghost_producer);

    for (unsigned int i = 0; i < pfp_h->size(); ++i)
    {


      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(i));
        if (!pfParticleMetadataList.empty())
        {
            const art::Ptr<recob::PFParticle> pParticle(pfp_h, i);
            for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
            {
                const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
                const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
                if (!pfParticlePropertiesMap.empty())

                for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
		  if(it->first=="TrackScore")
		    {
		      std::cout << " Found PFParticle " << pParticle->Self() << " with: " << std::endl;
                    std::cout << "  - " << it->first << " = " << it->second << std::endl;
		    ubxsec_event->pfp_trackscore.push_back(it->second);
		    }
            }
        }
    }

  ////
 // Get PID information
  art::FindManyP<anab::ParticleID> particleids_from_track (track_h, e, _particle_id_producer);
  if (!particleids_from_track.isValid()) {
    std::cout << "[UBXSec] anab::ParticleID is not valid." << std::endl;
  }
  // Fill a std::map Track->ParticleID
  std::map<art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> track_to_pid_map;
  for (auto track : track_p_v) {
    std::vector<art::Ptr<anab::ParticleID>> pids = particleids_from_track.at(track.key());
    if(pids.size() == 0)
      continue;
    for (auto pid : pids) {
      // Check that pid object contains plane-2 PIDA (because that's all we care about)
      std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pid->ParticleIDAlgScores();
      for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
	anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
	int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	if (AlgScore.fVariableType==anab::kPIDA && planenum==2){
	  track_to_pid_map[track] = pid;
	  continue;
	}
      }
    }
  }
   // Get MCSFitResult - Muon
   art::Handle<std::vector<recob::MCSFitResult> > mcsfitresult_mu_h;
   e.getByLabel(_mcsfitresult_mu_producer,mcsfitresult_mu_h);
   if(!mcsfitresult_mu_h.isValid()){
     std::cout << "[UBXSec] MCSFitResult product " << _mcsfitresult_mu_producer << " not found..." << std::endl;
   }
   std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresult_mu_v;
   art::fill_ptr_vector(mcsfitresult_mu_v, mcsfitresult_mu_h);

  // pandoraCosmic PFPs (for cosmic removal studies)
  art::Handle<std::vector<recob::PFParticle>> pfp_cosmic_h;
  e.getByLabel("pandoraCosmic",pfp_cosmic_h);
  if(pfp_cosmic_h.isValid()){
    std::vector<art::Ptr<recob::PFParticle>> pfp_cosmic_v;
    art::fill_ptr_vector(pfp_cosmic_v, pfp_cosmic_h);
    ubxsec_event->n_primary_cosmic_pfp = 0;
    for (auto p : pfp_cosmic_v) {
      if (!p->IsPrimary()) continue;
      ubxsec_event->n_primary_cosmic_pfp++;
    }
  } else {
    std::cout << "[UBXSec] pandoraCosmic PFP product not found..." << std::endl;
  }

  // Flashes                                                                                                                                                                                             
  ::art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);

  flash_time                                                 = 0.;
  flash_PEs                                                  = 0.;
  flash_z                                                    = 0.;
  flash_y                                                    = 0.;
  flash_z_width                                              = 0.;
  flash_y_width                                              = 0.;

  _nflashes_in_beamgate                                      = beamflash_h->size();
  _nflashes_in_beamspill                                     = 0;
  _nflashes_in_beamspill_window_passing_filter_PE_cut        = 0;

  if ( _nflashes_in_beamgate > 0 ) { 
    
    flash_PEs = -1.;

    for ( size_t flash_iter = 0; flash_iter < beamflash_h->size(); flash_iter++ ) {

      if ( beamflash_h->at( flash_iter ).Time() > _beam_spill_start && beamflash_h->at( flash_iter ).Time() < _beam_spill_end ) { 

	std::cout << "Flash Time = " << beamflash_h->at( flash_iter ).Time() << " us." << std::endl;
	std::cout << "Flash PEs = " << beamflash_h->at( flash_iter ).TotalPE() << " PEs." << std::endl;

	_nflashes_in_beamspill++;

	if ( beamflash_h->at( flash_iter ).Time() > flash_PEs ) {

	  flash_time         = beamflash_h->at( flash_iter ).Time();
	  flash_PEs          = beamflash_h->at( flash_iter ).TotalPE();
	  flash_z            = beamflash_h->at( flash_iter ).ZCenter();
	  flash_y            = beamflash_h->at( flash_iter ).YCenter();
	  flash_z_width      = beamflash_h->at( flash_iter ).ZWidth();
	  flash_y_width      = beamflash_h->at( flash_iter ).YWidth();

	}

	// Cut for the entire KDAR filter.
	if ( beamflash_h->at( flash_iter ).TotalPE() > 50.0 ) 
	  _nflashes_in_beamspill_window_passing_filter_PE_cut++;
    
      }
  
    }

  }

  ubxsec_event->nbeamfls = beamflash_h->size();
  ubxsec_event->beamfls_pe.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_time.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_z.resize(ubxsec_event->nbeamfls);
  ubxsec_event->beamfls_spec.resize(ubxsec_event->nbeamfls);

  for (size_t n = 0; n < beamflash_h->size(); n++) {
    auto const& flash                          = (*beamflash_h)[n];
    ubxsec_event->beamfls_pe[n]                = flash.TotalPE();
    ubxsec_event->beamfls_time[n]              = flash.Time();
    ubxsec_event->beamfls_spec[n].resize(32);
    ubxsec_event->candidate_flash_time         = 0.;
    ubxsec_event->candidate_flash_z            = 0.;
    double min_pe = -1;
    for (unsigned int i = 0; i < 32; i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      if (_do_opdet_swap && e.isRealData()) {
        opdet = _opdet_swap_map.at(opdet);
      }

      ubxsec_event->beamfls_spec[n][opdet] = flash.PE(i);

      if (ubxsec_event->beamfls_time[n] > _beam_spill_start && ubxsec_event->beamfls_time[n] < _beam_spill_end) {
        // Find largest flash above threshold
        if (flash.TotalPE() > _total_pe_cut && flash.TotalPE() > min_pe) {
          ubxsec_event->candidate_flash_time = flash.Time();
          ubxsec_event->candidate_flash_z    = flash.ZCenter();
          min_pe                             = flash.TotalPE();
        }
      }
    } // OpDet loop

    double Ycenter, Zcenter, Ywidth, Zwidth;
    GetFlashLocation(ubxsec_event->beamfls_spec[n], Ycenter, Zcenter, Ywidth, Zwidth);
    ubxsec_event->beamfls_z[n] = Zcenter;

  } // flash loop

  // Collecting GENIE particles
  if(_use_genie_info) {
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (e.getByLabel("generator",mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    this->PrintMC(mclist);
    ubxsec_event->fv            = 0;
    ubxsec_event->fv_sce        = 0;

    ubxsec_event->ccnc          = -1;
    ubxsec_event->nupdg         = -1;
    ubxsec_event->nu_e          = -1;
    ubxsec_event->lep_costheta  = -9999.;
    ubxsec_event->true_muon_mom = -9999.;

    ubxsec_event->ResizeGenieTruthVectors(mclist.size());
    for (size_t iList = 0; iList < mclist.size(); iList++) {

      // Check if the true neutrino vertex is in the FV
      double truth_nu_vtx[3] = {mclist[iList]->GetNeutrino().Nu().Vx(),
                                mclist[iList]->GetNeutrino().Nu().Vy(),
                                mclist[iList]->GetNeutrino().Nu().Vz()};

      if (_fiducial_volume.InFV(truth_nu_vtx)) {
        ubxsec_event->fv = 1;
      }

      // Save the vertex for all neutrinos
      ubxsec_event->tvtx_x.at(iList) = mclist[iList]->GetNeutrino().Nu().Vx();
      ubxsec_event->tvtx_y.at(iList) = mclist[iList]->GetNeutrino().Nu().Vy();
      ubxsec_event->tvtx_z.at(iList) = mclist[iList]->GetNeutrino().Nu().Vz();
 
      // Look at the space charge correction
      geo::Vector_t sce_corr    = _SCE->GetPosOffsets(geo::Point_t(mclist[iList]->GetNeutrino().Nu().Vx(),
								mclist[iList]->GetNeutrino().Nu().Vy(),
								mclist[iList]->GetNeutrino().Nu().Vz()));

      double g4Ticks            = _detector_clocks->TPCG4Time2Tick(mclist[iList]->GetNeutrino().Nu().T())
                                + _detector_properties->GetXTicksOffset(0,0,0)
                                - _detector_properties->TriggerOffset();


      double xOffset_mcc8       = _detector_properties->ConvertTicksToX(g4Ticks, 0, 0, 0) - sce_corr.X(); //(this was for mcc8)
      double xOffset            = sce_corr.X();
      double xOffset_mcc9       = _detector_properties->ConvertTicksToX(g4Ticks, 0, 0, 0);
      double yOffset            = sce_corr.Y();
      double zOffset            = sce_corr.Z();

      ubxsec_event->sce_corr_x  = xOffset;
      ubxsec_event->sce_corr_y  = yOffset;
      ubxsec_event->sce_corr_z  = zOffset;
      ubxsec_event->time_mcc9_x = xOffset_mcc9;
      ubxsec_event->timeminussce_mcc8_x = xOffset_mcc8;

      if (_fiducial_volume.InFV(mclist[iList]->GetNeutrino().Nu().Vx() + xOffset,
                                mclist[iList]->GetNeutrino().Nu().Vy() + yOffset,
                                mclist[iList]->GetNeutrino().Nu().Vz() + zOffset)) {
        ubxsec_event->fv_sce = 1;
      }

      int n_genie_particles         = 0;
      int n_genie_particles_charged = 0;
      for (int p = 0; p < mclist[iList]->NParticles(); p++) {
        const simb::MCParticle mc_par = mclist[iList]->GetParticle(p);
        if (mc_par.StatusCode() != 1) continue;
        n_genie_particles ++;
        const TParticlePDG* par_pdg   = _database_pdg->GetParticle(mc_par.PdgCode());
        if (!par_pdg) continue;
        if (par_pdg->Charge() == 0) continue;
        n_genie_particles_charged ++;
      }

      // Only save this if the true neutrino vertex is in the FV
      if (ubxsec_event->fv == 1) {
        ubxsec_event->ccnc            = mclist[iList]->GetNeutrino().CCNC();
        ubxsec_event->mode            = mclist[iList]->GetNeutrino().Mode();
        ubxsec_event->nupdg           = mclist[iList]->GetNeutrino().Nu().PdgCode();
        ubxsec_event->nu_e            = mclist[iList]->GetNeutrino().Nu().E();
        ubxsec_event->lep_costheta    = mclist[iList]->GetNeutrino().Lepton().Pz() / mclist[iList]->GetNeutrino().Lepton().P();
        ubxsec_event->lep_phi         = UBXSecHelper::GetPhi(mclist[iList]->GetNeutrino().Lepton().Px(),
                                                           mclist[iList]->GetNeutrino().Lepton().Py(),
                                                           mclist[iList]->GetNeutrino().Lepton().Pz());
        ubxsec_event->genie_mult      = n_genie_particles;
        ubxsec_event->genie_mult_ch   = n_genie_particles_charged;
      }

      ubxsec_event->nsignal = 0;
      if(ubxsec_event->nupdg==14 && ubxsec_event->ccnc==0 && ubxsec_event->fv==1) ubxsec_event->nsignal=1;

      // Also save muon momentum if is CC interaction
      ubxsec_event->true_muon_mom = -9999.;
      if (mclist[iList]->GetNeutrino().CCNC() == 0) {
        for (int p = 0; p < mclist[iList]->NParticles(); p++) {
          auto const & mcp = mclist[iList]->GetParticle(p);
          if (mcp.Mother() != 0) continue;
          if (mcp.PdgCode() != 13) continue;
          ubxsec_event->true_muon_mom = mcp.P();
        }
      }

    }
  }

  ubxsec_event->is_signal = false;
  if (ubxsec_event->ccnc == 0 && ubxsec_event->nupdg == 14 && ubxsec_event->fv == 1) {
    ubxsec_event->is_signal = true;
  }

  lar_pandora::PFParticlesToSpacePoints pfp_to_spacept;
  lar_pandora::SpacePointsToHits spacept_to_hits;

  lar_pandora::PFParticleVector temp2;
  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, temp2, pfp_to_spacept);

  lar_pandora::SpacePointVector temp3;
  lar_pandora::LArPandoraHelper::CollectSpacePoints(e, _pfp_producer, temp3, spacept_to_hits);

  art::Handle<std::vector<recob::OpHit>> ophit_h;
  e.getByLabel("ophitBeam", ophit_h);

  art::Handle<std::vector<recob::OpHit>> ophit_cosmic_h;
  e.getByLabel("ophitCosmic", ophit_cosmic_h);

  // Check if the muon is reconstructed
  for (auto p : pfp_v) {
    auto mcghosts = mcghost_from_pfp.at(p.key());
    if (mcghosts.size() > 0) {
      art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghosts.at(0).key()).at(0);
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
      if (mc_truth) {
        if (mc_truth->Origin() == simb::kBeamNeutrino
            && mcpar->PdgCode() == 13 && mcpar->Mother() == 0) {
          ubxsec_event->muon_is_reco = true;
        }
      }
    }
  }

  std::vector<lar_pandora::TrackVector     > track_v_v;
  std::vector<lar_pandora::ShowerVector    > shower_v_v;
  std::vector<lar_pandora::PFParticleVector> pfp_v_v;
  for (size_t slice = 0; slice < tpcobj_h->size(); slice++) {
    track_v_v.push_back(tpcobjToTrackAssns.at(slice));
    shower_v_v.push_back(tpcobjToShowerAssns.at(slice));
    pfp_v_v.push_back(tpcobjToPFPAssns.at(slice));
  }

  ubxsec_event->nslices = tpcobj_h->size();
  ubxsec_event->ResizeVectors(tpcobj_h->size());

  ubxsec_event->n_tpcobj_nu_origin     = 0;
  ubxsec_event->n_tpcobj_cosmic_origin = 0;

  std::vector<art::Ptr<recob::Track>> muon_candidate_track_per_slice_v;
  std::vector<art::Ptr<recob::PFParticle>> muon_candidate_pfparticle_per_slice_v;
  std::vector<art::Ptr<recob::Vertex>> neutrino_candidate_vertex_per_slice_v;
  muon_candidate_track_per_slice_v.resize(ubxsec_event->nslices);
  muon_candidate_pfparticle_per_slice_v.resize(ubxsec_event->nslices);
  neutrino_candidate_vertex_per_slice_v.resize(ubxsec_event->nslices);

  //
  // THIS IS THE MAIN LOOP OVER THE
  // TPCOBJECTS CANDIDATES IN THIS EVENT
  //

  for (unsigned int slice = 0; slice < tpcobj_h->size(); slice++){

    ubana::TPCObject tpcobj = (*tpcobj_h)[slice];

    ubxsec_event->slc_npfp[slice]    = tpcobj.GetNPFP();
    ubxsec_event->slc_ntrack[slice]  = tpcobj.GetNTracks();
    ubxsec_event->slc_nshower[slice] = tpcobj.GetNShowers();

    // Slice origin
    ubxsec_event->slc_origin[slice]  = tpcobj.GetOrigin();

    if (tpcobj.GetOrigin() == ubana::kBeamNeutrino || tpcobj.GetOrigin() == ubana::kMixed)
      ubxsec_event->n_tpcobj_nu_origin ++;
    else
      ubxsec_event->n_tpcobj_cosmic_origin ++;

    // Slice origin extra
    ubxsec_event->slc_origin_extra[slice] = tpcobj.GetOriginExtra();

    // Containment
    ubxsec_event->slc_iscontained[slice] = UBXSecHelper::TracksAreContained(tpcobj.GetTracks());

    // Cosmic tagging: stopping mu
    ubxsec_event->slc_stopmu_tagged[slice]=false;
    std::vector<art::Ptr<anab::CosmicTag>> stopmu_tags = tpcobjToStopMu.at(slice);
    if ( stopmu_tags.size() == 1 ) {
      auto smt = stopmu_tags.at(0);
      if (smt->CosmicType() == anab::CosmicTagID_t::kGeometry_Y){
	ubxsec_event->slc_stopmu_tagged[slice] = true;
      }
    }

    // Reco vertex
    double reco_nu_vtx_raw[3];
    recob::Vertex tpcobj_nu_vtx = tpcobj.GetVertex();
    //tpcobj_nu_vtx.XYZ(reco_nu_vtx_raw);
    std::vector<art::Ptr<recob::Vertex>> recob_vtx_v = tpcobjToVertexAssns.at(slice);
    if (recob_vtx_v.size() > 0) {
      recob_vtx_v.at(0)->XYZ(reco_nu_vtx_raw);
      neutrino_candidate_vertex_per_slice_v.at(slice) = recob_vtx_v.at(0);
    } else {
      reco_nu_vtx_raw[0] = reco_nu_vtx_raw[1] = reco_nu_vtx_raw[2] = -9999;
    }

    // X position correction (time offset)
    double reco_nu_vtx[3];
    UBXSecHelper::GetTimeCorrectedPoint(reco_nu_vtx_raw, reco_nu_vtx, ubxsec_event->candidate_flash_time, _drift_velocity);

    // Set these equal to the position of the vertex as I define it.
    ubxsec_event->slc_nuvtx_fv[slice] = (_fiducial_volume.InFV(reco_nu_vtx) ? 1 : 0);

    // Vertex resolution
    if (ubxsec_event->slc_origin[slice] == ubana::kBeamNeutrino) {
      ubxsec_event->vtx_resolution = sqrt( pow(ubxsec_event->slc_nuvtx_y[slice]-ubxsec_event->tvtx_y[0], 2) + pow(ubxsec_event->slc_nuvtx_z[slice]-ubxsec_event->tvtx_z[0], 2) );
    }

    // Multiplicity
    int p, t, s;
    tpcobj.GetMultiplicity(p, t, s);
    ubxsec_event->slc_mult_pfp[slice]             = p;
    ubxsec_event->slc_mult_track[slice]           = t;
    ubxsec_event->slc_mult_shower[slice]          = s;
    ubxsec_event->slc_mult_track_tolerance[slice] = tpcobj.GetNTracksCloseToVertex(_tolerance_track_multiplicity);

    // Candidate Consistency
    ubxsec_event->slc_consistency[slice]                    = true;
    std::vector<art::Ptr<anab::CosmicTag>> consistency_tags = tpcobjToConsistency.at(slice);
    auto ct = consistency_tags.at(0);
    if (ct->CosmicType() != anab::CosmicTagID_t::kNotTagged) {
      ubxsec_event->slc_consistency[slice] = false;
    }
    ubxsec_event->slc_consistency_score[slice] = ct->CosmicScore();

    // Neutrino Flash match: is this slice selected by the external flash matching as the neutrino slice?
    // Get PFPs from TPCObject
    auto pfps_from_tpcobj = tpcobjToPFPAssns.at(slice);
    // Get primary PFP (as in: the one that has no parent). Check PDG value of primary PFP. If this is the neutrino slice, it will have a neutrino PDG code (12, 14, or 16). If not, it will have a muon PDG code
    bool isnuslc = false;
    for (auto pfp : pfps_from_tpcobj){
      if (pfp->IsPrimary()){
        if (TMath::Abs(pfp->PdgCode())==12 || TMath::Abs(pfp->PdgCode()==14) || TMath::Abs(pfp->PdgCode()==16)){
          isnuslc = true;
          break;
        }
      }
    }
    ubxsec_event->slc_is_nu[slice]   = isnuslc;
  
    // Hits
    int nhits_u, nhits_v, nhits_w;
    UBXSecHelper::GetNumberOfHitsPerPlane(e, _pfp_producer, track_v_v[slice], nhits_u, nhits_v, nhits_w);
    ubxsec_event->slc_nhits_u[slice] = nhits_u;
    ubxsec_event->slc_nhits_v[slice] = nhits_v;
    ubxsec_event->slc_nhits_w[slice] = nhits_w;

    // Longest track and check boundary
    recob::Track lt;
    if (UBXSecHelper::GetLongestTrackFromTPCObj(track_v_v[slice], lt)){
      ubxsec_event->slc_longesttrack_length[slice]      = lt.Length();
      ubxsec_event->slc_longesttrack_phi[slice]         = UBXSecHelper::GetCorrectedPhi(lt, tpcobj_nu_vtx);
      ubxsec_event->slc_longesttrack_theta[slice]       = UBXSecHelper::GetCorrectedCosTheta(lt, tpcobj_nu_vtx);
      ubxsec_event->slc_longesttrack_iscontained[slice] = UBXSecHelper::TrackIsContained(lt);
      int vtx_ok;
      ubxsec_event->slc_crosses_top_boundary[slice]     = (UBXSecHelper::IsCrossingTopBoundary(lt, vtx_ok) ? 1 : 0);
    } else {
      ubxsec_event->slc_longesttrack_length[slice]      = -9999;
    }

    // Longest shower
    recob::Shower ls;
    if (UBXSecHelper::GetLongestShowerFromTPCObj(shower_v_v[slice], ls)) {
      ubxsec_event->slc_longestshower_length[slice]    = ls.Length();
      ubxsec_event->slc_longestshower_openangle[slice] = ls.OpenAngle();
      ubxsec_event->slc_longestshower_startx[slice]    = ls.ShowerStart().X();
      ubxsec_event->slc_longestshower_starty[slice]    = ls.ShowerStart().Y();
      ubxsec_event->slc_longestshower_startz[slice]    = ls.ShowerStart().Z();
      ubxsec_event->slc_longestshower_phi[slice]       = UBXSecHelper::GetPhi(ls.Direction());
      ubxsec_event->slc_longestshower_theta[slice]     = UBXSecHelper::GetCosTheta(ls.Direction());
    }

    // Track quality
    ubxsec_event->slc_kalman_chi2[slice] = -9999;

  bool goodTrack = true;
    for (auto trk : track_v_v[slice]) {
      if (deadRegionsFinder.NearDeadReg2P( (trk->Vertex()).Y(), (trk->Vertex()).Z(), _minimumDistDeadReg )  ||
          deadRegionsFinder.NearDeadReg2P( (trk->End()).Y(),    (trk->End()).Z(),    _minimumDistDeadReg )  ||
          deadRegionsFinder.NearDeadRegCollection(trk->Vertex().Z(), _minimumDistDeadReg) ||
          deadRegionsFinder.NearDeadRegCollection(trk->End().Z(),    _minimumDistDeadReg) ||
          !UBXSecHelper::TrackPassesHitRequirment(e, _pfp_producer, trk, _minimumHitRequirement) ) {
        goodTrack = false;
        break;
      }
    }

    if (goodTrack) ubxsec_event->slc_passed_min_track_quality[slice] = true;
    else ubxsec_event->slc_passed_min_track_quality[slice] = false;

    // Vertex quality
    recob::Vertex slice_vtx = tpcobj.GetVertex();
    double slice_vtx_xyz[3];
    slice_vtx.XYZ(slice_vtx_xyz);
    ubxsec_event->slc_passed_min_vertex_quality[slice] = true;
    if (deadRegionsFinder.NearDeadReg2P(slice_vtx_xyz[1], slice_vtx_xyz[2], _minimumDistDeadReg))
      ubxsec_event->slc_passed_min_vertex_quality[slice] = false;

    // Channel status
    ubxsec_event->slc_nuvtx_closetodeadregion_u[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 0) ? 1 : 0);
    ubxsec_event->slc_nuvtx_closetodeadregion_v[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 1) ? 1 : 0);
    ubxsec_event->slc_nuvtx_closetodeadregion_w[slice] = (UBXSecHelper::PointIsCloseToDeadRegion(reco_nu_vtx, 2) ? 1 : 0);

    // Vertex check
    ubxsec::VertexCheck vtxCheck(track_v_v[slice], slice_vtx);
    ubxsec_event->slc_vtxcheck_angle[slice] = vtxCheck.AngleBetweenLongestTracks();

    // OpHits                                                                                                                                                                                           
    std::vector<ubxsec::Hit3D_t> hit3d_v;
    hit3d_v.clear();
    for (auto pfp : pfp_v_v[slice]) {
      auto iter = pfp_to_spacept.find(pfp);
      if ( iter == pfp_to_spacept.end() )
	continue;
      for (auto sp_pt : (iter->second)) {
        auto iter2 = spacept_to_hits.find(sp_pt);
        if (iter2 == spacept_to_hits.end()) {
          continue;
        }
        // Save sp_pt position and hit charge for all the sp_pt you have                                                                                                                              
	auto hit = iter2->second;
	ubxsec::Hit3D_t thishit;
        thishit.x = sp_pt->XYZ()[0];
        thishit.y = sp_pt->XYZ()[1];
        thishit.z = sp_pt->XYZ()[2];
        thishit.q = hit->Integral();
        hit3d_v.emplace_back(thishit);
      }
    }

    // Now construct average position
    double sumx = 0, sumy = 0, sumz = 0;
    double totq = 0;
    for (auto hit3d : hit3d_v) {
      sumx += hit3d.q * hit3d.x;
      sumy += hit3d.q * hit3d.y;
      sumz += hit3d.q * hit3d.z;

      totq += hit3d.q;
    }
    double charge_center[3] = {sumx / totq, sumy / totq, sumz / totq};

    int this_opch = UBXSecHelper::GetClosestPMT(charge_center);

    // Look at the opHits from this pmt
    int n_intime_ophits = 0;
    double n_intime_pe  = 0;
    for (size_t oh = 0; oh < ophit_h->size(); oh++) {
      auto const & ophit = (*ophit_h)[oh];
      size_t opdet       = geo->OpDetFromOpChannel(ophit.OpChannel());
      if(_make_ophit_csv) _csvfile2 << oh << "," << opdet << "," << ophit.PeakTime() << "," << _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude()) << std::endl;
      if (ophit.OpChannel() != this_opch) continue;
      if (ophit.PeakTime() > _beam_spill_start && ophit.PeakTime() < _beam_spill_end) {
        n_intime_ophits ++;
        n_intime_pe     += _pecalib.BeamPE(opdet,ophit.Area(),ophit.Amplitude());
      }
    } // end loop ophit

    ubxsec_event->slc_n_intime_pe_closestpmt[slice] = n_intime_pe;

    // Muon Candidate
    _muon_finder.Reset();
    _muon_finder.SetTracks(track_v_v[slice]);
    _muon_finder.SetTrackToPIDMap(track_to_pid_map);
    art::Ptr<recob::Track> candidate_track;
    // cout<<candidate
    bool muon_cand_exists = _muon_finder.GetCandidateTrack(candidate_track);
    if (muon_cand_exists) {

      bool fully_contained = _fiducial_volume.InFV(candidate_track->Vertex<TVector3>(), candidate_track->End<TVector3>());
      
      ubxsec_event->slc_muoncandidate_exists[slice]    = true;
      ubxsec_event->slc_muoncandidate_contained[slice] = fully_contained;
      ubxsec_event->slc_muoncandidate_length[slice]    = candidate_track->Length();
      ubxsec_event->slc_muoncandidate_phi[slice]       = UBXSecHelper::GetCorrectedPhi((*candidate_track), tpcobj_nu_vtx);
      ubxsec_event->slc_muoncandidate_theta[slice]     = UBXSecHelper::GetCorrectedCosTheta((*candidate_track), tpcobj_nu_vtx);
      ubxsec_event->slc_muoncandidate_theta_xz[slice]  = UBXSecHelper::GetCorrectedCosThetaXZ((*candidate_track), tpcobj_nu_vtx);
      ubxsec_event->slc_muoncandidate_theta_yz[slice]  = UBXSecHelper::GetCorrectedCosThetaYZ((*candidate_track), tpcobj_nu_vtx);

      ubxsec_event->slc_muoncandidate_mom_range[slice] = _trk_mom_calculator.GetTrackMomentum(candidate_track->Length(), 13);
      ubxsec_event->slc_muoncandidate_mom_mcs[slice]   = _trk_mom_calculator.GetMomentumMultiScatterLLHD(candidate_track);
     
      // For MCS first check the track direction is rigth

      TVector3 temp(reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]);
      bool track_direction_correct = (candidate_track->Vertex<TVector3>() - temp).Mag() < (candidate_track->End<TVector3>() - temp).Mag();

      if (track_direction_correct) {
	ubxsec_event->slc_muoncandidate_mom_mcs[slice] = mcsfitresult_mu_v.at(candidate_track.key())->fwdMomentum();
      } else {
	ubxsec_event->slc_muoncandidate_mom_mcs[slice] = mcsfitresult_mu_v.at(candidate_track.key())->bwdMomentum();
      }

      // Also see if the track is recon going downwards (for cosmic studies)
      bool track_going_down = candidate_track->Vertex().Y() > candidate_track->End().Y();
      
      // Look at calorimetry for the muon candidate
      std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_track.at(candidate_track.key());
      
      ubxsec_event->slc_muoncandidate_dqdx_vector[slice]      = UBXSecHelper::GetDqDxVector(calos);
      ubxsec_event->slc_muoncandidate_dqdx_trunc[slice]       = UBXSecHelper::GetDqDxTruncatedMean(calos);
      ubxsec_event->slc_muoncandidate_dqdx_u_trunc[slice]     = UBXSecHelper::GetDqDxTruncatedMean(calos, 0);
      ubxsec_event->slc_muoncandidate_dqdx_v_trunc[slice]     = UBXSecHelper::GetDqDxTruncatedMean(calos, 1);

      ubxsec_event->slc_muoncandidate_res_range_y[slice]      = UBXSecHelper::GetResRange(calos,2);
      ubxsec_event->slc_muoncandidate_res_range_u[slice]      = UBXSecHelper::GetResRange(calos, 0);
      ubxsec_event->slc_muoncandidate_res_range_v[slice]      = UBXSecHelper::GetResRange(calos, 1);

      ubxsec_event->slc_muoncandidate_dEdx_y[slice]           = UBXSecHelper::GetdEdx(calos,2);
      ubxsec_event->slc_muoncandidate_dEdx_u[slice]           = UBXSecHelper::GetdEdx(calos, 0);
      ubxsec_event->slc_muoncandidate_dEdx_v[slice]           = UBXSecHelper::GetdEdx(calos, 1);

      ubxsec_event->slc_muoncandidate_dQdx_y[slice]           = UBXSecHelper::GetdQdx(calos,2);
      ubxsec_event->slc_muoncandidate_dQdx_u[slice]           = UBXSecHelper::GetdQdx(calos, 0);
      ubxsec_event->slc_muoncandidate_dQdx_v[slice]           = UBXSecHelper::GetdQdx(calos, 1);

      ubxsec_event->slc_muoncandidate_mip_consistency[slice]  = _muon_finder.MIPConsistency(ubxsec_event->slc_muoncandidate_dqdx_trunc[slice],
                                                                                           ubxsec_event->slc_muoncandidate_length[slice]);
      ubxsec_event->slc_muoncandidate_mip_consistency2[slice] = _muon_finder.SVMPredict(ubxsec_event->slc_muoncandidate_dqdx_trunc[slice],
                                                                                        ubxsec_event->slc_muoncandidate_length[slice]);

      // Get the related PFP
      art::Ptr<recob::PFParticle> candidate_pfp = pfp_from_track.at(candidate_track.key()).at(0);
      const auto mcghosts = mcghost_from_pfp.at(candidate_pfp.key());
      if (mcghosts.size() > 0) {
        art::Ptr<simb::MCParticle> mcpar = mcpar_from_mcghost.at(mcghost_from_pfp.at(candidate_pfp.key()).at(0).key()).at(0);
        const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpar->TrackId());
        ubxsec_event->slc_muoncandidate_truepdg[slice] = mcpar->PdgCode();
        if (mc_truth) {

          // Check the true origin of the candidate PFP
          if (mc_truth->Origin() == simb::kBeamNeutrino) {
            ubxsec_event->slc_muoncandidate_trueorigin[slice] = ubana::kBeamNeutrino;
          } else if (mc_truth->Origin() == simb::kCosmicRay) {
            ubxsec_event->slc_muoncandidate_trueorigin[slice] = ubana::kCosmicRay;
          }

          // Now make momentum distributions
          if (mc_truth->Origin() == simb::kBeamNeutrino && mcpar->PdgCode() == 13 && mcpar->Mother() == 0) {
            _h_mom_true_mcs->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            if (fully_contained) {
              _mom_true_contained          = mcpar->P();
              _mom_mcs_contained           = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
              _mom_range_contained         = ubxsec_event->slc_muoncandidate_mom_range[slice];
              _mom_tree_contained->Fill();
              _h_mom_true_mcs_contained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
              _h_mom_true_range_contained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_range[slice]);
              _h_mom_range_mcs_contained->Fill(ubxsec_event->slc_muoncandidate_mom_range[slice],
                                               ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            } else {
              _mom_true_uncontained = mcpar->P();
              _mom_mcs_uncontained = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
              _mom_tree_uncontained->Fill();
              _h_mom_true_mcs_uncontained->Fill(mcpar->P(), ubxsec_event->slc_muoncandidate_mom_mcs[slice]);
            }
          }
          if (mc_truth->Origin() == simb::kCosmicRay && (mcpar->PdgCode() == 13 || mcpar->PdgCode() == -13)) {
            _mom_cosmic_true = mcpar->P();
            _mom_cosmic_mcs = ubxsec_event->slc_muoncandidate_mom_mcs[slice];
            _mom_cosmic_mcs_downforced = track_going_down ?   mcsfitresult_mu_v.at(candidate_track.key())->fwdMomentum()
                                                            : mcsfitresult_mu_v.at(candidate_track.key())->bwdMomentum();
            _mom_cosmic_range = ubxsec_event->slc_muoncandidate_mom_range[slice];
            _mom_cosmic_down = track_going_down;
            _mom_cosmic_tree->Fill();
          }
        }
      }

      // 
      // Look at residuals
      //
      std::vector<TVector3> hit_v; // a vec of hits from coll plane
      std::vector<TVector3> track_v; // a vec of hits from coll plane

      // Collect hits
      auto iter = trackToHitsMap.find(candidate_track);
      if (iter != trackToHitsMap.end()) {
        std::vector<art::Ptr<recob::Hit>> hits = iter->second;
        for (auto hit : hits) {
          if (hit->View() == 2) {
            TVector3 h (hit->WireID().Wire, hit->PeakTime(), 0);
            hit_v.emplace_back(h);
          }
        }
      }

      // Collect track points
      for (size_t i = 0; i < candidate_track->NumberTrajectoryPoints(); i++) {
        try {
          if (candidate_track->HasValidPoint(i)) {
            TVector3 trk_pt = candidate_track->LocationAtPoint<TVector3>(i);
            double wire     = geo->NearestWire(trk_pt, 2);
            double time     = _detector_properties->ConvertXToTicks(trk_pt.X(), geo::PlaneID(0,0,2));
            TVector3 p (wire, time, 0.);
            track_v.emplace_back(p);
          }
        } catch (...) {
          continue;
        }
      }

      //if (_debug) std::cout << "[UBXSec] \t \t Hit points: " << hit_v.size() << ", track points: " << track_v.size() << std::endl;
      ubana::TrackQuality _track_quality;
      _track_quality.SetTrackPoints(track_v);
      _track_quality.SetHitCollection(hit_v);
      std::pair<double, double> residual_mean_std           = _track_quality.GetResiduals();
      std::pair<double, double> residual_truncated_mean_std = _track_quality.GetTruncatedResiduals();
      std::pair<double,int> dist_wire_pair                  = _track_quality.GetTrackGap();

      int start_wire = dist_wire_pair.second;
      int end_wire  = dist_wire_pair.second + dist_wire_pair.first;

      if (start_wire > end_wire)
        std::swap(start_wire, end_wire);

      // Create a channel to wire map
      std::map<int, int> wire_to_channel;
      for (unsigned int ch = 0; ch < 8256; ch++) {
        std::vector< geo::WireID > wire_v = geo->ChannelToWire(ch);
        wire_to_channel[wire_v[0].Wire] = ch;
      }

      int n_dead_wires = 0;

      for (int wire = start_wire; wire < end_wire; wire++) {

        int channel = wire_to_channel[wire];

        // Channel statuses: 1=dead, 3=noisy, 4=good
        if (chanFilt.Status(channel) < 4) {
          n_dead_wires++;
        }
      }

      double r = _track_quality.GetR();

      double n_hits_in_cluster = 0;
      auto it = recoParticlesToHits.find(candidate_pfp);
      if (it != recoParticlesToHits.end()) {
        for (auto h : it->second) {
          if (h->View() == 2) {
            n_hits_in_cluster++;
          }
        }
      }

      double ratio = (double)hit_v.size()/n_hits_in_cluster;

      // Also look at the scattering angle
      std::vector<TVector3> dir_v;
      for (size_t p = 0; p < candidate_track->NumberTrajectoryPoints(); p++) {
        if (!candidate_track->HasValidPoint(p)) continue;
        dir_v.push_back(candidate_track->DirectionAtPoint<TVector3>(p));
      }
      std::vector<double> angle_v;
      for (size_t p = 0; p < dir_v.size()-1; p++) {
        double angle = dir_v.at(p).Angle(dir_v.at(p+1));
        angle_v.push_back(angle);
      }
      std::sort(angle_v.begin(), angle_v.end());

      double max_angle = -1;
      if (angle_v.size() != 0)
        max_angle = angle_v.at(angle_v.size()-1) / TMath::Pi() * 180.;

      ubxsec_event->slc_muoncandidate_residuals_mean[slice]                 = residual_mean_std.first;
      ubxsec_event->slc_muoncandidate_residuals_std[slice]                  = residual_mean_std.second;
      ubxsec_event->slc_muoncandidate_residuals_truncatedmean[slice]        = residual_truncated_mean_std.first;
      ubxsec_event->slc_muoncandidate_residuals_truncatedstd[slice]         = residual_truncated_mean_std.second;
      ubxsec_event->slc_muoncandidate_wiregap[slice]                        = end_wire-start_wire;
      ubxsec_event->slc_muoncandidate_wiregap_dead[slice]                   = n_dead_wires;
      ubxsec_event->slc_muoncandidate_linearity[slice]                      = r;
      ubxsec_event->slc_muoncandidate_perc_used_hits_in_cluster[slice]      = ratio;
      ubxsec_event->slc_muoncandidate_maxscatteringangle[slice]             = max_angle;

      muon_candidate_track_per_slice_v.at(slice)                            = candidate_track;
      muon_candidate_pfparticle_per_slice_v.at(slice)                       = candidate_pfp;


    } else {

      ubxsec_event->slc_muoncandidate_exists[slice]    = false;
      ubxsec_event->slc_muoncandidate_length[slice]    = -9999;
      ubxsec_event->slc_muoncandidate_phi[slice]       = -9999;
      ubxsec_event->slc_muoncandidate_theta[slice]     = -9999;
      ubxsec_event->slc_muoncandidate_theta_xz[slice]  = -9999;
      ubxsec_event->slc_muoncandidate_theta_yz[slice]  = -9999;
      ubxsec_event->slc_muoncandidate_mom_range[slice] = -9999;
      ubxsec_event->slc_muoncandidate_mom_mcs[slice]   = -9999;

    }

    // Particle ID
    for (auto pfp : pfps_from_tpcobj){

      std::vector<art::Ptr<ubana::MCGhost>> mcghosts = mcghost_from_pfp.at(pfp.key());
      std::vector<art::Ptr<simb::MCParticle>> mcpars;
      int pdg = -1;
      if (mcghosts.size() == 0 || mcghosts.size() > 1 ) {
        continue;
      }

      mcpars = mcpar_from_mcghost.at(mcghosts[0].key());
      pdg = mcpars[0]->PdgCode();
      const auto mc_truth = UBXSecHelper::TrackIDToMCTruth(e, "largeant", mcpars[0]->TrackId());
      if (!mc_truth) {
        std::cerr << "[UBXSec] Problem with MCTruth pointer." << std::endl;
        continue;
      }
      if (mc_truth->Origin() == simb::kBeamNeutrino &&
          mcpars[0]->PdgCode() == 13 && mcpars[0]->Mother() == 0) {
        ubxsec_event->true_muon_mom_matched = mcpars[0]->P();

        if (_fiducial_volume.InFV(mcpars[0]->Vx(), mcpars[0]->Vy(), mcpars[0]->Vz()) &&
          _fiducial_volume.InFV(mcpars[0]->EndX(), mcpars[0]->EndY(), mcpars[0]->EndZ())) {
          ubxsec_event->mc_muon_contained = true;
        }
      }

      std::vector<art::Ptr<recob::Track>> tracks = tracks_from_pfp.at(pfp.key());
      for (auto track : tracks) {

        std::vector<art::Ptr<anab::ParticleID>> pids = particleids_from_track.at(track.key());
        if(pids.size() == 0){
	  continue;
	}

	std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[0]->ParticleIDAlgScores();

	// Loop though AlgScoresVec and find the variables we want
	double tmppida=-9999;
	for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
	  anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
	  int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
	  if (planenum<0 || planenum>2) continue;
	  if (AlgScore.fVariableType==anab::kPIDA){
	    if ( planenum == 2 ) {
	      tmppida=AlgScore.fValue;
	    }
	  }
	}

	// Don't fill things if PIDA was not found
	if (tmppida==-9999) continue;

	if (pdg == 13) {
	  _h_pida_muon->Fill(tmppida);
	  _h_pida_len_muon->Fill(tmppida, track->Length());
	  if( tmppida > 0 && tmppida < 50. && _make_pida_csv) _csvfile << tmppida << "," << track->Length() << "," << "1" << std::endl;
	} else if (pdg == 2212) {
	  _h_pida_proton->Fill(tmppida);
	  _h_pida_len_proton->Fill(tmppida, track->Length());
	  if( tmppida > 0 && tmppida < 50. && _make_pida_csv) _csvfile << tmppida << "," << track->Length() << "," << "0" << std::endl;
	} else if (pdg == 211) {
	  _h_pida_pion->Fill(tmppida);
	  _h_pida_len_pion->Fill(tmppida, track->Length());
	} else if (pdg == 321) {
	  _h_pida_kaon->Fill(tmppida);
	  _h_pida_len_kaon->Fill(tmppida, track->Length());
	}
      }
    }

    // MCS - Track direction, study cosmic direction
    // Take the muon candidate track for a TPCObject
    // with cosmic origin, and check the track direction
    // (up/down) from MCS result
    if (muon_cand_exists && tpcobj.GetOrigin() == ubana::kCosmicRay) {
      bool best_fwd = mcsfitresult_mu_v.at(candidate_track.key())->isBestFwd();
      bool down_track = candidate_track->Vertex().Y() > candidate_track->End().Y();
      if (down_track && best_fwd) {               // Track is reco going down and mcs agrees (true for cosmic)
        _h_mcs_cosmic_track_direction->Fill(0);
        _mcs_cosmic_track_direction = 0;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
      } else if (!down_track && best_fwd) {       // Track is reco going up and mcs agrees
        _h_mcs_cosmic_track_direction->Fill(1);
        _mcs_cosmic_track_direction = 1;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
      } else if (down_track && !best_fwd) {       // Track is reco going down and mcs disagrees
        _h_mcs_cosmic_track_direction->Fill(1);
        _mcs_cosmic_track_direction = 2;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
      } else if (!down_track && !best_fwd) {      // Track is reco going up and mcs disagrees (true for cosmic)
        _h_mcs_cosmic_track_direction->Fill(0);
        _mcs_cosmic_track_direction = 3;
        _mcs_cosmic_track_downll    = mcsfitresult_mu_v.at(candidate_track.key())->bwdLogLikelihood();
        _mcs_cosmic_track_upll      = mcsfitresult_mu_v.at(candidate_track.key())->fwdLogLikelihood();
      }
      _mcs_cosmic_track_direction_tree->Fill();
    }

    std::cout << "[UBXSec] --- SLICE INFORMATION SAVED" << std::endl;
  } // slice loop

  if (_is_mc) {

    // SW Trigger
    art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
    e.getByLabel("swtrigger", softwareTriggerHandle);

    if (softwareTriggerHandle.isValid()){
      if (softwareTriggerHandle->getNumberOfAlgorithms() == 1) {
        std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
        ubxsec_event->is_swtriggered = (softwareTriggerHandle->passedAlgo(algoNames[0]) ? 1 : 0);
      }
    }

    // MC Flash
    ::art::Handle<std::vector<recob::OpFlash> > nuMcflash_h;
    e.getByLabel("NeutrinoMCFlash",nuMcflash_h);
    if( !nuMcflash_h.isValid() || nuMcflash_h->empty() ) {
      std::cerr << "Don't have neutrino MC flashes." << std::endl;
    } else {

      auto const& flash = (*nuMcflash_h)[0];
      ubxsec_event->numc_flash_spec.resize(geo->NOpDets());
      for (unsigned int i = 0; i < geo->NOpDets(); i++) {
        unsigned int opdet = geo->OpDetFromOpChannel(i);
        ubxsec_event->numc_flash_spec[opdet] = flash.PE(i);
      }
    }

    // MCFlash vs op activity
    bool opActivityInBeamSpill = false;
    // Check if there are recon beam flashed in the beam spill window
    for (auto reco_fls_time : ubxsec_event->beamfls_time) {
      if (reco_fls_time > _beam_spill_start && reco_fls_time < _beam_spill_end) {
         opActivityInBeamSpill = true;
       }
    }
    if (nuMcflash_h->size() == 0) {
      if(opActivityInBeamSpill) {
        //std::cout << "No MCFlash but optical activity in the beam spill." << std::endl;
        ubxsec_event->no_mcflash_but_op_activity = true;
      }
    }

  }


  // POT
  art::Handle< sumdata::POTSummary > potsum_h;
  if(e.getByLabel(_potsum_producer, potsum_h))
    ubxsec_event->pot = potsum_h->totpot;
  else
    ubxsec_event->pot = 0.;

  // *********************
  // Event Selection
  // *********************

  _event_selection.SetEvent(ubxsec_event);

  size_t slice_index=-999;
  int slice_index_1=slice_index;
  std::string reason = "no_failure";
  std::map<std::string,bool> failure_map;

  failure_map.clear();

  failure_map.emplace(std::make_pair("a_crt_veto_cut", false));
  failure_map.emplace(std::make_pair("b_beam_disc_flashes", false));
  failure_map.emplace(std::make_pair("c_beam_spill_flash", false));
  failure_map.emplace(std::make_pair("d_has_slices", false));
  failure_map.emplace(std::make_pair("e_has_slice_tagged_as_neutrino", false));
  failure_map.emplace(std::make_pair("f_ntrack", false));
  failure_map.emplace(std::make_pair("g_track_length", false));
  failure_map.emplace(std::make_pair("h_residuals_std_up", false));
  failure_map.emplace(std::make_pair("i_perc_used_hits_in_cluster", false));

  bool is_selected = _event_selection.IsSelected(slice_index_1, failure_map);
  if (_debug) std::cout << "[UBXSec] >>>>>>>>>>>>>>>>>>>>>> Is Selected? " << (is_selected ? "YES" : "NO") << std::endl;

  for (auto iter : failure_map) {

    std::cout << "[UBXSec] Cut: " << iter.first << "  >>>  " << (iter.second ? "PASSED" : "NOT PASSED") << std::endl;
    if ( !iter.second ) {
      reason = "fail_" + iter.first;
      break;
    }

  }

  std::cout << "[UBXSec] Selection Failed at Cut: " << reason << std::endl;

  ::ubana::SelectionResult selection_result;
  selection_result.SetSelectionType("numu_cc_inclusive");
  selection_result.SetFailureReason(reason);
  selection_result.SetCutFlowStatus(failure_map);

  // See if one of the last three cases fails.
  failure_map.emplace(std::make_pair("j_containment_cut", true));
  failure_map.emplace(std::make_pair("k_fiducial_volume_cut", true));
  failure_map.emplace(std::make_pair("l_length_of_tracks_from_vertex_cut", true));
  failure_map.emplace(std::make_pair("m_new_two_track_length_requirement", true));
  failure_map.emplace(std::make_pair("n_number_of_tracks_in_TPCObject_requirement", true));
  failure_map.emplace(std::make_pair("o_number_of_tracks_originating_from_vertex_requirement", true));
  failure_map.emplace(std::make_pair("p_CRT_distance_cut", true));

  // *********************
  // Save Event Selection Output in the Event
  // *********************
  if ( !is_selected ) {

    total_num_events_failing++;

    total_num_events_failing_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;

    selection_result.SetSelectionStatus(false);

    ubxsec_event->is_selected = false;

    selectionResultVector->emplace_back(std::move(selection_result));
    
    if ( reason == "fail_b_beam_disc_flashes" ) {
      total_num_events_failed_beam_disc_flashes++;
      total_num_events_failed_beam_disc_flashes_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
      std::cout << "No beam discriminated flashes." << std::endl;
    }

    if ( reason == "fail_c_beam_spill_flash"  ) {
      total_num_events_failed_beam_spill_flash++;
      total_num_events_failed_beam_spill_flash_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
      std::cout << "No beam discriminated flash > 50 PEs in the beamspill window." << std::endl;
    }

    if ( reason == "fail_d_has_slices"  ) {
      total_num_events_failed_has_slices++;
      total_num_events_failed_has_slices_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
      std::cout << "There is no slice in this event." << std::endl;
    }

    if ( reason == "fail_f_ntrack"  ) {
      total_num_events_failed_ntrack++;
      total_num_events_failed_ntrack_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
      std::cout << "There are less than the required number of tracks in the TPC Object." << std::endl;
    }

    if ( reason == "fail_g_track_length"  ) {
      total_num_events_failed_track_length++;
      std::cout << "The longest muon candidate track in the event fails the track length cut." << std::endl;
    }

    if ( reason == "fail_h_residuals_std_up"  ) {
      total_num_events_failed_residuals_std_up++;
      total_num_events_failed_residuals_std_up_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
      std::cout << "This event has failed residuals std up cut." << std::endl;
    }

    if ( reason == "fail_i_perc_used_hits_in_cluster"  ) {
      total_num_events_failed_perc_used_hits_in_cluster++;
      total_num_events_failed_perc_used_hits_in_cluster_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
      std::cout << "This event has failed the percentage of hits used in the cluster." << std::endl;
    }

  } else if ( is_selected ) {

    // Save all of the information for the passing event.
    nslices_in_event = ubxsec_event->nslices;
    
    int nu_slc_idx = -1;

    // Find the slice that was tagged as a neutrino.
    for ( int slice_iter = 0; slice_iter < nslices_in_event; slice_iter++ ) {

      if ( ubxsec_event->slc_is_nu.at( slice_iter ) == true ) {
	nu_slc_idx = slice_iter;
	break;
      }

    }
 
    ubana::TPCObject tpcobj = (*tpcobj_h)[nu_slc_idx];

    auto tracks = tpcobjToTrackAssns.at( nu_slc_idx );

    number_of_tracks_in_TPCObject  = tracks.size();

    // Reco vertex                                                                                                                                                                                         
    double reconstructed_nu_vtx_raw[3];
    std::vector<art::Ptr<recob::Vertex>> reconstructed_vertex_v = tpcobjToVertexAssns.at(nu_slc_idx);
    if ( reconstructed_vertex_v.size() > 0 ) {
      reconstructed_vertex_v.at(0)->XYZ(reconstructed_nu_vtx_raw);
    } else {
      reconstructed_nu_vtx_raw[0] = reconstructed_nu_vtx_raw[1] = reconstructed_nu_vtx_raw[2] = -9999;
    }

    // X position correction (time offset)                                                                                                                                                                 
    double reconstructed_nu_vtx[3];

    reconstructed_nu_vtx[0] = ( reconstructed_nu_vtx_raw[0] - ( 0.1098 * ubxsec_event->candidate_flash_time ) );
    reconstructed_nu_vtx[1] = reconstructed_nu_vtx_raw[1];
    reconstructed_nu_vtx[2] = reconstructed_nu_vtx_raw[2];

    // Offset this point for space charge using the same technique as for the vertices below.
    geo::Vector_t tpcobj_vertex_offsets = SCE_data->GetCalPosOffsets(geo::Point_t{reconstructed_nu_vtx[0],reconstructed_nu_vtx[1],reconstructed_nu_vtx[2]});
    
    double reconstructed_nu_vtx_with_SCE_correction[3];

    reconstructed_nu_vtx_with_SCE_correction[0]  = reconstructed_nu_vtx[0] - tpcobj_vertex_offsets.X();
    reconstructed_nu_vtx_with_SCE_correction[1]  = reconstructed_nu_vtx[1] + tpcobj_vertex_offsets.Y();
    reconstructed_nu_vtx_with_SCE_correction[2]  = reconstructed_nu_vtx[2] + tpcobj_vertex_offsets.Z();

    pandora_vertex_x                             = reconstructed_nu_vtx_with_SCE_correction[0];
    pandora_vertex_y                             = reconstructed_nu_vtx_with_SCE_correction[1];
    pandora_vertex_z                             = reconstructed_nu_vtx_with_SCE_correction[2];

    std::cout << "Reconstructed Nu Vtx Before Correcting For SCE: x = " << reconstructed_nu_vtx[0] << " cm y = " << reconstructed_nu_vtx[1] << " cm z = " << reconstructed_nu_vtx[2] << " cm." << std::endl;
    std::cout << "Reconstructed Nu Vtx After Correcting For SCE: x = " << reconstructed_nu_vtx_with_SCE_correction[0] << " cm y = " << reconstructed_nu_vtx_with_SCE_correction[1] << " cm z = " << reconstructed_nu_vtx_with_SCE_correction[2] << " cm." << std::endl;

    distance_between_truth_vertex_and_pandora_vertex = TMath::Sqrt( ( nu_vtx_x_truth - reconstructed_nu_vtx_with_SCE_correction[0] ) * ( nu_vtx_x_truth - reconstructed_nu_vtx_with_SCE_correction[0] ) + ( nu_vtx_y_truth - reconstructed_nu_vtx_with_SCE_correction[1] ) * ( nu_vtx_y_truth - reconstructed_nu_vtx_with_SCE_correction[1] ) + ( nu_vtx_z_truth - reconstructed_nu_vtx_with_SCE_correction[2] ) * ( nu_vtx_z_truth - reconstructed_nu_vtx_with_SCE_correction[2] ) );

    std::vector< double > muon_candidate_track_IDs;
    std::vector< double > muon_candidate_track_lengths;

    muon_candidate_track_IDs.clear();
    muon_candidate_track_lengths.clear();

    // Declare vectors for the three variables going into the BDT.
    std::vector< double > Muon_Candidate_Muon_Track_Score_vector;
    std::vector< double > Muon_Candidate_Proton_PID_vector;
    std::vector< double > Muon_Candidate_NuScore_vector;
    std::vector< double > Muon_Candidate_Truth_Length;
    std::vector< int    > Muon_Candidate_PDG_Code;
    std::vector< double > muon_candidate_truth_start_x_vector;
    std::vector< double > muon_candidate_truth_start_y_vector;
    std::vector< double > muon_candidate_truth_start_z_vector;
    std::vector< double > muon_candidate_truth_end_x_vector;
    std::vector< double > muon_candidate_truth_end_y_vector;
    std::vector< double > muon_candidate_truth_end_z_vector;

    Muon_Candidate_Muon_Track_Score_vector.clear();
    Muon_Candidate_Proton_PID_vector.clear();
    Muon_Candidate_NuScore_vector.clear();
    Muon_Candidate_Truth_Length.clear();
    Muon_Candidate_PDG_Code.clear();
    muon_candidate_truth_start_x_vector.clear();
    muon_candidate_truth_start_y_vector.clear();
    muon_candidate_truth_start_z_vector.clear();
    muon_candidate_truth_end_x_vector.clear();
    muon_candidate_truth_end_y_vector.clear();
    muon_candidate_truth_end_z_vector.clear();

    // Declare variables for these variables for the longest track in the TPCObject.
    double Longest_Track_Muon_Track_Score_ = 0.;
    double Longest_Track_Proton_PID_       = 0.;
    double Longest_Track_NuScore_          = 0.;

    std::vector< double > Track_MCParticle_Starting_X;
    std::vector< double > Track_MCParticle_Starting_Y;
    std::vector< double > Track_MCParticle_Starting_Z;
    std::vector< double > Track_MCParticle_Ending_X;
    std::vector< double > Track_MCParticle_Ending_Y;
    std::vector< double > Track_MCParticle_Ending_Z;

    Track_MCParticle_Starting_X.clear();
    Track_MCParticle_Starting_Y.clear();
    Track_MCParticle_Starting_Z.clear();
    Track_MCParticle_Ending_X.clear();
    Track_MCParticle_Ending_Y.clear();
    Track_MCParticle_Ending_Z.clear();

    double longest_muon_candidate_length   = -1.;

    // Loop through the tracks to find out which MCParticle each of them corresponds to.                                                                                                                   
    for ( size_t tracks_iter = 0; tracks_iter < tracks.size(); tracks_iter++ ) {

      // Find the PID values corresponding to the track.                                                                                                                                                  
      for (auto const pfp : pfdaughters){

	if ( !pfp->IsPrimary() ){

	  if ( particlesToTracks.find(pfp) != particlesToTracks.end() ){

	    this_chiproton_    = -1.;
	    this_chipion_      = -1.;
	    this_chimuon_      = -1.;
	    muon_proton_ratio  = -1.;
	    muon_pion_ratio    = -1.;
	    Muon_Track_Score_  = -1.;
	    Proton_PID_        = -1.;
	    NuScore_           = -1.;

	    // Unpack the corresponding track.                                                                                                                                                            
	    const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track                                                                                                  

	    if ( tracks.at( tracks_iter )->ID() == this_track->ID() ) {

	      std::map<std::string, float> pid_map;

	      if( trackHelper.getPID(pid_map, this_track, trackPIDAssn ) ){
		this_chiproton_ = pid_map.at("chi2_proton");
		this_chipion_   = pid_map.at("chi2_pion");
		this_chimuon_   = pid_map.at("chi2_muon");

		muon_proton_ratio = ( this_chimuon_ / this_chiproton_ );
		muon_pion_ratio   = ( this_chimuon_ / this_chipion_ );

		// Push back the track into the muon candidate vector if it passes Thomas's cuts.                                                                                                         
		if ( muon_proton_ratio < 0.168 && muon_pion_ratio < 1.06 ) {
		  
		  muon_candidate_track_lengths.push_back( tracks.at( tracks_iter )->Length() );
		  muon_candidate_track_IDs.push_back( tracks.at( tracks_iter )->ID() );

		}

		// Track Score.                                                                                                                                                                           
		const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
		TrackScore_.push_back(pfp_properties.at("TrackScore"));

		Muon_Track_Score_ = pfp_properties.at("TrackScore");

		if ( muon_proton_ratio < 0.168 && muon_pion_ratio < 1.06 ) {

		  Muon_Candidate_Muon_Track_Score_vector.push_back( Muon_Track_Score_ );

		}

		// Neutrino Score.                                                                                                                                                                        
		lar_pandora::MetadataVector neutrino_metadata_vec;

		if ( pfneutrinos.size() > 0 ) {

		  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
		  neutrino_metadata_vec = particlesToMetadata.at(pfnu);

		}

		// Set 'NuScore_' to a default value.                                                                                                                                                     
		NuScore_ = -1.;

		if ( neutrino_metadata_vec.size() == 1 ) {

		  const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
		  NuScore_ = neutrino_properties.at("NuScore"); // nuscore $$                                                                                                                             

		  if ( muon_proton_ratio < 0.168 && muon_pion_ratio < 1.06 ) {

		    Muon_Candidate_NuScore_vector.push_back( NuScore_ );

		  }

		}
		else {

		  if ( muon_proton_ratio < 0.168 && muon_pion_ratio < 1.06 ) {

		    Muon_Candidate_NuScore_vector.push_back( -1. );

		  }

		}

		auto Trk_start_offset = SCE_data->GetCalPosOffsets(geo::Point_t(this_track->Start().X(), this_track->Start().Y(), this_track->Start().Z()));
		TrackStart_x_sce_.push_back(this_track->Start().X() - Trk_start_offset.X() );
		TrackStart_y_sce_.push_back(this_track->Start().Y() + Trk_start_offset.Y() );
		TrackStart_z_sce_.push_back(this_track->Start().Z() + Trk_start_offset.Z() );
		auto Trk_end_offset = SCE_data->GetCalPosOffsets(geo::Point_t(this_track->End().X(), this_track->End().Y(), this_track->End().Z()));
		TrackEnd_x_sce_.push_back(this_track->End().X() - Trk_end_offset.X() );
		TrackEnd_y_sce_.push_back(this_track->End().Y() + Trk_end_offset.Y() );
		TrackEnd_z_sce_.push_back(this_track->End().Z() + Trk_end_offset.Z() );

		TVector3 Trk_start_SCEcorr;
		TVector3 Trk_end_SCEcorr;
		Trk_start_SCEcorr.SetX(TrackStart_x_sce_.back());
		Trk_start_SCEcorr.SetY(TrackStart_y_sce_.back());
		Trk_start_SCEcorr.SetZ(TrackStart_z_sce_.back());
		Trk_end_SCEcorr.SetX(TrackEnd_x_sce_.back());
		Trk_end_SCEcorr.SetY(TrackEnd_y_sce_.back());
		Trk_end_SCEcorr.SetZ(TrackEnd_z_sce_.back());

		// Proton PID.                                                                                                                                                                            
		PID pid;
		pid.Chi2(trackPIDAssn,this_track, Trk_start_SCEcorr, Trk_end_SCEcorr);
		TrackPID_chiproton_.push_back(pid.PID_Chi2P_3pl);  // chi squqre cuts                                                                                                                     

		Proton_PID_ = pid.PID_Chi2P_3pl;

		if ( muon_proton_ratio < 0.168 && muon_pion_ratio < 1.06 ) {

		  Muon_Candidate_Proton_PID_vector.push_back( Proton_PID_ );

		}

		// Fill the variable for the longest track if this is (currently) the longest track in the TPCObject.                                                                                     
		if ( this_track->Length() > longest_muon_candidate_length ) {

		  Longest_Track_Muon_Track_Score_          = Muon_Track_Score_;
		  Longest_Track_NuScore_                   = NuScore_;
		  Longest_Track_Proton_PID_                = Proton_PID_;

		  longest_muon_candidate_length = this_track->Length();

		}

	      }

	    }

	  }

	}

      }

    }		  

    // 'nu_slc_idx' necessarily has to be 0 or greater.
    std::vector<lar_pandora::TrackVector     > track_v_v_for_event_tree;
    track_v_v_for_event_tree.push_back(tpcobjToTrackAssns.at(nu_slc_idx));

    // Find the longest_track.
    recob::Track lt_track;
    UBXSecHelper::GetLongestTrackFromTPCObj(track_v_v_for_event_tree[nu_slc_idx], lt_track );

    longest_TPCObject_track_length = lt_track.Length();

    // Find the track key corresponding to the true muon candidate.                                                                                                                                     
    double longest_track_length      = -1;
    int    muon_candidate_ID         = -1;

    Muon_Candidate_Muon_Track_Score_ = 0.;
    Muon_Candidate_Proton_PID_       = 0.;
    Muon_Candidate_NuScore_          = 0.;

    for ( size_t muon_candidate_iter = 0; muon_candidate_iter < muon_candidate_track_lengths.size(); muon_candidate_iter++ ) {

      if ( muon_candidate_track_lengths.at( muon_candidate_iter ) > longest_track_length ) {

        longest_track_length             = muon_candidate_track_lengths.at( muon_candidate_iter);
        muon_candidate_ID                = muon_candidate_track_IDs.at( muon_candidate_iter );

        Muon_Candidate_Muon_Track_Score_ = Muon_Candidate_Muon_Track_Score_vector.at( muon_candidate_iter );
        Muon_Candidate_Proton_PID_       = Muon_Candidate_Proton_PID_vector.at( muon_candidate_iter );
        Muon_Candidate_NuScore_          = Muon_Candidate_NuScore_vector.at( muon_candidate_iter );

      }

    }
    
    size_t muon_candidate_idx       = -1;
    art::Ptr<recob::Track> correct_muon_candidate_track;
    bool found_muon_candidate_track = false;

    for ( size_t pandora_track_iter = 0; pandora_track_iter < pandora_track_h->size(); pandora_track_iter++ ) {

      if ( pandora_track_h->at( pandora_track_iter ).ID() == muon_candidate_ID )
	muon_candidate_idx = pandora_track_iter;

    }

    for ( size_t pandora_track_pointers_iter = 0; pandora_track_pointers_iter < pandora_track_pointers_h.size(); pandora_track_pointers_iter++ ) {

      // Perform the same conditional for the pointer in my code.
      if ( pandora_track_pointers_h.at( pandora_track_pointers_iter )->ID() == muon_candidate_ID ) {
	correct_muon_candidate_track = pandora_track_pointers_h.at( pandora_track_pointers_iter );
	found_muon_candidate_track = true;
      }

    } // End of the loop over the pandora tracks to find the one that corresponds to the muon candidate.


    // Use the original method if 'muon_candidate_idx' is still equal to -1 (just take the longest track in the TPCObject as the muon candidate).                                                          
    if ( int(muon_candidate_idx) == -1 ) {

      // Set the 'Muon_Candidate_' variable equal to those for the longest track in the TPCObject.                                                                                                        
      Muon_Candidate_Muon_Track_Score_ = Longest_Track_Muon_Track_Score_;
      Muon_Candidate_Proton_PID_       = Longest_Track_Proton_PID_;
      Muon_Candidate_NuScore_          = Longest_Track_NuScore_;

      for ( size_t pandora_track_iter = 0; pandora_track_iter < pandora_track_h->size(); pandora_track_iter++ ) {

        if ( fabs( pandora_track_h->at( pandora_track_iter ).Length() - ubxsec_event->slc_longesttrack_length[nu_slc_idx] ) < 0.001 ) {
          muon_candidate_idx = pandora_track_iter;
          muon_candidate_ID  = pandora_track_h->at( pandora_track_iter ).ID();
        }

      } // End of the loop over the pandora tracks to find the one that corresponds to the muon candidate.                                                                                               

    }

    if ( found_muon_candidate_track == false ) {

      for ( size_t pandora_track_pointers_iter = 0; pandora_track_pointers_iter < pandora_track_pointers_h.size(); pandora_track_pointers_iter++ ) {

	if ( fabs( pandora_track_pointers_h.at( pandora_track_pointers_iter )->Length() - ubxsec_event->slc_longesttrack_length[nu_slc_idx] ) < 0.001 )
	  correct_muon_candidate_track = pandora_track_pointers_h.at( pandora_track_pointers_iter );

      } // End of the loop over the pandora tracks to find the one that corresponds to the muon candidate.                                                                                                

    }

    // Start Neutrino Energy Calculation Code.                                                                                                                                                          
    bool first_points_of_track_are_start = true;

    // Set a variable for the calorimety object of the muon candidate on the track.                                                                                                                      
    auto Calo_v                                        = trk_calo_assn_v.at( muon_candidate_idx );
    auto Calo_v_no_SCE_corrections                     = trk_calo_assn_v_no_SCE_corrections.at( muon_candidate_idx );

    auto vhit                                          = fmthm.at( muon_candidate_idx );
    auto vmeta                                         = fmthm.data( muon_candidate_idx );

    int  count                                         = 0;

    std::vector< float > xyz_v_x_coords;
    xyz_v_x_coords.clear();

    std::vector< float > xyz_v_y_coords;
    xyz_v_y_coords.clear();

    std::vector< float > xyz_v_z_coords;
    xyz_v_z_coords.clear();

    std::vector< float > xyz_v_no_SCE_corrections_x_coords;
    xyz_v_no_SCE_corrections_x_coords.clear();

    std::vector< float > xyz_v_no_SCE_corrections_y_coords;
    xyz_v_no_SCE_corrections_y_coords.clear();

    std::vector< float > xyz_v_no_SCE_corrections_z_coords;
    xyz_v_no_SCE_corrections_z_coords.clear();

    std::vector < float > dqdx_v;
    dqdx_v.clear();

    std::vector < float > dedx_v;
    dedx_v.clear();
    
    std::vector < float > dqdx_points_ordered;
    dqdx_points_ordered.clear();

    std::vector < float > dedx_points_ordered;
    dedx_points_ordered.clear();

    std::vector< float > dQdx_outliers_and_endpoints_removed_one_standard_deviation;
    dQdx_outliers_and_endpoints_removed_one_standard_deviation.clear();

    std::vector< float > dEdx_outliers_and_endpoints_removed;
    dEdx_outliers_and_endpoints_removed.clear();
    
    bool fill_vector                                   = true;

    for ( size_t pl = 0; pl < 3; pl++ ) {

      // Do not loop over planes with a meaningless index.
      if ( Calo_v[pl]->PlaneID().Plane != 0 && Calo_v[pl]->PlaneID().Plane != 1 && Calo_v[pl]->PlaneID().Plane != 2 ) 
	continue;

      dqdx_v.clear();
      dedx_v.clear();
      dqdx_points_ordered.clear();
      dedx_points_ordered.clear();
      dQdx_outliers_and_endpoints_removed_one_standard_deviation.clear();
      dEdx_outliers_and_endpoints_removed.clear();
      xyz_v_x_coords.clear();
      xyz_v_no_SCE_corrections_x_coords.clear();
      xyz_v_y_coords.clear();
      xyz_v_no_SCE_corrections_y_coords.clear();
      xyz_v_z_coords.clear();
      xyz_v_no_SCE_corrections_z_coords.clear();

      if ( Calo_v[pl]->dQdx().size() > 20 ) {

        if ( Calo_v[pl]->PlaneID().Plane == 0 )
	  u_plane_has_more_than_20_points = true;
      
	if ( Calo_v[pl]->PlaneID().Plane == 1 )
	  v_plane_has_more_than_20_points = true;

	if ( Calo_v[pl]->PlaneID().Plane == 2 )
	  y_plane_has_more_than_20_points = true;

      }
      
      // Loop through the vector and fill the relevant vectors.                                                                                                                                       
      for ( size_t hit_iter = 0; hit_iter < vmeta.size(); hit_iter++ ) {

	fill_vector = true;

	int ind = vmeta[hit_iter]->Index();

	// check that the traj point is in the calorimetry point                                                                     
	// and belongs to the plane we are interested in                                                                             
	if( pandora_track_h->at( muon_candidate_idx ).HasValidPoint(ind) && vhit[hit_iter]->WireID().Plane == pl && count < int( Calo_v[pl]->dQdx().size() ) ){

	  auto const& hit = vhit.at( hit_iter );

	  // Check to make sure the hit is not in one of the spikes.
	  float rms    = hit->RMS();
	  //float charge = hit->Integral();

	  //std::cout << "Plane #" << pl << ", Hit #" << hit_iter << " Width = " << rms << " Ticks, Charge = " << charge << " ADCs, dE/dx = " << Calo_v[pl]->dEdx()[count] << " MeV/cm." << std::endl;

	  // Filter on the peaks (we don't want to use these in those quantities).                                                  
	  if ( fabs( rms - 0.333 ) < 0.0001 || fabs( rms - 0.5 ) < 0.0001 || fabs( rms - 0.667 ) < 0.0001 || fabs( rms - 1.0 ) < 0.0001 || fabs( rms - 1.333 ) < 0.0001 || fabs( rms - 1.667 ) < 0.0001 || fabs( rms - 2.0 ) < 0.0001 || fabs( rms - 2.333 ) < 0.0001 || fabs( rms - 2.667 ) < 0.0001 || fabs( rms - 2.8572 ) < 0.0001 || fabs( rms - 3.0 ) < 0.0001 || fabs( rms - 3.333 ) < 0.0001 || fabs( rms - 3.667 ) < 0.0001 || fabs( rms - 4.0 ) < 0.0001 || fabs( rms - 4.333 ) < 0.0001 || fabs( rms - 4.667 ) < 0.0001 || fabs( rms - 5.0 ) < 0.0001 || fabs( rms - 5.333 ) < 0.0001 || fabs( rms - 5.667 ) < 0.0001 || fabs( rms - 6.0 ) < 0.0001 || fabs( rms - 6.333 ) < 0.0001 || fabs( rms - 6.5 ) < 0.0001 || fabs( rms - 7.0 ) < 0.0001 || fabs( rms - 7.333 ) < 0.0001 || fabs( rms - 7.5 ) < 0.0001 || fabs( rms - 7.667 ) < 0.0001 || fabs( rms - 8.0 ) < 0.0001 || fabs( rms - 8.5 ) < 0.0001 || fabs( rms - 9.0 ) < 0.0001 || fabs( rms - 9.5 ) < 0.0001 || fabs( rms - 10.0 ) < 0.0001 || fabs( rms - 10.5 ) < 0.0001 || fabs( rms - 11.0 ) < 0.0001 || fabs( rms - 11.5 ) < 0.0001 || fabs( rms - 12.0 ) < 0.0001 || fabs( rms - 12.5 ) < 0.0001 || fabs( rms - 13.0 ) < 0.0001 || fabs( rms - 13.5 ) < 0.0001 || fabs( rms - 14.0 ) < 0.0001 || fabs( rms - 14.5 ) < 0.0001 || fabs( rms - 15.0 ) < 0.0001 || fabs( rms - 15.5 ) < 0.0001 || fabs( rms - 16.0 ) < 0.0001 || fabs( rms - 16.5 ) < 0.0001 || fabs( rms - 17.0 ) < 0.0001 || fabs( rms - 17.5 ) < 0.0001 || fabs( rms - 18.0 ) < 0.0001 || fabs( rms - 18.5 ) < 0.0001 || fabs( rms - 19.0 ) < 0.0001 || fabs( rms - 19.5 ) < 0.0001 || fabs( rms - 20.0 ) < 0.0001 )
	    fill_vector = false;
	
	  if ( fill_vector == true ) {
	    dqdx_v.push_back( Calo_v[pl]->dQdx()[count] );
	    dedx_v.push_back( Calo_v[pl]->dEdx()[count] );
	    dqdx_points_ordered.push_back( Calo_v[pl]->dQdx()[count] );
	    dedx_points_ordered.push_back( Calo_v[pl]->dEdx()[count] );
	    xyz_v_x_coords.push_back( Calo_v[pl]->XYZ()[count].X() );
	    xyz_v_no_SCE_corrections_x_coords.push_back( Calo_v_no_SCE_corrections[pl]->XYZ()[count].X() );
	    xyz_v_y_coords.push_back( Calo_v[pl]->XYZ()[count].Y() );
	    xyz_v_no_SCE_corrections_y_coords.push_back( Calo_v_no_SCE_corrections[pl]->XYZ()[count].Y() );
	    xyz_v_z_coords.push_back( Calo_v[pl]->XYZ()[count].Z() );
	    xyz_v_no_SCE_corrections_z_coords.push_back( Calo_v_no_SCE_corrections[pl]->XYZ()[count].Z() );
	  }

	  count++;

	} // End of the conditional that the point actually exists along the track.

      } // End of the loop over the hits.

      if ( dedx_points_ordered.size() > 1 ) {

	// Now sort the vector (put it in order).
	std::sort( dqdx_points_ordered.begin(), dqdx_points_ordered.end() ); 
	std::sort( dedx_points_ordered.begin(), dedx_points_ordered.end() );
      	
	double dedx_median_value = dedx_points_ordered.at( int( dedx_points_ordered.size() / 2 ) );                                                                                                      
	
	// Find the standard deviation on this median.                                                                                                                                                    
	double rms_sum            = 0.;                                                                                                                                                                
	double difference_squared = 0.;                                                                                                                                                                  
	
	// Loop over the points to find the rms sum of the dedx points w.r.t. the median value.
	for ( size_t point_iter = 0; point_iter < dedx_v.size(); point_iter++ ) {                                                                                                                         
          
	  difference_squared = TMath::Power( ( dedx_median_value - dedx_v.at( point_iter ) ), 2);                                                                                                          
	  
	  rms_sum += difference_squared;                                                                                                                                                                   
                                                                                                                                                                                                          
	}                                                                                                                                                                                                  
                                                                                                                                                                                                         
	// Calculate the rms value.                                                                                                                                                                        
	double rms_value = TMath::Sqrt( ( rms_sum ) / double( dedx_v.size() ) );                                                                                                                           
      
	// Remove the points that are greater than two standard deviations from the median.                                                                                                                
	// New dE/dx vector.                                                                                                                                                                            
	std::vector< double > dEdx_outliers_and_endpoints_removed;                                                                                                                                         
	dEdx_outliers_and_endpoints_removed.clear();           
	
	for ( size_t energy_point_iter = 0; energy_point_iter < dedx_v.size(); energy_point_iter++ ) {        
	  
	  if ( fabs( dedx_median_value - dedx_v.at( energy_point_iter ) ) < ( 2.0 * rms_value ) )                                                                                                          
	    dEdx_outliers_and_endpoints_removed.push_back( dedx_v.at( energy_point_iter ) );                                                                                                              
	
	}
	
	// Find the direction of the track on this plane.
	// Fit this vector to a histogram in both directions.                                                                                                                                           
	TH1D* track_forward_hist    = new TH1D("track_forward_hist", "dE/dx vs. Point on Track", dEdx_outliers_and_endpoints_removed.size(), 0, dEdx_outliers_and_endpoints_removed.size() );          
	
	TH1D* track_backward_hist   = new TH1D("track_backward_hist", "dE/dx vs. Point on Track", dEdx_outliers_and_endpoints_removed.size(), 0, dEdx_outliers_and_endpoints_removed.size() );           
	
	// Loop through the points and fill both of the histograms.                                                                                                                                     
	for ( size_t histogram_iter = 0; histogram_iter < dEdx_outliers_and_endpoints_removed.size(); histogram_iter++ ) {                                                                               
	  
	  track_forward_hist->SetBinContent( histogram_iter, dEdx_outliers_and_endpoints_removed.at( histogram_iter ) );                                                                              
	  track_backward_hist->SetBinContent( ( dEdx_outliers_and_endpoints_removed.size() - 1 - histogram_iter ), dEdx_outliers_and_endpoints_removed.at( histogram_iter ) );                         
        
	}                                                                                                                                                                                               
                                                                       	
	double forward_slope  = 0.;                                                                                                                                                                    
	double backward_slope = 0.;                 
	
	// Declare the function and perform the fit.
	TF1 *st_y_function  = new TF1("st_y_function","[0]*x + [1]", 0, dEdx_outliers_and_endpoints_removed.size() );
		  
	// Forward fit.                                                                                                                                                                                 
	track_forward_hist->Fit(st_y_function,"R+");                                                                                                                                                     
	forward_slope = st_y_function->GetParameter(0);                                                                                                                                                 
	  
	// Backward fit.         
	track_backward_hist->Fit(st_y_function,"R+");                                                                                                                                                 
	backward_slope = st_y_function->GetParameter(0);                                                                                                                                             

	if ( backward_slope < forward_slope ) {
	
	  if ( Calo_v[pl]->PlaneID().Plane == 0 )
	    u_plane_has_track_forward = true;
	  
	  if ( Calo_v[pl]->PlaneID().Plane == 1 )
	    v_plane_has_track_forward = true;
	
	  if ( Calo_v[pl]->PlaneID().Plane == 2 )
	    y_plane_has_track_forward = true;
	
	}

      }
	
      count = 0;

    } // End of the loop over the planes.

    // Go through the nine cases of track orientation.                                                                                                                                                 
    if      ( u_plane_has_track_forward && u_plane_has_more_than_20_points && v_plane_has_track_forward && v_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else if (  u_plane_has_track_forward && u_plane_has_more_than_20_points && y_plane_has_track_forward && y_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else if ( v_plane_has_track_forward && v_plane_has_more_than_20_points && y_plane_has_track_forward && y_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else if ( !u_plane_has_more_than_20_points && !v_plane_has_more_than_20_points && y_plane_has_more_than_20_points && y_plane_has_track_forward )
      first_points_of_track_are_start = true;

    else if ( !u_plane_has_more_than_20_points && !y_plane_has_more_than_20_points && v_plane_has_more_than_20_points && v_plane_has_track_forward )
      first_points_of_track_are_start = true;

    else if ( !v_plane_has_more_than_20_points && !y_plane_has_more_than_20_points && u_plane_has_more_than_20_points && u_plane_has_track_forward )
      first_points_of_track_are_start = true;

    else if ( !u_plane_has_more_than_20_points && !v_plane_has_more_than_20_points && !y_plane_has_more_than_20_points )
      first_points_of_track_are_start = true;

    else
      first_points_of_track_are_start = false;

    std::cout << "Above the loop over the plane info for vertexing." << std::endl;

    for ( size_t plane_iter = 0; plane_iter < Calo_v.size(); plane_iter++ ) {

      auto xyz_v_coords = Calo_v[plane_iter]->XYZ();

      if ( Calo_v[plane_iter]->PlaneID().Plane != 0 && Calo_v[plane_iter]->PlaneID().Plane != 1 && Calo_v[plane_iter]->PlaneID().Plane != 2 ) 
	continue;

      vertex_location_x = -10000.;
      vertex_location_y = -10000.;
      vertex_location_z = -10000.;

      std::cout << "Above setting the vertex location." << std::endl;

      // Orient the track based on the verdict of track orientation.
      if ( xyz_v_coords.size() > 0 && xyz_v_no_SCE_corrections_x_coords.size() > 0 ) {

	if ( first_points_of_track_are_start == true ) {

	  vertex_location_x          = xyz_v_coords.at( 0 ).X() - ( flash_time * 0.1098 );                                                                                                              
	  vertex_location_y          = xyz_v_coords.at( 0 ).Y();                                                                                                                                          
	  vertex_location_z          = xyz_v_coords.at( 0 ).Z();                                                                                                                                           

	  other_end_location_x       = xyz_v_x_coords.at( xyz_v_x_coords.size() - 1 ) - ( flash_time * 0.1098 );                                                                                         
	  other_end_location_y       = xyz_v_y_coords.at( xyz_v_y_coords.size() - 1 );                                                                                                                    
	  other_end_location_z       = xyz_v_z_coords.at( xyz_v_z_coords.size() - 1 );                                                                                                                 
                                                                                                                                                                                                           
	}

	else {

	  vertex_location_x          = xyz_v_coords.at( xyz_v_coords.size() - 1 ).X() - ( flash_time * 0.1098 );                                                                                     
	  vertex_location_y          = xyz_v_coords.at( xyz_v_coords.size() - 1 ).Y();                                                                                                                 
	  vertex_location_z          = xyz_v_coords.at( xyz_v_coords.size() - 1 ).Z();                                                                                                                            
	  other_end_location_x       = xyz_v_x_coords.at( 0 ) - ( flash_time * 0.1098 );                                                                                                                 
	  other_end_location_y       = xyz_v_y_coords.at( 0 );                                                                                                                                            
	  other_end_location_z       = xyz_v_z_coords.at( 0 );                                                                                                                                                                                                                                                                                                                                                       
	}

      }

    }

    std::cout << "Below the loop over the plane info for vertexing." << std::endl;

    auto const& trktraj_repeat = pandora_track_h->at( muon_candidate_idx ).Trajectory();

    // Flip these for whatever reason.
    auto firstValid_again     = trktraj_repeat.LastValidPoint();
    auto lastValid_again      = trktraj_repeat.FirstValidPoint();
    
    if ( first_points_of_track_are_start == true ) {

      vertex_location_x          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_again ).X() - ( 0.1098 * flash_time );
      vertex_location_y          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_again ).Y();
      vertex_location_z          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_again ).Z();
      
      other_end_location_x       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid_again ).X() - ( 0.1098 * flash_time );
      other_end_location_y       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid_again ).Y();
      other_end_location_z       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid_again ).Z();

    }

    else {

      vertex_location_x          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid_again ).X() - ( 0.1098 * flash_time );
      vertex_location_y          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid_again ).Y();
      vertex_location_z          = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( lastValid_again ).Z();
      
      other_end_location_x       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_again ).X() - ( 0.1098 * flash_time );
      other_end_location_y       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_again ).Y();
      other_end_location_z       = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_again ).Z();
	
    }

    // Correct the vertex points for SCE in this case.
    geo::Vector_t track_vertex_offsets = SCE_data->GetCalPosOffsets(geo::Point_t{vertex_location_x, vertex_location_y, vertex_location_z});

    vertex_location_x    = vertex_location_x - track_vertex_offsets.X();
    vertex_location_y    = vertex_location_y + track_vertex_offsets.Y();
    vertex_location_z    = vertex_location_z + track_vertex_offsets.Z();

    // Do the same for the other end points for SCE.
    geo::Vector_t other_end_offsets   = SCE_data->GetCalPosOffsets(geo::Point_t{other_end_location_x, other_end_location_y, other_end_location_z});
  
    other_end_location_x = other_end_location_x - other_end_offsets.X();
    other_end_location_y = other_end_location_y + other_end_offsets.Y();
    other_end_location_z = other_end_location_z + other_end_offsets.Z();
	
    // Identify the points on the muon track for orientating the track in the 3D fitting algorithm.                                                                                                       
    auto const& trktraj = pandora_track_h->at( muon_candidate_idx ).Trajectory();

    auto firstValid     = trktraj.FirstValidPoint();
    auto lastValid      = trktraj.LastValidPoint();

    // Loop through the longest track and calculate the piecewise track distance.                                                                                                                         
    for ( size_t point_iter = firstValid; point_iter < ( pandora_track_h->at( muon_candidate_idx ).NumberTrajectoryPoints() - 1 ); point_iter++ )  {

      double point0_X = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter ).X();
      double point0_Y = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter ).Y();
      double point0_Z = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter ).Z();

      double point1_X = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter + 1 ).X();
      double point1_Y = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter + 1 ).Y();
      double point1_Z = track_h->at( muon_candidate_idx ).LocationAtPoint( point_iter + 1 ).Z();

      if ( point_iter == firstValid ) {

	muon_track_first_point_x = point0_X - ( 0.1098 * flash_time );
	muon_track_first_point_y = point0_Y;
	muon_track_first_point_z = point0_Z;

	// Find the offsets on each of the points.                                                                                                                                                        
	geo::Vector_t muon_vertex_offsets = SCE_data->GetCalPosOffsets(geo::Point_t{muon_track_first_point_x,muon_track_first_point_y,muon_track_first_point_z});

	muon_track_first_point_x  =  muon_track_first_point_x - muon_vertex_offsets.X();
	muon_track_first_point_y  =  muon_track_first_point_y + muon_vertex_offsets.Y();
	muon_track_first_point_z  =  muon_track_first_point_z + muon_vertex_offsets.Z();

      }

      if ( point_iter == ( lastValid - 1 ) ) {

	muon_track_last_point_x = point1_X - ( 0.1098 * flash_time );
	muon_track_last_point_y = point1_Y;
	muon_track_last_point_z = point1_Z;

	// Find the offsets on each of the points.                                                                                                                                                        
	geo::Vector_t muon_vertex_offsets = SCE_data->GetCalPosOffsets(geo::Point_t{muon_track_last_point_x,muon_track_last_point_y,muon_track_last_point_z});

	muon_track_last_point_x  =  muon_track_last_point_x - muon_vertex_offsets.X();
	muon_track_last_point_y  =  muon_track_last_point_y + muon_vertex_offsets.Y();
	muon_track_last_point_z  =  muon_track_last_point_z + muon_vertex_offsets.Z();

      }

    }
				
    muon_candidate_length               = pandora_track_h->at( muon_candidate_idx ).Length();
    muon_candidate_kinetic_energy_range = ( TMath::Sqrt( ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) * ( p_calculator_from_length.GetTrackMomentum( muon_candidate_length, 13 ) * 1000. ) + 105.7 * 105.7 ) - 105.7 ); 

    if ( muon_track_first_point_x > 20. && muon_track_first_point_x < 236.35 && muon_track_first_point_y > -96.5 && muon_track_first_point_y < 96.5 && muon_track_first_point_z > 20.0 && muon_track_first_point_z < 1016.8 )
      at_least_one_terminal_point_in_FV = 1;

    if ( muon_track_last_point_x > 20. && muon_track_last_point_x < 236.35 && muon_track_last_point_y > -96.5 && muon_track_last_point_y < 96.5 && muon_track_last_point_z > 20.0 && muon_track_last_point_z < 1016.8 )
      at_least_one_terminal_point_in_FV = 1;

    std::cout << "TRACK ORIENTATION PRINT STATEMENTS:" << std::endl;
    std::cout << "Location of the truth vertex: x = " << nu_vtx_x_truth << " cm y = " << nu_vtx_y_truth << " cm z = " << nu_vtx_z_truth << " cm." << std::endl;
    std::cout << "Location of the reco vertex BEFORE the changes: x = " << vertex_location_x  << " cm y = " << vertex_location_y << " cm z = " << vertex_location_z << " cm." << std::endl;
    std::cout << "Muon candidate length = " << muon_candidate_length << " cm." << std::endl;
    
    if ( fabs( muon_candidate_length - longest_TPCObject_track_length ) > 0.001 )
      std::cout << "Longest TPC Object Track Length = " << longest_TPCObject_track_length << " cm." << std::endl;

    // Flip the coordinates if they are not in the right order.
    double distance_track_start_truth_vertex = TMath::Sqrt( ( nu_vtx_x_truth - muon_track_first_point_x ) * ( nu_vtx_x_truth - muon_track_first_point_x ) + ( nu_vtx_y_truth - muon_track_first_point_y ) * ( nu_vtx_y_truth - muon_track_first_point_y ) + ( nu_vtx_z_truth - muon_track_first_point_z ) * ( nu_vtx_z_truth - muon_track_first_point_z ) );
    
    double distance_track_end_truth_vertex   = TMath::Sqrt( (nu_vtx_x_truth - muon_track_last_point_x ) * ( nu_vtx_x_truth - muon_track_last_point_x ) + ( nu_vtx_y_truth - muon_track_last_point_y )* ( nu_vtx_y_truth - muon_track_last_point_y ) + ( nu_vtx_z_truth - muon_track_last_point_z ) * ( nu_vtx_z_truth - muon_track_last_point_z ) );
    
    if ( distance_track_end_truth_vertex < distance_track_start_truth_vertex ) {
      
      double muon_track_first_x_holder = muon_track_first_point_x;
      double muon_track_first_y_holder = muon_track_first_point_y;
      double muon_track_first_z_holder = muon_track_first_point_z;
      
      double muon_track_last_x_holder  = muon_track_last_point_x;
      double muon_track_last_y_holder  = muon_track_last_point_y;
      double muon_track_last_z_holder  = muon_track_last_point_z;
      
      // Switch the coordinates.
      muon_track_first_point_x = muon_track_last_x_holder;
      muon_track_first_point_y = muon_track_last_y_holder;
      muon_track_first_point_z = muon_track_last_z_holder;
      
      muon_track_last_point_x  = muon_track_first_x_holder;
      muon_track_last_point_y  = muon_track_first_y_holder;
      muon_track_last_point_z  = muon_track_first_z_holder;
      
      distance_between_truth_vertex_and_closest_track_terminal_point = distance_track_end_truth_vertex;
      
      std::cout << "Distance between truth vertex and closest track terminal point = " << distance_track_end_truth_vertex << " cm." << std::endl;
      
    }

    else {

      distance_between_truth_vertex_and_closest_track_terminal_point = distance_track_start_truth_vertex;

      std::cout << "Distance between truth vertex and closest track terminal point = " << distance_track_start_truth_vertex << " cm." << std::endl;

    }
    
    muon_candidate_is_contained = 1;

      // Include print statements on the two terminal points of the muon candidate track.
      std::cout << "Muon track starting coordinates (w/ trigger offset): x = " << muon_track_first_point_x << " cm y = " << muon_track_first_point_y << " cm z = " << muon_track_first_point_z << " cm." << std::endl;

      if ( muon_track_first_point_x < 3. || muon_track_first_point_x > 253.35 || muon_track_first_point_y < -113.5 || muon_track_first_point_y > 113.5 || muon_track_first_point_z < 3.0 || muon_track_first_point_z > 1033.6 ) {
	std::cout << "The first track endpoint is not contained within the volume under consideration in the detector." << std::endl;
	muon_candidate_is_contained = 0;
      }

      std::cout << "Muon track ending coordinates (w/ trigger offset): x = " << muon_track_last_point_x << " cm y = " << muon_track_last_point_y << " cm z = " << muon_track_last_point_z << " cm." << std::endl;

      if ( muon_track_last_point_x < 3. || muon_track_last_point_x > 253.35 || muon_track_last_point_y < -113.5 || muon_track_last_point_y > 113.5 || muon_track_last_point_z < 3.0 || muon_track_last_point_z > 1033.6 ) {
	std::cout << "The last track endpoint is not contained within the volume under consideration in the detector."<< std::endl;
	 muon_candidate_is_contained = 0;
      }

      num_of_tracks_originating_from_vertex                          = 1; // Start with the muon.
      num_of_tracks_originating_from_end                             = 0;
      average_track_length_originating_from_vertex_other_than_muon   = 0.;
      average_track_length_originating_from_end                      = 0.;
      
      for ( int tpcobject_track_iter = 0; tpcobject_track_iter < int( tracks.size() ); tpcobject_track_iter++ ) {

	// Skip the longest track (the vertex is found with respect to that one).                                                                                                                      
	if (  tracks.at( tpcobject_track_iter )->ID() == muon_candidate_ID )
	  continue;

	// Find the starting and ending x, y, and z coordinates of the muon track and of the track that you are currently looking at.                                                                  
	auto track = tracks.at( tpcobject_track_iter );

	auto first_valid_point_idx = track->FirstValidPoint();
	auto last_valid_point_idx  = track->LastValidPoint();
	
	auto other_track_length    = track->Length();

	auto first_valid_point     = track->LocationAtPoint( first_valid_point_idx );
	auto last_valid_point      = track->LocationAtPoint( last_valid_point_idx );

	double other_track_x0      = first_valid_point.X() - ( 0.1098 * flash_time );
	double other_track_y0      = first_valid_point.Y();
	double other_track_z0      = first_valid_point.Z();
	double other_track_x1      = last_valid_point.X() - ( 0.1098 * flash_time );
	double other_track_y1      = last_valid_point.Y();
	double other_track_z1      = last_valid_point.Z();

	// Offset these track points for the space charge effect.                                                                                                                                    
	geo::Vector_t other_track_first_point_offsets  = SCE_data->GetCalPosOffsets(geo::Point_t{other_track_x0,other_track_y0,other_track_z0});
	
	double other_track_x0_with_offset = other_track_x0 - other_track_first_point_offsets.X();
	double other_track_y0_with_offset = other_track_y0 + other_track_first_point_offsets.Y();
	double other_track_z0_with_offset = other_track_z0 + other_track_first_point_offsets.Z();
	
	geo::Vector_t other_track_second_point_offsets = SCE_data->GetCalPosOffsets(geo::Point_t{other_track_x1,other_track_y1,other_track_z1});
	
	double other_track_x1_with_offset = other_track_x1 - other_track_second_point_offsets.X();
	double other_track_y1_with_offset = other_track_y1 + other_track_second_point_offsets.Y();
	double other_track_z1_with_offset = other_track_z1 + other_track_second_point_offsets.Z();

	// Find which point is closer to both the vertex and the end of the muon candidate.                                                                                                           
	double closest_distance_vertex = 0.;
	double closest_distance_end    = 0.;

	// Change this up to use the truth vertex.                                                                                                                                                         
        double distance_vertex_point0  = TMath::Sqrt( ( vertex_location_x - other_track_x0_with_offset ) * ( vertex_location_x - other_track_x0_with_offset ) + ( vertex_location_y - other_track_y0_with_offset ) * ( vertex_location_y - other_track_y0_with_offset ) + ( vertex_location_z - other_track_z0_with_offset ) * ( vertex_location_z - other_track_z0_with_offset ) );

        double distance_vertex_point1  = TMath::Sqrt( ( vertex_location_x - other_track_x1_with_offset ) * ( vertex_location_x - other_track_x1_with_offset ) + ( vertex_location_y - other_track_y1_with_offset ) * ( vertex_location_y - other_track_y1_with_offset ) + ( vertex_location_z - other_track_z1_with_offset ) * ( vertex_location_z - other_track_z1_with_offset ) );

        double distance_end_point0     = TMath::Sqrt( ( other_end_location_x - other_track_x0_with_offset ) * ( other_end_location_x - other_track_x0_with_offset ) + ( other_end_location_y - other_track_y0_with_offset ) * ( other_end_location_y - other_track_y0_with_offset ) + ( other_end_location_z - other_track_z0_with_offset ) * ( other_end_location_z - other_track_z0_with_offset ) );

        double distance_end_point1     = TMath::Sqrt( ( other_end_location_x - other_track_x1_with_offset ) * ( other_end_location_x - other_track_x1_with_offset ) + ( other_end_location_y - other_track_y1_with_offset ) * ( other_end_location_y - other_track_y1_with_offset ) + ( other_end_location_z - other_track_z1_with_offset ) * ( other_end_location_z - other_track_z1_with_offset ) );

	if ( distance_vertex_point0 < distance_vertex_point1 ) {

	  closest_distance_vertex = distance_vertex_point0;

	  std::cout << "The starting point of the truth track (other than the muon candidate) is closer; it's " << closest_distance_vertex << " away from the truth vertex." << std::endl;

	}

	else {

	  closest_distance_vertex = distance_vertex_point1;

	  std::cout << "The ending point of the truth track (other than the muon candidate) is closer; it's "<< closest_distance_vertex << " away from the truth vertex." <<std::endl;

	}
	
	
	if ( distance_end_point0 < distance_end_point1 ) {
	  closest_distance_end = distance_end_point0;
	 }
	
	else {
	  closest_distance_end = distance_end_point1;
	}

	// Now, determine if the track originates at the vertex or at the other end.                                                                                                                  
	if ( closest_distance_vertex < closest_distance_end ) { // && closest_distance_vertex < 3.0 ) { // The second cut is to ensure that this is a track that actually originates at the vertex.       

	  num_of_tracks_originating_from_vertex                        += 1;
	  average_track_length_originating_from_vertex_other_than_muon += other_track_length;
	

	}

	else {

	  num_of_tracks_originating_from_end                           += 1;
	  average_track_length_originating_from_end                    += other_track_length;
	  
	}

      } // End of the loop over the tracks in the event.

      if ( num_of_tracks_originating_from_vertex > 1 ) {

	average_track_length_originating_from_vertex_other_than_muon /= num_of_tracks_originating_from_vertex;

      }

      if ( num_of_tracks_originating_from_end > 0 ) {

	average_track_length_originating_from_end /= num_of_tracks_originating_from_end;

      }

    // Flip the orientation of the track based on how many tracks originate from the vertex vs. the other end of the track.
      if ( ( num_of_tracks_originating_from_end > ( num_of_tracks_originating_from_vertex - 1 ) || ( ( num_of_tracks_originating_from_vertex - 1 ) == num_of_tracks_originating_from_end && average_track_length_originating_from_end > average_track_length_originating_from_vertex_other_than_muon && number_of_tracks_in_TPCObject > 1 ) ) && !( muon_candidate_is_contained == 1 && muon_candidate_length < 150. ) ) {

      double vertex_location_x_temporary    = vertex_location_x;
      double vertex_location_y_temporary    = vertex_location_y;
      double vertex_location_z_temporary    = vertex_location_z;

      double other_end_location_x_temporary = other_end_location_x;
      double other_end_location_y_temporary = other_end_location_y;
      double other_end_location_z_temporary = other_end_location_z;

      vertex_location_x                     = other_end_location_x_temporary;
      vertex_location_y                     = other_end_location_y_temporary;
      vertex_location_z                     = other_end_location_z_temporary;

      other_end_location_x                  = vertex_location_x_temporary;
      other_end_location_y                  = vertex_location_y_temporary;
      other_end_location_z                  = vertex_location_z_temporary;
      
      }

    std::cout << "Number of tracks in the TPC Object = " << number_of_tracks_in_TPCObject << "." << std::endl;
    std::cout << "Number of tracks coming out of vertex of TPC Object (other than muon candidate) = " << ( num_of_tracks_originating_from_vertex - 1 ) << "." << std::endl;
    std::cout << "Number of tracks coming out of end of one of the tracks = " << num_of_tracks_originating_from_end << "." << std::endl;

    fails_fiducial_volume_req   = 0;

    if ( vertex_location_x < 20.0 || vertex_location_x > 236.35 || vertex_location_y < -96.5 || vertex_location_y > 96.5 || vertex_location_z < 20.0 || vertex_location_z > 1016.8 )
      fails_fiducial_volume_req = 1;

    bool fails_3cm_containment_cut = false;

    // Loop through the track in the TPC Object and find if they originate at the vertex.  Skip the longest track.                                                                                         
    for ( size_t tpcobject_track_iter = 0; tpcobject_track_iter < tracks.size(); tpcobject_track_iter++ ) {

      // Skip the longest track (the vertex is found with respect to that one).                                                                                                                           
      if (  tracks.at( tpcobject_track_iter )->ID() == muon_candidate_ID )
	continue;

      // Find the starting and ending x, y, and z coordinates of the muon track and of the track that you are currently looking at.                                                                       
      auto track                                     = *(tracks.at( tpcobject_track_iter ));
      auto first_valid_point_idx                     = track.FirstValidPoint();
      auto last_valid_point_idx                      = track.LastValidPoint();

      auto first_valid_point                         = track.LocationAtPoint( first_valid_point_idx );
      auto last_valid_point                          = track.LocationAtPoint( last_valid_point_idx );

      double other_track_x0                          = first_valid_point.X() - ( 0.1098 * flash_time );
      double other_track_y0                          = first_valid_point.Y();
      double other_track_z0                          = first_valid_point.Z();
      double other_track_x1                          = last_valid_point.X() - ( 0.1098 * flash_time );
      double other_track_y1                          = last_valid_point.Y();
      double other_track_z1                          = last_valid_point.Z();

      // Offset these track points for the space charge effect.                                                                                                                                         
      geo::Vector_t other_track_first_point_offsets  = sce->GetCalPosOffsets(geo::Point_t{other_track_x0,other_track_y0,other_track_z0});

      double other_track_x0_with_offset              = other_track_x0 - other_track_first_point_offsets.X();
      double other_track_y0_with_offset              = other_track_y0 + other_track_first_point_offsets.Y();
      double other_track_z0_with_offset              = other_track_z0 + other_track_first_point_offsets.Z();

      geo::Vector_t other_track_second_point_offsets = sce->GetCalPosOffsets(geo::Point_t{other_track_x1,other_track_y1,other_track_z1});

      double other_track_x1_with_offset              = other_track_x1 - other_track_second_point_offsets.X();
      double other_track_y1_with_offset              = other_track_y1 + other_track_second_point_offsets.Y();
      double other_track_z1_with_offset              = other_track_z1 + other_track_second_point_offsets.Z();

      // Test all of the points for their position in the TPC.                                                                                                                                            
      // 3 cm containment requirement.                                                                                                                                                                    
      if ( other_track_x0_with_offset < 3.0 || other_track_x0_with_offset > 253.35 )
        fails_3cm_containment_cut = true;

      if ( other_track_y0_with_offset < -113.5 || other_track_y0_with_offset > 113.5 )
        fails_3cm_containment_cut = true;

      if ( other_track_z0_with_offset < 3.0 || other_track_z0_with_offset > 1033.8 )
        fails_3cm_containment_cut = true;

      if ( other_track_x1_with_offset < 3.0 || other_track_x1_with_offset > 253.35 )
	fails_3cm_containment_cut = true;

      if ( other_track_y1_with_offset < -113.5 || other_track_y1_with_offset > 113.5 )
        fails_3cm_containment_cut = true;

      if ( other_track_z1_with_offset < 3.0 || other_track_z1_with_offset > 1033.8 )
        fails_3cm_containment_cut = true;

    }

    bool vertex_point_fails_containment = false;

    // Do this because the track endpoints are not SCE-corrected.
    auto const& trktraj_in_loop   = pandora_track_h->at( muon_candidate_idx ).Trajectory();

    // Flip these for whatever reason.                                                                                                                                                                  
    auto firstValid_in_loop              = trktraj_in_loop.LastValidPoint();

    double muon_vertex_x                 = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_in_loop ).X() - ( 0.1098 * flash_time );
    double muon_vertex_y                 = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_in_loop ).Y();
    double muon_vertex_z                 = pandora_track_h->at( muon_candidate_idx ).LocationAtPoint( firstValid_in_loop ).Z();

    double vtx_dist_to_track_start = TMath::Sqrt( ( muon_vertex_x - vertex_location_x ) * ( muon_vertex_x - vertex_location_x ) + ( muon_vertex_y - vertex_location_y ) * ( muon_vertex_y - vertex_location_y ) + ( muon_vertex_z - vertex_location_z ) * ( muon_vertex_z - vertex_location_z ) );

    double end_dist_to_track_start = TMath::Sqrt( ( muon_vertex_x - other_end_location_x ) * ( muon_vertex_x - other_end_location_x ) + ( muon_vertex_y - other_end_location_y ) * ( muon_vertex_y - other_end_location_y ) + ( muon_vertex_z - other_end_location_z ) * ( muon_vertex_z - other_end_location_z ) );

    if ( vtx_dist_to_track_start < end_dist_to_track_start ) {

      if ( vertex_location_x < 3.0 || vertex_location_x > 253.35 || vertex_location_y < -113.5 || vertex_location_x > 113.5 || vertex_location_z < 3.0 || vertex_location_z > 1033.8 )
	vertex_point_fails_containment = true;

    }

    else {

      if ( other_end_location_x < 3.0 || other_end_location_x > 253.35 || other_end_location_y < -113.5 || other_end_location_x > 113.5 || other_end_location_z < 3.0 || other_end_location_z > 1033.8 )
	vertex_point_fails_containment = true;

    }   

    // Calculate the energy of the track based on this information.
      if ( muon_candidate_is_contained == 0 ) {

	// Use the backwards assumption for the MCS here.
	if ( vertex_point_fails_containment == true ) {

	  muon_candidate_kinetic_energy_mcs = ( TMath::Sqrt( ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->bwdMomentum() * 1000. ) * ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->bwdMomentum() * 1000. ) + 105.7 * 105.7 ) - 105.7 );

	  std::cout << "Non-contained Muon backwards momentum from MCS = " << pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->bwdMomentum() << " GeV." << std::endl;
	  std::cout << "Length of the muon candidate = " << muon_candidate_length << " cm." << std::endl;

	}

	// Use the forwards assumption for the MCS here.
	else {

	  muon_candidate_kinetic_energy_mcs = ( TMath::Sqrt( ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->fwdMomentum() * 1000. ) * ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->fwdMomentum() * 1000. ) + 105.7 * 105.7 ) - 105.7 );

	  std::cout << "Non-contained Muon forwards momentum from MCS = " << pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->fwdMomentum() << " GeV." << std::endl;
	  std::cout << "Length of the muon candidate = " << muon_candidate_length << " cm." << std::endl;

	}

      }

      else {

	double dist_pandora_vtx_to_my_vtx = TMath::Sqrt( ( vertex_location_x - pandora_vertex_x ) * ( vertex_location_x - pandora_vertex_x ) + ( vertex_location_y - pandora_vertex_y ) * ( vertex_location_y - pandora_vertex_y ) + ( vertex_location_z - pandora_vertex_z ) * ( vertex_location_z - pandora_vertex_z ) );
	
	std::cout << "This event is contained." << std::endl;
	std::cout << "The distance from my vertex to the pandora vertex = " << dist_pandora_vtx_to_my_vtx << " cm." << std::endl;

	if ( dist_pandora_vtx_to_my_vtx < ( muon_candidate_length / 2.0 ) ) {

	  muon_candidate_kinetic_energy_mcs = ( TMath::Sqrt( ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->fwdMomentum() * 1000. ) * ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->fwdMomentum() * 1000. ) + 105.7 * 105.7 ) - 105.7 );
	  std::cout << "Using forward momentum." << std::endl;

	}

	else {

	  muon_candidate_kinetic_energy_mcs = ( TMath::Sqrt( ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->bwdMomentum() * 1000. ) * ( pandora_MCSResult_pointers_h.at( correct_muon_candidate_track.key() )->bwdMomentum() * 1000. ) + 105.7 * 105.7 ) - 105.7 );
	  std::cout << "Using backward momentum." << std::endl;

	}

      }

      std::cout << "Truth muon kinetic energy = " << truth_muon_kinetic_energy << " MeV." << std::endl;
      std::cout << "Muon candidate kinetic energy MCS = " << muon_candidate_kinetic_energy_mcs << " MeV." << std::endl;
      std::cout << "Muon candidate kinetic energy range = " << muon_candidate_kinetic_energy_range << " MeV." << std::endl;

      std::cout << "Abs Difference Between MCS KE and Truth = " << fabs( muon_candidate_kinetic_energy_mcs - truth_muon_kinetic_energy ) << " MeV." << std::endl;
      std::cout << "Is muon candidate contained? " << muon_candidate_is_contained << "." << std::endl;

      std::cout << "Below the containment conditionals." << std::endl;
    
      sum_of_TPCObject_track_lengths = 0.;
  
      // See how many of the tracks in the TPC object have a length greater than 5 cm.
      for ( size_t tpcobj_track_iter = 0; tpcobj_track_iter < tracks.size(); tpcobj_track_iter++ ) {
      
	sum_of_TPCObject_track_lengths += tracks.at( tpcobj_track_iter )->Length();

      }

      // Include a variable for the vertex position before it is t0-corrected.
      ubxsec_muon_phi                = ubxsec_event->slc_longesttrack_phi[nu_slc_idx];
      ubxsec_muon_cos_theta           = ubxsec_event->slc_longesttrack_theta[nu_slc_idx];

      // Calculate the distance between the truth vertex and the reconstructed vertex.
      truth_reco_vtx_distance          = TMath::Sqrt( ( nu_vtx_x_truth - vertex_location_x ) * ( nu_vtx_x_truth - vertex_location_x )  + ( nu_vtx_y_truth - vertex_location_y ) * ( nu_vtx_y_truth - vertex_location_y ) + ( nu_vtx_z_truth - vertex_location_z ) * ( nu_vtx_z_truth - vertex_location_z ) );
      
      fails_residuals_req                         = 0;
      muon_candidate_residuals                    = ubxsec_event->slc_muoncandidate_residuals_std.at(nu_slc_idx);
    
      if ( muon_candidate_residuals > 2.5 )
	fails_residuals_req    = 1;

      fails_perc_used_hits_in_cluster             = 0;
      muon_candidate_percent_used_hits_in_cluster = ubxsec_event->slc_muoncandidate_perc_used_hits_in_cluster.at(nu_slc_idx);
      
      if ( muon_candidate_percent_used_hits_in_cluster < 0.7 )
	fails_perc_used_hits_in_cluster           = 1;
      
      if ( fails_fiducial_volume_req == 1 ) {

	std::cout << "Failed fiducial volume." << std::endl;
	failure_map["k_fiducial_volume_cut"]                   = false;
	reason                                                 = "fail_fiducial_volume";
	total_num_events_failed_fiducial_volume++;
	
	total_num_events_failing_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
	total_num_events_failed_fiducial_volume_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
	

      }

      std::cout << "In the passing section." << std::endl;

      total_num_events_passing_weighted += spline_fix_mcweight * central_value_mcweight * rootino_fix_mcweight * ppfx_value_mcweight;
        
      total_num_events_passing++;
    
      entire_event_is_contained = 1;

      if ( fails_3cm_containment_cut == true || muon_candidate_is_contained == 0 )  {
	std::cout << "Entire event is not contained." << std::endl;
	entire_event_is_contained = 0;
      }
      
      // Fill both trees.                                                                                                                                                                   
      _passing_events_tree->Fill();
      
      std::cout << "Filled other trees." << std::endl;
      
      selection_result.SetSelectionStatus(true);
    
      ubxsec_event->is_selected = true;
      
      std::cout << "Set event as selected." << std::endl;
      
      // Grab the selected TPCObject
      std::vector<art::Ptr<ubana::TPCObject>> tpcobj_v;
      art::fill_ptr_vector(tpcobj_v, tpcobj_h);
      art::Ptr<ubana::TPCObject> new_tpcobj = tpcobj_v.at(slice_index_1);
      
      std::cout << "Grabbed the selected TPC Object." << std::endl;
    
      if (_debug) std::cout << "[UBXSec] >>>>>>>>>>>>>>>>>>>>>> Selected TPCObject with index " << slice_index_1 << std::endl;
      
      // Prepare the tpcobj output
      std::vector<art::Ptr<ubana::TPCObject>>  out_tpcobj_v;
      out_tpcobj_v.resize(1);
      out_tpcobj_v.at(0) = new_tpcobj;
    
      std::cout << "Prepared the TPC Object output." << std::endl;
      
      selectionResultVector->emplace_back(std::move(selection_result));
      util::CreateAssn(*this, e, *selectionResultVector, out_tpcobj_v, *assnOutSelectionResultTPCObject);
    
      std::cout << "Moved the selection result." << std::endl;
    
      // For the TPCNeutrinoID Filter
      util::CreateAssn(*this, e, muon_candidate_track_per_slice_v.at(slice_index_1), neutrino_candidate_vertex_per_slice_v.at(slice_index_1), *vertexTrackAssociations);
      util::CreateAssn(*this, e, muon_candidate_pfparticle_per_slice_v.at(slice_index_1), neutrino_candidate_vertex_per_slice_v.at(slice_index_1), *vertexPFParticleAssociations);

      std::cout << "At the bottom of the passing section." << std::endl;
    
  } // End of the case in which the event is selected.
  
  else {
   
    std::cout << "This event failed because of " << reason << "."<< std::endl;
   
  }
  
  // Total # of Events            
  std::cout << "The total number of events = "                                                                         << total_num_events                                                        << "." << std::endl;
  std::cout << "The total number of events, weighted = "                                                               << total_num_events_weighted                                               << "." << std::endl;
  std::cout << "The total number of events passing = "                                                                 << total_num_events_passing                                                << "." << std::endl;
  std::cout << "The total number of events passing, weighted = "                                                       << total_num_events_passing_weighted                                       << "." << std::endl; 
  std::cout << "The total number of events failing = "                                                                 << total_num_events_failing                                                << "." << std::endl;
  std::cout << "The total number of events failing, weighted = "                                                       << total_num_events_failing_weighted                                       << "." << std::endl;
  
  // Cut Breakdown                                                              
  std::cout << "The total number of events with no 'simpleFlashBeam' flashes = "                                       << total_num_events_failed_beam_disc_flashes                               << "." << std::endl;
  std::cout << "The total number of events with no 'simpleFlashBeam' flashes, weighted = "                             << total_num_events_failed_beam_disc_flashes_weighted                      << "." << std::endl;
  std::cout << "The total number of events without a 'simpleFlashBeam' flash of > 50 PEs in the beamspill window = "   << total_num_events_failed_beam_spill_flash                                << "." << std::endl;
  std::cout << "The total number of events without a 'simpleFlashBeam' flash of > 50 PEs in the beamspill window, weighted = "   << total_num_events_failed_beam_spill_flash_weighted             << "." << std::endl;
  std::cout << "The total number of events with no TPC Objects = "                                                     << total_num_events_failed_has_slices                                      << "." << std::endl;
  std::cout << "The total number of events with no TPC Objects, weighted = "                                           << total_num_events_failed_has_slices_weighted                             << "." << std::endl;
  std::cout << "The total number of events with no TPC Object that is tagged as the neutrino = "                       << total_num_events_failed_has_slice_tagged_as_neutrino                    << "." << std::endl;
  std::cout << "The total number of events with a TPC Object failing the > 40 cm track length cut = "                  << total_num_events_failed_track_length                                    << "." << std::endl;
  std::cout << "The total number of events failing the fiducial volume cut = "                                         << total_num_events_failed_fiducial_volume                                 << "." << std::endl;
  std::cout << "The total number of events failing the fiducial volume cut, weighted = "                               << total_num_events_failed_fiducial_volume_weighted                        << "." << std::endl;
  std::cout << "The total number of events failing the number of tracks cut = "                                        << total_num_events_failed_ntrack                                          << "." << std::endl;
  std::cout << "The total number of events failing the number of tracks cut, weighted = "                              << total_num_events_failed_ntrack_weighted                                 << "." << std::endl;
  std::cout << "The total number of events failing the residuals std up cut = "                                        << total_num_events_failed_residuals_std_up                                << "." << std::endl;
  std::cout << "The total number of events failing the residuals std up cut, weighted = "                              << total_num_events_failed_residuals_std_up_weighted                       << "." << std::endl;
  std::cout << "The total number of events failing the perc hits used in cluster = "                                   << total_num_events_failed_perc_used_hits_in_cluster                       << "." << std::endl;
  std::cout << "The total number of events failing the perc hits used in cluster, weighted = "                         << total_num_events_failed_perc_used_hits_in_cluster_weighted              << "." << std::endl;
  
  if(_debug) std::cout << "[UBXSec] Filling tree now." << std::endl;
  _tree1->Fill();
  
  e.put(std::move(selectionResultVector));
  e.put(std::move(assnOutSelectionResultTPCObject));
  
  e.put(std::move(vertexTrackAssociations));
  e.put(std::move(vertexPFParticleAssociations));

  if(_debug) std::cout << "********** UBXSec ends" << std::endl;

  return;

} // End of the UBXSec function.



void UBXSec::endSubRun(art::SubRun& sr) {

  if (_debug) std::cout << "[UBXSec::endSubRun] Starts" << std::endl;

  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> potsum_h;

  // MC
  if (_is_mc) {
    if (_debug) std::cout << "[UBXSec::endSubRun] Getting POT for MC" << std::endl;
    if(sr.getByLabel(_potsum_producer, potsum_h)) {
      if (_debug) std::cout << "[UBXSec::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot;
    }
    else
      _sr_pot = 0.;
  }

  // Data - Use Zarko's script instead
  if (_is_data) {
 //   if (_debug) std::cout << "[UBXSec::endSubRun] Getting POT for DATA, producer " << _potsum_producer << ", instance " << _potsum_instance << std::endl;
    if (sr.getByLabel(_potsum_producer,  potsum_h)){
      if (_debug) std::cout << "[UBXSec::endSubRun] POT are valid" << std::endl;
      _sr_pot = potsum_h->totpot;
    }
    else
      _sr_pot = 0;
  }

  _sr_tree->Fill();

  if (_debug) std::cout << "[UBXSec::endSubRun] Ends" << std::endl;
}




void UBXSec::PrintMC(std::vector<art::Ptr<simb::MCTruth>> mclist) {

  std::cout << "[UBXSec] ================= MC Information ================= [UBXSec]" << std::endl;

  int iList = 0;
  std::cout << " NEUTRINO:" << std::endl;
  if (mclist[iList]->NeutrinoSet()) {
    std::cout << "\tPDG      " << mclist[iList]->GetNeutrino().Nu().PdgCode() << std::endl;
    std::cout << "\tCC/NC?   " << (mclist[iList]->GetNeutrino().CCNC() == 0 ? "CC" : "NC") << std::endl;
    std::cout << "\tMode     " << mclist[iList]->GetNeutrino().Mode() << std::endl;
    std::cout << "\tQSqr     " << mclist[iList]->GetNeutrino().QSqr() << std::endl;
    std::cout << "\tW        " << mclist[iList]->GetNeutrino().W() << std::endl;
    std::cout << "\tX        " << mclist[iList]->GetNeutrino().X() << std::endl;
    std::cout << "\tY        " << mclist[iList]->GetNeutrino().Y() << std::endl;
    std::cout << "\tHitNuc   " << mclist[iList]->GetNeutrino().HitNuc() << std::endl;
    std::cout << "\tE        " << mclist[iList]->GetNeutrino().Nu().E() << std::endl;
    std::cout << "\tVx       " << mclist[iList]->GetNeutrino().Nu().Vx() << std::endl;
    std::cout << "\tVy       " << mclist[iList]->GetNeutrino().Nu().Vy() << std::endl;
    std::cout << "\tVz       " << mclist[iList]->GetNeutrino().Nu().Vz() << std::endl;

  } else
    std::cout << "\t---No Neutrino information---" << std::endl;

  std::cout << std::endl;
  std::cout << " PRIMARIES (only with status code==1):" << std::endl;
  for (int p = 0; p < mclist[0]->NParticles(); p++) {
    const simb::MCParticle mc_par = mclist[0]->GetParticle(p);
    if (mc_par.StatusCode() != 1) continue;
    std::cout << "\tPDG           " << mc_par.PdgCode() << std::endl;
    std::cout << "\tStart process " << mc_par.Process() << std::endl;
    std::cout << "\tEnd process   " << mc_par.EndProcess() << std::endl;
    std::cout << "\tEnergy        " << mc_par.E() << std::endl;
    std::cout << "\tMomentum      " << mc_par.P() << std::endl;
    std::cout << "\tVertex        " << mc_par.Vx() << ", " << mc_par.Vy() << ", " << mc_par.Vz() << std::endl;
    std::cout << "\tStatus Code   " << mc_par.StatusCode() << std::endl << std::endl;
  }

  std::cout << "[UBXSec] ================= MC Information ================= [UBXSec]" << std::endl;
}


//_______________________________________________________________________________________
float UBXSec::GetTrackShowerScore(art::Ptr<recob::PFParticle> pfParticle, const lar_pandora::PFParticlesToMetadata pfParticleToMetadata)
{
  // Find the PFParticle in the metadata map
   auto itParticle = pfParticleToMetadata.find(pfParticle);
  // Check the PFParticle was found
  if (itParticle == pfParticleToMetadata.end())
    throw cet::exception("WorkshopTrackShowerHelperComplete") << "PFParticle has no metadata" << std::endl;
  // There should only be one metadata for each PFParticle
  if (itParticle->second.size() != 1)
    throw cet::exception("WorkshopTrackShowerHelperComplete") << "PFParticle has mutiple metadata" << std::endl;
  // The metadata vector has size one as required, so just take the first element
  const auto metadata = itParticle->second.front();
  // Now get the properties map - this is a map from a string (the property name) to a float (the property value)
  const auto propertiesMap = metadata->GetPropertiesMap();
  // Look for the track shower ID score
  const auto itScore = propertiesMap.find("TrackScore");
  // Check the track score was available
  if (itScore == propertiesMap.end())
    throw cet::exception("WorkshopTrackShowerHelperComplete") << "PFParticle has no track score" << std::endl;
  return itScore->second;
}
void UBXSec::GetFlashLocation(std::vector<double> pePerOpDet,
                              double& Ycenter,
                              double& Zcenter,
                              double& Ywidth,
                              double& Zwidth)
{

  // Reset variables
  Ycenter = Zcenter = 0.;
  Ywidth  = Zwidth  = -999.;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  for (unsigned int opdet = 0; opdet < pePerOpDet.size(); opdet++) {

    // Get physical detector location for this opChannel
    double PMTxyz[3];
    ::art::ServiceHandle<geo::Geometry> geo;
    geo->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);

    // Add up the position, weighting with PEs
    sumy    += pePerOpDet[opdet]*PMTxyz[1];
    sumy2   += pePerOpDet[opdet]*PMTxyz[1]*PMTxyz[1];
    sumz    += pePerOpDet[opdet]*PMTxyz[2];
    sumz2   += pePerOpDet[opdet]*PMTxyz[2]*PMTxyz[2];

    totalPE += pePerOpDet[opdet];
  }

  Ycenter = sumy/totalPE;
  Zcenter = sumz/totalPE;

  // This is just sqrt(<x^2> - <x>^2)
  if ( (sumy2*totalPE - sumy*sumy) > 0. )
    Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

  if ( (sumz2*totalPE - sumz*sumz) > 0. )
    Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
}
//____________________________________________________________________________________________
void UBXSec::GetTaggedPFP(art::Event const & e, std::string cosmictag_producer, double score_cut, lar_pandora::PFParticleVector & pfpTaggedOut,std::vector<int> & tagid_v){

  pfpTaggedOut.clear();
  tagid_v.clear();

  if (_debug_cr) std::cout << "Getting cosmic tags from " << cosmictag_producer << std::endl;

  // Get the CosmicTag from the ART event
  art::Handle<std::vector<anab::CosmicTag>> cosmicTagHandle;
  e.getByLabel(cosmictag_producer, cosmicTagHandle);

  if (!cosmicTagHandle.isValid() || cosmicTagHandle->empty()){
    std::cerr << "Cosmic tag " << cosmictag_producer << " is not valid or empty." << std::endl;
    std::cout <<"valid? "<<cosmicTagHandle.isValid()<<"  empty?? "<<cosmicTagHandle->empty()<<std::endl;
    return;
  }

  // Look up the associations to PFPs
  art::FindManyP<recob::PFParticle> cosmicPFPAssns(cosmicTagHandle, e, cosmictag_producer);

  if (_debug_cr) std::cout << " cosmicPFPAssns.size(): " << cosmicPFPAssns.size() << std::endl;

  // Loop over the cosmic tags
  for (unsigned int ct = 0; ct < cosmicPFPAssns.size(); ct++) {

    // Get the cosmic tag
    art::Ptr<anab::CosmicTag> cosmicTag(cosmicTagHandle, ct);
    //  if(_debug_cr) std::cout << "This cosmic tag (" << ct << ") has type: " << cosmicTag->CosmicType() << " and score: " << cosmicTag->CosmicScore() << std::endl;

    // Get the PFP associated with this CT
    std::vector<art::Ptr<recob::PFParticle>> cosmicTagToPFP_v = cosmicPFPAssns.at(cosmicTag.key());
    //if(_debug) std::cout << "Number of PFP associated with this Cosmic Tag: " << cosmicTagToPFP_v.size() << std::endl;

    if (score_cut < 0) {
      pfpTaggedOut.emplace_back(cosmicTagToPFP_v.at(0));
      tagid_v.emplace_back(cosmicTag->CosmicType());
    } else {
      if (cosmicTag->CosmicScore() > score_cut) {
        pfpTaggedOut.emplace_back(cosmicTagToPFP_v.at(0));
        tagid_v.emplace_back(cosmicTag->CosmicType());
      }
    }
  }

}
//____________________________________________________________________________________________
int UBXSec::PFPInCommon(lar_pandora::PFParticleVector first, lar_pandora::PFParticleVector second){

  int nInCommon = 0;


  for (unsigned int f = 0; f < first.size(); f++){

    for (unsigned int s = 0; s < second.size(); s++){

      if(first.at(f) == second.at(s)) {

        nInCommon++;

        //std::cout << "In common found, flash is  " << flash << " and geo is " << geo << std::endl;

      }
    }

  }
  return nInCommon;

}



DEFINE_ART_MODULE(UBXSec)
