#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/L1TMuonPhase2/interface/Phase2L1EMTFHit.h"
#include "DataFormats/L1TMuonPhase2/interface/Phase2L1EMTFTrack.h"
#include "L1Trigger/Phase2L1EMTF/interface/Toolbox.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "L1TMuonSimulations/NtupleTools/interface/SubsystemMCTruth.h"
#include "L1TMuonSimulations/NtupleTools/interface/UnaryOperation.h"
#include "L1TMuonSimulations/NtupleTools/interface/UnaryPredicate.h"

namespace detail {
  // Like std::remove_pointer_t, but works with std::unique_ptr<T>.
  template <typename T>
  using remove_pointer_t = typename std::pointer_traits<T>::element_type;

  // If a variable is a container, clear the container. Otherwise, reset the variable to zero.
  template <typename T, typename = void>
  struct is_clearable : std::false_type {};

  template <typename T>
  struct is_clearable<T, std::void_t<decltype(std::declval<T>().clear())> > : std::true_type {};

  template <typename T>
  inline constexpr bool is_clearable_v = is_clearable<T>::value;

  template <typename T>
  typename std::enable_if_t<is_clearable_v<T> > reset_to_zero(T* t) {
    t->clear();
  }

  template <typename T>
  typename std::enable_if_t<!is_clearable_v<T> > reset_to_zero(T* t) {
    *t = T{};
  }
}  // namespace detail

// Mixin that encapsulates the TTree and its associated vectors
class TreeMixin {
protected:
  void initTree();
  void fillTree();
  void writeTree();

  // TTree
  TTree* tree;

  // clang-format off
  // Hits
  std::unique_ptr<std::vector<int16_t> >   vh_endcap;          //
  std::unique_ptr<std::vector<int16_t> >   vh_station;         //
  std::unique_ptr<std::vector<int16_t> >   vh_ring;            //
  std::unique_ptr<std::vector<int16_t> >   vh_sector;          //
  std::unique_ptr<std::vector<int16_t> >   vh_subsector;       //
  std::unique_ptr<std::vector<int16_t> >   vh_chamber;         //
  std::unique_ptr<std::vector<int16_t> >   vh_cscid;           //
  std::unique_ptr<std::vector<int16_t> >   vh_bx;              //
  std::unique_ptr<std::vector<int16_t> >   vh_subsystem;       // kDT, kCSC, kRPC, kGEM, kME0
  std::unique_ptr<std::vector<int16_t> >   vh_neighbor;        //
  std::unique_ptr<std::vector<int16_t> >   vh_strip;           //
  std::unique_ptr<std::vector<int16_t> >   vh_wire;            //
  std::unique_ptr<std::vector<int16_t> >   vh_roll;            //
  std::unique_ptr<std::vector<int16_t> >   vh_quality;         //
  std::unique_ptr<std::vector<int16_t> >   vh_pattern;         //
  std::unique_ptr<std::vector<int16_t> >   vh_bend;            //
  std::unique_ptr<std::vector<int16_t> >   vh_time;            //
  std::unique_ptr<std::vector<int16_t> >   vh_fr;              //
  std::unique_ptr<std::vector<int32_t> >   vh_emtf_phi;        //
  std::unique_ptr<std::vector<int32_t> >   vh_emtf_theta;      //
  std::unique_ptr<std::vector<float> >     vh_sim_phi;         // in degrees
  std::unique_ptr<std::vector<float> >     vh_sim_theta;       // in degrees
  std::unique_ptr<std::vector<float> >     vh_sim_r;           // in cm
  std::unique_ptr<std::vector<float> >     vh_sim_z;           // in cm
  std::unique_ptr<std::vector<int32_t> >   vh_sim_tp1;         //
  std::unique_ptr<std::vector<int32_t> >   vh_sim_tp2;         //
  std::unique_ptr<int32_t>                 vh_size;            //

  // SimHits
  std::unique_ptr<std::vector<int16_t> >   vc_subsystem;       // kDT, kCSC, kRPC, kGEM, kME0
  std::unique_ptr<std::vector<int16_t> >   vc_endcap;          //
  std::unique_ptr<std::vector<int16_t> >   vc_station;         //
  std::unique_ptr<std::vector<int16_t> >   vc_ring;            //
  std::unique_ptr<std::vector<int16_t> >   vc_chamber;         //
  std::unique_ptr<std::vector<int16_t> >   vc_layer;           //
  std::unique_ptr<std::vector<float> >     vc_phi;             // in degrees
  std::unique_ptr<std::vector<float> >     vc_theta;           // in degrees
  std::unique_ptr<std::vector<float> >     vc_perp;            // in cm
  std::unique_ptr<std::vector<float> >     vc_z;               // in cm
  std::unique_ptr<std::vector<float> >     vc_local_x;         //
  std::unique_ptr<std::vector<float> >     vc_local_y;         //
  std::unique_ptr<std::vector<int32_t> >   vc_sim_tp;          //
  std::unique_ptr<std::vector<int32_t> >   vc_pdgid;           // particleType()
  std::unique_ptr<std::vector<int16_t> >   vc_process;         // processType()
  std::unique_ptr<std::vector<float> >     vc_mom_p;           // pabs()
  std::unique_ptr<std::vector<float> >     vc_mom_phi;         // phiAtEntry()
  std::unique_ptr<std::vector<float> >     vc_mom_theta;       // thetaAtEntry()
  std::unique_ptr<std::vector<float> >     vc_tof;             // timeOfFlight()
  std::unique_ptr<int32_t>                 vc_size;            //

  // Tracks
  std::unique_ptr<std::vector<float> >     vt_pt;              //
  std::unique_ptr<std::vector<float> >     vt_xml_pt;          //
  std::unique_ptr<std::vector<float> >     vt_pt_dxy;          //
  std::unique_ptr<std::vector<float> >     vt_dxy;             //
  std::unique_ptr<std::vector<float> >     vt_invpt;           //
  std::unique_ptr<std::vector<float> >     vt_invpt_dxy;       //
  std::unique_ptr<std::vector<float> >     vt_phi;             // in degrees
  std::unique_ptr<std::vector<float> >     vt_theta;           // in degrees
  std::unique_ptr<std::vector<float> >     vt_eta;             //
  std::unique_ptr<std::vector<int16_t> >   vt_q;               // charge
  std::unique_ptr<std::vector<uint64_t> >  vt_address;         //
  std::unique_ptr<std::vector<int16_t> >   vt_mode;            //
  std::unique_ptr<std::vector<int16_t> >   vt_endcap;          //
  std::unique_ptr<std::vector<int16_t> >   vt_sector;          //
  std::unique_ptr<std::vector<int16_t> >   vt_bx;              //
  std::unique_ptr<std::vector<int16_t> >   vt_nhits;           //
  std::unique_ptr<std::vector<int32_t> >   vt_hitref1;         //
  std::unique_ptr<std::vector<int32_t> >   vt_hitref2;         //
  std::unique_ptr<std::vector<int32_t> >   vt_hitref3;         //
  std::unique_ptr<std::vector<int32_t> >   vt_hitref4;         //
  std::unique_ptr<int32_t>                 vt_size;            //

  // L1TrackTrigger tracks
  std::unique_ptr<std::vector<float> >     vd_pt;              //
  std::unique_ptr<std::vector<float> >     vd_phi;             // in radians
  std::unique_ptr<std::vector<float> >     vd_theta;           // in radians
  std::unique_ptr<std::vector<float> >     vd_eta;             //
  std::unique_ptr<std::vector<float> >     vd_vx;              // POCA, in cm
  std::unique_ptr<std::vector<float> >     vd_vy;              // POCA, in cm
  std::unique_ptr<std::vector<float> >     vd_vz;              // POCA, in cm
  std::unique_ptr<std::vector<int16_t> >   vd_q;               // charge
  std::unique_ptr<std::vector<float> >     vd_rinv;            //
  std::unique_ptr<std::vector<float> >     vd_chi2;            //
  std::unique_ptr<std::vector<int16_t> >   vd_ndof;            //
  std::unique_ptr<std::vector<int16_t> >   vd_phisector;       //
  std::unique_ptr<std::vector<int16_t> >   vd_etasector;       //
  std::unique_ptr<std::vector<int16_t> >   vd_hitpattern;      //
  std::unique_ptr<std::vector<float> >     vd_sim_pt;          //
  std::unique_ptr<std::vector<float> >     vd_sim_phi;         //
  std::unique_ptr<std::vector<float> >     vd_sim_eta;         //
  std::unique_ptr<std::vector<int32_t> >   vd_sim_tp;          //
  std::unique_ptr<std::vector<int32_t> >   vd_pdgid;           //
  std::unique_ptr<std::vector<int16_t> >   vd_genuine;         // isUnknown, isCombinatoric, isLooselyGenuine, isGenuine
  std::unique_ptr<int32_t>                 vd_size;            //

  // TrackingParticles
  std::unique_ptr<std::vector<float> >     vp_pt;              //
  std::unique_ptr<std::vector<float> >     vp_phi;             // in radians
  std::unique_ptr<std::vector<float> >     vp_theta;           // in radians
  std::unique_ptr<std::vector<float> >     vp_eta;             //
  std::unique_ptr<std::vector<float> >     vp_vx;              // in cm
  std::unique_ptr<std::vector<float> >     vp_vy;              // in cm
  std::unique_ptr<std::vector<float> >     vp_vz;              // in cm
  std::unique_ptr<std::vector<float> >     vp_invpt;           //
  std::unique_ptr<std::vector<float> >     vp_d0;              // in cm
  std::unique_ptr<std::vector<float> >     vp_beta;            //
  std::unique_ptr<std::vector<float> >     vp_mass;            //
  std::unique_ptr<std::vector<int16_t> >   vp_q;               // charge
  std::unique_ptr<std::vector<int16_t> >   vp_bx;              //
  std::unique_ptr<std::vector<int16_t> >   vp_event;           //
  std::unique_ptr<std::vector<int32_t> >   vp_pdgid;           //
  std::unique_ptr<std::vector<int16_t> >   vp_status;          //
  std::unique_ptr<std::vector<int16_t> >   vp_decay;           //
  std::unique_ptr<std::vector<int32_t> >   vp_genp;            //
  std::unique_ptr<int32_t>                 vp_size;            //

  // PileupSummaryInfo
  std::unique_ptr<std::vector<uint64_t> >  ve_event;           //
  std::unique_ptr<std::vector<uint32_t> >  ve_run;             //
  std::unique_ptr<std::vector<uint32_t> >  ve_lumi;            //
  std::unique_ptr<std::vector<float> >     ve_npv;             // getTrueNumInteractions()
  std::unique_ptr<std::vector<int16_t> >   ve_nvertices;       // getPU_NumInteractions()
  std::unique_ptr<int32_t>                 ve_size;            //
  // clang-format on
};

// _____________________________________________________________________________
void TreeMixin::initTree() {
  // Create tree using TFileService
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "tree");

  // Define a macro to call std::make_unique<T> with the right type, then call tree->Branch()
  // by providing the pointer.
  // E.g.
  //     vh_endcap = std::make_unique<std::vector<int16_t> >();
  //     tree->Branch("vh_endcap", vh_endcap.get());
#define MY_MAKE_UNIQUE_AND_BRANCH(NAME)                                   \
  (NAME) = std::make_unique<detail::remove_pointer_t<decltype(NAME)> >(); \
  tree->Branch(#NAME, (NAME).get());

  // Hits
  MY_MAKE_UNIQUE_AND_BRANCH(vh_endcap)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_station)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_ring)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_sector)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_subsector)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_chamber)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_cscid)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_bx)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_subsystem)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_neighbor)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_strip)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_wire)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_roll)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_quality)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_pattern)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_bend)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_time)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_fr)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_emtf_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_emtf_theta)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_sim_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_sim_theta)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_sim_r)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_sim_z)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_sim_tp1)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_sim_tp2)
  MY_MAKE_UNIQUE_AND_BRANCH(vh_size)

  // SimHits
  MY_MAKE_UNIQUE_AND_BRANCH(vc_subsystem)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_endcap)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_station)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_ring)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_chamber)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_layer)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_theta)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_perp)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_z)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_local_x)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_local_y)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_sim_tp)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_pdgid)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_process)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_mom_p)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_mom_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_mom_theta)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_tof)
  MY_MAKE_UNIQUE_AND_BRANCH(vc_size)

  // Tracks
  MY_MAKE_UNIQUE_AND_BRANCH(vt_pt)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_xml_pt)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_pt_dxy)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_dxy)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_invpt)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_invpt_dxy)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_theta)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_eta)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_q)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_address)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_mode)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_endcap)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_sector)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_bx)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_nhits)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_hitref1)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_hitref2)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_hitref3)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_hitref4)
  MY_MAKE_UNIQUE_AND_BRANCH(vt_size)

  // L1TrackTrigger tracks
  MY_MAKE_UNIQUE_AND_BRANCH(vd_pt)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_theta)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_eta)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_vx)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_vy)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_vz)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_q)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_rinv)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_chi2)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_ndof)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_phisector)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_etasector)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_hitpattern)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_sim_pt)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_sim_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_sim_eta)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_sim_tp)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_pdgid)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_genuine)
  MY_MAKE_UNIQUE_AND_BRANCH(vd_size)

  // TrackingParticles
  MY_MAKE_UNIQUE_AND_BRANCH(vp_pt)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_phi)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_theta)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_eta)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_vx)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_vy)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_vz)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_invpt)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_d0)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_beta)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_mass)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_q)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_bx)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_event)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_pdgid)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_status)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_decay)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_genp)
  MY_MAKE_UNIQUE_AND_BRANCH(vp_size)

  // PileupSummaryInfo
  MY_MAKE_UNIQUE_AND_BRANCH(ve_event)
  MY_MAKE_UNIQUE_AND_BRANCH(ve_run)
  MY_MAKE_UNIQUE_AND_BRANCH(ve_lumi)
  MY_MAKE_UNIQUE_AND_BRANCH(ve_npv)
  MY_MAKE_UNIQUE_AND_BRANCH(ve_nvertices)
  MY_MAKE_UNIQUE_AND_BRANCH(ve_size)

#undef MY_MAKE_UNIQUE_AND_BRANCH
}

void TreeMixin::fillTree() {
  // Fill tree
  tree->Fill();

  // Define a macro to either clear a container or reset a variable to zero.
  // E.g.
  //     vh_endcap->clear();
  // Or
  //     (*vh_size) = 0;
#define MY_RESET_TO_ZERO(NAME) detail::reset_to_zero((NAME).get());

  // Hits
  MY_RESET_TO_ZERO(vh_endcap)
  MY_RESET_TO_ZERO(vh_station)
  MY_RESET_TO_ZERO(vh_ring)
  MY_RESET_TO_ZERO(vh_sector)
  MY_RESET_TO_ZERO(vh_subsector)
  MY_RESET_TO_ZERO(vh_chamber)
  MY_RESET_TO_ZERO(vh_cscid)
  MY_RESET_TO_ZERO(vh_bx)
  MY_RESET_TO_ZERO(vh_subsystem)
  MY_RESET_TO_ZERO(vh_neighbor)
  MY_RESET_TO_ZERO(vh_strip)
  MY_RESET_TO_ZERO(vh_wire)
  MY_RESET_TO_ZERO(vh_roll)
  MY_RESET_TO_ZERO(vh_quality)
  MY_RESET_TO_ZERO(vh_pattern)
  MY_RESET_TO_ZERO(vh_bend)
  MY_RESET_TO_ZERO(vh_time)
  MY_RESET_TO_ZERO(vh_fr)
  MY_RESET_TO_ZERO(vh_emtf_phi)
  MY_RESET_TO_ZERO(vh_emtf_theta)
  MY_RESET_TO_ZERO(vh_sim_phi)
  MY_RESET_TO_ZERO(vh_sim_theta)
  MY_RESET_TO_ZERO(vh_sim_r)
  MY_RESET_TO_ZERO(vh_sim_z)
  MY_RESET_TO_ZERO(vh_sim_tp1)
  MY_RESET_TO_ZERO(vh_sim_tp2)
  MY_RESET_TO_ZERO(vh_size)

  // SimHits
  MY_RESET_TO_ZERO(vc_subsystem)
  MY_RESET_TO_ZERO(vc_endcap)
  MY_RESET_TO_ZERO(vc_station)
  MY_RESET_TO_ZERO(vc_ring)
  MY_RESET_TO_ZERO(vc_chamber)
  MY_RESET_TO_ZERO(vc_layer)
  MY_RESET_TO_ZERO(vc_phi)
  MY_RESET_TO_ZERO(vc_theta)
  MY_RESET_TO_ZERO(vc_perp)
  MY_RESET_TO_ZERO(vc_z)
  MY_RESET_TO_ZERO(vc_local_x)
  MY_RESET_TO_ZERO(vc_local_y)
  MY_RESET_TO_ZERO(vc_sim_tp)
  MY_RESET_TO_ZERO(vc_pdgid)
  MY_RESET_TO_ZERO(vc_process)
  MY_RESET_TO_ZERO(vc_mom_p)
  MY_RESET_TO_ZERO(vc_mom_phi)
  MY_RESET_TO_ZERO(vc_mom_theta)
  MY_RESET_TO_ZERO(vc_tof)
  MY_RESET_TO_ZERO(vc_size)

  // Tracks
  MY_RESET_TO_ZERO(vt_pt)
  MY_RESET_TO_ZERO(vt_xml_pt)
  MY_RESET_TO_ZERO(vt_pt_dxy)
  MY_RESET_TO_ZERO(vt_dxy)
  MY_RESET_TO_ZERO(vt_invpt)
  MY_RESET_TO_ZERO(vt_invpt_dxy)
  MY_RESET_TO_ZERO(vt_phi)
  MY_RESET_TO_ZERO(vt_theta)
  MY_RESET_TO_ZERO(vt_eta)
  MY_RESET_TO_ZERO(vt_q)
  MY_RESET_TO_ZERO(vt_address)
  MY_RESET_TO_ZERO(vt_mode)
  MY_RESET_TO_ZERO(vt_endcap)
  MY_RESET_TO_ZERO(vt_sector)
  MY_RESET_TO_ZERO(vt_bx)
  MY_RESET_TO_ZERO(vt_nhits)
  MY_RESET_TO_ZERO(vt_hitref1)
  MY_RESET_TO_ZERO(vt_hitref2)
  MY_RESET_TO_ZERO(vt_hitref3)
  MY_RESET_TO_ZERO(vt_hitref4)
  MY_RESET_TO_ZERO(vt_size)

  // L1TrackTrigger tracks
  MY_RESET_TO_ZERO(vd_pt)
  MY_RESET_TO_ZERO(vd_phi)
  MY_RESET_TO_ZERO(vd_theta)
  MY_RESET_TO_ZERO(vd_eta)
  MY_RESET_TO_ZERO(vd_vx)
  MY_RESET_TO_ZERO(vd_vy)
  MY_RESET_TO_ZERO(vd_vz)
  MY_RESET_TO_ZERO(vd_q)
  MY_RESET_TO_ZERO(vd_rinv)
  MY_RESET_TO_ZERO(vd_chi2)
  MY_RESET_TO_ZERO(vd_ndof)
  MY_RESET_TO_ZERO(vd_phisector)
  MY_RESET_TO_ZERO(vd_etasector)
  MY_RESET_TO_ZERO(vd_hitpattern)
  MY_RESET_TO_ZERO(vd_sim_pt)
  MY_RESET_TO_ZERO(vd_sim_phi)
  MY_RESET_TO_ZERO(vd_sim_eta)
  MY_RESET_TO_ZERO(vd_sim_tp)
  MY_RESET_TO_ZERO(vd_pdgid)
  MY_RESET_TO_ZERO(vd_genuine)
  MY_RESET_TO_ZERO(vd_size)

  // TrackingParticles
  MY_RESET_TO_ZERO(vp_pt)
  MY_RESET_TO_ZERO(vp_phi)
  MY_RESET_TO_ZERO(vp_theta)
  MY_RESET_TO_ZERO(vp_eta)
  MY_RESET_TO_ZERO(vp_vx)
  MY_RESET_TO_ZERO(vp_vy)
  MY_RESET_TO_ZERO(vp_vz)
  MY_RESET_TO_ZERO(vp_invpt)
  MY_RESET_TO_ZERO(vp_d0)
  MY_RESET_TO_ZERO(vp_beta)
  MY_RESET_TO_ZERO(vp_mass)
  MY_RESET_TO_ZERO(vp_q)
  MY_RESET_TO_ZERO(vp_bx)
  MY_RESET_TO_ZERO(vp_event)
  MY_RESET_TO_ZERO(vp_pdgid)
  MY_RESET_TO_ZERO(vp_status)
  MY_RESET_TO_ZERO(vp_decay)
  MY_RESET_TO_ZERO(vp_genp)
  MY_RESET_TO_ZERO(vp_size)

  // PileupSummaryInfo
  MY_RESET_TO_ZERO(ve_event)
  MY_RESET_TO_ZERO(ve_run)
  MY_RESET_TO_ZERO(ve_lumi)
  MY_RESET_TO_ZERO(ve_npv)
  MY_RESET_TO_ZERO(ve_nvertices)
  MY_RESET_TO_ZERO(ve_size)

#undef MY_RESET_TO_ZERO
}

void TreeMixin::writeTree() {
  // Writing is handled by TFileService
}

// _____________________________________________________________________________
// EDAnalyzer
class NtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>, public TreeMixin {
public:
  explicit NtupleMaker(const edm::ParameterSet& iConfig);
  ~NtupleMaker() override;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void endRun(const edm::Run&, const edm::EventSetup&) override;

  // Main function
  void process(const edm::Event&, const edm::EventSetup&);

  // Aux functions
  void getHandles(const edm::Event&, const edm::EventSetup&);

  // Typedefs
  typedef l1t::Phase2L1EMTFHit EMTFHit;
  typedef l1t::Phase2L1EMTFHitCollection EMTFHitCollection;
  typedef l1t::Phase2L1EMTFTrack EMTFTrack;
  typedef l1t::Phase2L1EMTFTrackCollection EMTFTrackCollection;

  typedef TTTrack<Ref_Phase2TrackerDigi_> TTTFTrack;
  typedef std::vector<TTTFTrack> TTTFTrackCollection;
  typedef TTTrackAssociationMap<Ref_Phase2TrackerDigi_> TTTFTrackAssociator;

  // Helper objects
  std::unique_ptr<emtf::phase2::SubsystemMCTruth> truth_;

  // InputTags
  const edm::InputTag genPartTag_, simTrackTag_, trkPartTag_, pileupInfoTag_;
  const edm::InputTag emtfHitTag_, emtfTrackTag_, tttfTrackTag_, tttfTrackAssocTag_;

  // Output filename
  const std::string fileName_;

  // Verbosity level
  int verbose_;

  // Tokens
  edm::EDGetTokenT<reco::GenParticleCollection> genPartToken_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trkPartToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
  edm::EDGetTokenT<EMTFHitCollection> emtfHitToken_;
  edm::EDGetTokenT<EMTFTrackCollection> emtfTrackToken_;
  edm::EDGetTokenT<TTTFTrackCollection> tttfTrackToken_;
  edm::EDGetTokenT<TTTFTrackAssociator> tttfTrackAssocToken_;

  // Products
  const reco::GenParticleCollection* genParts_;
  const edm::SimTrackContainer* simTracks_;
  const TrackingParticleCollection* trkParts_;
  const std::vector<PileupSummaryInfo>* pileupInfo_;
  const EMTFHitCollection* emtfHits_;
  const EMTFTrackCollection* emtfTracks_;
  const TTTFTrackCollection* tttfTracks_;
  const TTTFTrackAssociator* tttfTrackAssoc_;
};

// _____________________________________________________________________________
NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig)
    : truth_(std::make_unique<emtf::phase2::SubsystemMCTruth>(iConfig, consumesCollector())),
      genPartTag_(iConfig.getParameter<edm::InputTag>("genPartTag")),
      simTrackTag_(iConfig.getParameter<edm::InputTag>("simTrackTag")),
      trkPartTag_(iConfig.getParameter<edm::InputTag>("trkPartTag")),
      pileupInfoTag_(iConfig.getParameter<edm::InputTag>("pileupInfoTag")),
      emtfHitTag_(iConfig.getParameter<edm::InputTag>("emtfHitTag")),
      emtfTrackTag_(iConfig.getParameter<edm::InputTag>("emtfTrackTag")),
      tttfTrackTag_(iConfig.getParameter<edm::InputTag>("tttfTrackTag")),
      tttfTrackAssocTag_(iConfig.getParameter<edm::InputTag>("tttfTrackAssocTag")),
      fileName_(iConfig.getParameter<std::string>("fileName")),
      verbose_(iConfig.getUntrackedParameter<int>("verbosity", 0)) {
  // SharedResources
  usesResource("TFileService");

  genPartToken_ = consumes<reco::GenParticleCollection>(genPartTag_);
  simTrackToken_ = consumes<edm::SimTrackContainer>(simTrackTag_);
  trkPartToken_ = consumes<TrackingParticleCollection>(trkPartTag_);
  pileupInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(pileupInfoTag_);
  emtfHitToken_ = consumes<EMTFHitCollection>(emtfHitTag_);
  emtfTrackToken_ = consumes<EMTFTrackCollection>(emtfTrackTag_);
  tttfTrackToken_ = consumes<TTTFTrackCollection>(tttfTrackTag_);
  tttfTrackAssocToken_ = consumes<TTTFTrackAssociator>(tttfTrackAssocTag_);

  genParts_ = nullptr;
  simTracks_ = nullptr;
  trkParts_ = nullptr;
  pileupInfo_ = nullptr;
  emtfHits_ = nullptr;
  emtfTracks_ = nullptr;
  tttfTracks_ = nullptr;
  tttfTrackAssoc_ = nullptr;
}

NtupleMaker::~NtupleMaker() {}

void NtupleMaker::beginJob() { initTree(); }

void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  getHandles(iEvent, iSetup);
  process(iEvent, iSetup);
}

void NtupleMaker::endJob() { writeTree(); }

void NtupleMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

void NtupleMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

// _____________________________________________________________________________
void NtupleMaker::getHandles(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // MC truth for EMTFHit
  truth_->initEvent(iEvent, iSetup);

  if (!iEvent.isRealData()) {
    auto genParts_handle = iEvent.getHandle(genPartToken_);
    auto simTracks_handle = iEvent.getHandle(simTrackToken_);
    auto trkParts_handle = iEvent.getHandle(trkPartToken_);
    auto pileupInfo_handle = iEvent.getHandle(pileupInfoToken_);

    genParts_ = genParts_handle.product();
    simTracks_ = simTracks_handle.product();
    trkParts_ = trkParts_handle.product();
    pileupInfo_ = pileupInfo_handle.product();
  }

  auto emtfHits_handle = iEvent.getHandle(emtfHitToken_);
  auto emtfTracks_handle = iEvent.getHandle(emtfTrackToken_);
  auto tttfTracks_handle = iEvent.getHandle(tttfTrackToken_);
  auto tttfTrackAssoc_handle = iEvent.getHandle(tttfTrackAssocToken_);

  emtfHits_ = emtfHits_handle.product();
  emtfTracks_ = emtfTracks_handle.product();
  tttfTracks_ = tttfTracks_handle.product();
  tttfTrackAssoc_ = tttfTrackAssoc_handle.product();
}

void NtupleMaker::process(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Preprocessing
  EMTFHitCollection emtf_hits;
  std::copy_if(emtfHits_->begin(), emtfHits_->end(), std::back_inserter(emtf_hits), emtf::phase2::EMTFHitPred{});

  EMTFTrackCollection emtf_tracks;
  std::copy_if(
      emtfTracks_->begin(), emtfTracks_->end(), std::back_inserter(emtf_tracks), emtf::phase2::EMTFTrackPred{});

  TrackingParticleCollection trk_particles;
  std::copy_if(
      trkParts_->begin(), trkParts_->end(), std::back_inserter(trk_particles), emtf::phase2::TrackingParticlePred{});

  // MC truth for EMTFHit
  truth_->buildTrackingParticleLinks(trk_particles);

  using SimHitInfoCollection = emtf::phase2::SubsystemMCTruth::SimHitInfoCollection;
  const SimHitInfoCollection& sim_hits = truth_->findSimHits();

  using emtf::phase2::toolbox::deg_to_rad;
  using emtf::phase2::toolbox::rad_to_deg;

  // Verbose
  if (verbose_ > 0) {
    std::cout << "[DEBUG] # hits: " << emtf_hits.size() << " # tracks: " << emtf_tracks.size()
              << " # particles: " << trk_particles.size() << std::endl;
  }

  // ___________________________________________________________________________
  // Hits
  for (const auto& hit : emtf_hits) {
    //TODO: implement this
  }
  (*vh_size) = emtf_hits.size();

  // ___________________________________________________________________________
  // SimHits
  for (const auto& hit : sim_hits) {
    vc_subsystem->push_back(hit.subsystem);
    vc_endcap->push_back(hit.endcap);
    vc_station->push_back(hit.station);
    vc_ring->push_back(hit.ring);
    vc_chamber->push_back(hit.chamber);
    vc_layer->push_back(hit.layer);
    vc_phi->push_back(rad_to_deg(hit.globalPosition.phi().value()));
    vc_theta->push_back(rad_to_deg(hit.globalPosition.theta().value()));
    vc_perp->push_back(hit.globalPosition.perp());
    vc_z->push_back(hit.globalPosition.z());
    vc_local_x->push_back(hit.localPosition.x());
    vc_local_y->push_back(hit.localPosition.y());
    vc_sim_tp->push_back(hit.sim_tp);
    vc_pdgid->push_back(hit.pSimHit.particleType());
    vc_process->push_back(hit.pSimHit.processType());
    vc_mom_phi->push_back(hit.pSimHit.pabs());
    vc_mom_phi->push_back(rad_to_deg(hit.pSimHit.phiAtEntry().value()));
    vc_mom_theta->push_back(rad_to_deg(hit.pSimHit.thetaAtEntry().value()));
    vc_tof->push_back(hit.pSimHit.timeOfFlight());
  }
  (*vc_size) = sim_hits.size();

  // ___________________________________________________________________________
  // Tracks
  for (const auto& trk : emtf_tracks) {
    //TODO: implement this
  }
  (*vt_size) = emtf_tracks.size();

  // ___________________________________________________________________________
  // L1TrackTrigger tracks
  if ((tttfTracks_ != nullptr) and (tttfTrackAssoc_ != nullptr)) {
    int itrk = 0;

    auto tttfTracks_handle = iEvent.getHandle(tttfTrackToken_);

    for (const auto& trk : *tttfTracks_) {
      edm::Ptr<TTTFTrack> trkPtr(tttfTracks_handle, itrk++);
      edm::Ptr<TrackingParticle> tpPtr = tttfTrackAssoc_->findTrackingParticlePtr(trkPtr);

      float sim_pt = 0.;
      float sim_phi = 0.;
      float sim_eta = 0.;
      int sim_tp = -1;
      int sim_pdgid = 0;
      int sim_genuine = 0;
      if (tpPtr.isNonnull()) {
        sim_pt = tpPtr->pt();
        sim_eta = tpPtr->eta();
        sim_phi = tpPtr->phi();
        sim_tp = tpPtr.key();
        sim_pdgid = tpPtr->pdgId();
      }
      if (tttfTrackAssoc_->isGenuine(trkPtr)) {
        sim_genuine = 3;
      } else if (tttfTrackAssoc_->isLooselyGenuine(trkPtr)) {
        sim_genuine = 2;
      } else if (tttfTrackAssoc_->isCombinatoric(trkPtr)) {
        sim_genuine = 1;
      } else if (tttfTrackAssoc_->isUnknown(trkPtr)) {
        sim_genuine = 0;
      }

      const GlobalVector& momentum = trk.momentum();
      const GlobalPoint& poca = trk.POCA();
      const double rinv = trk.rInv();

      vd_pt->push_back(momentum.perp());
      vd_phi->push_back(momentum.phi());
      vd_theta->push_back(momentum.theta());
      vd_eta->push_back(momentum.eta());
      vd_vx->push_back(poca.x());
      vd_vy->push_back(poca.y());
      vd_vz->push_back(poca.z());
      vd_q->push_back(rinv >= 0 ? 1 : -1);
      vd_rinv->push_back(rinv);
      vd_chi2->push_back(trk.chi2());
      vd_ndof->push_back(trk.getStubRefs().size() * 2 - trk.nFitPars());
      vd_phisector->push_back(trk.phiSector());
      vd_etasector->push_back(trk.etaSector());
      vd_hitpattern->push_back(trk.hitPattern());
      vd_sim_pt->push_back(sim_pt);
      vd_sim_phi->push_back(sim_phi);
      vd_sim_eta->push_back(sim_eta);
      vd_sim_tp->push_back(sim_tp);
      vd_pdgid->push_back(sim_pdgid);
      vd_genuine->push_back(sim_genuine);
    }
    (*vd_size) = tttfTracks_->size();
  }

  // ___________________________________________________________________________
  // TrackingParticles
  for (const auto& part : trk_particles) {
    int igenPart = -1;
    if (!part.genParticles().empty()) {
      igenPart = (part.genParticles().begin())->key();
    }

    auto calc_invpt_fn = emtf::phase2::InvptOp{};
    auto calc_d0_fn = emtf::phase2::DzeroOp{};

    vp_pt->push_back(part.pt());
    vp_phi->push_back(part.phi());
    vp_theta->push_back(part.theta());
    vp_eta->push_back(part.eta());
    vp_vx->push_back(part.vx());
    vp_vy->push_back(part.vy());
    vp_vz->push_back(part.vz());
    vp_invpt->push_back(calc_invpt_fn(part.charge(), part.pt()));
    vp_d0->push_back(calc_d0_fn(calc_invpt_fn(part.charge(), part.pt()), part.phi(), part.vx(), part.vy()));
    vp_beta->push_back(part.p() / part.energy());
    vp_mass->push_back(part.mass());
    vp_q->push_back(part.charge());
    vp_bx->push_back(part.eventId().bunchCrossing());
    vp_event->push_back(part.eventId().event());
    vp_pdgid->push_back(part.pdgId());
    vp_status->push_back(part.status());
    vp_decay->push_back(part.decayVertices().size());
    vp_genp->push_back(igenPart);
  }
  (*vp_size) = trk_particles.size();

  // ___________________________________________________________________________
  // PileupSummaryInfo
  if (pileupInfo_ != nullptr) {
    float npv = 0.;
    int nvertices = 0;
    for (const auto& pui : *pileupInfo_) {
      if (pui.getBunchCrossing() == 0) {  // BX=0
        npv = pui.getTrueNumInteractions();
        nvertices = pui.getPU_NumInteractions();
        break;
      }
    }

    {
      ve_event->push_back(iEvent.id().event());
      ve_run->push_back(iEvent.id().run());
      ve_lumi->push_back(iEvent.id().luminosityBlock());
      ve_npv->push_back(npv);
      ve_nvertices->push_back(nvertices);
    }
    (*ve_size) = 1;
  }

  // ___________________________________________________________________________
  // Fill
  fillTree();
}

// _____________________________________________________________________________
void NtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  {
    // used by class SubsystemMCTruth
    desc.add<edm::InputTag>("cscSimHitsTag", edm::InputTag("g4SimHits", "MuonCSCHits"));
    desc.add<edm::InputTag>("cscSimHitsXFTag", edm::InputTag("mix", "g4SimHitsMuonCSCHits"));
    desc.add<edm::InputTag>("cscStripSimLinksTag", edm::InputTag("simMuonCSCDigis", "MuonCSCStripDigiSimLinks"));
    desc.add<edm::InputTag>("cscWireSimLinksTag", edm::InputTag("simMuonCSCDigis", "MuonCSCWireDigiSimLinks"));
    desc.add<edm::InputTag>("rpcSimHitsTag", edm::InputTag("g4SimHits", "MuonRPCHits"));
    desc.add<edm::InputTag>("rpcSimHitsXFTag", edm::InputTag("mix", "g4SimHitsMuonRPCHits"));
    desc.add<edm::InputTag>("rpcDigiSimLinksTag", edm::InputTag("simMuonRPCDigis", "RPCDigiSimLink"));
    desc.add<edm::InputTag>("gemSimHitsTag", edm::InputTag("g4SimHits", "MuonGEMHits"));
    desc.add<edm::InputTag>("gemSimHitsXFTag", edm::InputTag("mix", "g4SimHitsMuonGEMHits"));
    desc.add<edm::InputTag>("gemDigiSimLinksTag", edm::InputTag("simMuonGEMDigis", "GEM"));
    desc.add<edm::InputTag>("gemStripSimLinksTag", edm::InputTag("simMuonGEMDigis", "GEM"));
    desc.add<edm::InputTag>("me0SimHitsTag", edm::InputTag("g4SimHits", "MuonME0Hits"));
    desc.add<edm::InputTag>("me0SimHitsXFTag", edm::InputTag("mix", "g4SimHitsMuonME0Hits"));
    desc.add<edm::InputTag>("me0DigiSimLinksTag", edm::InputTag("simMuonME0Digis", "ME0"));
    desc.add<edm::InputTag>("me0StripSimLinksTag", edm::InputTag("simMuonME0Digis", "ME0"));
    desc.add<edm::InputTag>("dtSimHitsTag", edm::InputTag("g4SimHits", "MuonDTHits"));
    desc.add<edm::InputTag>("dtSimHitsXFTag", edm::InputTag("mix", "g4SimHitsMuonDTHits"));
    desc.add<edm::InputTag>("dtDigiSimLinksTag", edm::InputTag("simMuonDTDigis"));
    desc.add<bool>("crossingFrame", false);
  }
  desc.add<edm::InputTag>("genPartTag", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("simTrackTag", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("trkPartTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("pileupInfoTag", edm::InputTag("addPileupInfo"));
  desc.add<edm::InputTag>("emtfHitTag", edm::InputTag("phase2L1EMTFProducer"));
  desc.add<edm::InputTag>("emtfTrackTag", edm::InputTag("phase2L1EMTFProducer"));
  desc.add<edm::InputTag>("tttfTrackTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<edm::InputTag>("tttfTrackAssocTag", edm::InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"));
  desc.add<std::string>("fileName", "ntuple.root");
  desc.addUntracked<int>("verbosity", 0);
  descriptions.add("ntupler", desc);

  //edm::ParameterSetDescription default_desc;
  //default_desc.setUnknown();
  //descriptions.addDefault(default_desc);
}

// define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(NtupleMaker);
