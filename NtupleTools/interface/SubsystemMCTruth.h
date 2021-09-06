#ifndef L1TMuonSimulations_NtupleTools_SubsystemMCTruth_h
#define L1TMuonSimulations_NtupleTools_SubsystemMCTruth_h

// References:
//     SimMuon/MCTruth/interface/MuonTruth.h
//     SimMuon/MCTruth/interface/CSCHitAssociator.h
//     SimMuon/MCTruth/interface/RPCHitAssociator.h
//     SimMuon/MCTruth/interface/GEMHitAssociator.h

#include <cassert>
#include <cstdint>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"
#include "SimDataFormats/GEMDigiSimLink/interface/ME0DigiSimLink.h"
#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLink.h"
#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLinkCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"

//#include "SimDataFormats/Track/interface/SimTrack.h"
//#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

// Forward declarations
namespace l1t {
  class Phase2L1EMTFHit;
}  // namespace l1t

namespace emtf {

  namespace phase2 {

    class SubsystemMCTruth {
    public:
      // Typedefs
      // From SimMuon/MCTruth/interface/MuonTruth.h
      typedef edm::DetSetVector<StripDigiSimLink> StripDigiSimLinks;
      typedef edm::DetSetVector<StripDigiSimLink> WireDigiSimLinks;
      typedef edm::DetSet<StripDigiSimLink> LayerLinks;
      typedef edm::DetSetVector<RPCDigiSimLink> RPCDigiSimLinks;
      typedef edm::DetSet<RPCDigiSimLink> RPCLayerLinks;
      typedef edm::DetSetVector<GEMDigiSimLink> GEMDigiSimLinks;
      typedef edm::DetSet<GEMDigiSimLink> GEMLayerLinks;
      typedef edm::DetSetVector<ME0DigiSimLink> ME0DigiSimLinks;
      typedef edm::DetSet<ME0DigiSimLink> ME0LayerLinks;
      typedef DTDigiSimLinkCollection DTDigiSimLinks;

      typedef std::pair<uint32_t, EncodedEventId> SimHitIdpr;

      struct SimHitInfo {
        PSimHit pSimHit;
        LocalPoint localPosition;
        GlobalPoint globalPosition;
        int subsystem;
        int endcap;
        int station;
        int ring;
        int chamber;
        int layer;
        int sim_tp;
      };

      typedef std::vector<SimHitInfo> SimHitInfoCollection;

      typedef l1t::Phase2L1EMTFHit EMTFHit;

      explicit SubsystemMCTruth(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iConsumes);

      ~SubsystemMCTruth();

      // Set up handles from Event and EventSetup
      void initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      // Build sim track -> tracking particle links
      void build(const TrackingParticleCollection& trkPartColl);

      // Find the matching tracking particle and return its index
      std::pair<int, int> findTrackingParticle(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;

      // Find sim hits
      SimHitInfoCollection findSimHits() const;

    private:
      int findTrackingParticleCSCStrip(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
      int findTrackingParticleCSCWire(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
      int findTrackingParticleRPC(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
      int findTrackingParticleGEM(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
      int findTrackingParticleME0(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;
      int findTrackingParticleDT(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const;

      int findTrackingParticleFromMatches(const std::map<SimHitIdpr, float>& matches,
                                          const TrackingParticleCollection& trkPartColl) const;

      int findTrackingParticleFromPSimHit(const PSimHit& pSimHit) const;

      // InputTags
      const edm::InputTag cscSimHitsTag_, cscSimHitsXFTag_, cscStripSimLinksTag_, cscWireSimLinksTag_;
      const edm::InputTag rpcSimHitsTag_, rpcSimHitsXFTag_, rpcDigiSimLinksTag_;
      const edm::InputTag gemSimHitsTag_, gemSimHitsXFTag_, gemDigiSimLinksTag_, gemStripSimLinksTag_;
      const edm::InputTag me0SimHitsTag_, me0SimHitsXFTag_, me0DigiSimLinksTag_, me0StripSimLinksTag_;
      const edm::InputTag dtSimHitsTag_, dtSimHitsXFTag_, dtDigiSimLinksTag_;

      // Use PSimHits or PSimHits-with-crossing-frame
      const bool crossingFrame_;

      // PSimHits and SimLinks
      edm::EDGetTokenT<edm::PSimHitContainer> cscSimHitsToken_;
      edm::EDGetTokenT<CrossingFrame<PSimHit> > cscSimHitsXFToken_;
      edm::EDGetTokenT<StripDigiSimLinks> cscStripSimLinksToken_;
      edm::EDGetTokenT<WireDigiSimLinks> cscWireSimLinksToken_;
      //
      edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHitsToken_;
      edm::EDGetTokenT<CrossingFrame<PSimHit> > rpcSimHitsXFToken_;
      edm::EDGetTokenT<RPCDigiSimLinks> rpcDigiSimLinksToken_;
      //
      edm::EDGetTokenT<edm::PSimHitContainer> gemSimHitsToken_;
      edm::EDGetTokenT<CrossingFrame<PSimHit> > gemSimHitsXFToken_;
      edm::EDGetTokenT<GEMDigiSimLinks> gemDigiSimLinksToken_;
      edm::EDGetTokenT<StripDigiSimLinks> gemStripSimLinksToken_;
      //
      edm::EDGetTokenT<edm::PSimHitContainer> me0SimHitsToken_;
      edm::EDGetTokenT<CrossingFrame<PSimHit> > me0SimHitsXFToken_;
      edm::EDGetTokenT<ME0DigiSimLinks> me0DigiSimLinksToken_;
      edm::EDGetTokenT<StripDigiSimLinks> me0StripSimLinksToken_;
      //
      edm::EDGetTokenT<edm::PSimHitContainer> dtSimHitsToken_;
      edm::EDGetTokenT<CrossingFrame<PSimHit> > dtSimHitsXFToken_;
      edm::EDGetTokenT<DTDigiSimLinks> dtDigiSimLinksToken_;

      const edm::PSimHitContainer* cscSimHits_;
      const CrossingFrame<PSimHit>* cscSimHitsXF_;
      const StripDigiSimLinks* cscStripSimLinks_;
      const WireDigiSimLinks* cscWireSimLinks_;
      //
      const edm::PSimHitContainer* rpcSimHits_;
      const CrossingFrame<PSimHit>* rpcSimHitsXF_;
      const RPCDigiSimLinks* rpcDigiSimLinks_;
      //
      const edm::PSimHitContainer* gemSimHits_;
      const CrossingFrame<PSimHit>* gemSimHitsXF_;
      const GEMDigiSimLinks* gemDigiSimLinks_;
      const StripDigiSimLinks* gemStripSimLinks_;
      //
      const edm::PSimHitContainer* me0SimHits_;
      const CrossingFrame<PSimHit>* me0SimHitsXF_;
      const ME0DigiSimLinks* me0DigiSimLinks_;
      const StripDigiSimLinks* me0StripSimLinks_;
      //
      const edm::PSimHitContainer* dtSimHits_;
      const CrossingFrame<PSimHit>* dtSimHitsXF_;
      const DTDigiSimLinks* dtDigiSimLinks_;

      // Geometry
      edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomToken_;
      edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
      edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;
      edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
      edm::ESGetToken<ME0Geometry, MuonGeometryRecord> me0GeomToken_;

      const DTGeometry* dtGeom_;
      const CSCGeometry* cscGeom_;
      const RPCGeometry* rpcGeom_;
      const GEMGeometry* gemGeom_;
      const ME0Geometry* me0Geom_;

      std::map<SimHitIdpr, int> trackingParticleLinks_;
    };

  }  // namespace phase2

}  // namespace emtf

#endif  // L1TMuonSimulations_NtupleTools_SubsystemMCTruth_h not defined
