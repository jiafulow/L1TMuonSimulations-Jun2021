#include "L1TMuonSimulations/NtupleTools/interface/SubsystemMCTruth.h"

#include <algorithm>
#include <iostream>

#include "DataFormats/MuonDetId/interface/DTWireId.h"

#include "DataFormats/L1TMuonPhase2/interface/Phase2L1EMTFHit.h"

using namespace emtf::phase2;

SubsystemMCTruth::SubsystemMCTruth(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iConsumes)
    : cscSimHitsTag_(iConfig.getParameter<edm::InputTag>("cscSimHitsTag")),
      cscSimHitsXFTag_(iConfig.getParameter<edm::InputTag>("cscSimHitsXFTag")),
      cscStripSimLinksTag_(iConfig.getParameter<edm::InputTag>("cscStripSimLinksTag")),
      cscWireSimLinksTag_(iConfig.getParameter<edm::InputTag>("cscWireSimLinksTag")),
      rpcSimHitsTag_(iConfig.getParameter<edm::InputTag>("rpcSimHitsTag")),
      rpcSimHitsXFTag_(iConfig.getParameter<edm::InputTag>("rpcSimHitsXFTag")),
      rpcDigiSimLinksTag_(iConfig.getParameter<edm::InputTag>("rpcDigiSimLinksTag")),
      gemSimHitsTag_(iConfig.getParameter<edm::InputTag>("gemSimHitsTag")),
      gemSimHitsXFTag_(iConfig.getParameter<edm::InputTag>("gemSimHitsXFTag")),
      gemDigiSimLinksTag_(iConfig.getParameter<edm::InputTag>("gemDigiSimLinksTag")),
      gemStripSimLinksTag_(iConfig.getParameter<edm::InputTag>("gemStripSimLinksTag")),
      me0SimHitsTag_(iConfig.getParameter<edm::InputTag>("me0SimHitsTag")),
      me0SimHitsXFTag_(iConfig.getParameter<edm::InputTag>("me0SimHitsXFTag")),
      me0DigiSimLinksTag_(iConfig.getParameter<edm::InputTag>("me0DigiSimLinksTag")),
      me0StripSimLinksTag_(iConfig.getParameter<edm::InputTag>("me0StripSimLinksTag")),
      dtSimHitsTag_(iConfig.getParameter<edm::InputTag>("dtSimHitsTag")),
      dtSimHitsXFTag_(iConfig.getParameter<edm::InputTag>("dtSimHitsXFTag")),
      dtDigiSimLinksTag_(iConfig.getParameter<edm::InputTag>("dtDigiSimLinksTag")),
      crossingFrame_(iConfig.getParameter<bool>("crossingFrame")) {
  // PSimHits and SimLinks
  if (!crossingFrame_) {
    cscSimHitsToken_ = iConsumes.consumes<edm::PSimHitContainer>(cscSimHitsTag_);
  } else {
    cscSimHitsXFToken_ = iConsumes.consumes<CrossingFrame<PSimHit> >(cscSimHitsXFTag_);
  }
  cscStripSimLinksToken_ = iConsumes.consumes<StripDigiSimLinks>(cscStripSimLinksTag_);
  cscWireSimLinksToken_ = iConsumes.consumes<WireDigiSimLinks>(cscWireSimLinksTag_);
  //
  if (!crossingFrame_) {
    rpcSimHitsToken_ = iConsumes.consumes<edm::PSimHitContainer>(rpcSimHitsTag_);
  } else {
    rpcSimHitsXFToken_ = iConsumes.consumes<CrossingFrame<PSimHit> >(rpcSimHitsXFTag_);
  }
  rpcDigiSimLinksToken_ = iConsumes.consumes<RPCDigiSimLinks>(rpcDigiSimLinksTag_);
  //
  if (!crossingFrame_) {
    gemSimHitsToken_ = iConsumes.consumes<edm::PSimHitContainer>(gemSimHitsTag_);
  } else {
    gemSimHitsXFToken_ = iConsumes.consumes<CrossingFrame<PSimHit> >(gemSimHitsXFTag_);
  }
  gemDigiSimLinksToken_ = iConsumes.consumes<GEMDigiSimLinks>(gemDigiSimLinksTag_);
  gemStripSimLinksToken_ = iConsumes.consumes<StripDigiSimLinks>(gemStripSimLinksTag_);
  //
  if (!crossingFrame_) {
    me0SimHitsToken_ = iConsumes.consumes<edm::PSimHitContainer>(me0SimHitsTag_);
  } else {
    me0SimHitsXFToken_ = iConsumes.consumes<CrossingFrame<PSimHit> >(me0SimHitsXFTag_);
  }
  me0DigiSimLinksToken_ = iConsumes.consumes<ME0DigiSimLinks>(me0DigiSimLinksTag_);
  me0StripSimLinksToken_ = iConsumes.consumes<StripDigiSimLinks>(me0StripSimLinksTag_);
  //
  if (!crossingFrame_) {
    dtSimHitsToken_ = iConsumes.consumes<edm::PSimHitContainer>(dtSimHitsTag_);
  } else {
    dtSimHitsXFToken_ = iConsumes.consumes<CrossingFrame<PSimHit> >(dtSimHitsXFTag_);
  }
  dtDigiSimLinksToken_ = iConsumes.consumes<DTDigiSimLinks>(dtDigiSimLinksTag_);

  cscSimHits_ = nullptr;
  cscSimHitsXF_ = nullptr;
  cscStripSimLinks_ = nullptr;
  cscWireSimLinks_ = nullptr;
  //
  rpcSimHits_ = nullptr;
  rpcSimHitsXF_ = nullptr;
  rpcDigiSimLinks_ = nullptr;
  //
  gemSimHits_ = nullptr;
  gemSimHitsXF_ = nullptr;
  gemDigiSimLinks_ = nullptr;
  gemStripSimLinks_ = nullptr;
  //
  me0SimHits_ = nullptr;
  me0SimHitsXF_ = nullptr;
  me0DigiSimLinks_ = nullptr;
  me0StripSimLinks_ = nullptr;
  //
  dtSimHits_ = nullptr;
  dtSimHitsXF_ = nullptr;
  dtDigiSimLinks_ = nullptr;

  // Geometry
  dtGeomToken_ = iConsumes.esConsumes<DTGeometry, MuonGeometryRecord>();
  cscGeomToken_ = iConsumes.esConsumes<CSCGeometry, MuonGeometryRecord>();
  rpcGeomToken_ = iConsumes.esConsumes<RPCGeometry, MuonGeometryRecord>();
  gemGeomToken_ = iConsumes.esConsumes<GEMGeometry, MuonGeometryRecord>();
  me0GeomToken_ = iConsumes.esConsumes<ME0Geometry, MuonGeometryRecord>();

  dtGeom_ = nullptr;
  cscGeom_ = nullptr;
  rpcGeom_ = nullptr;
  gemGeom_ = nullptr;
  me0Geom_ = nullptr;
}

SubsystemMCTruth::~SubsystemMCTruth() {}

void SubsystemMCTruth::initEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // PSimHits
  if (!iEvent.isRealData()) {
    if (!crossingFrame_) {
      auto cscSimHits_handle = iEvent.getHandle(cscSimHitsToken_);
      auto rpcSimHits_handle = iEvent.getHandle(rpcSimHitsToken_);
      auto gemSimHits_handle = iEvent.getHandle(gemSimHitsToken_);
      auto me0SimHits_handle = iEvent.getHandle(me0SimHitsToken_);
      auto dtSimHits_handle = iEvent.getHandle(dtSimHitsToken_);

      cscSimHits_ = cscSimHits_handle.product();
      rpcSimHits_ = rpcSimHits_handle.product();
      gemSimHits_ = gemSimHits_handle.product();
      me0SimHits_ = me0SimHits_handle.product();
      dtSimHits_ = dtSimHits_handle.product();
    }
  }

  // PSimHits (with crossing frame)
  if (!iEvent.isRealData()) {
    if (crossingFrame_) {
      auto cscSimHitsXF_handle = iEvent.getHandle(cscSimHitsXFToken_);
      auto rpcSimHitsXF_handle = iEvent.getHandle(rpcSimHitsXFToken_);
      auto gemSimHitsXF_handle = iEvent.getHandle(gemSimHitsXFToken_);
      auto me0SimHitsXF_handle = iEvent.getHandle(me0SimHitsXFToken_);
      auto dtSimHitsXF_handle = iEvent.getHandle(dtSimHitsXFToken_);

      cscSimHitsXF_ = cscSimHitsXF_handle.product();
      rpcSimHitsXF_ = rpcSimHitsXF_handle.product();
      gemSimHitsXF_ = gemSimHitsXF_handle.product();
      me0SimHitsXF_ = me0SimHitsXF_handle.product();
      dtSimHitsXF_ = dtSimHitsXF_handle.product();
    }
  }

  // SimLinks
  if (!iEvent.isRealData()) {
    auto cscStripSimLinks_handle = iEvent.getHandle(cscStripSimLinksToken_);
    auto cscWireSimLinks_handle = iEvent.getHandle(cscWireSimLinksToken_);
    auto rpcDigiSimLinks_handle = iEvent.getHandle(rpcDigiSimLinksToken_);
    auto gemDigiSimLinks_handle = iEvent.getHandle(gemDigiSimLinksToken_);
    auto me0DigiSimLinks_handle = iEvent.getHandle(me0DigiSimLinksToken_);
    auto dtDigiSimLinks_handle = iEvent.getHandle(dtDigiSimLinksToken_);

    cscStripSimLinks_ = cscStripSimLinks_handle.product();
    cscWireSimLinks_ = cscWireSimLinks_handle.product();
    rpcDigiSimLinks_ = rpcDigiSimLinks_handle.product();
    gemDigiSimLinks_ = gemDigiSimLinks_handle.product();
    me0DigiSimLinks_ = me0DigiSimLinks_handle.product();
    dtDigiSimLinks_ = dtDigiSimLinks_handle.product();
  }

  // Geometry
  auto dtGeom_handle = iSetup.getHandle(dtGeomToken_);
  auto cscGeom_handle = iSetup.getHandle(cscGeomToken_);
  auto rpcGeom_handle = iSetup.getHandle(rpcGeomToken_);
  auto gemGeom_handle = iSetup.getHandle(gemGeomToken_);
  auto me0Geom_handle = iSetup.getHandle(me0GeomToken_);

  dtGeom_ = dtGeom_handle.product();
  cscGeom_ = cscGeom_handle.product();
  rpcGeom_ = rpcGeom_handle.product();
  gemGeom_ = gemGeom_handle.product();
  me0Geom_ = me0Geom_handle.product();
}

void SubsystemMCTruth::buildTrackingParticleLinks(const TrackingParticleCollection& trkPartColl) {
  trackingParticleLinks_.clear();

  for (auto it_trkpart = trkPartColl.begin(); it_trkpart != trkPartColl.end(); ++it_trkpart) {
    for (auto it_simtrk = it_trkpart->g4Track_begin(); it_simtrk != it_trkpart->g4Track_end(); ++it_simtrk) {
      unsigned int simTrackId = it_simtrk->trackId();
      EncodedEventId eventId = it_simtrk->eventId();
      SimHitIdpr matchId(simTrackId, eventId);
      //assert(trackingParticleLinks_.find(matchId) == trackingParticleLinks_.end());  // matchId should not already exist
      trackingParticleLinks_[matchId] = std::distance(trkPartColl.begin(), it_trkpart);
    }
  }
}

std::pair<int, int> SubsystemMCTruth::findTrackingParticle(const EMTFHit& hit,
                                                           const TrackingParticleCollection& trkPartColl) const {
  int sim_tp1 = -1;
  int sim_tp2 = -1;
  if (hit.subsystem() == L1TMuon::kCSC) {
    sim_tp1 = findTrackingParticleCSCStrip(hit, trkPartColl);
    sim_tp2 = findTrackingParticleCSCWire(hit, trkPartColl);
  } else if (hit.subsystem() == L1TMuon::kRPC) {
    sim_tp1 = findTrackingParticleRPC(hit, trkPartColl);
    sim_tp2 = sim_tp1;
  } else if (hit.subsystem() == L1TMuon::kGEM) {
    sim_tp1 = findTrackingParticleGEM(hit, trkPartColl);
    sim_tp2 = sim_tp1;
  } else if (hit.subsystem() == L1TMuon::kME0) {
    sim_tp1 = findTrackingParticleME0(hit, trkPartColl);
    sim_tp2 = sim_tp1;
  } else if (hit.subsystem() == L1TMuon::kDT) {
    sim_tp1 = findTrackingParticleDT(hit, trkPartColl);
    sim_tp2 = sim_tp1;
  }
  return std::make_pair(sim_tp1, sim_tp2);
}

int SubsystemMCTruth::findTrackingParticleCSCStrip(const EMTFHit& hit,
                                                   const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  int strip0 = hit.strip();
  int strip1 = (strip0 - 1) / 2 + 1;  // different convention used in CSC StripDigiSimLink (fullstrip)

  // Check all 6 CSC layers
  for (unsigned ilayer = 0; ilayer < 6; ++ilayer) {
    const CSCDetId detid0(hit.rawDetId());
    const CSCDetId detid1(detid0.endcap(), detid0.station(), detid0.ring(), detid0.chamber(), ilayer + 1);

    StripDigiSimLinks::const_iterator cscStripLayerLinks = cscStripSimLinks_->find(detid1);
    if (cscStripLayerLinks != cscStripSimLinks_->end()) {
      for (LayerLinks::const_iterator linkItr = cscStripLayerLinks->begin(); linkItr != cscStripLayerLinks->end();
           ++linkItr) {
        unsigned int channel = linkItr->channel();
        unsigned int simTrackId = linkItr->SimTrackId();
        //unsigned int xfPosition = linkItr->CFposition();
        //unsigned int tofBin = linkItr->TofBin();
        EncodedEventId eventId = linkItr->eventId();
        float fraction = linkItr->fraction();

        if (std::abs(strip1 - static_cast<int>(channel)) <= 3) {  // allow +/-3
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())
            matches[matchId] = 0.;
          matches[matchId] += fraction;
        }
      }
    }
  }
  return findTrackingParticleFromMatches(matches, trkPartColl);
}

int SubsystemMCTruth::findTrackingParticleCSCWire(const EMTFHit& hit,
                                                  const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  int wire0 = hit.wire1();
  int wire1 = (wire0 + 100) + 1;  // different convention used in CSC StripDigiSimLink

  // Check all 6 CSC layers
  for (unsigned ilayer = 0; ilayer < 6; ++ilayer) {
    const CSCDetId detid0(hit.rawDetId());
    const CSCDetId detid1(detid0.endcap(), detid0.station(), detid0.ring(), detid0.chamber(), ilayer + 1);

    WireDigiSimLinks::const_iterator cscWireLayerLinks = cscWireSimLinks_->find(detid1);
    if (cscWireLayerLinks != cscWireSimLinks_->end()) {
      for (LayerLinks::const_iterator linkItr = cscWireLayerLinks->begin(); linkItr != cscWireLayerLinks->end();
           ++linkItr) {
        unsigned int channel = linkItr->channel();
        unsigned int simTrackId = linkItr->SimTrackId();
        //unsigned int xfPosition = linkItr->CFposition();
        //unsigned int tofBin = linkItr->TofBin();
        EncodedEventId eventId = linkItr->eventId();
        float fraction = linkItr->fraction();

        if (std::abs(wire1 - static_cast<int>(channel)) <= 3) {  // allow +/-3
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())
            matches[matchId] = 0.;
          matches[matchId] += fraction;
        }
      }
    }
  }
  return findTrackingParticleFromMatches(matches, trkPartColl);
}

int SubsystemMCTruth::findTrackingParticleRPC(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  // Check all strips in the RPC cluster
  const RPCDetId detid(hit.rawDetId());
  int stripA = hit.stripLo();
  int stripB = hit.stripHi();
  int bx = hit.bx();

#if 0
  // This is giving me an error: Assertion `std::distance(p.first, p.second) == 1' failed.
  RPCDigiSimLinks::const_iterator rpcDigiLayerLinks = rpcDigiSimLinks_->find(detid);
  if (rpcDigiLayerLinks != rpcDigiSimLinks_->end()) {
    for (RPCLayerLinks::const_iterator linkItr = rpcDigiLayerLinks->begin(); linkItr != rpcDigiLayerLinks->end();
         ++linkItr) {
      unsigned int simStrip = linkItr->getStrip();
      unsigned int simBX = linkItr->getBx();
      unsigned int simTrackId = linkItr->getTrackId();
      EncodedEventId eventId = linkItr->getEventId();

      for (int strip0 = stripA; strip0 < stripB + 1; ++strip0) {
        if ((std::abs(strip0 - static_cast<int>(simStrip)) <= 1) && (std::abs(bx - static_cast<int>(simBX)) <= 1)) {  // allow +/-1
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())
            matches[matchId] = 0.;
          matches[matchId] += 1.0;
        }
      }
    }
  }
#else
  // My temporary fix
  for (RPCDigiSimLinks::const_iterator linkItr1 = rpcDigiSimLinks_->begin(); linkItr1 != rpcDigiSimLinks_->end();
       ++linkItr1) {
    for (RPCLayerLinks::const_iterator linkItr = linkItr1->begin(); linkItr != linkItr1->end(); ++linkItr) {
      unsigned int detUnitId = linkItr->getDetUnitId();
      unsigned int simStrip = linkItr->getStrip();
      unsigned int simBX = linkItr->getBx();
      unsigned int simTrackId = linkItr->getTrackId();
      EncodedEventId eventId = linkItr->getEventId();

      if (detUnitId == detid.rawId()) {
        for (int strip0 = stripA; strip0 < stripB + 1; ++strip0) {
          if ((std::abs(strip0 - static_cast<int>(simStrip)) <= 1) && (std::abs(bx - static_cast<int>(simBX)) <= 1)) {  // allow +/-1
            SimHitIdpr matchId(simTrackId, eventId);
            if (matches.find(matchId) == matches.end())
              matches[matchId] = 0.;
            matches[matchId] += 1.0;

            //// Debug
            //std::cout << "Dump RPCDigiSimLink - strip " << linkItr->getStrip() << " bx " << linkItr->getBx()
            //          << " entry " << linkItr->getEntryPoint() << " p4 " << linkItr->getMomentumAtEntry()
            //          << " tof " << linkItr->getTimeOfFlight() << " eloss " << linkItr->getEnergyLoss()
            //          << " pdg " << linkItr->getParticleType() << " process " << linkItr->getProcessType()
            //          << " trackId " << linkItr->getTrackId() << " eventId " << linkItr->getEventId().bunchCrossing()
            //          << ", " << linkItr->getEventId().event() << std::endl;
          }
        }
      }
    }
  }
#endif

  return findTrackingParticleFromMatches(matches, trkPartColl);
}

int SubsystemMCTruth::findTrackingParticleGEM(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  // Check all strips in the GEM cluster
  int stripA = hit.stripLo();
  int stripB = hit.stripHi();
  int bx = hit.bx();

  // Check both GEM layers
  for (unsigned ilayer = 0; ilayer < 2; ++ilayer) {
    const GEMDetId detid0(hit.rawDetId());
    const GEMDetId detid1(
        detid0.region(), detid0.ring(), detid0.station(), ilayer + 1, detid0.chamber(), detid0.roll());

    GEMDigiSimLinks::const_iterator gemDigiLayerLinks = gemDigiSimLinks_->find(detid1);
    if (gemDigiLayerLinks != gemDigiSimLinks_->end()) {
      for (GEMLayerLinks::const_iterator linkItr = gemDigiLayerLinks->begin(); linkItr != gemDigiLayerLinks->end();
           ++linkItr) {
        unsigned int simStrip = linkItr->getStrip();
        unsigned int simBX = linkItr->getBx();
        unsigned int simTrackId = linkItr->getTrackId();
        EncodedEventId eventId = linkItr->getEventId();

        //int simPad = 1 + static_cast<int>(p->padOfStrip(simStrip));
        unsigned int simPad = (hit.station() == 1) ? ((simStrip + 1) / 2) : ((simStrip + 1) / 2);

        for (int strip0 = stripA; strip0 < stripB + 1; ++strip0) {
          if ((std::abs(strip0 - static_cast<int>(simPad)) <= 3) && (std::abs(bx - static_cast<int>(simBX)) <= 1)) {  // allow +/-3
            SimHitIdpr matchId(simTrackId, eventId);
            if (matches.find(matchId) == matches.end())
              matches[matchId] = 0.;
            matches[matchId] += 1.0;

            //// Debug
            //std::cout << "Dump GEMDigiSimLink - strip " << linkItr->getStrip() << " bx " << linkItr->getBx()
            //          << " entry " << linkItr->getEntryPoint() << " p4 " << linkItr->getMomentumAtEntry()
            //          << " tof " << linkItr->getTimeOfFlight() << " eloss " << linkItr->getEnergyLoss()
            //          << " pdg " << linkItr->getParticleType() << " process " << linkItr->getProcessType()
            //          << " trackId " << linkItr->getTrackId() << " eventId " << linkItr->getEventId().bunchCrossing()
            //          << ", " << linkItr->getEventId().event() << std::endl;
          }
        }
      }
    }
  }
  return findTrackingParticleFromMatches(matches, trkPartColl);
}

int SubsystemMCTruth::findTrackingParticleME0(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  int strip0 = hit.strip();  // in half-strip unit
  int strip1 = (strip0 >> 1);
  int partition = hit.wire1();  // in half-roll unit
  int iroll = (partition >> 1) + 1;
  int bx = hit.bx();

  // Check all 6 ME0 layers.
  for (unsigned ilayer = 0; ilayer < 6; ++ilayer) {
    // Check neighbor rolls
    unsigned iroll_first = (iroll == 1) ? iroll : iroll - 1;
    unsigned iroll_last = (iroll == 8) ? iroll : iroll + 1;
    for (unsigned iroll = iroll_first; iroll <= iroll_last; ++iroll) {
      ME0DetId detid0(hit.rawDetId());
      ME0DetId detid1(detid0.region(), ilayer + 1, detid0.chamber(), iroll);

      ME0DigiSimLinks::const_iterator me0DigiLayerLinks = me0DigiSimLinks_->find(detid1);
      if (me0DigiLayerLinks != me0DigiSimLinks_->end()) {
        for (ME0LayerLinks::const_iterator linkItr = me0DigiLayerLinks->begin(); linkItr != me0DigiLayerLinks->end();
             ++linkItr) {
          unsigned int simStrip = linkItr->getStrip();
          unsigned int simBX = linkItr->getBx();
          unsigned int simTrackId = linkItr->getTrackId();
          EncodedEventId eventId = linkItr->getEventId();

          if ((std::abs(strip1 - static_cast<int>(simStrip)) <= 3) && (std::abs(bx - static_cast<int>(simBX)) <= 1)) {  // allow +/-3
            SimHitIdpr matchId(simTrackId, eventId);
            if (matches.find(matchId) == matches.end())
              matches[matchId] = 0.;
            matches[matchId] += 1.0;
          }
        }
      }
    }
  }
  return findTrackingParticleFromMatches(matches, trkPartColl);
}

int SubsystemMCTruth::findTrackingParticleDT(const EMTFHit& hit, const TrackingParticleCollection& trkPartColl) const {
  std::map<SimHitIdpr, float> matches;

  DTChamberId detid(hit.rawDetId());
  int bx = hit.bx();

  for (DTDigiSimLinkCollection::DigiRangeIterator detUnit = dtDigiSimLinks_->begin(); detUnit != dtDigiSimLinks_->end();
       ++detUnit) {
    const DTLayerId& layerid = (*detUnit).first;
    const DTDigiSimLinkCollection::Range& range = (*detUnit).second;

    if (detid == layerid.chamberId()) {
      for (DTDigiSimLinkCollection::const_iterator linkItr = range.first; linkItr != range.second; ++linkItr) {
        // Unfortunately lost all these info in the L1 data format L1MuDTChambPhDigi
        //float digitime = linkItr->time();
        //int wire = linkItr->wire();
        //int digiNumber = linkItr->number();
        unsigned int simTrackId = linkItr->SimTrackId();
        EncodedEventId eventId = linkItr->eventId();

        if (bx == eventId.bunchCrossing()) {
          SimHitIdpr matchId(simTrackId, eventId);
          if (matches.find(matchId) == matches.end())
            matches[matchId] = 0.;
          matches[matchId] += 1.0;
        }
      }
    }
  }
  return findTrackingParticleFromMatches(matches, trkPartColl);
}

int SubsystemMCTruth::findTrackingParticleFromMatches(const std::map<SimHitIdpr, float>& matches,
                                                      const TrackingParticleCollection& trkPartColl) const {
  // Find the matched TP with the highest weight
  int best_sim_tp = -1;  // index of the tp
  float max_weight = 0.;
  for (auto it_match = matches.begin(); it_match != matches.end(); ++it_match) {
    auto found = trackingParticleLinks_.find(it_match->first);
    if (found != trackingParticleLinks_.end()) {
      int sim_tp = found->second;
      //auto it_trkpart = std::next(trkPartColl.begin(), sim_tp);
      if (max_weight < it_match->second) {
        max_weight = it_match->second;
        best_sim_tp = sim_tp;
      }
    }
  }
  return best_sim_tp;
}

int SubsystemMCTruth::findTrackingParticleFromPSimHit(const PSimHit& pSimHit) const {
  int sim_tp = -1;  // index of the tp
  unsigned int simTrackId = pSimHit.trackId();
  EncodedEventId eventId = pSimHit.eventId();
  SimHitIdpr matchId(simTrackId, eventId);
  auto found = trackingParticleLinks_.find(matchId);
  if (found != trackingParticleLinks_.end()) {
    sim_tp = found->second;
  }
  return sim_tp;
}

SubsystemMCTruth::SimHitInfoCollection SubsystemMCTruth::findSimHits() const {
  SimHitInfoCollection result;

  // CSC
  std::transform(cscSimHits_->begin(), cscSimHits_->end(), std::back_inserter(result), [this](const PSimHit& pSimHit) {
    unsigned int detUnitId = pSimHit.detUnitId();
    const CSCLayer* layer = cscGeom_->layer(detUnitId);
    const CSCDetId& detid = layer->id();
    const LocalPoint& lp = pSimHit.localPosition();
    const GlobalPoint& gp = layer->toGlobal(lp);
    int subsystem = L1TMuon::kCSC;
    int endcap = (detid.endcap() == 2) ? -1 : detid.endcap();
    int sim_tp = findTrackingParticleFromPSimHit(pSimHit);
    SimHitInfo out{pSimHit, lp, gp, subsystem, endcap, detid.station(), detid.ring(), detid.chamber(), detid.layer(), sim_tp};
    return out;
  });

  // RPC
  edm::PSimHitContainer rpcSimHitsNoBarrel;
  std::copy_if(
      rpcSimHits_->begin(), rpcSimHits_->end(), std::back_inserter(rpcSimHitsNoBarrel), [this](const PSimHit& pSimHit) {
        unsigned int detUnitId = pSimHit.detUnitId();
        const RPCRoll* roll = rpcGeom_->roll(detUnitId);
        const RPCDetId& detid = roll->id();

        // Identifier for barrel RPC
        const bool is_barrel = (detid.region() == 0);
        return (not is_barrel);
      });

  std::transform(
      rpcSimHitsNoBarrel.begin(), rpcSimHitsNoBarrel.end(), std::back_inserter(result), [this](const PSimHit& pSimHit) {
        unsigned int detUnitId = pSimHit.detUnitId();
        const RPCRoll* roll = rpcGeom_->roll(detUnitId);
        const RPCDetId& detid = roll->id();
        const LocalPoint& lp = pSimHit.localPosition();
        const GlobalPoint& gp = roll->toGlobal(lp);
        int subsystem = L1TMuon::kRPC;
        int endcap = detid.region();
        int sim_tp = findTrackingParticleFromPSimHit(pSimHit);

        // Identifier for barrel RPC
        const bool is_barrel = (detid.region() == 0);
        // Identifier for iRPC (RE3/1, RE4/1)
        const bool is_irpc = ((not is_barrel) and (detid.station() >= 3) and (detid.ring() == 1));
        const int detid_chamber = ((detid.sector() - 1) * (is_irpc ? 3 : 6)) + detid.subsector();

        SimHitInfo out{pSimHit, lp, gp, subsystem, endcap, detid.station(), detid.ring(), detid_chamber, detid.layer(), sim_tp};
        return out;
      });

  // GEM
  std::transform(gemSimHits_->begin(), gemSimHits_->end(), std::back_inserter(result), [this](const PSimHit& pSimHit) {
    unsigned int detUnitId = pSimHit.detUnitId();
    const GEMEtaPartition* roll = gemGeom_->etaPartition(detUnitId);
    const GEMDetId& detid = roll->id();
    const LocalPoint& lp = pSimHit.localPosition();
    const GlobalPoint& gp = roll->toGlobal(lp);
    int subsystem = L1TMuon::kGEM;
    int endcap = detid.region();
    int sim_tp = findTrackingParticleFromPSimHit(pSimHit);
    SimHitInfo out{pSimHit, lp, gp, subsystem, endcap, detid.station(), detid.ring(), detid.chamber(), detid.layer(), sim_tp};
    return out;
  });

  // ME0
  std::transform(me0SimHits_->begin(), me0SimHits_->end(), std::back_inserter(result), [this](const PSimHit& pSimHit) {
    unsigned int detUnitId = pSimHit.detUnitId();
    const ME0EtaPartition* roll = me0Geom_->etaPartition(detUnitId);
    const ME0DetId& detid = roll->id();
    const LocalPoint& lp = pSimHit.localPosition();
    const GlobalPoint& gp = roll->toGlobal(lp);
    int subsystem = L1TMuon::kME0;
    int endcap = detid.region();
    int sim_tp = findTrackingParticleFromPSimHit(pSimHit);
    SimHitInfo out{pSimHit, lp, gp, subsystem, endcap, detid.station(), 4, detid.chamber(), detid.layer(), sim_tp};
    return out;
  });

#if 0
  // DT
  std::transform(dtSimHits_->begin(), dtSimHits_->end(), std::back_inserter(result), [this](const PSimHit& pSimHit) {
    unsigned int detUnitId = pSimHit.detUnitId();
    const DTWireId wireid(detUnitId);
    const DTLayer* layer = dtGeom_->layer(wireid.layerId());
    const DTLayerId& detid = layer->id();
    const LocalPoint& lp = pSimHit.localPosition();
    const GlobalPoint& gp = layer->toGlobal(lp);
    int subsystem = L1TMuon::kDT;
    int endcap = 0;
    int sim_tp = findTrackingParticleFromPSimHit(pSimHit);
    SimHitInfo out{pSimHit, lp, gp, subsystem, endcap, detid.station(), 1, detid.sector(), detid.layer(), sim_tp};
    return out;
  });
#endif

  //// Debug
  //for (auto&& simhit : result) {
  //  std::cout << "type, lay, cmb, phi, theta, eta, r, z: " << simhit.subsystem << " " << simhit.layer << " "
  //            << simhit.chamber << " " << simhit.globalPosition.phi() << " " << simhit.globalPosition.theta() << " "
  //            << simhit.globalPosition.eta() << " " << simhit.globalPosition.perp() << " " << simhit.globalPosition.z()
  //            << std::endl;
  //}
  return result;
}
