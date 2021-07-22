#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>
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

class NtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {
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

  void initTree();
  void fillTree();
  void writeTree();

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

  // TTree
  TTree* tree_;
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

  // ___________________________________________________________________________
  // SimHits
  for (const auto& hit : sim_hits) {
    //TODO: implement this
  }

  // ___________________________________________________________________________
  // Tracks
  for (const auto& trk : emtf_tracks) {
    //TODO: implement this
  }

  // ___________________________________________________________________________
  // L1TrackTrigger tracks
  if ((tttfTracks_ != nullptr) and (tttfTrackAssoc_ != nullptr)) {
    //TODO: implement this
  }

  // ___________________________________________________________________________
  // TrackingParticles
  for (const auto& part : trk_particles) {
    //TODO: implement this
  }

  // ___________________________________________________________________________
  // PileupSummaryInfo
  if (pileupInfo_ != nullptr) {
    //TODO: implement this
  }

  // ___________________________________________________________________________
  // Fill
  fillTree();
}

// _____________________________________________________________________________
void NtupleMaker::initTree() {
  // Create tree using TFileService
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");
}

void NtupleMaker::fillTree() { tree_->Fill(); }

void NtupleMaker::writeTree() {
  // Writing is handled by TFileService
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
