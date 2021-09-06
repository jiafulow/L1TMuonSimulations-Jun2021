#ifndef L1TMuonSimulations_NtupleTools_UnaryPredicate_h
#define L1TMuonSimulations_NtupleTools_UnaryPredicate_h

#include <cmath>

namespace emtf {

  namespace phase2 {

    struct EMTFHitPred {
      template <typename T>
      constexpr bool operator()(const T& hit) const {
        bool passed = ((-2 <= hit.bx()) and (hit.bx() <= 2));  // only BX=[-2..+2]
        return passed;
      }
    };

    struct EMTFTrackPred {
      template <typename T>
      constexpr bool operator()(const T& trk) const {
        bool passed = (trk.bx() == 0);  // only BX=0
        return passed;
      }
    };

    struct TrackingParticlePred {
      template <typename T>
      constexpr bool operator()(const T& part) const {
        // Keep leptons, W/Z/H, charged particles, and high-pt gamma.
        int pdgid = std::abs(part.pdgId());
        bool initial_passed = (pdgid == 13);
        initial_passed |= ((pdgid == 23) or (pdgid == 24) or (pdgid == 25) or (pdgid == 37));
        initial_passed |= (((pdgid == 11) or (pdgid == 15)) and (part.pt() > 1.));
        initial_passed |= (((pdgid == 211) or (pdgid == 321) or (pdgid == 2212)) and (part.pt() > 1.));
        initial_passed |= ((part.charge() != 0) and (part.status() == 1) and (part.pt() > 2.));  // charged above 2 GeV
        initial_passed |= ((pdgid == 22) and (part.status() == 1) and (part.pt() > 10.));        // gamma above 10 GeV
        initial_passed |= (((1000001 <= pdgid) and (pdgid <= 1000039)) or
                           ((2000001 <= pdgid) and (pdgid <= 2000015)));  // SUSY fiction particles

        // Signal event
        //bool signal = ((part.eventId().bunchCrossing() == 0) and (part.eventId().event() == 0));
        // In time bunch-crossing
        //bool intime = (part.eventId().bunchCrossing() == 0);
        // In time + out of time bunch-crossing
        bool kindof_intime = ((-2 <= part.eventId().bunchCrossing()) and (part.eventId().bunchCrossing() <= 2));

        // Primary: pt > 0.1 GeV, |eta| < 3.0, |rho0| < 0.5 cm, |z0| < 30 cm
        //bool primary = ((part.pt() > 0.1) and (std::abs(part.eta()) <= 3.0) and
        //                (std::hypot(part.vx(), part.vy()) < 0.5) and (std::abs(part.vz()) < 30.0));
        // Primary + secondary: pt > 0.1 GeV, |eta| < 5.0, |x0| < 600 cm, |y0| < 600 cm, |z0| < 1200 cm
        bool kindof_primary = ((part.pt() > 0.1) and (std::abs(part.eta()) <= 5.0) and (std::abs(part.vx()) <= 600.) and
                               (std::abs(part.vy()) <= 600.) and (std::abs(part.vz()) <= 1200.));
        // Keep trk particles coming from gen particles
        bool from_genp = (not part.genParticles().empty());
        // Remove trk particles without muon hits
        bool killed = ((part.status() == -99) and (part.genParticles().empty()) and
                       ((part.numberOfHits() == 0) or (part.numberOfHits() == part.numberOfTrackerHits())));
        // Not decayed
        //bool nodecay = (part.decayVertices().empty());
        bool passed = (initial_passed and kindof_intime and ((kindof_primary or from_genp) and (not killed)));
        return passed;
      }
    };

  }  // namespace phase2

}  // namespace emtf

#endif  // L1TMuonSimulations_NtupleTools_UnaryPredicate_h not defined
