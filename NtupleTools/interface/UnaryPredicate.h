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
        // Signal event
        //bool signal = (part.eventId().event() == 0);
        // In time bunch-crossing
        //bool intime = (part.eventId().bunchCrossing() == 0);
        // In time + out-of-time bunch-crossing (-2 <= BX <= +2)
        bool kindof_intime = ((-2 <= part.eventId().bunchCrossing()) and (part.eventId().bunchCrossing() <= 2));
        // Primary: charged, pt > 0.2 GeV, |eta| < 3.0, |rho0| < 0.5 cm, |z0| < 30 cm
        //bool primary = ((part.charge() != 0) and (part.pt() > 0.2) and (std::abs(part.eta()) <= 3.0) and
        //                (std::hypot(part.vx(), part.vy()) < 0.5) and (std::abs(part.vz()) < 30.0));
        // Primary+secondary: charged, pt > 0.5 GeV, |eta| < 3.0, |x0| < 300 cm, |y0| < 300 cm, |z0| < 500 cm
        bool kindof_primary =
            ((part.charge() != 0) and (part.pt() > 0.5) and (std::abs(part.eta()) <= 3.0) and
             (std::abs(part.vx()) <= 300.) and (std::abs(part.vy()) <= 300.) and (std::abs(part.vz()) <= 500.));
        // Not decayed
        //bool nodecay = (part.decayVertices().empty());

        bool passed = kindof_intime and kindof_primary;
        return passed;
      }
    };

  }  // namespace phase2

}  // namespace emtf

#endif  // L1TMuonSimulations_NtupleTools_UnaryPredicate_h not defined
