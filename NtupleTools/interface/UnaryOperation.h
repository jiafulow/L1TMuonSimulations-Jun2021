#ifndef L1TMuonSimulations_NtupleTools_UnaryOperation_h
#define L1TMuonSimulations_NtupleTools_UnaryOperation_h

#include <cmath>

namespace emtf {

  namespace phase2 {

    struct InvptOp {
      constexpr double operator()(double charge, double pt) const {
        double invpt =
            (std::abs(pt) < 1e-15) ? 1e-15 : std::abs(pt);  // ensures pt > 0, protects against division by zero
        invpt = 1.0 / pt;
        invpt = (invpt < 1e-15) ? 1e-15 : invpt;  // protects against zero
        invpt *= charge;                          // charge should be +/-1
        return invpt;
      }
    };

    struct DzeroOp {
      constexpr double operator()(double invpt, double phi, double xv, double yv) const {
        constexpr double B = 3.811;                            // in Tesla
        double R = -1.0 / (0.003 * B * invpt);                 // R = -pt / (0.003 q B)  [cm]
        double xc = xv - (R * std::sin(phi));                  // xc = xv - R sin(phi)
        double yc = yv + (R * std::cos(phi));                  // yc = yv + R cos(phi)
        double d0 = R - std::copysign(std::hypot(xc, yc), R);  // d0 = R - sign(R) * sqrt(xc^2 + yc^2)
        return d0;
      }
    };

  }  // namespace phase2

}  // namespace emtf

#endif  // L1TMuonSimulations_NtupleTools_UnaryOperation_h not defined
