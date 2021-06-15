#ifndef FlatRandomPtGunProducer2_H
#define FlatRandomPtGunProducer2_H

/** \class FlatRandomPtGunProducer2
 *
 * Generates single particle gun in HepMC format
 * Julia Yarba 12/2005 
 ***************************************/

#include "L1TMuonSimulations/ParticleGuns/interface/BaseFlatGunProducer2.h"

namespace edm {

  class FlatRandomPtGunProducer2 : public BaseFlatGunProducer2 {
  public:
    explicit FlatRandomPtGunProducer2(const ParameterSet&);
    ~FlatRandomPtGunProducer2() override;

    void produce(Event&, const EventSetup&) override;

  private:
    double fMinPt;
    double fMaxPt;
    double fMinOneOverPt;
    double fMaxOneOverPt;
    bool fRandomCharge;
    std::string fPtSpectrum;
  };
}  // namespace edm

#endif
