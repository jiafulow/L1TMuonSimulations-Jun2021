#ifndef L1TMuonSimulations_NtupleTools_GenParticleDescendentFinder_h
#define L1TMuonSimulations_NtupleTools_GenParticleDescendentFinder_h

// References:
//     PhysicsTools/HepMCCandAlgos/plugins/GenParticlePruner.cc
//     PhysicsTools/HepMCCandAlgos/src/GenParticlesHelper.cc

#include <cassert>
#include <vector>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

namespace emtf {

  namespace phase2 {

    class GenParticleDescendentFinder {
    public:
      typedef reco::GenParticleRefVector::key_type size_type;  // should be unsigned int

      void build(const reco::GenParticleCollection& src);

      // Check if GenParticle (given by its index) descends from W/Z/H
      bool isImportant(size_type index) const { return flags_.at(index); }

    private:
      void recursiveFlagDaughters(size_type index,
                                  const reco::GenParticleCollection& src,
                                  bool keepOrDrop,
                                  std::vector<size_type>& allIndices);

      std::vector<bool> flags_;
    };

  }  // namespace phase2

}  // namespace emtf

#endif  // L1TMuonSimulations_NtupleTools_GenParticleDescendentFinder_h not defined
