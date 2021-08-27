#include "L1TMuonSimulations/NtupleTools/interface/GenParticleDescendentFinder.h"

#include <algorithm>
#include <iostream>

using namespace emtf::phase2;

void GenParticleDescendentFinder::build(const reco::GenParticleCollection& src) {
  bool keepOrDropAll = false;
  bool keepOrDrop = true;
  flags_.clear();
  flags_.resize(src.size(), keepOrDropAll);

  for (reco::GenParticleCollection::const_iterator i = src.begin(); i != src.end(); ++i) {
    int pdgid = std::abs(i->pdgId());
    // Mark W/Z/H descendents
    if ((pdgid == 23) or (pdgid == 24) or (pdgid == 25) or (pdgid == 37)) {
      size_type index = std::distance(src.begin(), i);
      std::vector<size_type> allIndices;
      flags_.at(index) = keepOrDrop;
      recursiveFlagDaughters(index, src, keepOrDrop, allIndices);
    }
  }
}

void GenParticleDescendentFinder::recursiveFlagDaughters(size_type index,
                                                         const reco::GenParticleCollection& src,
                                                         bool keepOrDrop,
                                                         std::vector<size_type>& allIndices) {
  const reco::GenParticleRefVector& daughters = src[index].daughterRefVector();
  // avoid infinite recursion if the daughters are set to "this" particle.
  size_type cachedIndex = index;

  for (reco::GenParticleRefVector::const_iterator i = daughters.begin(); i != daughters.end(); ++i) {
    index = i->key();
    // To also avoid infinite recursion if a "loop" is found in the daughter list,
    // check to make sure the index hasn't already been added.
    if (std::find(allIndices.begin(), allIndices.end(), index) == allIndices.end()) {
      allIndices.push_back(index);
      if (cachedIndex != index) {
        flags_.at(index) = keepOrDrop;
        recursiveFlagDaughters(index, src, keepOrDrop, allIndices);
      }
    }
  }
}
