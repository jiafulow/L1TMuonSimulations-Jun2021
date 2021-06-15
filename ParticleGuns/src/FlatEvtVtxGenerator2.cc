

#include "L1TMuonSimulations/ParticleGuns/interface/FlatEvtVtxGenerator2.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
//#include "CLHEP/Vector/ThreeVector.h"
#include "HepMC/SimpleVector.h"

FlatEvtVtxGenerator2::FlatEvtVtxGenerator2(const edm::ParameterSet& p) : BaseEvtVtxGenerator2(p) {
  fMinX = p.getParameter<double>("MinX") * cm;
  fMinY = p.getParameter<double>("MinY") * cm;
  fMinZ = p.getParameter<double>("MinZ") * cm;
  fMaxX = p.getParameter<double>("MaxX") * cm;
  fMaxY = p.getParameter<double>("MaxY") * cm;
  fMaxZ = p.getParameter<double>("MaxZ") * cm;
  fMinT = p.getParameter<double>("MinT") * ns * c_light;
  fMaxT = p.getParameter<double>("MaxT") * ns * c_light;

  if (fMinX > fMaxX) {
    throw cms::Exception("Configuration") << "Error in FlatEvtVtxGenerator2: "
                                          << "MinX is greater than MaxX";
  }
  if (fMinY > fMaxY) {
    throw cms::Exception("Configuration") << "Error in FlatEvtVtxGenerator2: "
                                          << "MinY is greater than MaxY";
  }
  if (fMinZ > fMaxZ) {
    throw cms::Exception("Configuration") << "Error in FlatEvtVtxGenerator2: "
                                          << "MinZ is greater than MaxZ";
  }
  if (fMinT > fMaxT) {
    throw cms::Exception("Configuration") << "Error in FlatEvtVtxGenerator2: "
                                          << "MinT is greater than MaxT";
  }
}

FlatEvtVtxGenerator2::~FlatEvtVtxGenerator2() {}

//Hep3Vector * FlatEvtVtxGenerator2::newVertex() {
HepMC::FourVector FlatEvtVtxGenerator2::newVertex(CLHEP::HepRandomEngine* engine) const {
  double aX, aY, aZ, aT;
  aX = CLHEP::RandFlat::shoot(engine, fMinX, fMaxX);
  aY = CLHEP::RandFlat::shoot(engine, fMinY, fMaxY);
  aZ = CLHEP::RandFlat::shoot(engine, fMinZ, fMaxZ);
  aT = CLHEP::RandFlat::shoot(engine, fMinT, fMaxT);

  return HepMC::FourVector(aX, aY, aZ, aT);
}

void FlatEvtVtxGenerator2::minX(double min) { fMinX = min; }

void FlatEvtVtxGenerator2::minY(double min) { fMinY = min; }

void FlatEvtVtxGenerator2::minZ(double min) { fMinZ = min; }

void FlatEvtVtxGenerator2::maxX(double max) { fMaxX = max; }

void FlatEvtVtxGenerator2::maxY(double max) { fMaxY = max; }

void FlatEvtVtxGenerator2::maxZ(double max) { fMaxZ = max; }
