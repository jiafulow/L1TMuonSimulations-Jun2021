/*
 *  \author Julia Yarba
 */

#include "L1TMuonSimulations/ParticleGuns/interface/BaseFlatGunProducer2.h"

#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

using namespace edm;
using namespace CLHEP;

BaseFlatGunProducer2::BaseFlatGunProducer2(const ParameterSet& pset)
    : fEvt(nullptr)
// fPDGTable( new DefaultConfig::ParticleDataTable("PDG Table") )
{
  Service<RandomNumberGenerator> rng;
  if (!rng.isAvailable()) {
    throw cms::Exception("Configuration")
        << "The RandomNumberProducer module requires the RandomNumberGeneratorService\n"
           "which appears to be absent.  Please add that service to your configuration\n"
           "or remove the modules that require it.";
  }

  const ParameterSet& pgun_params = pset.getParameter<ParameterSet>("PGunParameters");

  // although there's the method ParameterSet::empty(),
  // it looks like it's NOT even necessary to check if it is,
  // before trying to extract parameters - if it is empty,
  // the default values seem to be taken
  fPartIDs = pgun_params.getParameter<std::vector<int> >("PartID");
  fMinEta = pgun_params.getParameter<double>("MinEta");
  fMaxEta = pgun_params.getParameter<double>("MaxEta");
  fMinPhi = pgun_params.getParameter<double>("MinPhi");
  fMaxPhi = pgun_params.getParameter<double>("MaxPhi");

  //
  //fPDGTablePath = "/afs/cern.ch/sw/lcg/external/clhep/1.9.2.1/slc3_ia32_gcc323/data/HepPDT/" ;
  /*
  string HepPDTBase( std::getenv("HEPPDT_PARAM_PATH") ) ; 
  fPDGTablePath = HepPDTBase + "/data/" ;
  fPDGTableName = "PDG_mass_width_2004.mc"; // should it be 2004 table ?

  string TableFullName = fPDGTablePath + fPDGTableName ;
  std::ifstream PDFile( TableFullName.c_str() ) ;
  if( !PDFile ) 
  {
      throw cms::Exception("FileNotFound", "BaseFlatGunProducer2::BaseFlatGunProducer2()")
	<< "File " << TableFullName << " cannot be opened.\n";
  }

  HepPDT::TableBuilder tb(*fPDGTable) ;
  if ( !addPDGParticles( PDFile, tb ) ) { std::cout << " Error reading PDG !" << std::endl; }
  // the tb dtor fills fPDGTable
*/

  fVerbosity = pset.getUntrackedParameter<int>("Verbosity", 0);

  fAddAntiParticle = pset.getParameter<bool>("AddAntiParticle");

  produces<GenRunInfoProduct, Transition::EndRun>();
}

BaseFlatGunProducer2::~BaseFlatGunProducer2() {
  // no need to cleanup GenEvent memory - done in HepMCProduct
  // if (fEvt != NULL) delete fEvt ; // double check
  // delete fPDGTable;
}

void BaseFlatGunProducer2::beginRun(const Run& run, const EventSetup& es) {
  es.getData(fPDGTable);
  return;
}
void BaseFlatGunProducer2::endRun(const Run& run, const EventSetup& es) {}

void BaseFlatGunProducer2::endRunProduce(Run& run, const EventSetup& es) {
  // just create an empty product
  // to keep the EventContent definitions happy
  // later on we might put the info into the run info that this is a PGun
  std::unique_ptr<GenRunInfoProduct> genRunInfo(new GenRunInfoProduct());
  run.put(std::move(genRunInfo));
}
