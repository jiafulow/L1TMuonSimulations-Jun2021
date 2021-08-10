import FWCore.ParameterSet.Config as cms

def customise_ntupler_step(process):
    # Add phase2L1EMTFSequence
    process.load('L1Trigger.Phase2L1EMTF.rpcRecHitsForEMTF_cfi')
    process.load('L1Trigger.Phase2L1EMTF.phase2L1EMTFProducer_cfi')
    process.phase2L1EMTFTask = cms.Task(process.rpcRecHitsForEMTF, process.phase2L1EMTFProducer)
    process.phase2L1EMTFSequence = cms.Sequence(process.phase2L1EMTFTask)
    #process.phase2L1EMTFSequence = cms.Sequence(process.rpcRecHitsForEMTF+process.phase2L1EMTFProducer)  # without cms.Task

    # Add ntuplerSequence, ntupler_step
    process.load('L1TMuonSimulations.NtupleTools.ntupler_cfi')
    process.TFileService = cms.Service('TFileService', fileName = process.ntupler.fileName)
    process.ntuplerSequence = cms.Sequence(process.ntupler)
    try:
      process.ntupler_step = cms.Path(process.generator+process.phase2L1EMTFSequence+process.ntuplerSequence)
    except AttributeError:
      process.ntupler_step = cms.Path(process.phase2L1EMTFSequence+process.ntuplerSequence)
    process.schedule.extend([process.ntupler_step])

    # Remove cms.EndPath instances from schedule
    paths_in_schedule = [path for path in process.schedule if not isinstance(path, cms.EndPath)]
    process.schedule = cms.Schedule(*paths_in_schedule)
    return process