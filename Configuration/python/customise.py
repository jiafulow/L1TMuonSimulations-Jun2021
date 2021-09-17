import FWCore.ParameterSet.Config as cms

def customise_ntupler_step(process):
    # Add phase2L1EMTFSequence
    process.load('L1Trigger.Phase2L1EMTF.simCscTriggerPrimitiveDigisForEMTF_cfi')
    process.load('L1Trigger.Phase2L1EMTF.rpcRecHitsForEMTF_cfi')
    process.load('L1Trigger.Phase2L1EMTF.phase2L1EMTFProducer_cfi')
    process.phase2L1EMTFTask = cms.Task(
        process.simCscTriggerPrimitiveDigisForEMTF,
        process.rpcRecHitsForEMTF,
        process.phase2L1EMTFProducer
    )
    process.phase2L1EMTFSequence = cms.Sequence(process.phase2L1EMTFTask)
    ## Alternatively, do it without cms.Task
    #process.phase2L1EMTFSequence = cms.Sequence(
    #    process.simCscTriggerPrimitiveDigisForEMTF +
    #    process.rpcRecHitsForEMTF +
    #    process.phase2L1EMTFProducer
    #)

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

def remove_L1TrackTrigger_step(process):
    # Remove L1TrackTrigger_step from schedule
    if process.L1TrackTrigger_step in process.schedule:
        process.schedule.remove(process.L1TrackTrigger_step)
    return process

def remove_L1simulation_step(process):
    # Remove L1simulation_step from schedule
    if process.L1simulation_step in process.schedule:
        process.schedule.remove(process.L1simulation_step)
    return process
