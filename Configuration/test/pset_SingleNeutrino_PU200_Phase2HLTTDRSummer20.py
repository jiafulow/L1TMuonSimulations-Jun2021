# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step1 -n 100 -s L1TrackTrigger,L1 --conditions auto:phase2_realistic --era Phase2C9 --no_exec --fileout file:step1.root --eventcontent FEVTDEBUGHLT --filein /store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/031CF6BB-92EA-D340-A3B9-8B6CD1B068EE.root --datatier GEN-SIM-DIGI-RAW --python_filename pset_SingleNeutrino_PU200_Phase2HLTTDRSummer20.py --customise SimGeneral/MixingModule/customiseStoredTPConfig.higherPtTP,SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110,L1TMuonSimulations/Configuration/customise.customise_ntupler_step,L1TMuonSimulations/Configuration/customise.remove_L1TrackTrigger_step,L1TMuonSimulations/Configuration/customise.remove_L1simulation_step --geometry Extended2026D49 --mc
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

process = cms.Process('L1',Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/031CF6BB-92EA-D340-A3B9-8B6CD1B068EE.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/04996FF3-DE23-904D-9D48-6A45C0877945.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/0DB38DFA-3D11-8D48-A2B3-2CA9ACEB9B33.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/11B78060-E5B6-4E46-A017-2F7F628C826A.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/1AA8F787-C9DB-7A4D-A3E1-736FEDBABF7A.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/2533D662-0B27-3840-8EC0-C5A324855128.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/2B735638-793C-5D4E-869F-65F53210AED6.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/3600A30B-EDD1-AA45-B95C-7EFB1DA9C3FC.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/453F4996-C90F-8D43-914C-FF93AD1D42B8.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/516C1162-3ABA-2346-9F2B-303F3E4CD80C.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/6393392D-7560-BE48-888D-5AD477409053.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/6BD701F5-70A1-8E4B-9293-F3E69A5D0081.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/7121B81D-D065-9041-B84D-22D38E178ECB.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/7541EE82-3EE5-E142-B4E1-9AD87FF356CC.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/787E8DF4-131E-0341-9824-57810B201F34.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/85A993AF-6C06-BB4D-903E-AEA8D5968307.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/8892C18F-A84E-904A-BA1B-3F93ADFABD22.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/90DFD05E-62CD-4944-A77E-56B013A69850.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/97E6B34C-71BB-9E4C-B692-9315C3BFF64F.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/A83F6422-B237-0349-821B-5CD6A3AFB634.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/AABDB026-9E38-E843-A141-840AB3317E2A.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/B05B2365-8797-0B4A-9A81-DA79CAEE2975.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/B363CA03-2B83-CD4D-B16B-B2A534B8D28B.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/B4FD847D-9DE7-C14C-9CB9-407EAE509ADA.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/B63E1704-CBFC-7943-AFE5-15825E5B3A2E.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/B9B97EC7-75B7-C244-AE3F-99E476C2A428.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/CA4A2DF9-0CED-7B4A-85F6-6922538F7E28.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/D32D341E-D0B7-004E-9B02-7355903E99F7.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/D4F32B1D-9FFC-CC44-AE5C-BCE1DB97D4F7.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/E450C0FC-9741-2D42-BBD3-888E4842E058.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110002/EEDFBB9A-B6C7-2544-BC19-564B55D26646.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/03A80EA6-8BCE-8B40-B51A-7845A27B5F12.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/0C57E928-98D1-DD48-84C1-53978A69D47F.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/1F58CAA6-2458-2346-BD3F-C7B4F151464E.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/2441067F-E6C1-8242-9018-9289CE59880D.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/3976F74B-314E-954A-89D3-41469507E3DC.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/436490B4-3BBF-AD4C-B0AD-D4F5B6170FCB.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/5311B3D0-187C-A34E-B914-79FF1E30FFD1.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/6929D0AE-DBA4-AA41-91FB-5E0B7473270B.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/6A09DDC7-A0CA-3541-992B-03E02A2F81AA.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/9300CADC-4881-5847-ABC4-41A485F65B09.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/986DCE90-6D7F-094B-9A2C-91591D7BDDC2.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/A4C69CDB-79F6-C241-ACD8-330422A3A620.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/B44EE0E4-77A8-A043-9B25-4CD01D83CC35.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/CF5CC650-A488-1845-906D-1003150A4039.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/D37CA82A-EE42-F646-A63B-BD14BF4BE2BF.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/DD7D7AA3-2391-8542-950F-57EF66AA2C46.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/E45EEE77-8B85-1B49-BC9B-E2AEF81D8C9C.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/110003/E5917D7F-5A65-2A43-90A1-BBA3697DFD3C.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/01FBD78D-2AFF-1E4A-9BEA-F6AAC5BFFCB1.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/0651BB29-1704-3C44-A616-320E46BF69A0.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/090FA07A-A311-DE43-ACE3-D67555DF8521.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/0E0D3FFE-64E0-C64D-AD64-FB020F038F6F.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/20C86B94-03BC-8E4F-B2F5-9718F4A1A762.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/2361F9B0-6CCF-CE48-83B3-5DAD529912CD.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/42A4749A-4562-7E41-9DCB-EEA13D42F7D6.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/52E976AC-0110-6647-A5B9-4A6FB09FD3C1.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/6EF2E6B0-299E-F94A-A08F-CE0EBC9EBDD9.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/79E1CA0F-E849-464B-B072-068D993AF561.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/C5A1C347-53FB-B24A-BFF2-0055CDF1A502.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/CCCD8F33-0146-D742-91BA-19573816279F.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/D949CA34-0EDE-BB4A-9011-31C7FCE0DE16.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/DD962499-0B34-8947-A26D-79EA8D4F9581.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/ECEBB680-A175-DA4F-A505-512A860EFB99.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/EDB8BFE0-E714-6544-BE76-7DD28DD16902.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/F0D93AB0-289A-A440-80C8-B827A3885568.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/F559D91C-1D9E-FA4A-9F0D-326956029FBC.root',
        '/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v1/270000/F5CA488E-68D5-714B-9CBE-CDA1EE401C70.root',
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Path and EndPath definitions
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.pL1TkPrimaryVertex = cms.Path(process.L1TkPrimaryVertex)
process.pL1TkPhotonsCrystal = cms.Path(process.L1TkPhotonsCrystal)
process.pL1TkIsoElectronsCrystal = cms.Path(process.L1TkIsoElectronsCrystal)
process.pL1TkElectronsLooseCrystal = cms.Path(process.L1TkElectronsLooseCrystal)
process.pL1TkElectronsHGC = cms.Path(process.L1TkElectronsHGC)
process.pL1TkMuon = cms.Path(process.L1TkMuons+process.L1TkMuonsTP)
process.pL1TkElectronsLooseHGC = cms.Path(process.L1TkElectronsLooseHGC)
process.pL1TkElectronsEllipticMatchHGC = cms.Path(process.L1TkElectronsEllipticMatchHGC)
process.pL1TkElectronsCrystal = cms.Path(process.L1TkElectronsCrystal)
process.pL1TkPhotonsHGC = cms.Path(process.L1TkPhotonsHGC)
process.pL1TkIsoElectronsHGC = cms.Path(process.L1TkIsoElectronsHGC)
process.pL1TkElectronsEllipticMatchCrystal = cms.Path(process.L1TkElectronsEllipticMatchCrystal)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.customiseStoredTPConfig
from SimGeneral.MixingModule.customiseStoredTPConfig import higherPtTP

# call to customisation function higherPtTP imported from SimGeneral.MixingModule.customiseStoredTPConfig
process = higherPtTP(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000

# call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customisePhase2TTOn110
from L1Trigger.Configuration.customisePhase2TTOn110 import customisePhase2TTOn110

# call to customisation function customisePhase2TTOn110 imported from L1Trigger.Configuration.customisePhase2TTOn110
process = customisePhase2TTOn110(process)

# Automatic addition of the customisation function from L1TMuonSimulations.Configuration.customise
from L1TMuonSimulations.Configuration.customise import customise_ntupler_step,remove_L1TrackTrigger_step,remove_L1simulation_step

# call to customisation function customise_ntupler_step imported from L1TMuonSimulations.Configuration.customise
process = customise_ntupler_step(process)

# call to customisation function remove_L1TrackTrigger_step imported from L1TMuonSimulations.Configuration.customise
process = remove_L1TrackTrigger_step(process)

# call to customisation function remove_L1simulation_step imported from L1TMuonSimulations.Configuration.customise
process = remove_L1simulation_step(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
