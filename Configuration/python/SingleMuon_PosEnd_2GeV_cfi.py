import FWCore.ParameterSet.Config as cms

# Reference:
#     https://github.com/cms-sw/genproductions/blob/master/genfragments/EightTeV/SingleMuMinusFlatPt0p2To100_cff.py

generator = cms.EDProducer("FlatRandomPtGunProducer2",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(7000.0),
        MinPt = cms.double(2.0),
        PartID = cms.vint32(-13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.141592653589793),
        MinEta = cms.double(1.2),
        MinPhi = cms.double(-3.141592653589793),
        RandomCharge = cms.bool(True),
        PtSpectrum = cms.string('flatOneOverPt'),
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single muon+/- pt 2 to 7000 flat in 1/pt positive endcap'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
