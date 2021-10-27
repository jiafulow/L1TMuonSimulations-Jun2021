import FWCore.ParameterSet.Config as cms

# Reference:
#     https://github.com/cms-sw/cmssw/blob/master/IOMC/EventVertexGenerators/python/VtxSmearedFlat_cfi.py

from IOMC.EventVertexGenerators.VtxSmearedParameters_cfi import FlatVtxSmearingParameters,VtxSmearedCommon
VtxSmeared = cms.EDProducer("FlatEvtVtxGenerator2",
    FlatVtxSmearingParameters,
    VtxSmearedCommon
)

# Important note: flat independent distributions in Z and T are not correct for physics production
# In reality, if two flat beams interact the real distribution will not be flat with independent Z and T
# but Z and T will be correlated, as example in GaussEvtVtxGenerator.
# Can restore correlation via MinT += (MinZ - MaxZ)/2 and MaxT += (MaxZ - MinZ)/2
# in [ns] units (recall c_light = 29.98cm/ns)
VtxSmeared.MaxX = cms.double(150)
VtxSmeared.MaxY = cms.double(150)
VtxSmeared.MaxZ = cms.double(30)
VtxSmeared.MaxT = cms.double(30 / 29.9792458)
VtxSmeared.MinX = cms.double(-1 * VtxSmeared.MaxX.value())
VtxSmeared.MinY = cms.double(-1 * VtxSmeared.MaxY.value())
VtxSmeared.MinZ = cms.double(-1 * VtxSmeared.MaxZ.value())
VtxSmeared.MinT = cms.double(-1 * VtxSmeared.MaxT.value())
