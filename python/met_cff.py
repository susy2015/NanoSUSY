import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def setupMETSig(process):
    process.metSigTable = cms.EDProducer("METSignificanceInputProducer",
        srcPFCandidates = cms.InputTag('packedPFCandidates'),
        srcLeptons      = cms.VInputTag('slimmedElectrons','slimmedMuons','slimmedPhotons'), # photons added here because of good resolution
        srcPfJets       = cms.InputTag('slimmedJets'),#instead of slimmed
        METSigParams    = cms.VPSet(
            cms.PSet( name = cms.string("dRMatch"), dRMatch = cms.double(0.4), jetThreshold = cms.double(15.) ),
            ),
    )
