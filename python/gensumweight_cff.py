import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *



##################### User floats producers, selectors ##########################


def setupCustomizedGenWeight(process, runOnMC=False, path=None):
    ##################### Tables for final output and docs ##########################
    process.GenWeightTable = cms.EDProducer("SignalSumGenWeight",
        genFilterInfoTag = cms.InputTag("genFilterEfficiencyProducer"),
        pvSrc      = cms.InputTag("offlineSlimmedPrimaryVertices"),
        pvName     = cms.string("SigWTab"),
        goodPvCut  = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
        svSrc      = cms.InputTag("slimmedSecondaryVertices"),
        Jetsrc     = cms.InputTag("linkedObjects","jets"),
        svCut      = cms.string("pt<20 && numberOfDaughters >=3"),
        svName     = cms.string("SB"),
        svDoc      = cms.string("secondary vertices from IVF algorithm"),
    )

    #after cross linkining
    process.GenWeightTables = cms.Task( process.GenWeightTable)
                                          #process.softbCandidateTable)
    if path is None:
        process.schedule.associate(process.GenWeightTables)
    else:
        getattr(process, path).associate(process.GenWeightTables)
    return process
