import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *



##################### User floats producers, selectors ##########################


def setupCustomizedSB(process, runOnMC=False, path=None):
    ##################### Tables for final output and docs ##########################
    process.softbTable = cms.EDProducer("SoftbTableProducer",
        pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
        goodPvCut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
        svSrc = cms.InputTag("slimmedSecondaryVertices"),
        Jetsrc = cms.InputTag("updatedJets"),
        svCut = cms.string("pt<20 && numberOfDaughters >=3"),
        dlenMin = cms.double(0),
        dlenSigMin = cms.double(3),
        svName = cms.string("SB"),
        svDoc  = cms.string("secondary vertices from IVF algorithm"),
    )

    process.softbCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("softbTable"),
        cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
        name = cms.string("SB"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(True),
        variables = cms.PSet( P4Vars,
            ntracks= Var("numberOfDaughters()", int, doc = "number of tracks"),
            ndof   = Var("vertexNdof()", float, doc = "number of degrees of freedom",precision=8),
            chi2   = Var("vertexNormalizedChi2()", float, doc = "reduced chi2, i.e. chi/ndof",precision=8),
        ),
    )
    process.softbCandidateTable.variables.pt.precision=10
    process.softbCandidateTable.variables.phi.precision=12


    #after cross linkining
    process.softbTables = cms.Task( process.softbTable,
                                          process.softbCandidateTable)
    if path is None:
        process.schedule.associate(process.softbTables)
    else:
        getattr(process, path).associate(process.softbTables)

