import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def setupTauMVAVariables(process, runOnMC=False, path=None):
    process.MVATable = cms.EDProducer("TauMVAProducer",
                                      pfCandSrc  = cms.InputTag("packedPFCandidates"),
                                      jetsSrc    = cms.InputTag("linkedObjects","jets"),
                                      minCandPt  = cms.double(10.0),
                                      maxCandEta = cms.double(2.6),
                                      tauName    = cms.string("PFcand"),
                                      )

    process.TauMVATable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("MVATable"),
        cut = cms.string(""),
        name = cms.string("PFcand"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(True),
        variables = cms.PSet( P4Vars,
            dz    = Var("dz()", float, doc = "pfcand info dz", precision=8),
            fromPV= Var("fromPV()", float, doc = "pfcand info from Primary Vertex", precision=8),
        ),
    )

    process.TauMVATable.variables.pt.precision=10
    process.TauMVATable.variables.eta.precision=12
    process.TauMVATable.variables.phi.precision=10
    process.TauMVATable.variables.mass.precision=10
    process.TauMVATables = cms.Task(process.MVATable,
                                    process.TauMVATable)

    if path is None:
        process.schedule.associate(process.TauMVATables)
    else:
        getattr(process, path).associate(process.TauMVATables)

    return process
