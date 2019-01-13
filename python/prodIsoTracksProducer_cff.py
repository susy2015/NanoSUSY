import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def setupprodIsoTracksVariables(process, runOnMC=False, path=None):
    process.prodIsoTracksTable = cms.EDProducer("prodIsoTracksProducer",
                                      vtxSrc            = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      metSrc            = cms.InputTag('slimmedMETs'),
                                      forVetoIsoTrkSrc  = cms.InputTag("trackIsolation"),
                                      loose_isoTrkSrc   = cms.InputTag("packedPFCandidates"),
                                      exclPdgIdVec      = cms.vint32(),
                                      dR_ConeSize       = cms.double(0.3),
                                      dz_CutValue       = cms.double(0.1),
                                      minPt_PFCandidate = cms.double(5.0),
                                      isoCut            = cms.double(0.5),
                                      debug             = cms.bool(False),
                                      )

    process.prodIsoTracksTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("prodIsoTracksTable"),
        cut = cms.string(""),
        name = cms.string("IsoTracks"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(True),
        variables = cms.PSet( P4Vars,
            q     = Var("charge()", float, doc = "isoTrack info charge", precision=8),
            dz    = Var("dz()", float, doc = "isoTrack info dz", precision=8),
            pdgId = Var("pdgId()", int, doc = "isoTrack info pdgId"),
            fromPV= Var("fromPV()", float, doc = "isoTrack info from Primary Vertex", precision=8),
        ),
    )

    process.prodIsoTracksTable.variables.pt.precision=10
    process.prodIsoTracksTable.variables.eta.precision=12
    process.prodIsoTracksTable.variables.phi.precision=10
    process.prodIsoTracksTable.variables.mass.precision=10
    process.prodIsoTracksTables = cms.Task(process.prodIsoTracksTable,
                                    process.prodIsoTracksTable)

    if path is None:
        process.schedule.associate(process.prodIsoTracksTables)
    else:
        getattr(process, path).associate(process.prodIsoTracksTables)

    return process
