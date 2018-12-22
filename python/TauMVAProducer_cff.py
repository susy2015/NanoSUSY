import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def setupTauMVAVariables(process):

	process.TauMVAProducer = cms.EDProducer(
		"TauMVAProducer",	
		pfCandSrc 		= cms.InputTag("packedPFCandidates"),
		jetsSrc 		= cms.InputTag("slimmedJets"),
		minCandPt       	= cms.double(10.0),
		maxCandEta      	= cms.double(2.6),
	)

	process.TauMVATable = cms.EDProducer("SimpleCandidateFlatTableProducer",
		src = cms.InputTag("packedPFCandidates"),
		cut = cms.string(""),
		name = cms.string("pfcands"),
		singleton = cms.bool(True),
		extention = cms.bool(True),
		variables = cms.PSet(
			P4Vars,
		)
	)

	process.TauMVATable.variables.pt.precision=10
	process.TauMVATable.variables.eta.precision=10
	process.TauMVATable.variables.phi.precision=10
	process.TauMVATable.variables.mass.precision=10

        process.TauMVATask = cms.Task(process.TauMVAProducer, process.TauMVATable)

        process.schedule.associate(process.TauMVATask)

        return process
