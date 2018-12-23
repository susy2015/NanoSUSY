import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def setupTauMVAVariables(process, runOnMC=False, path=None):

	process.TauMVA = cms.EDProducer("TauMVAProducer",	
		pfCandSrc 		= cms.InputTag("packedPFCandidates"),
		jetsSrc 		= cms.InputTag("linkedObjects","jets"),
		minCandPt       	= cms.double(10.0),
		maxCandEta      	= cms.double(2.6),
		tauName			= cms.string("pfcands"),
	)

	process.TauMVATable = cms.EDProducer("SimpleCandidateFlatTableProducer",
		src = cms.InputTag("TauMVA"),
		cut = cms.string(""),
		name = cms.string("pfcands"),
		singleton = cms.bool(False),
		extention = cms.bool(True),
		variables = cms.PSet( P4Vars,
		),
	)

	process.TauMVATable.variables.pt.precision=10
	process.TauMVATable.variables.eta.precision=10
	process.TauMVATable.variables.phi.precision=10
	process.TauMVATable.variables.mass.precision=10

        process.TauMVATask = cms.Task(process.TauMVA, process.TauMVATable)

	if path is None:
		process.schedule.associate(process.TauMVATask)
	else:
		getattr(process, path).associate(process.TauMVATask)

	return process
