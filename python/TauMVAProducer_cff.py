import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def setupTauMVAVariables(process):

	process.TauMVAProducer = cms.EDProducer(
		"TauMVAProducer",	
		pfCandSrc 		= cms.InputTag("packedPFCandidates"),
		jetsSrc 		= cms.InputTag("slimmedJets"),
		minCandPt       	= cms.double(10.0),
		maxCandEta      	= cms.double(2.6),
		variables		= cms.PSet(
			P4Vars,
		),
	)

	process.TauMVAProducer.variables.pt.precision=10
	process.TauMVAProducer.variables.eta.precision=10
	process.TauMVAProducer.variables.phi.precision=10
	process.TauMVAProducer.variables.mass.precision=10

        process.TauMVATask = cms.Task(process.TauMVAProducer)

        process.schedule.associate(process.TauMVATask)

        return process
