import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoSUSY.softb_cff import setupCustomizedSB
from PhysicsTools.NanoSUSY.TauMVAProducer_cff import setupTauMVAVariables
from PhysicsTools.NanoSUSY.prodIsoTracksProducer_cff import setupprodIsoTracksVariables

def nanoSUSY_customizeCommon(process):
    setupCustomizedSB(process)
    setupTauMVAVariables(process)
    process.particleLevelSequence.remove(process.genParticles2HepMCHiggsVtx);
    process.particleLevelSequence.remove(process.rivetProducerHTXS);
    process.particleLevelTables.remove(process.HTXSCategoryTable)
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
    ## Needed to avoid segfault erros in the output moule when producing flat ntuple
    ## From https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
    process.add_(cms.Service("InitRootHandlers", EnableIMT = cms.untracked.bool(False)))
    return process

def nanoSUSY_customize80XLegacy(process):
    process = nanoSUSY_customizeCommon(process)
    ## Adding prodIsoTrack for 80Xlegacy
    setupprodIsoTracksVariables(process)
    return process

def nanoSUSY_customizeData(process):
    process = nanoSUSY_customizeCommon(process, False)
    process.NANOAODoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process


def nanoSUSY_customizeData_METMuEGClean(process):
    process = nanoSUSY_customizeCommon(process, False)

    from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
    corMETFromMuonAndEG(process,
                        pfCandCollection="",  # not needed
                        electronCollection="slimmedElectronsBeforeGSFix",
                        photonCollection="slimmedPhotonsBeforeGSFix",
                        corElectronCollection="slimmedElectrons",
                        corPhotonCollection="slimmedPhotons",
                        allMETEGCorrected=True,
                        muCorrection=False,
                        eGCorrection=True,
                        runOnMiniAOD=True,
                        postfix="MuEGClean"
                        )
    process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
    process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
    process.slimmedMETsMuEGClean.rawVariation = cms.InputTag("patPFMetRawMuEGClean")
    process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
    del process.slimmedMETsMuEGClean.caloMET
    process.metTable.src = cms.InputTag('slimmedMETsMuEGClean')

    process.NANOAODoutput.fakeNameForCrab = cms.untracked.bool(True)  # hack for crab publication
    return process


def nanoSUSY_customizeMC(process):
    process = nanoSUSY_customizeCommon(process, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
