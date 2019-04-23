# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Commented_config -s NANO (type) --mc(Data)(Fast) --eventcontent NANOAODSIM(Needs to be this and different from central production) --datatier NANOAODSIM --filein gotten from DAS --no_exec --conditions auto:phase1_2017_realistic TAG -n 100events --era Run2_2017,run2_nanoAOD_94XMiniAODv1 --customise TopTagger/TopTagger/resolvedTagger_cff.customizeResolvedTaggerAllCanidiatesAndVariables Things we wrote --customise PhysicsTools/NanoSUSY/nanoSUSY_cff.nanoSUSY_customizeCommon
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2017,eras.run2_nanoAOD_94XMiniAODv1)
#NANO is the type of out put you want and era's can be found here: https://github.com/cms-sw/cmssw/tree/02d4198c0b6615287fd88e9a8ff650aea994412e/Configuration/Eras/python or you can check here with the campaign that was run: https://cms-pdmv.cern.ch/mcm/campaigns?prepid=RunIIFall17NanoAOD&page=0&shown=131135 in the sequence

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#These are general list of things run
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('Root file you want run can be found here: https://cmsweb.cern.ch/das/request?instance=prod/global&input=file+dataset%3D%2FST_tWnunu_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8%2FRunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3%2FMINIAODSIM Example: /store/mc/RunIIAutumn18MiniAOD/ST_tWnunu_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/260000/F9FC4C58-35D6-0C41-9465-7534C1141E2B.root '),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('This is the -p file.py nevts:-n number of events'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('Output name.root'),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)#This changes nanoSequenceFS or Data depending on what you are running over
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
#This also changes depending on if you are runninf FS or Data:
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

# Automatic addition of the customisation function from TopTagger.TopTagger.resolvedTagger_cff
from TopTagger.TopTagger.resolvedTagger_cff import customizeResolvedTaggerAllCanidiatesAndVariables 
#This is the code we wrote ourselves to add varaibles
#call to customisation function customizeResolvedTaggerAllCanidiatesAndVariables imported from TopTagger.TopTagger.resolvedTagger_cff
process = customizeResolvedTaggerAllCanidiatesAndVariables(process)

# Automatic addition of the customisation function from PhysicsTools.NanoSUSY.nanoSUSY_cff
from PhysicsTools.NanoSUSY.nanoSUSY_cff import nanoSUSY_customizeCommon 

#call to customisation function nanoSUSY_customizeCommon imported from PhysicsTools.NanoSUSY.nanoSUSY_cff
process = nanoSUSY_customizeCommon(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
