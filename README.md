# NanoSUSY

This is an attempt for a modified NanoAOD, suited for the analysis need from
SUSY, mainly focus on the stop all-hadroinc searches.

*The code follows contributions from various sources as below:*
* NanoHRT from huilin: https://github.com/hqucms/NanoHRT
* METSig from Daniel: https://github.com/danbarto/cmssw/blob/sumPtForMETSig/PhysicsTools/NanoAOD/plugins/METSignificanceInputProducer.cc


### Set up CMSSW

```tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc700
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src/
cmsenv
```

### Get customized NanoAOD producers

```tcsh
git cms-merge-topic -u pastika:AddAxis1_1026
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone git@github.com:susy2015/TopTagger.git
git clone git@github.com:susy2015/NanoSUSY.git PhysicsTools/NanoSUSY
scram build -j 8
```

### Testing

```bash
cd PhysicsTools/NanoSUSY/crab
mkdir -p ${CMSSW_BASE}/src/TopTagger/TopTagger/data
getTaggerCfg.sh -o -n -t DeepResolved_DeepCSV_GR_noDisc_Release_v1.0.0 -d $CMSSW_BASE/src/TopTagger/TopTagger/data
```

##### 2016 Sample

For 2016 80X MC, '/\*/RunIISummer16MiniAODv2\*/MINIAODSIM'
``` 
cmsDriver.py prod2016MCv2 -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --filein /store/mc/RunIISummer16MiniAODv2/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/621E65CB-98BE-E611-BC32-00259029E7FC.root --no_exec --conditions 80X_mcRun2_asymptotic_2016_TrancheIV_v6 -n 100 --era Run2_2016,run2_miniAOD_80XLegacy --customise TopTagger/TopTagger/resolvedTagger_cff.customizeResolvedTaggerAllCanidiatesAndVariables --customise PhysicsTools/NanoSUSY/nanoSUSY_cff.nanoSUSY_customizeCommon
cmsRun prod2016MCv2_NANO.py
```

For 2016 80X Fastsim, '/\*/RunIISpring16MiniAODv\*/MINIAODSIM'
``` 
cmsDriver.py prod2016Fastv2 -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --filein /store/mc/RunIISpring16MiniAODv2/SMS-T2tt_mStop-150to250_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/10000/1C8E082A-9140-E611-8C60-008CFA1974D8.root --no_exec --conditions 80X_mcRun2_asymptotic_2016_miniAODv2_v0 -n 100 --era Run2_2016,run2_miniAOD_80XLegacy --customise TopTagger/TopTagger/resolvedTagger_cff.customizeResolvedTaggerAllCanidiatesAndVariables --customise PhysicsTools/NanoSUSY/nanoSUSY_cff.nanoSUSY_customizeCommon --fast
cmsRun prod2016Fastv2_NANO.py
```

For 2016 94X MC, '/\*/RunIISummer16MiniAODv3\*/MINIAODSIM'
``` 
cmsDriver.py prod2016MC -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --filein /store/mc/RunIISummer16MiniAODv3/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/1E4DDC30-2AEB-E811-BBB1-0CC47A5FA3BD.root --no_exec --conditions 94X_mcRun2_asymptotic_v3 -n 100 --era Run2_2016,run2_nanoAOD_94X2016 --customise TopTagger/TopTagger/resolvedTagger_cff.customizeResolvedTaggerAllCanidiatesAndVariables --customise PhysicsTools/NanoSUSY/nanoSUSY_cff.nanoSUSY_customizeCommon
cmsRun prod2016MC_NANO.py
```


##### 2017 Sample

For 2017 MC, '/\*/RunIIFall17\*12Apr2018\*/MINIAODSIM'
``` 
cmsDriver.py prod2017MC -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --filein /store/user/benwu/Stop18/NtupleSyncMiniAOD/00257B91-1808-E811-BD39-0242AC130002.root --no_exec --conditions auto:phase1_2017_realistic -n 100 --era Run2_2017,run2_nanoAOD_94XMiniAODv1 --customise TopTagger/TopTagger/resolvedTagger_cff.customizeResolvedTaggerAllCanidiatesAndVariables --customise PhysicsTools/NanoSUSY/nanoSUSY_cff.nanoSUSY_customizeCommon
cmsRun prod2017MC_NANO.py
```

### Submit crab jobs
```tcsh
source /cvmfs/cms.cern.ch/crab3/crab.csh
grid-proxy-init -debug -verify
python crab.py  -p prod2016MC.py -i example_16MC -o /store/user/benwu/StopStudy/TestNanoSUSY -t test
```

<details> <summary> To be updated: </summary>

Data:

```bash
cmsDriver.py test_nanoHRT_data -n 1000 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 94X_dataRun2_v4 --step NANO --nThreads 4 --era Run2_2016,run2_miniAOD_80XLegacy --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeData_METMuEGClean --filein /store/data/Run2016G/JetHT/MINIAOD/03Feb2017-v1/100000/006E7AF2-AEEC-E611-A88D-7845C4FC3B00.root --fileout file:nano_data.root >& test_data.log &

less +F test_data.log
```


### Production

**Step 0**: switch to the crab production directory and set up grid proxy, CRAB environment, etc.

```bash
cd $CMSSW_BASE/PhysicsTools/NanoHRT/crab
# set up grid proxy
voms-proxy-init -rfc -voms cms --valid 168:00
# set up CRAB env (must be done after cmsenv)
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

**Step 1**: generate the python config file with `cmsDriver.py` with the following commands:

MC (80X, MiniAODv2):

```bash
cmsDriver.py mc -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 94X_mcRun2_asymptotic_v2 --step NANO --nThreads 4 --era Run2_2016,run2_miniAOD_80XLegacy --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --filein file:step-1.root --fileout file:nano.root --no_exec
```

Data (`23Sep2016` ReReco):

```bash
cmsDriver.py data -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 94X_dataRun2_v4 --step NANO --nThreads 4 --era Run2_2016,run2_miniAOD_80XLegacy --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeData_METMuEGClean --filein file:step-1.root --fileout file:nano.root --no_exec
```

**Step 2**: use the `crab.py` script to submit the CRAB jobs:

For MC:

`python crab.py -p mc_NANO.py -o /store/group/lpcjme/noreplica/NanoHRT/mc/[version] -t NanoTuples-[version] -i mc_[ABC].txt --num-cores 4 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_mc_[ABC] --dryrun`

For data:

`python crab.py -p data_NANO.py -o /store/group/lpcjme/noreplica/NanoHRT/data/[version] -t NanoTuples-[version] -i data.txt --num-cores 4 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_data --dryrun`

A JSON file can be applied for data samples with the `-j` options. By default, we use the golden JSON for 2016:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
```

These command will perform a "dryrun" to print out the CRAB configuration files. Please check everything is correct (e.g., the output path, version number, requested number of cores, etc.) before submitting the actual jobs. To actually submit the jobs to CRAB, just remove the `--dryrun` option at the end.

**Step 3**: check job status

The status of the CRAB jobs can be checked with:

```bash
./crab.py --status --work-area crab_projects_[ABC]
```

Note that this will also resubmit failed jobs automatically.

The crab dashboard can also be used to get a quick overview of the job status:
`https://dashb-cms-job.cern.ch/dashboard/templates/task-analysis`

More options of this `crab.py` script can be found with:

```bash
./crab.py -h
```
</details>
