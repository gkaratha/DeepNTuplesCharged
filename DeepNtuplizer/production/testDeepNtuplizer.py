
import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
### parsing job options 
import sys

options = VarParsing.VarParsing()

options.register('inputScript','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"input Script")
options.register('outputFile','output',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('maxEvents', 50001,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "total jobs")
options.register('reportEvery', 1000, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "report every")
options.register('gluonReduction', 0.0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "gluon reduction")
options.register('selectJets', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "select jets with good gen match")
options.register('phase2', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "apply jet selection for phase 2. Currently sets JetEtaMax to 3.0 and picks slimmedJetsPuppi as jet collection.")
options.register('puppi', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use puppi jets")
options.register('eta', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use eta up to 5.0")


import os
release=os.environ['CMSSW_VERSION'][6:]
print("Using release "+release)


options.register(
	'inputFiles','',
	VarParsing.VarParsing.multiplicity.list,
	VarParsing.VarParsing.varType.string,
	"input files (default is the tt RelVal)"
	)

if hasattr(sys, "argv"):
    options.parseArguments()

if options.puppi:
    usePuppi = True
else:
    usePuppi = False
process = cms.Process("DNNFiller")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#'auto:run2_mc'
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_postBPix_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),  
   wantSummary=cms.untracked.bool(False)
)

process.load('DeepNTuples.DeepNtuplizer.samples.TTJetsPhase1_cfg') #default input

if options.inputFiles:
	process.source.fileNames = options.inputFiles

if options.inputScript != '' and options.inputScript != 'DeepNTuples.DeepNtuplizer.samples.TTJetsPhase1_cfg':
    process.load(options.inputScript)

numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job


process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
if options.nJobs > 1:
    print ("running over these files:")
    print (process.source.fileNames)

process.source.fileNames = [ 'file:/eos/cms/store/cmst3/group/softJets/gkaratha/chain_m70_dm20_cfgRun24_133X_Run2024_test_10172024/Mini/chain_m70_dm20_'+str(i)+'_step5_mini.root' for i in range(1,50)]

process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents  = cms.untracked.PSet( 
    input = cms.untracked.int32 (options.maxEvents) 
)
releases = release.split("_")

bTagInfos = ['pfDeepFlavourTagInfos',
             'pfImpactParameterTagInfos',
             'pfInclusiveSecondaryVertexFinderTagInfos',
             'pfParticleNetAK4TagInfos',] #['pfParticleTransformerAK4TagInfos',]

from RecoBTag.ONNXRuntime.pfParticleNetAK4_cff import _pfParticleNetAK4JetTagsAll as pfParticleNetAK4JetTagsAll
from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import _pfParticleNetFromMiniAODAK4PuppiCentralJetTagsProbs
from RecoBTag.ONNXRuntime.pfUnifiedParticleTransformerAK4_cff import _pfUnifiedParticleTransformerAK4JetTagsAll

if (int(releases[0])>8) or ( (int(releases[0])==8) and (int(releases[1]) >= 4) ) :
 bTagDiscriminators = [
     'pfDeepCSVJetTags:probudsg', #to be fixed with new names
     'pfDeepCSVJetTags:probb',
     'pfDeepCSVJetTags:probc',
     'pfDeepCSVJetTags:probbb',
     'pfDeepFlavourJetTags:probb',
     'pfDeepFlavourJetTags:probbb',
     'pfDeepFlavourJetTags:problepb',
     'pfDeepFlavourJetTags:probc',
     'pfDeepFlavourJetTags:probuds',
     'pfDeepFlavourJetTags:probg',
     'pfParticleTransformerAK4JetTags:probb',
     'pfParticleTransformerAK4JetTags:probbb',
     'pfParticleTransformerAK4JetTags:problepb',
     'pfParticleTransformerAK4JetTags:probc',
     'pfParticleTransformerAK4JetTags:probuds',
     'pfParticleTransformerAK4JetTags:probg',
 ] + _pfParticleNetFromMiniAODAK4PuppiCentralJetTagsProbs + pfParticleNetAK4JetTagsAll + _pfUnifiedParticleTransformerAK4JetTagsAll
else :
  bTagDiscriminators = [
      'pfDeepCSVJetTags:probudsg', #to be fixed with new names
      'pfDeepCSVJetTags:probb',
      'pfDeepCSVJetTags:probc',
      'pfDeepCSVJetTags:probbb',
      'pfDeepCSVJetTags:probcc',
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
      'pfDeepFlavourJetTags:probc',
      'pfDeepFlavourJetTags:probuds',
      'pfDeepFlavourJetTags:probg',
      'pfParticleTransformerAK4JetTags:probb',
      'pfParticleTransformerAK4JetTags:probbb',
      'pfParticleTransformerAK4JetTags:problepb',
      'pfParticleTransformerAK4JetTags:probc',
      'pfParticleTransformerAK4JetTags:probuds',
      'pfParticleTransformerAK4JetTags:probg',
 ] + _pfParticleNetFromMiniAODAK4PuppiCentralJetTagsProbs + pfParticleNetAK4JetTagsAll + _pfUnifiedParticleTransformerAK4JetTagsAll

jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if options.phase2:
    usePuppi = True

if usePuppi:
    jetCorrectionsAK4 = ('AK4PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    
###########################################################################
#
# Setup puppi modules and set them to recalculate weights
#
###########################################################################
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD(process, False)
process.puppi.useExistingWeights = False
process.puppiNoLep.useExistingWeights = False

###########################################################################
#
# Make function wrapper around PatAlgos helper functions
#
###########################################################################
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
def addProcessAndTask(proc, label, module):
  task = getPatAlgosToolsTask(proc)
  addToProcessAndTask(label, module, proc, task)

###########################################################################
#
# Recluster AK4 Puppi jets
#
###########################################################################
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJetsPuppi
jetCollectionRecluster = "ak4PFJetsPuppiRecluster"
addProcessAndTask(process, jetCollectionRecluster, ak4PFJetsPuppi.clone(
      src = "packedPFCandidates",
      srcWeights = "puppi"
    )
)

###########################################################################
#
# Patify reclustered AK4 Puppi jets
#
###########################################################################
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(
  process,
  postfix            = "Recluster",
  labelName          = "AK4Puppi",
  jetSource          = cms.InputTag(jetCollectionRecluster),
  algo               = "ak", #name of algo must be in this format
  rParam             = 0.4,
  pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
  pfCandidates       = cms.InputTag("packedPFCandidates"),
  svSource           = cms.InputTag("slimmedSecondaryVertices"),
  muSource           = cms.InputTag("slimmedMuons"),
  elSource           = cms.InputTag("slimmedElectrons"),
  genJetCollection   = cms.InputTag("ak4GenJetsRecluster"), # This is setup below
  genParticles       = cms.InputTag("prunedGenParticles"),
  jetCorrections     = jetCorrectionsAK4,
)

process.patJetsAK4PuppiRecluster.getJetMCFlavour = True
getattr(process, "patJetFlavourAssociationAK4PuppiRecluster").weights = cms.InputTag("puppi")

if usePuppi:
    jet_collection = 'patJetsAK4PuppiRecluster'
else:
    jet_collection = 'slimmedJets'

updateJetCollection(
        process,
        labelName = "DeepFlavour",
        jetSource = cms.InputTag(jet_collection),  # 'ak4Jets'
        jetCorrections = jetCorrectionsAK4,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
)

if hasattr(process,'updatedPatJetsTransientCorrectedDeepFlavour'):
  process.updatedPatJetsTransientCorrectedDeepFlavour.addTagInfos = cms.bool(True) 
  process.updatedPatJetsTransientCorrectedDeepFlavour.addBTagInfo = cms.bool(True)
#else:
#  raise ValueError('I could not find updatedPatJetsTransientCorrectedDeepFlavour to embed the tagInfos, please check the cfg')


# Very Loose IVF SV collection
from PhysicsTools.PatAlgos.tools.helpers import loadWithPrefix
loadWithPrefix(process, 'RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff', "looseIVF")
process.looseIVFinclusiveCandidateVertexFinder.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFinclusiveCandidateVertexFinder.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLen2DSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.vertexMinDLenSig = cms.double(0.)
process.looseIVFinclusiveCandidateVertexFinder.fitterSigmacut = 20

process.looseIVFcandidateVertexArbitrator.primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.looseIVFcandidateVertexArbitrator.tracks = cms.InputTag("packedPFCandidates")
process.looseIVFcandidateVertexArbitrator.secondaryVertices = cms.InputTag("looseIVFcandidateVertexMerger")
process.looseIVFcandidateVertexArbitrator.fitterSigmacut = 20


outFileName = options.outputFile + '_' + str(options.job) +  '.root'
print ('Using output file ' + outFileName)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(outFileName))

process.load("DeepNTuples.DeepNtuplizer.GenBanalysis_cfi")
process.genbanalizer.pfChargedJets = cms.InputTag('selectedUpdatedPatJetsDeepFlavour')
process.genbanalizer.jets = cms.InputTag("slimmedJets")

#1631
process.ProfilerService = cms.Service (
      "ProfilerService",
       firstEvent = cms.untracked.int32(1631),
       lastEvent = cms.untracked.int32(1641),
       paths = cms.untracked.vstring('p') 
)

#Trick to make it work in 9_1_X
process.tsk = cms.Task()
for mod in process.producers_().values(): #.itervalues():
    process.tsk.add(mod)
for mod in process.filters_().values(): #.itervalues():
    process.tsk.add(mod)

process.patAlgosToolsTask = getPatAlgosToolsTask(process)

process.p = cms.Path(
    process.tsk, 
    process.patAlgosToolsTask 
)

process.ep = cms.EndPath(process.genbanalizer)
#process.schedule = cms.Schedule(process.p, process.ep)
