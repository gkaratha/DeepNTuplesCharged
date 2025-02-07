
import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
### parsing job options 
import sys

options = VarParsing.VarParsing()

options.register('inputScript','',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"input Script")
options.register('outputFile','outputSelectedPat',VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string,"output File (w/o .root)")
options.register('maxEvents', 1001,VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int,"maximum events")
options.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "total jobs")
options.register('reportEvery', 1000, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, "report every")
options.register('gluonReduction', 0.0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.float, "gluon reduction")
options.register('selectJets', True, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "select jets with good gen match")
options.register('phase2', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "apply jet selection for phase 2. Currently sets JetEtaMax to 3.0 and picks slimmedJetsPuppi as jet collection.")
options.register('puppi', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use puppi jets")
options.register('eta', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "use eta up to 5.0")


import os
release=os.environ['CMSSW_VERSION'][6:]
print("Using release "+release)


#inputFiles =[ 'file:/eos/cms/store/cmst3/group/softJets/gkaratha/chain_m70_dm20_cfgRun24_133X_Run2024_test_10172024/Mini/chain_m70_dm20_'+str(i)+'_step5_mini.root' for i in range(1,250)]

inputFiles = ['/store/mc/RunIII2024Summer24MiniAOD/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v2/100000/001b1bf1-e811-4ec1-94ff-fa8a45fb629f.root']

if hasattr(sys, "argv"):
    options.parseArguments()

UsePuppiReclusterForStdJet=False
UsePFReclusterForStdJet=False
UseSlimmedForStdJet=True
UsePuppiForTrkJet=False
UsePFRForTrkJet=False


if (not UsePuppiReclusterForStdJet) and (not UsePFReclusterForStdJet) and (not UseSlimmedForStdJet):
   print("provide std jet")
   exit()

if (UsePuppiReclusterForStdJet + UsePFReclusterForStdJet + UseSlimmedForStdJet)>1:
   print("too many std jet")
   exit()

process = cms.Process("DNNFiller")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_postBPix_v2', '')


process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),  
   wantSummary=cms.untracked.bool(False)
)


process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(inputFiles), secondaryFileNames = cms.untracked.vstring())
process.source.skipEvents = cms.untracked.uint32(options.skipEvents)
process.maxEvents  = cms.untracked.PSet( 
    input = cms.untracked.int32 (options.maxEvents) 
)
releases = release.split("_")


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
    
###########################################################################
#
# Setup puppi modules and set them to recalculate weights
#
###########################################################################
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD(process, False)
process.puppi.useExistingWeights = False
process.puppiNoLep.useExistingWeights = False

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection


###########################################################################
#
# Make function wrapper around PatAlgos helper functions
#
###########################################################################
def addProcessAndTask(proc, label, module):
  task = getPatAlgosToolsTask(proc)
  addToProcessAndTask(label, module, proc, task)

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJetsPuppi
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets




###############################################################################
####################### Reclustered Standard AK4 jets ##########################
###############################################################################

############################## Puppi recluster jets ###########################
jetCorrectionsAK4 = ('AK4PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
if UsePuppiReclusterForStdJet:
  addProcessAndTask(process, "ak4PFJetsPuppiRecluster", ak4PFJetsPuppi.clone(
        src = "packedPFCandidates",
        srcWeights = "puppi",
        doAreaFastjet = True,
        jetPtMin=5
     )
  )

  addJetCollection(
      process,
      postfix            = "Recluster",
      labelName          = "AK4Puppi",
      jetSource          = cms.InputTag("ak4PFJetsPuppiRecluster"),
      algo               = "AK", #name of algo must be in this format
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
  #process.selectedPatJetsAK4PuppiRecluster.cut="pt > 10"
  updateJetCollection(
        process,
        labelName = "AK4PuppiR",
        jetSource = cms.InputTag("selectedPatJetsAK4PuppiRecluster"),  # 'ak4Jets'
        jetCorrections = jetCorrectionsAK4,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = None,
        btagDiscriminators = None,
        explicitJTA = False
  )


################################# PF recluster jets #########################
if UsePFReclusterForStdJet:
  addProcessAndTask(process, "ak4PFReclusterJets", ak4PFJets.clone(
          src = "packedPFCandidates",
          jetPtMin=5,
          doAreaFastjet = True,
        )
  )
  
  addJetCollection(
        process,
        labelName          = "AK4PFRecluster",
        jetSource          = cms.InputTag("ak4PFReclusterJets"),
        algo               = "ak", #name of algo must be in this format
        rParam             = 0.4,
        pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
        pfCandidates       = cms.InputTag("packedPFCandidates"),
        svSource           = cms.InputTag("slimmedSecondaryVertices"),
        muSource           = cms.InputTag("slimmedMuons"),
        elSource           = cms.InputTag("slimmedElectrons"),
        genJetCollection   = cms.InputTag("ak4GenJetsRecluster"),
        genParticles       = cms.InputTag("prunedGenParticles"),
        jetCorrections     = None ,
  )
  
  updateJetCollection(
        process,
        labelName = "AK4PFR",
        jetSource = cms.InputTag("patJetsAK4PFRecluster"),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        jetCorrections = None,
        btagDiscriminators = None,
        btagInfos = None,
        explicitJTA = False
  )


#############################################################################
############################## trackjet #####################################
#############################################################################
#### filter charged PF cands
addProcessAndTask(process, "packedPFCandidatesChg",cms.EDFilter("CandPtrSelector",
     src = cms.InputTag("packedPFCandidates"),
     cut = cms.string("charge != 0")
     )
)

############################## PF trkjet ###################################
if UsePFRForTrkJet:
   addProcessAndTask(process, "ak4PFChgJets", ak4PFJets.clone(
            src = "packedPFCandidatesChg",
            jetPtMin=5,
            doAreaFastjet = True
            )
   )

   addJetCollection(
      process,
      labelName          = "AK4PFChg",
      jetSource          = cms.InputTag("ak4PFChgJets"),
      algo               = "ak", #name of algo must be in this format
      rParam             = 0.4,
      pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidates       = cms.InputTag("packedPFCandidates"),
      svSource           = cms.InputTag("slimmedSecondaryVertices"),
      muSource           = cms.InputTag("slimmedMuons"),
      elSource           = cms.InputTag("slimmedElectrons"),
      genJetCollection   = cms.InputTag("ak4GenJetsRecluster"),
      genParticles       = cms.InputTag("prunedGenParticles"),
      jetCorrections     = None ,
   )

   updateJetCollection(
        process,
        labelName = "AK4PFChgR",
        jetSource = cms.InputTag("patJetsAK4PFChg"),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        jetCorrections = None,
        btagDiscriminators = None,
        btagInfos = None,
        explicitJTA = False
   )

################################# Puppi trkjet ###############################
if UsePuppiForTrkJet:
  addProcessAndTask(process, "ak4PuppiChgJets", ak4PFJetsPuppi.clone(
        src = "packedPFCandidatesChg",
        srcWeights = "puppi",
        doAreaFastjet = True,
        jetPtMin=5
     )
  )

  addJetCollection(
      process,
      labelName          = "AK4PuppiChgJets",
      jetSource          = cms.InputTag("ak4PFJetsPuppiChgJets"),
      algo               = "AK", #name of algo must be in this format
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
  process.selectedPatJetsAK4PuppiChgJets.getJetMCFlavour = True
  getattr(process, "selectedPatJetFlavourAssociationAK4PuppiChgJets").weights = cms.InputTag("puppi")

  updateJetCollection(
        process,
        labelName = "AK4PuppiR",
        jetSource = cms.InputTag("patJetsAK4PuppiChgJets"),  # 'ak4Jets'
        jetCorrections = jetCorrectionsAK4,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = None,
        btagDiscriminators = None,
        explicitJTA = False
  )


   addJetCollection(
      process,
      labelName          = "AK4PFChg",
      jetSource          = cms.InputTag("ak4PFChgJets"),
      algo               = "ak", #name of algo must be in this format
      rParam             = 0.4,
      pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidates       = cms.InputTag("packedPFCandidates"),
      svSource           = cms.InputTag("slimmedSecondaryVertices"),
      muSource           = cms.InputTag("slimmedMuons"),
      elSource           = cms.InputTag("slimmedElectrons"),
      genJetCollection   = cms.InputTag("ak4GenJetsRecluster"),
      genParticles       = cms.InputTag("prunedGenParticles"),
      jetCorrections     = None ,
   )

   updateJetCollection(
        process,
        labelName = "AK4PFChgR",
        jetSource = cms.InputTag("patJetsAK4PFChg"),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        jetCorrections = None,
        btagDiscriminators = None,
        btagInfos = None,
        explicitJTA = False
   )



if UsePFRForTrkJet:
   options.outputFile+="PFCvs"

if UseSlimmedForStdJet:
   std_jet_collection = "slimmedJets"
   options.outputFile+="Slimmed"
if UsePFReclusterForStdJet:
   std_jet_collection = "selectedUpdatedPatJetsAK4PFR"
   options.outputFile+="PFRecluster"
if UsePuppiReclusterForStdJet:
   std_jet_collection = "selectedUpdatedPatJetsAK4PuppiR"
   options.outputFile+="PuppiRecluster"


from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsWithNu = ak4GenJets.clone(src ='packedGenParticles')

 ## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16 && charge!=0"))
 ## Define GenJets
process.ak4GenJetsRecluster = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')


process.genJetReclusterTask = cms.Task(process.packedGenParticlesForJetsNoNu,process.ak4GenJetsWithNu,process.ak4GenJetsRecluster)

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

# GenBanalyzer
process.load("DeepNTuples.DeepNtuplizer.GenBanalysis_cfi")
process.genbanalizer.pfChargedJets = cms.InputTag('selectedUpdatedPatJetsAK4PFChgR')
process.genbanalizer.jets = cms.InputTag(std_jet_collection)
process.genbanalizer.writeJetPart = cms.bool(True)
process.genbanalizer.writePFcands = cms.bool(True)


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
    process.patAlgosToolsTask,
)
process.ep = cms.EndPath(process.genbanalizer)

