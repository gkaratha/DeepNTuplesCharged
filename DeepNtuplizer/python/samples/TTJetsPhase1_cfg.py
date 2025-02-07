import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/mc/PhaseIFall16MiniAOD/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PhaseIFall16PUFlat20to50_PhaseIFall16_81X_upgrade2017_realistic_v26-v1/90000/EC1E43BF-EFE2-E611-A9FF-901B0E542974.root' ] );


secFiles.extend( [
               ] )
