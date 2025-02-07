import FWCore.ParameterSet.Config as cms

genbanalizer = cms.EDAnalyzer('GenBNtuplizer',
                                pfChargedJets       = cms.InputTag("slimmedJetsPuppi"),
                                jets       = cms.InputTag("slimmedJets"),
                                secVertices = cms.InputTag("slimmedSecondaryVertices"),
                                genparts = cms.InputTag("prunedGenParticles"),
                                pfcands = cms.InputTag("packedPFCandidates"),
                                writeJetPart = cms.bool(True),
                                writePFcands = cms.bool(True),
                        
                                )
