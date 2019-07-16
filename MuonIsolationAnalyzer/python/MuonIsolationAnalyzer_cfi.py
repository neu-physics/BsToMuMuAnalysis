import FWCore.ParameterSet.Config as cms

MuonIsolationAnalyzer = cms.EDAnalyzer('MuonIsolationAnalyzer',
                                       genTag = cms.untracked.InputTag("genParticles"),
                                       #simHitsBTLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsBarrel"),
                                       #recHitsBTLTag = cms.untracked.InputTag("mtdRecHits:FTLBarrel"),
                                       #uncal_recHitsBTLTag = cms.untracked.InputTag("mtdUncalibratedRecHits:FTLBarrel"),
                                       #clustersBTLTag = cms.untracked.InputTag("mtdClusters:FTLBarrel"),
                                       #simHitsETLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsEndcap"),
                                       #recHitsETLTag = cms.untracked.InputTag("mtdRecHits:FTLEndcap"),
                                       #clustersETLTag = cms.untracked.InputTag("mtdClusters:FTLEndcap"),
                                       tracksTag = cms.untracked.InputTag("generalTracks"),
                                       #tracksLengthTag = cms.untracked.InputTag("trackExtenderWithMTD:pathLength"),
                                       #trackstmtdTag = cms.untracked.InputTag("trackExtenderWithMTD:tmtd"),
                                       #genVtxTag = cms.untracked.InputTag("g4SimHits"),
                                       #crysLayout = cms.untracked.int32(0),
                                       #track_hit_DRMax = cms.double(0.05),
                                       #track_hit_distMax = cms.double(99999.),
                                       #treeName = cms.untracked.string("DumpHits"),
                                       #verbosity = cms.bool(False),
                                       isZmumuSignal = cms.bool(False),
                                       #dumpRecHits = cms.bool(False)
                                       muonsTag = cms.untracked.InputTag("muons"),
                                       vertexTag = cms.untracked.InputTag("offlinePrimaryVertices")
                                   )
