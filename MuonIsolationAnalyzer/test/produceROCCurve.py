import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonIsolationAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")


#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_4_0_mtd5/RelValTTbar_Tauola_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/E7A00E21-D1DA-1F48-A14C-96D121FC2209.root'
                                #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
                            )
                        )

#process.muonIsoAnalyzer = cms.EDAnalyzer('MuonIsolationAnalyzer')
process.load('BsToMuMuAnalysis.MuonIsolationAnalyzer.MuonIsolationAnalyzer_cfi')
muonIsoAnalyzer = process.MuonIsolationAnalyzer


process.runseq = cms.Sequence()
process.runseq += muonIsoAnalyzer
process.path = cms.Path(process.runseq)
process.schedule = cms.Schedule(process.path)

#process.p = cms.Path(process.muonIsoAnalyzer)
