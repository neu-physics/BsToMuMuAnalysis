import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonIsolationAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

isZmumu = True   # boolean for signal/background
#isZmumu = False  # boolean for signal/background

is200PU = True
#is200PU = False

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15) )


import FWCore.Utilities.FileUtils as FileUtils
sample = ''
nPU = ''
if isZmumu:
    sample = 'zmumu'
    nPU = '_200PU' if is200PU else '_NoPU'
else:
    sample = 'ttbar'
    nPU = '_200PU' if is200PU else '_NoPU'


#inputDatafileList = FileUtils.loadListFromFile ('MuonIsolationAnalyzer/data/filelists/{0}_files_xrd{1}.txt'.format(sample, nPU) ) 
    
process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                              'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19MiniAOD/BsToMuMu_TuneCP5_14TeV-pythia8/MINIAODSIM/NoPU_106X_upgrade2023_realistic_v3-v1/50000/99E994AB-739C-D947-9DDB-3BEB4A622120.root'
                              #'root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v2/90001/6F484AE2-E28D-3848-8A3E-F1367E69F4BC.root'
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_4_0_mtd5/RelValTTbar_Tauola_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/C95B6330-0D18-2749-95E5-6DA8A137C097.root'
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_4_0_mtd5/RelValDYToLL_M_50_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/F280AA95-BF6E-7145-9166-BA967B3067FE.root' # z->mumu
        #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_4_0_mtd5/RelValTTbar_Tauola_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/E7A00E21-D1DA-1F48-A14C-96D121FC2209.root' # ttbar -> stuff
                                #cms.untracked.vstring( *inputDatafileList)
                            )
                        )

#process.muonIsoAnalyzer = cms.EDAnalyzer('MuonIsolationAnalyzer')
process.load('BsToMuMuAnalysis.MuonIsolationAnalyzer.MuonIsolationAnalyzer_cfi')
muonIsoAnalyzer = process.MuonIsolationAnalyzer
muonIsoAnalyzer.isZmumuSignal = isZmumu

process.TFileService = cms.Service("TFileService", fileName = cms.string("muonIsolation_output_{0}{1}.root".format(sample, nPU)) )

process.runseq = cms.Sequence()
process.runseq += muonIsoAnalyzer
process.path = cms.Path(process.runseq)
process.schedule = cms.Schedule(process.path)

#process.p = cms.Path(process.muonIsoAnalyzer)
