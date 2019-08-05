import FWCore.ParameterSet.Config as cms
import sys, os


process = cms.Process("MuonIsolationAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

# opt 1: isZmumu, opt2: isPU, opt3: input file
if len(sys.argv) != 4:
    print('INPUTS INCORRECT. RECIEVED {0} . EXITING\n'.format( sys.argv ) )
    exit()

# boolean for signal/background
isZmumu = True if sys.argv[1] == 'isZmumu' else False  
# boolean for 0/200PU
is200PU = True if sys.argv[2] == 'is200PU' else False 

#make output file
os.system('touch tempInputFile.txt')
os.system('echo {0} tempInputFile.txt'.format(sys.argv[3]) )
filename = sys.argv[3]
lastFourFileID = filename[len(filename)-9:len(filename)-5]

#####################################################################################
#####                     Now start normal CMSSW stuff                          #####
#####################################################################################

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15) )


import FWCore.Utilities.FileUtils as FileUtils
sample = ''
nPU = ''
if isZmumu:
    sample = 'zmumu'
    nPU = '_PU' if is200PU else ''
else:
    sample = 'ttbar'
    nPU = '_PU' if is200PU else ''


#inputDatafileList = FileUtils.loadListFromFile ('MuonIsolationAnalyzer/data/filelists/{0}_files_xrd{1}.txt'.format(sample, nPU) ) 
inputDatafileList = FileUtils.loadListFromFile ('tempInputFile.txt' ) 
    
process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_4_0_mtd5/RelValDYToLL_M_50_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/F280AA95-BF6E-7145-9166-BA967B3067FE.root' # z->mumu
                                #'root://cms-xrd-global.cern.ch//store/relval/CMSSW_10_4_0_mtd5/RelValTTbar_Tauola_14TeV/GEN-SIM-RECO/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/E7A00E21-D1DA-1F48-A14C-96D121FC2209.root' # ttbar -> stuff
                                cms.untracked.vstring( *inputDatafileList)
                            )
                        )

#process.muonIsoAnalyzer = cms.EDAnalyzer('MuonIsolationAnalyzer')
process.load('BsToMuMuAnalysis.MuonIsolationAnalyzer.MuonIsolationAnalyzer_cfi')
muonIsoAnalyzer = process.MuonIsolationAnalyzer
muonIsoAnalyzer.isZmumuSignal = isZmumu

process.TFileService = cms.Service("TFileService", fileName = cms.string("muonIsolation_output_{0}{1}_{2}.root".format(sample, nPU, lastFourFileID)) )

process.runseq = cms.Sequence()
process.runseq += muonIsoAnalyzer
process.path = cms.Path(process.runseq)
process.schedule = cms.Schedule(process.path)

#process.p = cms.Path(process.muonIsoAnalyzer)
