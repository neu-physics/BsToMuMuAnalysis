# Studies to show improvements in Bs->mumu studies at CMS using MTD information

Code written in CMSSW_10_4_0_mtd5 using an SL7 architecture. You can also check out and compile code in SL6.

## Samples

For Z->mumu muons, /RelValDYToLL_M_50_14TeV/CMSSW_10_4_0_mtd5-103X_upgrade2023_realistic_v2_2023D35noPU-v1/GEN-SIM-RECO 

For ttbar "fake" muons, /RelValTTbar_Tauola_14TeV/CMSSW_10_4_0_mtd5-103X_upgrade2023_realistic_v2_2023D35noPU-v1/GEN-SIM-RECO 

## Running muon isolation studies

>$ cmsRun BsToMuMuAnalysis/MuonIsolationAnalyzer/test/produceROCCurve.py

## Using Condor Batch on LPC

>$ source submitSampleToCondor.py <date>_<version number>