python submitSampleToCondor.py --outputDir zmumu_0PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/zmumu_files_xrd.txt
sleep 15

python submitSampleToCondor.py --outputDir ttbar_0PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/ttbar_files_xrd.txt
sleep 15

python submitSampleToCondor.py --outputDir zmumu_200PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/zmumu_files_xrd_PU.txt
sleep 15

python submitSampleToCondor.py --outputDir ttbar_200PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/ttbar_files_xrd_PU.txt
