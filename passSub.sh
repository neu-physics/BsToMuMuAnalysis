python submitSampleToCondor.py --outputDir zmumu_0PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/zmumu_files_xrd_NoPU.txt
sleep 15
python submitSampleToCondor.py --outputDir ttbar_0PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/ttbar_files_xrd_NoPU.txt
sleep 15

python submitSampleToCondor.py --outputDir zmumu_200PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/zmumu_files_xrd_200PU.txt --nFilesPerJob 10
sleep 15
python submitSampleToCondor.py --outputDir ttbar_200PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/ttbar_files_xrd_200PU.txt --nFilesPerJob 10
sleep 15

python submitSampleToCondor.py --outputDir zmumu_MTDTDR_200PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/zmumu_files_xrd_MTDTDR_200PU.txt --nFilesPerJob 25
sleep 15
python submitSampleToCondor.py --outputDir ttbar_MTDTDR_200PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/ttbar_files_xrd_MTDTDR_200PU.txt --nFilesPerJob 25

python submitSampleToCondor.py --outputDir zmumu_MTDTDR_0PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/zmumu_files_xrd_MTDTDR_NoPU.txt --nFilesPerJob 25
sleep 15
python submitSampleToCondor.py --outputDir ttbar_MTDTDR_0PU_$1 --inputTXTfile MuonIsolationAnalyzer/data/filelists/ttbar_files_xrd_MTDTDR_NoPU.txt --nFilesPerJob 10
