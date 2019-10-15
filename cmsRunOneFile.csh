#!/bin/bash

echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
export SCRAM_ARCH=slc7_amd64_gcc700
eval "scramv1 project CMSSW CMSSW_10_4_0_mtd5"
cd CMSSW_10_4_0_mtd5/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers           
echo "CMSSW: "$CMSSW_BASE                                               


xrdcp -s root://cmseos.fnal.gov//store/user/ali/$5/$4 .
tar -xf $4
rm $4

scram b clean
scram b
cd BsToMuMuAnalysis

#run actual command
cmsRun MuonIsolationAnalyzer/test/singleFileMuonIsolation.py $1 $2 $3 # isZmumu, is200PU, inputFile
       
### Now that the run is over, there is one or more root files created
echo "List all root files = "
#ls $1/*.root
ls ./*.root
echo "List all files"
ls 
echo "*******************************************"
OUTDIR=root://cmseos.fnal.gov//store/user/ali/$5/
echo "xrdcp output for condor"
for FILE in ./*.root
do
  echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
  fi
  rm ${FILE}
done

cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_4_0_mtd5

