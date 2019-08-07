#!/bin/tcsh

echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.csh  ## if a bash script, use .sh instead of .csh

xrdcp -s root://cmseos.fnal.gov//store/user/benjtann/$5/$4 .
tar -xf $4
rm $4

#setenv SCRAM_ARCH slc7_amd64_gcc700
setenv SCRAM_ARCH slc6_amd64_gcc700
cd CMSSW_10_4_0_mtd5/src/
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers

#source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
#cd ${_CONDOR_SCRATCH_DIR}
#xrdcp -s root://cmseos.fnal.gov//store/user/benjtann/condor_tarballs/$4 .
#tar -xf $4
#rm $4
#root -l -q -b 'trigEffStudy_2017data.C("'${1}'","'${2}'","'${3}'")'

cmsRun MuonIsolationAnalyzer/test/singleFileMuonIsolation.py $1 $2 $3 # isZmumu, is200PU, inputFile

### Now that the run is over, there is one or more root files created
echo "List all root files = "
#ls $1/*.root
ls ./*.root
echo "List all files"
ls 
echo "*******************************************"
#OUTDIR=root://cmseos.fnal.gov//store/user/benjtann/$1/
OUTDIR=root://cmseos.fnal.gov//store/user/benjtann/$5/
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

