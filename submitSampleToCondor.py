# /usr/bin/python

#Author: Ben Tannenwald
#Date: August 5th, 2019
#Purpose: Script to submit condor jobs for all files in sample

import os,sys, argparse

# *** 0. setup parser for command line
parser = argparse.ArgumentParser()
parser.add_argument("--outputDir", help="output directory for processed histograms and roofiles")
parser.add_argument("--inputTXTfile", help=".txt file containing list of input files for a given sample")
args = parser.parse_args()

if (len(vars(args)) != 2): # 2 --> three for options
    os.system('python submitSampleToCondor.py -h')
    quit()

# ** A. Test output directory existence and create if DNE
if(args.outputDir is None):
    print "#### Need to specify output directory using --outputDir <desired output directory> ####\nEXITING"
    quit()
else:
    if ( not os.path.exists(args.outputDir) ):
        print "Specified input file ({0}) DNE.\nCREATING NOW".format(args.inputTXTfile)
        os.system("mkdir {0}".format(args.outputDir))

    if ( not os.path.exists( (args.outputDir + '/condor_logs/') ) ):
        os.system("mkdir {0}".format( (args.outputDir + '/condor_logs/')) )
    if ( not os.path.exists( (args.outputDir + '/condor_err/') ) ):
        os.system("mkdir {0}".format( (args.outputDir + '/condor_err/')) )
    if ( not os.path.exists( (args.outputDir + '/condor_out/') ) ):
        os.system("mkdir {0}".format( (args.outputDir + '/condor_out/')) )

    print '-- Setting outputDir = {0}'.format(args.outputDir)

# ** B. Test input .txt file and exit if DNE
if(args.inputTXTfile is None):
    print "#### Need to specify input .txt file using --inputTXTfile <address to .txt file> ####\nEXITING"
    quit()
else:
    if ( not os.path.exists(args.inputTXTfile) ):
        print "#### Specified input file ({0}) DNE ####.\nEXITING".format(args.inputTXTfile)
        quit()
    else:
        print '-- Setting inputTXTfile = {0}'.format(args.inputTXTfile)

# FIXME, 05-08-19 BBT
## ** C. Exit if no grid proxy
#if ( not os.path.exists(os.path.expandvars("$X509_USER_PROXY")) ):
#    print "#### No GRID PROXY detected. Please do voms-proxy-init -voms cms before submitting Condor jobs ####.\nEXITING"
#    quit()


# *** 1. Create .tar of directory and store in personal EOS
print "##########     Tarring workdir     ##########"
tarball_name = "{0}.tar.gz".format(args.outputDir)
os.system("cd /uscms_data/d2/benjtann/MTD_DPG/CMSSW_10_4_0_mtd5/src/; tar -cvzf {0} BsToMuMuAnalysis --exclude '.git' --exclude 'test_*' --exclude 'submitOneFile_' --exclude '*.tar.gz' --exclude '*-19_*' --exclude '*2019' --exclude 'pass*' --exclude '.SCRAM*' --exclude 'tmp' --exclude 'lib' --exclude 'config' --exclude 'external' --exclude '*.root'; cd /uscms_data/d2/benjtann/MTD_DPG/CMSSW_10_4_0_mtd5/src/BsToMuMuAnalysis ; mv /uscms_data/d2/benjtann/MTD_DPG/CMSSW_10_4_0_mtd5/src/{0} .".format(tarball_name))
if ( not os.path.exists("/eos/uscms/store/user/benjtann/{0}/".format(args.outputDir)) ):
    os.system("mkdir /eos/uscms/store/user/benjtann/{0}/".format(args.outputDir))

os.system("xrdcp {0} root://cmseos.fnal.gov//store/user/benjtann/{1}/".format(tarball_name, args.outputDir))

# *** 2. Create temporary .pdl file for condor submission
print "\n##########     Submitting Condor jobs     ##########\n"
txtfile = open(args.inputTXTfile, 'r')

for line in txtfile:
    infile = line.split('\n')[0]
    lastFourFileID = infile[len(infile)-9:len(infile)-5] # get last four digits of filename
    jdl_filename = "submitOneFile_{0}_{1}.jdl".format(args.outputDir, lastFourFileID)
    isZmumu = 'isZmumu' if 'zmumu' in args.inputTXTfile else 'isTTbar'
    is200PU = 'is200PU' if 'PU' in args.inputTXTfile else 'is0PU'

    os.system("touch {0}".format(jdl_filename))
    os.system("echo universe = vanilla > {0}".format(jdl_filename))
    os.system("echo Executable = cmsRunOneFile.csh >> {0}".format(jdl_filename))
    os.system("echo Should_Transfer_Files = YES >> {0}".format(jdl_filename))
    os.system("echo WhenToTransferOutput = ON_EXIT >> {0}".format(jdl_filename))
    os.system("echo Transfer_Input_Files = cmsRunOneFile.csh, MuonIsolationAnalyzer/test/singleFileMuonIsolation.py, {0} >> {1}".format(tarball_name, jdl_filename))
    #os.system("echo notify_user = benjamin.tannenwald@CERN.CH >> {0}".format(jdl_filename))
    #os.system("notify_user = benjtann@FNAL.GOV >> {0}".format(jdl_filename))
    os.system("echo Output = {0}/condor_out/outfile_{1}.out  >> {2}".format(args.outputDir, lastFourFileID, jdl_filename))
    os.system("echo Error = {0}/condor_err/outfile_{1}.err >> {2}".format(args.outputDir, lastFourFileID, jdl_filename))
    os.system("echo Log = {0}/condor_logs/outfile_{1}.log >> {2}".format(args.outputDir, lastFourFileID, jdl_filename))
    os.system("echo x509userproxy = ${{X509_USER_PROXY}} >> {0}".format(jdl_filename))
    os.system("echo Arguments = {0} {1} {2} {3} {4} >> {5}".format(isZmumu, is200PU, infile, tarball_name, args.outputDir, jdl_filename))
    os.system("echo Queue 1 >> {0}".format(jdl_filename))       
    os.system("condor_submit {0}".format(jdl_filename))


# *** 3. Cleanup submission directory
print "\n##########     Cleanup submission directory     ##########\n"
os.system("rm *.jdl")
