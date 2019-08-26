# /usr/bin/python

#Author: Ben Tannenwald
#Date: August 5th, 2019
#Purpose: Script to submit condor jobs for all files in sample

import os,sys, argparse

# *** 0. setup parser for command line
parser = argparse.ArgumentParser()
parser.add_argument("--outputDir", help="output directory for processed histograms and roofiles")
parser.add_argument("--inputTXTfile", help=".txt file containing list of input files for a given sample")
parser.add_argument("--nFilesPerJob", help="number of files to run per job... use if lots of files in a dataset e.g. MTDTDR samples", default='1')
args = parser.parse_args()

if ( not (len(vars(args)) != 2 or len(vars(args)) != 3) ): # 2 OR 3 --> three for options
    os.system('python submitSampleToCondor.py -h')
    quit()

# ** A. Test output directory existence and create if DNE
if(args.outputDir is None):
    print "#### Need to specify output directory using --outputDir <desired output directory> ####\nEXITING"
    quit()
else:
    if ( not os.path.exists(args.outputDir) ):
        print "Specified output directory ({0}) DNE.\nCREATING NOW".format(args.inputTXTfile)
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
    print "#### Need to specify input .txt file using --inputTXTfile <address to .txt file> ####\nEXITING\n"
    quit()
else:
    if ( not os.path.exists(args.inputTXTfile) ):
        print "#### Specified input file ({0}) DNE ####.\nEXITING\n".format(args.inputTXTfile)
        quit()
    else:
        print '-- Setting inputTXTfile = {0}'.format(args.inputTXTfile)

# ** C. Test nFilesPerJob flag and exit if not sensible
if (args.nFilesPerJob).isdigit() is True:
    print "-- Running with {0} files per job ####\n".format(args.nFilesPerJob)
else:
    print "#### WARNING: Passed option of {0} files per job makes no sense. DNE ####\nEXITING\n".format(args.nFilesPerJob)
    quit()

# FIXME, 05-08-19 BBT
## ** D. Exit if no grid proxy
#if ( not os.path.exists(os.path.expandvars("$X509_USER_PROXY")) ):
#    print "#### No GRID PROXY detected. Please do voms-proxy-init -voms cms before submitting Condor jobs ####.\nEXITING"
#    quit()


# *** 1. Create .tar of directory and store in personal EOS
print "##########     Tarring workdir     ##########"
tarball_name = "{0}.tar.gz".format(args.outputDir)
os.system("cd /uscms_data/d2/benjtann/MTD_DPG/CMSSW_10_4_0_mtd5/src/; tar -cvzf {0} BsToMuMuAnalysis --exclude '.git' --exclude 'test_*' --exclude 'submitOneFile_' --exclude '*.tar.gz' --exclude '*-19_*' --exclude '*2019' --exclude 'pass*' --exclude '.SCRAM*' --exclude 'tmp' --exclude 'lib' --exclude 'config' --exclude 'external' --exclude '*notes_*.txt*' --exclude '08-*19' --exclude '*.root'; cd /uscms_data/d2/benjtann/MTD_DPG/CMSSW_10_4_0_mtd5/src/BsToMuMuAnalysis ; mv /uscms_data/d2/benjtann/MTD_DPG/CMSSW_10_4_0_mtd5/src/{0} .".format(tarball_name))
if ( not os.path.exists("/eos/uscms/store/user/benjtann/{0}/".format(args.outputDir)) ):
    os.system("mkdir /eos/uscms/store/user/benjtann/{0}/".format(args.outputDir))

os.system("xrdcp {0} root://cmseos.fnal.gov//store/user/benjtann/{1}/".format(tarball_name, args.outputDir))

# *** 2. Create temporary .pdl file for condor submission
print "\n##########     Submitting Condor jobs     ##########\n"

with open(args.inputTXTfile, 'r') as txtfile:
    line1 = txtfile.readline().split('\n')[0]
    while(line1):
        filesForJob = [line1]
        # * i. Get list of correct number of files
        while len(filesForJob) < int(args.nFilesPerJob):
            lineN = txtfile.readline().split('\n')[0]
            filesForJob.append(lineN)
        line1 = txtfile.readline().split('\n')[0]

        # * ii. Remove blanks
        filesForJob = [x for x in filesForJob if x != '']

        # * iii. Submit a job
        lastFourFileIDs = '-'.join( [lastFour[len(lastFour)-9:len(lastFour)-5] for lastFour in filesForJob] ) # get and concatenate last four digits of filenames
        filelist = ','.join( filesForJob )
        jdl_filename = "submitOneFile_{0}_{1}.jdl".format(args.outputDir, lastFourFileIDs)
        isZmumu = 'isZmumu' if 'zmumu' in args.inputTXTfile else 'isTTbar'
        is200PU = 'is200PU' if '200PU' in args.inputTXTfile else 'is0PU'
        

        os.system("touch {0}".format(jdl_filename))
        os.system("echo universe = vanilla > {0}".format(jdl_filename))
        os.system("echo Executable = cmsRunOneFile.csh >> {0}".format(jdl_filename))
        os.system("echo Should_Transfer_Files = YES >> {0}".format(jdl_filename))
        os.system("echo WhenToTransferOutput = ON_EXIT >> {0}".format(jdl_filename))
        os.system("echo Transfer_Input_Files = cmsRunOneFile.csh, MuonIsolationAnalyzer/test/singleFileMuonIsolation.py, {0} >> {1}".format(tarball_name, jdl_filename))
        os.system("echo Output = {0}/condor_out/outfile_{1}.out  >> {2}".format(args.outputDir, lastFourFileIDs, jdl_filename))
        os.system("echo Error = {0}/condor_err/outfile_{1}.err >> {2}".format(args.outputDir, lastFourFileIDs, jdl_filename))
        os.system("echo Log = {0}/condor_logs/outfile_{1}.log >> {2}".format(args.outputDir, lastFourFileIDs, jdl_filename))
        os.system("echo x509userproxy = ${{X509_USER_PROXY}} >> {0}".format(jdl_filename))
        os.system("echo Arguments = {0} {1} {2} {3} {4} >> {5}".format(isZmumu, is200PU, filelist, tarball_name, args.outputDir, jdl_filename))
        os.system("echo Queue 1 >> {0}".format(jdl_filename))       
        os.system("condor_submit {0}".format(jdl_filename))


# *** 3. Cleanup submission directory
print "\n##########     Cleanup submission directory     ##########\n"
os.system("rm *.jdl")
