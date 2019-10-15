import ROOT as R
import math as M
import argparse
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--inputDir',dest='inputDir')
parser.add_argument('--pattern',dest='pattern')
parser.add_argument('--output',dest="output")

args = parser.parse_args()

if (args.inputDir != ""):
  print args.inputDir
  chain=R.TChain("MuonIsolationAnalyzer/tree_isolation")
  files = []
  if args.inputDir[-1] != '/':
    args.inputDir += '/'
  print ('>> Creating list of files from : \n'+args.inputDir)
  command = '/bin/find '+args.inputDir+' -type f | grep root | grep -v failed | grep '+args.pattern
  str_files = subprocess.check_output(command,shell=True).splitlines()
  files.extend(['file:'+ifile for ifile in str_files])
  for file in files:
    print ">>Adding "+file
    chain.Add(file)

else:
  print('!!! no inputDir!!!')

Iso_range = [4,10,20,27,38,50]

histos = {}

histos["muon_iso"]=R.TH1F("muon_iso","muon_iso",5000,0,5.0)
histos["muon_iso_dt"]=R.TH1F("muon_iso_dt","muon_iso_dt",5000,0,5.0)
histos["muon_pt"]=R.TH1F("muon_pt","muon_pt",1000,0,100)

for r in Iso_range:
  h = 'muon_iso_{0}'.format(r)
  ht = 'muon_iso_dt_{0}'.format(r)
  histos[h]=R.TH1F(h,h,5000,0,5.0)
  histos[ht]=R.TH1F(ht,ht,5000,0,5.0)

E_threshold = 4.0

#loop over all events
for ievt,evt in enumerate(chain):
  if(ievt%100000==0): print ('analyzing event {0}'.format(ievt))
  for imu in range(0,len(evt.muon_pt)):
    if (evt.muon_pt[imu]<E_threshold):
      continue
    nevt = evt.nevt[imu]
    eta = evt.muon_eta[imu]
    pt = evt.muon_pt[imu]
    iso = evt.muon_pfCand_noDxy[imu]
    isodt = evt.muon_pfCand_noDxy_dt[imu]
    histos["muon_pt"].Fill(evt.muon_pt[imu])
    if (abs(evt.muon_eta[imu])<1.5):
      histos["muon_iso"].Fill(evt.muon_pfCand_noDxy[imu])
      histos["muon_iso_dt"].Fill(evt.muon_pfCand_noDxy_dt[imu])

      for i in range (0,len(Iso_range)):
        name = 'muon_iso_{0}'.format(Iso_range[i])
        namet = 'muon_iso_dt_{0}'.format(Iso_range[i])
        if (i!=(len(Iso_range)-1)):
          if ( evt.muon_pt[imu]>=Iso_range[i] and evt.muon_pt[imu]<Iso_range[i+1] ):
            histos[name].Fill(evt.muon_pfCand_noDxy[imu])
            histos[namet].Fill(evt.muon_pfCand_noDxy_dt[imu])
        else:
          if (evt.muon_pt[imu]>=Iso_range[i]):
            histos[name].Fill(evt.muon_pfCand_noDxy[imu])
            histos[namet].Fill(evt.muon_pfCand_noDxy_dt[imu])

    #print('nevt: {0} eta: {1} pt: {2} iso: {3} iso_dt: {4}'.format(nevt, eta, pt, iso, isodt))


fOut=R.TFile(args.output,"RECREATE")
for hn, histo in histos.iteritems():
  histo.Write()
fOut.Close()
print "Saved histos in "+args.output
