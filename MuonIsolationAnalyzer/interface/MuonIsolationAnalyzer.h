
// -*- C++ -*-
//
// Package:    BsToMuMuAnalysis/MuonIsolationAnalyzer
// Class:      MuonIsolationAnalyzer
//
/**\class MuonIsolationAnalyzer MuonIsolationAnalyzer.cc BsToMuMuAnalysis/MuonIsolationAnalyzer/plugins/MuonIsolationAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Tannenwald
//         Created:  Fri, 12 Jul 2019 15:17:35 GMT
//
//


// system include files
#include <memory>

#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "BsToMuMuAnalysis/MuonIsolationAnalyzer/interface/Utils.h"

#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace std;

struct eventInfo
{
  vector <int> nevt;
  vector <double> muon_eta;
  vector <double> muon_pt;
  vector <double> muon_pfCand;
  vector <double> muon_pfCand_dt;
  vector <double> muon_pfCand_noDxy;
  vector <double> muon_pfCand_noDxy_dt;
};

//using reco::TrackCollection;
//using reco::MuonCollection;

class MuonIsolationAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonIsolationAnalyzer(const edm::ParameterSet&);
      ~MuonIsolationAnalyzer();

      //typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
      //typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      float getMuonPFRelIso(const reco::Muon&) const;
      int getMuonPFCandIso(const reco::Muon&, edm::Handle<std::vector<reco::PFCandidate> >& pfCandHandle, const SimVertex *genPV, std::vector<float> &isolations, edm::Handle<ValueMap<float> > trackFastSimTimeValueMap, edm::Handle<ValueMap<float> > trackFastSimTimeErrValueMap, double timeResolution=-1, int n_evt=0 ) const;
      bool  isGoodMuon(const reco::Muon&,  edm::Handle<std::vector<reco::GenParticle> >& genHandle, edm::Handle<std::vector<reco::GenJet> >& genJetHandle, const SimVertex *genPV) const;
      int   ttbarDecayMode(edm::Handle<std::vector<reco::GenParticle> >& genHandle);
      void  getPromptTruthMuons(edm::Handle<std::vector<reco::GenParticle> >& genHandle);
      const reco::Candidate * GetObjectJustBeforeDecay( const reco::Candidate * particle );
      bool WBosonDecaysProperly( const reco::Candidate *W );
      int WBosonDecayMode( const reco::Candidate *W );
      bool passDeltaR( float coneSize, float eta1, float phi1, float eta2, float phi2) const;

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      void initEventStructure();

      // ----------member data ---------------------------
      //edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      //edm::Handle<reco::TrackCollection > tracksHandle_;
      
      //edm::EDGetTokenT<reco::MuonCollection > muonsToken_;
      //edm::Handle<reco::MuonCollection > muonsHandle_;
      edm::EDGetTokenT< std::vector<reco::Track> > tracksToken_;
      edm::Handle< std::vector<reco::Track> > tracksHandle_;
      edm::EDGetTokenT< std::vector<reco::Muon> > muonsToken_;
      edm::Handle< std::vector<reco::Muon> > muonsHandle_;
      edm::EDGetTokenT< std::vector<reco::Vertex> > vertex3DToken_;
      edm::Handle< std::vector<reco::Vertex> > vertex3DHandle_;
      edm::EDGetTokenT< std::vector<reco::Vertex> > vertex4DToken_;
      edm::Handle< std::vector<reco::Vertex> > vertex4DHandle_;
      edm::EDGetTokenT< std::vector<reco::GenParticle> > genToken_;
      edm::Handle< std::vector<reco::GenParticle> > genHandle_;
      edm::EDGetTokenT< std::vector<reco::GenJet> > genJetToken_;
      edm::Handle< std::vector<reco::GenJet> > genJetHandle_;
      edm::EDGetTokenT< std::vector<reco::PFCandidate> > pfCandToken_;
      edm::Handle< std::vector<reco::PFCandidate> > pfCandHandle_;
      edm::EDGetTokenT<vector<SimVertex> >                 genVertexToken_;
      edm::Handle<vector<SimVertex> >                      genVertexHandle_;    
      //edm::EDGetTokenT<genXYZ>                             genXYZToken_;
      edm::EDGetTokenT< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > genXYZToken_;
      edm::Handle< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > genXYZHandle_;
      edm::EDGetTokenT<float>                              genT0Token_;    
      edm::Handle<float>                                   genT0Handle_;
      edm::EDGetTokenT<ValueMap<float> > trackTimeToken_;
      edm::Handle<ValueMap<float> > trackTimeValueMap;
      edm::EDGetTokenT<ValueMap<float> > trackTimeErrToken_;
      edm::Handle<ValueMap<float> > trackTimeErrValueMap;
      edm::EDGetTokenT<ValueMap<float> > trackFastSimTimeToken_;
      edm::Handle<ValueMap<float> > trackFastSimTimeValueMap;
      edm::EDGetTokenT<ValueMap<float> > trackFastSimTimeErrToken_;
      edm::Handle<ValueMap<float> > trackFastSimTimeErrValueMap;

      reco::Vertex vertex3D;
      reco::Vertex vertex4D;
      //SimVertex genPV;

      bool isZmumuSignal_;
      TRandom *gRandom;
      TRandom *gRandom2;
      TRandom *gRandom3;
      //string lastFourFileID_;
      vector<const reco::Candidate*> promptMuonTruthCandidates_; 
      float coneSize_muonToGetJet;
      float coneSize_muonToPromptTruth;

      float pfCandIdentityCut;
      double nPromptMuons;
      double nNonPromptMuons;

      float dxy_muonVertex;
      float dz_muonVertex;
      float dxy_pfCandVertex;
      float dz_pfCandVertex;

      double btlEfficiency;
      double etlEfficiency;

      //---outputs
      TTree *eventTree;
      eventInfo *evInfo;
      TH1* h_event_cutflow_;    
      TH1* h_event_nPromptMuons_;    
      TH1* h_event_nNonPromptMuons_;    
      TH1* h_muon_pfRelIso03_BTL_;    
      TH1* h_muon_pfRelIso03_ETL_;   

      // isolation plots for ROC curves
      TH1* h_muon_pfCandIso03_BTL_;    
      TH1* h_muon_pfCandIso03_ETL_;    
      TH1* h_muon_pfCandIso03_noDxy_BTL_;    
      TH1* h_muon_pfCandIso03_noDxy_ETL_;    

      // isolation plots for ROC curves, add 3 sigma (sigma=30 ps) timing cut
      TH1* h_muon_pfCandIso03_dt30_BTL_;    
      TH1* h_muon_pfCandIso03_dt30_ETL_;    
      TH1* h_muon_pfCandIso03_dt30_noDxy_BTL_;    
      TH1* h_muon_pfCandIso03_dt30_noDxy_ETL_;    

      TH1* h_muon_pfCandIso03_dxy_BTL_;    
      TH1* h_muon_pfCandIso03_dz_BTL_;    
      TH1* h_muon_cutflow_;    
      TH1* h_muon_dxy_BTL_;    
      TH1* h_muon_dz_BTL_;    
      TH1* h_ttbarDecayMode_;    
      TH2* h_ttbarDecayMode_vs_nPromptMuons_;    
      TH1* h_muon_pfCandpT;
      TH1* h_muon_pT;
      TH1* h_muon_numCand;
      TH1* h_pfCandidate_cutflow_;    

};
