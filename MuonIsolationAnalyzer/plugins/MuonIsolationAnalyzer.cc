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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TH1.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace std;

//using reco::TrackCollection;
//using reco::MuonCollection;

class MuonIsolationAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonIsolationAnalyzer(const edm::ParameterSet&);
      ~MuonIsolationAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      float getMuonPFRelIso(const reco::Muon&) const;
      bool  isGoodMuon(const reco::Muon&) const;

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      //edm::Handle<reco::TrackCollection > tracksHandle_;
      
      //edm::EDGetTokenT<reco::MuonCollection > muonsToken_;
      //edm::Handle<reco::MuonCollection > muonsHandle_;
      edm::EDGetTokenT< std::vector<reco::Track> > tracksToken_;
      edm::Handle< std::vector<reco::Track> > tracksHandle_;
      edm::EDGetTokenT< std::vector<reco::Muon> > muonsToken_;
      edm::Handle< std::vector<reco::Muon> > muonsHandle_;
      edm::EDGetTokenT< std::vector<reco::Vertex> > vertexToken_;
      edm::Handle< std::vector<reco::Vertex> > vertexHandle_;
      edm::EDGetTokenT< std::vector<reco::GenParticle> > genToken_;
      edm::Handle< std::vector<reco::GenParticle> > genHandle_;

      reco::Vertex vertex;
      bool isZmumuSignal_;


      //---outputs
      TH1* h_event_cutflow_;    
      TH1* h_muon_pfRelIso03_;    
      TH1* h_muon_cutflow_;    

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonIsolationAnalyzer::MuonIsolationAnalyzer(const edm::ParameterSet& iConfig):
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  //muonsToken_(consumes<MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("muonsTag")))
  tracksToken_(consumes<std::vector<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  muonsToken_(consumes<std::vector<reco::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonsTag"))),
  vertexToken_(consumes<std::vector<reco::Vertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag"))),
  genToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genTag"))),
  isZmumuSignal_(iConfig.getParameter<bool>("isZmumuSignal"))
{
   //now do what ever initialization is needed

  //---outputs
  //edm::Service<TFileService> fs_; 

}


MuonIsolationAnalyzer::~MuonIsolationAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonIsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Handle<TrackCollection> tracks;
   //iEvent.getByToken(tracksToken_, tracks);

   // *** 0. Start cutflow
   h_event_cutflow_->Fill("All Events", 1);

   // *** 1. Load vertices
   iEvent.getByToken(vertexToken_, vertexHandle_);
   auto vertices = *vertexHandle_.product();

   int numpv=0; int iPV=0;
   bool firstGoodPV = false;
   reco::Vertex vertex;
   if( vertexHandle_.isValid() ){
     for(unsigned iVertex = 0; iVertex < vertices.size(); ++iVertex)   {
       auto vtx = vertices.at(iVertex);
       //for( reco::VertexCollection::const_iterator vtx = vtxs.begin(); vtx!=vtxs.end(); ++vtx ){
      
       iPV++;
       bool isGood = ( !(vtx.isFake()) &&
		       (vtx.ndof() >= 4.0) &&
		       (fabs(vtx.z()) <= 24.0) &&
		       (fabs(vtx.position().Rho()) <= 2.0) 
		       );
       
       if( !isGood ) continue;
       
       if( iPV==1 ){
	 firstGoodPV = true;
	 //vertex = (*vtx);
	 vertex = vtx;
       }
       numpv++;
     }
   }
   if (firstGoodPV == false) // no good PV found
     return;
   h_event_cutflow_->Fill("Good Vertex", 1);

   // *** 2. Load tracks
   iEvent.getByToken(tracksToken_,tracksHandle_);
   auto tracks = *tracksHandle_.product();

   for(unsigned iTrack = 0; iTrack < tracks.size(); ++iTrack)   {
     auto track = tracks.at(iTrack);
     //std::cout << "Track number " << iTrack << " has pT = " << track.pt() << std::endl;
   }

   // *** 3. Load generator MC particles
   iEvent.getByToken(genToken_,genHandle_);
   auto genParticles = *genHandle_.product();

   for(unsigned iGenParticle = 0; iGenParticle < genParticles.size(); ++iGenParticle)   {
     auto genParticle = genParticles.at(iGenParticle);
     //std::cout << "Track number " << iTrack << " has pT = " << track.pt() << std::endl;
   }


   // *** 4. Load muons
   iEvent.getByToken(muonsToken_, muonsHandle_);
   auto muons = *muonsHandle_.product();
   bool hasGoodMuon = false;

   for(unsigned iMuon = 0; iMuon < muons.size(); ++iMuon)   {
     auto muon = muons.at(iMuon);
     if ( !isGoodMuon( muon )) continue;

     std::cout << "Muon number " << iMuon << " has pT = " << muon.pt() << std::endl;
     std::cout << "Muon number " << iMuon << " has sum track pT (dR=0.3) isolation = " << muon.isolationR03().sumPt << std::endl;
     std::cout << "Muon number " << iMuon << " has sum non-PV track pT (dR=0.3) isolation = " << muon.pfIsolationR03().sumPUPt << std::endl;
     float pfRelIso03 = getMuonPFRelIso( muon );
     std::cout << "Muon number " << iMuon << " has PFRelIso (R=0.3) = " << pfRelIso03 << std::endl;

     // Fill information about muon PF relIso (R=0.3)
     h_muon_pfRelIso03_->Fill( pfRelIso03 );
     if (hasGoodMuon == false)
       hasGoodMuon = true;
   }

   if (hasGoodMuon)
     h_event_cutflow_->Fill(">= 1 Good Muon", 1);

   /*
   for(TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }
   */

   /*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   */
}


// ------------ method called once each job just before starting event loop  ------------
void
MuonIsolationAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");

  h_event_cutflow_ = fileService->make<TH1D>("h_event_cutflow", "h_event_cutflow", 7, 0, 7);

  h_muon_pfRelIso03_ = fileService->make<TH1F>("h_muon_pfRelIso03", "h_muon_pfRelIso03", 250, 0, 5.0);
  h_muon_cutflow_ = fileService->make<TH1D>("h_muon_cutflow", "h_muon_cutflow", 7, 0, 7);

}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonIsolationAnalyzer::endJob()
{
}

float MuonIsolationAnalyzer::getMuonPFRelIso(const reco::Muon& iMuon) const
{
  float result = 9999; 

  double pfIsoCharged = iMuon.pfIsolationR03().sumChargedHadronPt;
  double pfIsoNeutral = iMuon.pfIsolationR03().sumNeutralHadronEt + iMuon.pfIsolationR03().sumPhotonEt;

  double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - 0.5*iMuon.pfIsolationR03().sumPUPt );

  result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
  
  return result;
}


bool MuonIsolationAnalyzer::isGoodMuon(const reco::Muon& iMuon) const
{
  
  h_muon_cutflow_->Fill("All Muons", 1);
  // *** 0. pT > 10 GeV ---> may need to change this for Bs->mumu studies, but keep for now to reproduce TDR results
  if (iMuon.pt() < 10.)
    return false;
  h_muon_cutflow_->Fill("pT > 10", 1);

  // *** 1. |eta| < 2.4
  if ( fabs(iMuon.eta()) > 2.4 )
    return false;
  h_muon_cutflow_->Fill("|eta| < 2.4", 1);
  
  // *** 2. Loose ID
  if ( iMuon.passed('CutBasedIdLoose'))//userFloat("CutBasedIdLoose") ) // CutBasedIdLoose, MvaLoose, others?
    return false;
  h_muon_cutflow_->Fill("Loose ID", 1);

  // *** just check to make sure muon track available
  if ( !(iMuon.muonBestTrack().isAvailable()) )
    return false;

  // *** 3. d0 cut
  if ( fabs(iMuon.muonBestTrack()->dxy(vertex.position())) > 0.2 )
    return false;
  h_muon_cutflow_->Fill("d0 < 0.2 cm", 1);

  // *** 4. z0 cut
  if ( fabs(iMuon.muonBestTrack()->dz(vertex.position())) > 0.5 )
    return false;
  h_muon_cutflow_->Fill("z0 < 0.5 cm", 1);
  

  return true;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonIsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsolationAnalyzer);
