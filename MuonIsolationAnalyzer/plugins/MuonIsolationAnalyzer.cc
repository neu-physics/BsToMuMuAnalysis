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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TH1.h"
#include "TH2.h"

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

      //typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
      //typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      float getMuonPFRelIso(const reco::Muon&) const;
      float getMuonPFCandIso(const reco::Muon&, edm::Handle<std::vector<reco::PFCandidate> >& pfCandHandle) const;
      bool  isGoodMuon(const reco::Muon&) const;
      int   ttbarDecayMode(edm::Handle<std::vector<reco::GenParticle> >& genHandle);
      void  getPromptMuons(edm::Handle<std::vector<reco::GenParticle> >& genHandle);
      const reco::Candidate * GetObjectJustBeforeDecay( const reco::Candidate * particle );
      bool WBosonDecaysProperly( const reco::Candidate *W );
      int WBosonDecayMode( const reco::Candidate *W );

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
      edm::EDGetTokenT< std::vector<reco::PFCandidate> > pfCandToken_;
      edm::Handle< std::vector<reco::PFCandidate> > pfCandHandle_;
      edm::EDGetTokenT<vector<SimVertex> >                 genVertexToken_;
      edm::Handle<vector<SimVertex> >                      genVertexHandle_;    
  //edm::EDGetTokenT<genXYZ>                             genXYZToken_;
  edm::EDGetTokenT< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > genXYZToken_;
  edm::Handle< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > genXYZHandle_;
      edm::EDGetTokenT<float>                              genT0Token_;    
      edm::Handle<float>                                   genT0Handle_;

      reco::Vertex vertex;
      bool isZmumuSignal_;
      vector<const reco::Candidate*> promptMuonTruthCandidates_;

      //---outputs
      TH1* h_event_cutflow_;    
      TH1* h_muon_pfRelIso03_BTL_;    
      TH1* h_muon_pfRelIso03_ETL_;   
      TH1* h_muon_pfCandIso03_BTL_;    
      TH1* h_muon_pfCandIso03_ETL_;    
      TH1* h_muon_cutflow_;    
      TH1* h_ttbarDecayMode_;    
      TH2* h_ttbarDecayMode_vs_nPromptMuons_;    
      TH1* h_muon_pfCandpT;
      TH1* h_muon_pT;

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
  pfCandToken_(consumes<std::vector<reco::PFCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("pfCandTag"))),
  genVertexToken_(consumes<vector<SimVertex> >(iConfig.getUntrackedParameter<edm::InputTag>("genVtxTag"))),
  //genXYZToken_(consumes<genXYZ>(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),    
  genXYZToken_(consumes< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> >(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
  genT0Token_(consumes<float>(iConfig.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
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

   iEvent.getByToken(genXYZToken_, genXYZHandle_);
   iEvent.getByToken(genT0Token_, genT0Handle_);
   iEvent.getByToken(genVertexToken_, genVertexHandle_);    

   //---get truth PV
   const SimVertex *genPV = NULL;// = SimVertex();
   if(genVertexHandle_.isValid())
     //genPV = genVertexHandle_.product()->at(0);
     genPV = &(genVertexHandle_.product()->at(0));
   else
     {
       //auto xyz = genXYZHandle_.product();
       //auto t = *genT0Handle_.product();
       //auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
       std::cout << "oh actually sometimes we don't find a gen pv" << std::endl;
       return;
       //genPV = *(SimVertex(v, t));
     }
   
   int vtx_index = -1;
   
   //---find the vertex this muon is best associated to..
   double min_dz = std::numeric_limits<double>::max();
   //double min_dzdt = std::numeric_limits<double>::max();
   for( unsigned i = 0; i < vertexHandle_->size(); ++i ) {
       const auto& vtx = (*vertexHandle_)[i];
       const float dz = std::abs(vtx.z() - genPV->position().z());
       if( dz < min_dz )
	 {
	   min_dz = dz;
	   vtx_index = i;
	 }
       // if( dz < 0.1 )
       // {
       //     if(isTimingSample_)
       //     {
       //         const double dzdt = pow((vtx.z() - genPV->position().z())/vtx.zError(), 2) +
       //             pow((vtx.t()-genPV->position().t())/vtx.tError(), 2);                            
       //         if( dzdt < min_dzdt )
       //         {
       //             min_dzdt = dzdt;
       //             vtx_index = i;
       //         }
       //     }
       //     else if( dz < min_dz )
       //     {
       //             min_dz = dz;
       //             vtx_index = i;
       //     }
       // }
   }
   if (vtx_index == -1) // no good PV found
     return;
   vertex = (*vertexHandle_)[vtx_index];
   h_event_cutflow_->Fill("Good Vertex", 1);

   /*
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
   */

   // *** 2. Load tracks
   iEvent.getByToken(tracksToken_,tracksHandle_);
   auto tracks = *tracksHandle_.product();

   for(unsigned iTrack = 0; iTrack < tracks.size(); ++iTrack)   {
     auto track = tracks.at(iTrack);
   }


   // *** 3. Load generator MC particles
   int decayMode = -1;
   iEvent.getByToken(genToken_,genHandle_);

   if (!isZmumuSignal_) // background ttbar
     decayMode = ttbarDecayMode( genHandle_ );
   else // signal Z->mumu
     getPromptMuons( genHandle_ );

   h_ttbarDecayMode_->Fill( decayMode );
   
   
   // *** 4. Load reco muons
   iEvent.getByToken(muonsToken_, muonsHandle_);
   iEvent.getByToken(pfCandToken_, pfCandHandle_);
   auto muons = *muonsHandle_.product();
   bool hasGoodMuon = false;

   for(unsigned iMuon = 0; iMuon < muons.size(); ++iMuon)   {
     auto muon = muons.at(iMuon);
     //h_muon_pT->Fill(muon.pt());
     if ( !isGoodMuon( muon )) continue;

     // ** A. "Analysis-like" isolation borrowed from ttH
     float pfRelIso03 = getMuonPFRelIso( muon );
     //std::cout << "Muon number " << iMuon << " has pT = " << muon.pt() << " and PFRelIso (R=0.3) = " << pfRelIso03 << std::endl;
     // ** B. "TDR-like" isolation using PFCandidates
     float pfCandIso03 = getMuonPFCandIso( muon, pfCandHandle_ );

     h_muon_pT->Fill(muon.pt());     

     // Fill information about muon PF relIso (R=0.3)
     if ( fabs(muon.eta()) < 1.5){
       h_muon_pfRelIso03_BTL_->Fill( pfRelIso03 );
       h_muon_pfCandIso03_BTL_->Fill( pfCandIso03 );
     }
     else if ( fabs(muon.eta()) > 1.5 && fabs(muon.eta()) < 2.8){
       h_muon_pfRelIso03_ETL_->Fill( pfRelIso03 );
       h_muon_pfCandIso03_ETL_->Fill( pfCandIso03 );
     }
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
  h_ttbarDecayMode_ = fileService->make<TH1D>("h_ttbarDecayMode", "h_ttbarDecayMode", 7, 0, 7);
  h_ttbarDecayMode_vs_nPromptMuons_ = fileService->make<TH2D>("h_ttbarDecayMode_vs_nPromptMuons", "h_ttbarDecayMode_vs_nPromptMuons", 7, 0, 7, 4, 0, 4);

  h_muon_cutflow_ = fileService->make<TH1D>("h_muon_cutflow", "h_muon_cutflow", 7, 0, 7);
  h_muon_pfRelIso03_BTL_ = fileService->make<TH1F>("h_muon_pfRelIso03_BTL", "h_muon_pfRelIso03_BTL", 250, 0, 5.0);
  h_muon_pfRelIso03_ETL_ = fileService->make<TH1F>("h_muon_pfRelIso03_ETL", "h_muon_pfRelIso03_ETL", 250, 0, 5.0);
  h_muon_pfCandIso03_BTL_ = fileService->make<TH1F>("h_muon_pfCandIso03_BTL", "h_muon_pfCandIso03_BTL", 250, 0, 5.0);
  h_muon_pfCandIso03_ETL_ = fileService->make<TH1F>("h_muon_pfCandIso03_ETL", "h_muon_pfCandIso03_ETL", 250, 0, 5.0);
  h_muon_pfCandpT = fileService->make<TH1F>("h_muon_pfCandpT", "h_muon_pfCandpT", 500, 0, 100.0);
  h_muon_pT = fileService->make<TH1F>("h_muon_pT", "h_muon_pT", 500, 0, 100.0);

}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonIsolationAnalyzer::endJob()
{
}

float MuonIsolationAnalyzer::getMuonPFCandIso(const reco::Muon& iMuon, edm::Handle<std::vector<reco::PFCandidate> >& pfCandHandle) const
{
  float result = 0; 
  float isoCone = 0.3;
  float identityCone = 0.05;
  int numberOfAssociatedPFCandidates = 0;

   auto pfCandidates = *pfCandHandle.product();

   for(unsigned iPFCand = 0; iPFCand < pfCandidates.size(); ++iPFCand)   {
     auto pfCandidate = pfCandidates.at(iPFCand);
     h_muon_pfCandpT->Fill(pfCandidate.pt());
     // skip neutrals
     if (pfCandidate.charge()==0) 
       continue;
     reco::TrackRef track = pfCandidate.trackRef();
     if(track.isNull() || !track->quality(reco::TrackBase::highPurity))
       continue;

     // kinematic cuts on PF candidate
     if ( fabs(track->eta())<1.5){ // BTL acceptance
       //std::cout << "BTL  " << fabs(track->dz(vertex.position())) << std::endl;
       if (track->pt() < 0.7)
	 continue;
       if ( fabs(track->dz(vertex.position())) > 0.1 )
	 continue;
    }
     else if ( fabs(track->eta())>1.5 && fabs(track->eta())<2.8) { // ETL acceptance
       //std::cout << "ETL  " << fabs(track->dz(vertex.position())) << std::endl;
       if (track->pt() < 0.4)
	 continue;
       if ( fabs(track->dz(vertex.position())) > 0.2 )
         continue;
     }
     else
       continue;
     
     // sum candidates in cone
     float Deta = iMuon.eta() - pfCandidate.eta();
     float Dphi = deltaPhi( iMuon.phi(), pfCandidate.phi());
     float DR = sqrt(Deta*Deta+Dphi*Dphi);
     
     if (DR < isoCone && DR > identityCone){
       numberOfAssociatedPFCandidates++;
       result += pfCandidate.pt();
     }
   }
   //std::cout << "Isolation = " << result << " with " << numberOfAssociatedPFCandidates << " associated to muon" << std::endl;
   
   return result;
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

bool MuonIsolationAnalyzer::WBosonDecaysProperly( const reco::Candidate *W ){ 

  bool up_quark = false ; 
  bool down_quark = false ; 

  bool up_lepton = false ; 
  bool down_lepton = false ; 

  for ( unsigned int i = 0 ; i <  W -> numberOfDaughters(); i++ ){

    if( abs( W -> daughter( i ) -> pdgId() ) == 6
	||  abs( W -> daughter( i ) -> pdgId() ) == 4
	||  abs( W -> daughter( i ) -> pdgId() ) == 2 ) up_quark = true ; 
      
    if( abs( W -> daughter( i ) -> pdgId() ) == 12
	||  abs( W -> daughter( i ) -> pdgId() ) == 14
	||  abs( W -> daughter( i ) -> pdgId() ) == 16 ) up_lepton = true ; 

    if( abs( W -> daughter( i ) -> pdgId() ) == 5
	||  abs( W -> daughter( i ) -> pdgId() ) == 3
	||  abs( W -> daughter( i ) -> pdgId() ) == 1 ) down_quark = true ;

    if( abs( W -> daughter( i ) -> pdgId() ) == 11
	||  abs( W -> daughter( i ) -> pdgId() ) == 13
	||  abs( W -> daughter( i ) -> pdgId() ) == 15 ) down_lepton = true ; 

  }// end for-loop
  
  if( up_quark && down_quark ) return true; 
  if( up_lepton && down_lepton ) return true ; 

  return false; 
  
}

int MuonIsolationAnalyzer::WBosonDecayMode( const reco::Candidate *W ){

  int upDaughter_pdgId = 0;
  int downDaughter_pdgId = 0;

  // =================================================
  // (case-1) In some MC, W boson decays but stays (example : W -> u+d+W)
  // (case-2) In some MC, W boson decays after some step (example : W->W->u+d)
  // (case-3) In other MC, W boson obtains extra particles (example : W->Wee (+/-))
  //  In order to handle those cases,
  //   - check if the W boson has a pair of expected decay products, otherwise keep following the W boson.
  
  //   bool W_boson_decays  = false ;
  //   for ( unsigned int i = 0 ; i <  W -> numberOfDaughters(); i++ ){
  //     if(       abs( W -> daughter( i ) -> pdgId() ) == 6
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 4
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 2
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 12
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 14
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 16
  //       || abs( W -> daughter( i ) -> pdgId() ) == 5
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 3
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 1
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 11
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 13
  //       ||  abs( W -> daughter( i ) -> pdgId() ) == 15 ){
  //       W_boson_decays = true ;
  //       std::cout <<"W boson decays true. pdgid = " << ( W -> daughter( i ) -> pdgId()  ) << std::endl ;
  //     } // end if
  //   }// end for-loop
  //   std::cout <<"test"<< std::endl ;
  //   if( ! W_boson_decays ){
  //     W = GetObjectJustBeforeDecay( W ) ;
  //     std::cout <<"W boson : GetObjectJustBeforeDecay  is called." << std::endl ; 
  //     for ( unsigned int i = 0 ; i <  W -> numberOfDaughters(); i++ ){
  //       std::cout << "in the loop, the daughters are "<< W -> daughter( i ) -> pdgId() << std::endl ; 
  //     }// end for-loop                     
  //   }
  
  while ( ! WBosonDecaysProperly( W ) ){
    
    bool find_next_candidate = false ; 
    for ( unsigned int i = 0 ; (! find_next_candidate) && ( i < W -> numberOfDaughters()  ) ; i++ ){
      
      if( W -> daughter( i ) -> pdgId() == W -> pdgId() ){
	find_next_candidate = true ; 
	W = W -> daughter( i ) ; 
      } 
      
    }// end for-loop
  }
  
  for ( unsigned int i = 0 ; i <  W -> numberOfDaughters(); i++ ){
    
    if(       abs( W -> daughter( i ) -> pdgId() ) == 6
	      ||  abs( W -> daughter( i ) -> pdgId() ) == 4
	      ||  abs( W -> daughter( i ) -> pdgId() ) == 2
	      ||  abs( W -> daughter( i ) -> pdgId() ) == 12
	      ||  abs( W -> daughter( i ) -> pdgId() ) == 14
	      ||  abs( W -> daughter( i ) -> pdgId() ) == 16 ){
      // --- up type
      upDaughter_pdgId  =  fabs(W->daughter( i )->pdgId());
    }else if ( abs( W -> daughter( i ) -> pdgId() ) == 5
	       ||  abs( W -> daughter( i ) -> pdgId() ) == 3
	       ||  abs( W -> daughter( i ) -> pdgId() ) == 1
	       ||  abs( W -> daughter( i ) -> pdgId() ) == 11
	       ||  abs( W -> daughter( i ) -> pdgId() ) == 13
	       ||  abs( W -> daughter( i ) -> pdgId() ) == 15 ){
      // --- down type
      downDaughter_pdgId  =  fabs(W->daughter( i )->pdgId());
    }
    // =================================================
    
  }  
  if (0>1)
    std::cout << "top W has daughter0 pdgId = " << upDaughter_pdgId << " , daughter1 pdgId = " << downDaughter_pdgId << std::endl;
  
  return downDaughter_pdgId;

}

const reco::Candidate * MuonIsolationAnalyzer::GetObjectJustBeforeDecay( const reco::Candidate * particle ){

  for ( unsigned int i = 0 ; i <  particle -> numberOfDaughters(); i++ ){
    if( particle -> daughter( i ) -> pdgId()  ==  particle -> pdgId() ){

      return GetObjectJustBeforeDecay( particle -> daughter (i) );

    } // end if
  } // end for

  return particle ;

}

void MuonIsolationAnalyzer::getPromptMuons(edm::Handle<std::vector<reco::GenParticle> >& genHandle)
{
  unsigned int numberOfPromptMuons = 0;
  promptMuonTruthCandidates_ = {};

  // *** 1. Loop over MC particles
  for(size_t iGenParticle=0; iGenParticle<genHandle->size();iGenParticle++){
    
    // ** A. Prompt muons
    if ( fabs((*genHandle)[iGenParticle].pdgId()) == 13 && (*genHandle)[iGenParticle].isPromptFinalState()){
      numberOfPromptMuons++;
      const reco::Candidate * promptMuon = &(*genHandle)[iGenParticle] ;
      promptMuonTruthCandidates_.push_back( promptMuon );
    }
  }
  if (numberOfPromptMuons > 2)
    std::cout << "!!!!!!!!!!!!!!!! PANIC !!!!!! in z->mumu event, # of prompt final state muons = " << numberOfPromptMuons << std::endl;
  
  return;
}

int MuonIsolationAnalyzer::ttbarDecayMode(edm::Handle<std::vector<reco::GenParticle> >& genHandle)
{
  int decayMode = -1; // 0 == SL electron, 1 == SL muon, 2 == SL tau, 3 == DL, 4 == all had
  bool processedTopQuark = false;
  bool processedAntiTopQuark = false;
  int topWDecayMode = -1;
  int antitopWDecayMode = -1;
  unsigned int numberOfPromptMuons = 0;
  promptMuonTruthCandidates_ = {};

  // *** 1. Loop over MC particles
  for(size_t iGenParticle=0; iGenParticle<genHandle->size();iGenParticle++){

    // ** A. Prompt muons
    if ( fabs((*genHandle)[iGenParticle].pdgId()) == 13 && (*genHandle)[iGenParticle].isPromptFinalState()){
      numberOfPromptMuons++;
      const reco::Candidate * promptMuon = &(*genHandle)[iGenParticle] ;
      promptMuonTruthCandidates_.push_back( promptMuon );
    }

    // ** B. Top quark loop
    if ( (*genHandle)[iGenParticle].pdgId() == 6 && processedTopQuark == false) { // check top quark
      const reco::Candidate * genCandidate = &(*genHandle)[iGenParticle] ;
      genCandidate = GetObjectJustBeforeDecay( genCandidate );
      processedTopQuark = true;
      if (genCandidate->numberOfDaughters() == 2) {
	for (unsigned int iDaughter = 0; iDaughter < genCandidate->numberOfDaughters(); iDaughter++){

	  if ( fabs(genCandidate->daughter( iDaughter )->pdgId()) == 24) {// W boson
	    const reco::Candidate * W = genCandidate->daughter(iDaughter);
	    topWDecayMode = WBosonDecayMode( W );
	  } // end W boson loop

	} // end loop on top quark daughters
      } // if statement requiring two daughters of last top quark
      else
	std::cout << " top quark DOES NOT HAVE 2 DAUGHTERS BUT IS END OF CHAIN. PANNIC!!!!!!!!!!!!!!!!!" << std::endl;
    } // end top loop

    // ** C. Anti-top quark loop
    if ( (*genHandle)[iGenParticle].pdgId() == -6 && processedAntiTopQuark == false) { // check anti-top quark
      const reco::Candidate * genCandidate = &(*genHandle)[iGenParticle] ;
      genCandidate = GetObjectJustBeforeDecay( genCandidate );
      processedAntiTopQuark = true;

      if (genCandidate->numberOfDaughters() == 2) {
	for (unsigned int iDaughter = 0; iDaughter < genCandidate->numberOfDaughters(); iDaughter++){

	  if ( fabs(genCandidate->daughter( iDaughter )->pdgId()) == 24) {// W boson
	    const reco::Candidate * W = genCandidate->daughter(iDaughter);
	    antitopWDecayMode = WBosonDecayMode( W );
	  } // end W boson loop

	} // end loop on antitop quark daughters
      } // if statement requiring two daughters of last antitop quark
      else
	std::cout << " anti-top quark DOES NOT HAVE 2 DAUGHTERS BUT IS END OF CHAIN. PANNIC!!!!!!!!!!!!!!!!!" << std::endl;
    } // end anti-top loop
  } // end loop over mc particles

  // *** 2. Classify ttbar decay mode
  if (antitopWDecayMode == -1 || topWDecayMode == -1){
    std::cout << "One or more W bosons from top quark were not properly classified" << std::endl;
    decayMode = -11;
  }
  else if ( (antitopWDecayMode == 11 && topWDecayMode < 10) || (topWDecayMode == 11 && antitopWDecayMode < 10) ) // W->enu, W->had
    decayMode = 0;
  else if ( (antitopWDecayMode == 13 && topWDecayMode < 10) || (topWDecayMode == 13 && antitopWDecayMode < 10) ) // W->munu, W->had
    decayMode = 1;
  else if ( (antitopWDecayMode == 15 && topWDecayMode < 10) || (topWDecayMode == 15 && antitopWDecayMode < 10) ) // W->taunu, W->had
    decayMode = 2;
  else if ( antitopWDecayMode > 10 && topWDecayMode > 10 ) // DL
    decayMode = 3;
  else if ( antitopWDecayMode < 10 && topWDecayMode < 10 ) // all-had
    decayMode = 4;
  else{
    std::cout << "no idea what happened in W boson classification. both non -1 but no total mode identified" << std::endl;
    decayMode = -1;
  }
  h_ttbarDecayMode_vs_nPromptMuons_->Fill( decayMode, numberOfPromptMuons );

  //std::cout << "ttbar decay mode = " << decayMode << " , number of prompt muons = " << numberOfPromptMuons << " , size of truth indices vector = " << promptMuonTruthIndices_.size() << std::endl;
  if ( numberOfPromptMuons != promptMuonTruthCandidates_.size())
    std::cout << "NUMBER OF PROMPT MUONS NOT EQUAL TO NUMBER OF INDICES SAVED.... MISCOUNT SOMEWHERE !!!!!!!!!!!!!!!!!!! " << std::endl;

  return decayMode;
}

bool MuonIsolationAnalyzer::isGoodMuon(const reco::Muon& iMuon) const
{
  
  h_muon_cutflow_->Fill("All Muons", 1);
  // *** 0. pT > 20 GeV ---> may need to change this for Bs->mumu studies, but keep for now to reproduce TDR results
  if (iMuon.pt() < 20.)
    return false;
  h_muon_cutflow_->Fill("pT > 20", 1);

  // *** 1. |eta| < 2.4
  if ( fabs(iMuon.eta()) > 2.4 )
    return false;
  h_muon_cutflow_->Fill("|eta| < 2.4", 1);
  
  // *** 2. Loose ID
  //if ( iMuon.passed('CutBasedIdLoose'))//userFloat("CutBasedIdLoose") ) // CutBasedIdLoose, MvaLoose, others?
  if ( !(muon::isLooseMuon(iMuon)) )
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
  

  bool recoMuonMatchedToPromptTruth = false;
  // *** 5A. accept only prompt muons if Z->mumu signal
  if (isZmumuSignal_) {
    for( unsigned int iTruthMuon = 0; iTruthMuon < promptMuonTruthCandidates_.size(); iTruthMuon++){
      const reco::Candidate * promptTruthMuon = ( promptMuonTruthCandidates_[iTruthMuon]);
      float Deta = iMuon.eta() - promptTruthMuon->eta();
      float Dphi = deltaPhi( iMuon.phi(), promptTruthMuon->phi());
      float DR = sqrt(Deta*Deta+Dphi*Dphi);
      
      if (DR < 0.05)
	recoMuonMatchedToPromptTruth = true;
    }
    if (!recoMuonMatchedToPromptTruth)
      return false;
    h_muon_cutflow_->Fill("Signal Prompt Muon", 1);
  }
  // *** 5B. reject prompt muons if ttbar background
  else if (!isZmumuSignal_) {
    for( unsigned int iTruthMuon = 0; iTruthMuon < promptMuonTruthCandidates_.size(); iTruthMuon++){
      const reco::Candidate * promptTruthMuon = ( promptMuonTruthCandidates_[iTruthMuon]);
      float Deta = iMuon.eta() - promptTruthMuon->eta();
      float Dphi = deltaPhi( iMuon.phi(), promptTruthMuon->phi());
      float DR = sqrt(Deta*Deta+Dphi*Dphi);
      
      if (DR < 0.05)
	recoMuonMatchedToPromptTruth = true;
    }
    if (recoMuonMatchedToPromptTruth)
      return false;
    h_muon_cutflow_->Fill("Non-prompt Bkg Muon", 1);
    
  }

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
