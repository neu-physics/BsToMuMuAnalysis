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
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "BsToMuMuAnalysis/MuonIsolationAnalyzer/interface/Utils.h"

#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"

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
      int getMuonPFCandIso(const reco::Muon&, edm::Handle<std::vector<reco::PFCandidate> >& pfCandHandle, const SimVertex *genPV, std::vector<float> &isolations, edm::Handle<ValueMap<float> > trackFastSimTimeValueMap, edm::Handle<ValueMap<float> > trackFastSimTimeErrValueMap, double timeResolution=-1 ) const;
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
  vertex3DToken_(consumes<std::vector<reco::Vertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertex3DTag"))),
  vertex4DToken_(consumes<std::vector<reco::Vertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertex4DTag"))),
  genToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genTag"))),
  genJetToken_(consumes<std::vector<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genJetTag"))),
  pfCandToken_(consumes<std::vector<reco::PFCandidate> >(iConfig.getUntrackedParameter<edm::InputTag>("pfCandTag"))),
  genVertexToken_(consumes<vector<SimVertex> >(iConfig.getUntrackedParameter<edm::InputTag>("genVtxTag"))),
  //genXYZToken_(consumes<genXYZ>(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),    
  genXYZToken_(consumes< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> >(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
  genT0Token_(consumes<float>(iConfig.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
  trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) ),
  trackTimeErrToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeErrValueMapTag" ) ) ),
  trackFastSimTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackFastSimTimeValueMapTag" ) ) ),
  trackFastSimTimeErrToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackFastSimTimeErrValueMapTag" ) ) ),
  isZmumuSignal_(iConfig.getParameter<bool>("isZmumuSignal"))
  //lastFourFileID_(iConfig.getParameter<string>("lastFourFileID"))
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
   coneSize_muonToGetJet = 0.3;
   coneSize_muonToPromptTruth = 0.2;

   // Set constants
   dxy_muonVertex   = 0.2;
   dz_muonVertex    = 0.5;
   dxy_pfCandVertex = 0.02;
   dz_pfCandVertex  = 0.2; //for BTL
   //dz_pfCandVertex  = 0.1;

   btlEfficiency = 0.85;
   etlEfficiency = 0.90;
   
   //Handle<TrackCollection> tracks;
   //iEvent.getByToken(tracksToken_, tracks);

   // *** 0. Start cutflow
   h_event_cutflow_->Fill("All Events", 1);

   // *** 1. Load vertices
   iEvent.getByToken(vertex3DToken_, vertex3DHandle_);
   auto vertices3D = *vertex3DHandle_.product();

   iEvent.getByToken(vertex4DToken_, vertex4DHandle_);
   auto vertices4D = *vertex4DHandle_.product();

   iEvent.getByToken(genJetToken_, genJetHandle_);

   iEvent.getByToken(genXYZToken_, genXYZHandle_);
   iEvent.getByToken(genT0Token_, genT0Handle_);
   iEvent.getByToken(genVertexToken_, genVertexHandle_);    
  
   iEvent.getByToken( trackTimeToken_, trackTimeValueMap );
   iEvent.getByToken( trackTimeErrToken_, trackTimeErrValueMap );
   iEvent.getByToken( trackFastSimTimeToken_, trackFastSimTimeValueMap );
   iEvent.getByToken( trackFastSimTimeErrToken_, trackFastSimTimeErrValueMap );


   //---get truth PV
   const SimVertex *genPV = NULL;// = SimVertex();
   
   if(genVertexHandle_.isValid()){
     const vector<SimVertex>& genVertices = *genVertexHandle_;   
     genPV = &(genVertices.at(0));   
   }
   else {
     std::cout << "oh actually sometimes we don't find a gen pv" << std::endl;
     return;
     /*auto xyz = genXYZHandle_.product();
     auto t = *genT0Handle_.product();
     auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
     //genPV = const_cast<SimVertex&>(v, t);
     genPV = SimVertex(v, t);
     */
   }

      
   int vtx_index_3D = -1;
   int vtx_index_4D = -1;
   /* 
   //---find the vertex this muon is best associated to..
   double min_dz = std::numeric_limits<double>::max();
   //double min_dzdt = std::numeric_limits<double>::max();
   for( unsigned i = 0; i < vertex3DHandle_->size(); ++i ) {
       const auto& vtx = (*vertex3DHandle_)[i];
       const float dz = std::abs(vtx.z() - genPV->position().z());
       if( dz < min_dz )
	     {
	       min_dz = dz;
	       vtx_index_3D = i;
	     }
   }

   min_dz = std::numeric_limits<double>::max();
   //float mindist = std::numeric_limits<float>::max();
   for( unsigned i = 0; i < vertex4DHandle_->size(); ++i ) {
       const auto& vtx = (*vertex4DHandle_)[i];
       float dz = std::abs(vtx.z() - genPV->position().z());
       //const float dz = std::abs(vtx.z() - genPV->position().z())/vtx.zError();
       //if ( vtx.tError() <=0 ) continue;
       //float dt = std::abs(vtx.t() - genPV->position().t()*1000000000.) / vtx.tError();
       //float dist = sqrt(dz*dz+dt*dt);
       //if( dist < mindist )
	     if( dz < min_dz )
       {
	       //mindist = dist;
	       min_dz = dz;
         vtx_index_4D = i;
	     }
   }*/
       // if( dz < 0.1 )
       // {
       //     if(isTimingSample_)
       //     {
       //         const double dzdt = pow((vtx.z() - genPV.position().z())/vtx.zError(), 2) +
       //             pow((vtx.t()-genPV.position().t())/vtx.tError(), 2);                            
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
   
   
   if (vtx_index_3D == -1) // no good PV found
     vtx_index_3D = 0;
   if (vtx_index_4D == -1)
     vtx_index_4D = 0;
   vertex3D = (*vertex3DHandle_)[vtx_index_3D];
   vertex4D = (*vertex4DHandle_)[vtx_index_4D];
           
   //vertex = (*genPV);
   //vertex = (genVertexHandle_.product()->at(0));
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
   /*
   iEvent.getByToken(tracksToken_,tracksHandle_);
   auto tracks = *tracksHandle_.product();

   for(unsigned iTrack = 0; iTrack < tracks.size(); ++iTrack)   {
     auto track = tracks.at(iTrack);
   }

   */


   // *** 3. Load generator MC particles
   int decayMode = -1;
   iEvent.getByToken(genToken_,genHandle_);
   getPromptTruthMuons( genHandle_ );

   if (!isZmumuSignal_) // background ttbar
     decayMode = ttbarDecayMode( genHandle_ );
   h_ttbarDecayMode_->Fill( decayMode );
   
   
   // *** 4. Load reco muons
   iEvent.getByToken(muonsToken_, muonsHandle_);
   iEvent.getByToken(pfCandToken_, pfCandHandle_);
   auto muons = *muonsHandle_.product();
   bool hasGoodMuon = false;
   nPromptMuons = 0;
   nNonPromptMuons = 0;

   for(unsigned iMuon = 0; iMuon < muons.size(); ++iMuon)   {
     auto muon = muons.at(iMuon);
     
     if ( !isGoodMuon( muon, genHandle_, genJetHandle_, genPV )) continue;
     h_event_nPromptMuons_->Fill( nPromptMuons ) ;
     h_event_nNonPromptMuons_->Fill( nNonPromptMuons ) ;
     h_muon_pT->Fill(muon.pt());     

     // ** A. "Analysis-like" isolation borrowed from ttH
     float pfRelIso03 = getMuonPFRelIso( muon );
     // ** B. "TDR-like" isolation using PFCandidates
     std::vector<float> pfCandIso03 = {0, 0, 0, 0}; // 0th: nominal isolation, 1st: no dxy cut, 2nd: nominal + 3sigma timing cut, 3rd: no dxy + 3sigma timing cut
     int nPFCandidatesInCone = getMuonPFCandIso( muon, pfCandHandle_, genPV, pfCandIso03, trackFastSimTimeValueMap, trackFastSimTimeErrValueMap, 0.040 );
     //if (nPFCandidatesInCone == 0) continue;      // FIXME: come up with way to skip filling if pfCand sum is 0 --> this may be what causes slightly higher inefficiency


     // Fill information about muon PF relIso (R=0.3)
     if ( fabs(muon.eta()) < 1.5){
       h_muon_pfRelIso03_BTL_->Fill( pfRelIso03 );
       //if ( pfCandIso03>=0 )
       h_muon_pfCandIso03_BTL_->Fill( pfCandIso03.at(0) );
       h_muon_pfCandIso03_noDxy_BTL_->Fill( pfCandIso03.at(1) );
       h_muon_pfCandIso03_dt30_BTL_->Fill( pfCandIso03.at(2) );
       h_muon_pfCandIso03_dt30_noDxy_BTL_->Fill( pfCandIso03.at(3) );
     }
     else if ( fabs(muon.eta()) > 1.5 && fabs(muon.eta()) < 2.8){
       h_muon_pfRelIso03_ETL_->Fill( pfRelIso03 );
       //if ( pfCandIso03>=0 )
       h_muon_pfCandIso03_ETL_->Fill( pfCandIso03.at(0) );
       h_muon_pfCandIso03_noDxy_ETL_->Fill( pfCandIso03.at(1) );
       h_muon_pfCandIso03_dt30_ETL_->Fill( pfCandIso03.at(2) );
       h_muon_pfCandIso03_dt30_noDxy_ETL_->Fill( pfCandIso03.at(3) );
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
  h_event_nPromptMuons_ = fileService->make<TH1D>("h_event_nPromptMuons", "h_event_nPromptMuons", 7, 0, 7);
  h_event_nNonPromptMuons_ = fileService->make<TH1D>("h_event_nNonPromptMuons", "h_event_nNonPromptMuons", 7, 0, 7);

  h_ttbarDecayMode_ = fileService->make<TH1D>("h_ttbarDecayMode", "h_ttbarDecayMode", 7, 0, 7);
  h_ttbarDecayMode_vs_nPromptMuons_ = fileService->make<TH2D>("h_ttbarDecayMode_vs_nPromptMuons", "h_ttbarDecayMode_vs_nPromptMuons", 7, 0, 7, 4, 0, 4);

  h_muon_cutflow_ = fileService->make<TH1D>("h_muon_cutflow", "h_muon_cutflow", 7, 0, 7);
  h_muon_dxy_BTL_ = fileService->make<TH1F>("h_muon_dxy_BTL", "h_muon_dxy_BTL", 1000, 0, 1.0);
  h_muon_dz_BTL_  = fileService->make<TH1F>("h_muon_dz_BTL", "h_muon_dz_BTL", 1000, 0, 1.0);
  h_muon_pfRelIso03_BTL_ = fileService->make<TH1F>("h_muon_pfRelIso03_BTL", "h_muon_pfRelIso03_BTL", 250, 0, 2.0);
  h_muon_pfRelIso03_ETL_ = fileService->make<TH1F>("h_muon_pfRelIso03_ETL", "h_muon_pfRelIso03_ETL", 250, 0, 2.0);
  h_muon_pfCandIso03_BTL_ = fileService->make<TH1F>("h_muon_pfCandIso03_BTL", "h_muon_pfCandIso03_BTL", 250, 0, 2.0);
  h_muon_pfCandIso03_ETL_ = fileService->make<TH1F>("h_muon_pfCandIso03_ETL", "h_muon_pfCandIso03_ETL", 250, 0, 2.0);
  h_muon_pfCandIso03_noDxy_BTL_ = fileService->make<TH1F>("h_muon_pfCandIso03_noDxy_BTL", "h_muon_pfCandIso03_noDxy_BTL", 5000, 0, 5.0);
  h_muon_pfCandIso03_noDxy_ETL_ = fileService->make<TH1F>("h_muon_pfCandIso03_noDxy_ETL", "h_muon_pfCandIso03_noDxy_ETL", 5000, 0, 5.0);

  h_muon_pfCandIso03_dt30_BTL_ = fileService->make<TH1F>("h_muon_pfCandIso03_dt30_BTL", "h_muon_pfCandIso03_dt30_BTL", 250, 0, 2.0);
  h_muon_pfCandIso03_dt30_ETL_ = fileService->make<TH1F>("h_muon_pfCandIso03_dt30_ETL", "h_muon_pfCandIso03_dt30_ETL", 250, 0, 2.0);
  h_muon_pfCandIso03_dt30_noDxy_BTL_ = fileService->make<TH1F>("h_muon_pfCandIso03_dt30_noDxy_BTL", "h_muon_pfCandIso03_dt30_noDxy_BTL", 5000, 0, 5.0);
  h_muon_pfCandIso03_dt30_noDxy_ETL_ = fileService->make<TH1F>("h_muon_pfCandIso03_dt30_noDxy_ETL", "h_muon_pfCandIso03_dt30_noDxy_ETL", 5000, 0, 5.0);

  h_muon_pfCandIso03_dxy_BTL_ = fileService->make<TH1F>("h_muon_pfCandIso03_dxy_BTL", "h_muon_pfCandIso03_dxy_BTL", 500, 0, 0.5);
  h_muon_pfCandIso03_dz_BTL_  = fileService->make<TH1F>("h_muon_pfCandIso03_dz_BTL", "h_muon_pfCandIso03_dz_BTL", 500, 0, 0.5);

  h_muon_pfCandpT = fileService->make<TH1F>("h_muon_pfCandpT", "h_muon_pfCandpT", 500, 0, 100.0);
  h_muon_pT = fileService->make<TH1F>("h_muon_pT", "h_muon_pT", 500, 0, 100.0);
  h_muon_numCand = fileService->make<TH1F>("h_muon_numCand", "h_muon_numCand", 50, 0, 50.0);

  h_pfCandidate_cutflow_ = fileService->make<TH1D>("h_pfCandidate_cutflow", "h_pfCandidate_cutflow", 7, 0, 7);

}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonIsolationAnalyzer::endJob()
{
}

int MuonIsolationAnalyzer::getMuonPFCandIso(const reco::Muon& iMuon, edm::Handle<std::vector<reco::PFCandidate> >& pfCandHandle, const SimVertex *genPV, std::vector<float> &isolations, edm::Handle<ValueMap<float> > trackFastSimTimeValueMap, edm::Handle<ValueMap<float> > trackFastSimTimeErrValueMap, double timeResolution) const
{
  //float sumPFCandPtInCone = 0.0; 
  float isoCone = 0.3;
  int numberOfAssociatedPFCandidates = 0;
  bool thisCandPassesDxy = false;
  
  auto pfCandidates = *pfCandHandle.product();

   for(unsigned iPFCand = 0; iPFCand < pfCandidates.size(); ++iPFCand)   {

     auto pfCandidate = pfCandidates.at(iPFCand);
     h_muon_pfCandpT->Fill(pfCandidate.pt());
     h_pfCandidate_cutflow_->Fill("PF Candidate", 1);
     thisCandPassesDxy = false;

     // skip neutrals
     if (pfCandidate.charge()==0) 
       continue;
     h_pfCandidate_cutflow_->Fill("Charge != 0", 1);

     reco::TrackRef pfTrack = pfCandidate.trackRef();
     
     // keep PFCands with good tracks
     if(pfTrack.isNull())
       continue;
     if(!(pfTrack->quality(reco::TrackBase::highPurity)))
       continue;
     h_pfCandidate_cutflow_->Fill("HighPurity Track", 1);

     // reject PFCands matching muon track
     if ( pfTrack == iMuon.track() )
       continue;
     h_pfCandidate_cutflow_->Fill("no muonRef match", 1);

     // calculate dxy/dz 
     float dz_sim  = std::abs( pfTrack->vz() - genPV->position().z() ); 
     float dxy_sim = sqrt ( pow(pfTrack->vx() - genPV->position().x(),2) + pow(pfTrack->vy() - genPV->position().y(),2) ); 
     float dz4D = std::abs( pfTrack->dz(vertex4D.position()) );
     float dz3D = std::abs( pfTrack->dz(vertex3D.position()) );
     float dzmu = std::abs( pfTrack->dz(vertex4D.position()) - iMuon.track()->dz(vertex4D.position()) );
    
     
     if ( dz_sim > dz_pfCandVertex )
       continue;
     h_pfCandidate_cutflow_->Fill("dz < 0.1", 1);
     
     if ( dxy_sim < dxy_pfCandVertex ){
       thisCandPassesDxy = true;
       h_pfCandidate_cutflow_->Fill("dxy < 0.02", 1);
     }
     else
       thisCandPassesDxy = false;
     
     // kinematic cuts on PF candidate
     if ( fabs(pfCandidate.eta())<1.5){ // BTL acceptance
     //if ( fabs(pfCandidate.eta())<1.48){ // BTL acceptance
       h_pfCandidate_cutflow_->Fill("In BTL Volume", 1);
     
       if (pfTrack->pt() < 0.7)
	       continue;
       h_pfCandidate_cutflow_->Fill("pT > 0.7 GeV", 1);
       
       h_muon_pfCandIso03_dxy_BTL_->Fill( dxy_sim );
       h_muon_pfCandIso03_dz_BTL_->Fill( dz_sim );
       
     }
     //else if ( fabs(track->eta())>1.5 && fabs(track->eta())<2.8) { // ETL acceptance
     else if ( fabs(pfCandidate.eta())>1.5 && fabs(pfCandidate.eta())<2.8) { // ETL acceptance
     //else if ( fabs(pfCandidate.eta())>1.48 && fabs(pfCandidate.eta())<2.8) { // ETL acceptance
       if (pfTrack->pt() < 0.4)
	       continue;       
     }
     else
       continue;
     
     // sum candidates in cone
     if ( passDeltaR( isoCone, iMuon.eta(), iMuon.phi(), pfCandidate.eta(), pfCandidate.phi()) ) {
	 numberOfAssociatedPFCandidates++;

       if ( timeResolution!=-1) { // timeResolution either passed by user (!=-1) or defaul (==-1, don't use)
   double pfcandtimeFastSim = (*trackFastSimTimeValueMap)[pfTrack];
	 double pfcandtimeErrFastSim = (*trackFastSimTimeErrValueMap)[pfTrack];
   //double pfcandtimeFastSim = (*trackTimeValueMap)[pfTrack];
	 //double pfcandtimeErrFastSim = (*trackTimeErrValueMap)[pfTrack];
   TRandom *gRandom = new TRandom();
   TRandom *gRandom2 = new TRandom();
   TRandom *gRandom3 = new TRandom();
   double targetTimeResol = timeResolution;
	 //double defaultTimeResol  = 0.;
   double defaultTimeResolFastSim  = pfcandtimeErrFastSim;
	 //if ( pfCandidate.isTimeValid() ) 
	 //  defaultTimeResol = 0.035; // only using mtd5 samples I think..?

	 //double extra_resol = 0.;
	 double extra_resol_FastSim = 0.;
   //if ( targetTimeResol > defaultTimeResol) 
	 //  extra_resol = sqrt(targetTimeResol*targetTimeResol - defaultTimeResol*defaultTimeResol); 

   if ( targetTimeResol > defaultTimeResolFastSim)
     extra_resol_FastSim = sqrt(targetTimeResol*targetTimeResol - defaultTimeResolFastSim*defaultTimeResolFastSim);


	 double dtsim = 0.;
	 double pfCandidateTime = -999.;
	 if(1) {
   //if ( pfCandidate.isTimeValid() && !isnan( pfCandidate.time() )) {
	   // -- emulate BTL and ETL efficiency
	   bool keepTrack = true;

	   // introduce inefficiency loss to mirror more realistic detector behaviour, [BBT 08-23-19: COMMENT OUT FOR NOW]
	   double rndEff = gRandom2->Uniform(0.,1.);
     //std::cout << "eff: " << rndEff << std::endl;
     //std::cout << "KeepTrack: " << (int) keepTrack << std::endl;
	   if ( std::abs(pfCandidate.eta()) < 1.5 && rndEff > btlEfficiency ) keepTrack = false; 
	   if ( std::abs(pfCandidate.eta()) > 1.5 && rndEff > etlEfficiency ) keepTrack = false; 

     if ( pfcandtimeErrFastSim !=-1 ){
       double rndFastSim = gRandom->Gaus(0., extra_resol_FastSim);
       pfCandidateTime = pfcandtimeFastSim + rndFastSim;
       //std::cout << "rndFastSim = " << rndFastSim;
     }
	   if ((keepTrack) && (pfCandidateTime!=-999) ) {
	     // -- extra smearing to emulate different time resolution
	     
       //double rnd   = gRandom->Gaus(0., extra_resol);
	     //pfCandidateTime = pfCandidate.time() + rnd;
	     // extra extra smearing 
       double extra_smearing = sqrt(std::abs(targetTimeResol*targetTimeResol - 0.035*0.035));
       pfCandidateTime = pfCandidateTime + gRandom3->Gaus(0,extra_smearing);
       dtsim = std::abs(pfCandidateTime - genPV->position().t()*1000000000.);
	     //cout << "  target time resol = "<< targetTimeResol << "  extra_resol = "<< extra_resol_FastSim << "  extra rnd = " << extra_smearing << "  pfCandidate time = " << pfCandidateTime << "  dtsim = " << dtsim << endl;
	   }
	   else
	     dtsim = 0.;
	 }
	 else
	   dtsim = 0.;

	 // *** A. Fill all candidates
	 // ** 0. Nominal
	 if(thisCandPassesDxy)
	   isolations.at(0) += pfCandidate.pt();
	 // ** 1. No Dxy
	 isolations.at(1) += pfCandidate.pt(); 
	 // *** B. Fill candidate times if within 3 sigma
	 if (dtsim < 3.*targetTimeResol){
	   // ** 2. Nominal + 3sigma cut on timing
	   if(thisCandPassesDxy)
	     isolations.at(2) += pfCandidate.pt(); // nominal + 3sigma timing
	   // ** 3. No Dxy + 3sigma cut on timing
	   isolations.at(3) += pfCandidate.pt();  // no Dxy + 3sigma timing
	 }

       } // end timing
       else { // no check on timing
	 // ** 0. Nominal
	 if(thisCandPassesDxy)
	   isolations.at(0) += pfCandidate.pt();
	 // ** 1. No Dxy
	 isolations.at(1) += pfCandidate.pt(); 
	 
       }
     }// end of if condition for within isoCone
     
   } // PF Candidate loop

   //std::cout << "Isolation = " << result << " with " << numberOfAssociatedPFCandidates << " associated to muon" << std::endl;
   h_muon_numCand->Fill(numberOfAssociatedPFCandidates);
   
   // divide all sum PFCand in cone pT by muon pT
   for(unsigned int i=0; i < isolations.size(); i++)
     isolations.at(i) = isolations.at(i) / iMuon.pt();

   return numberOfAssociatedPFCandidates;
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

void MuonIsolationAnalyzer::getPromptTruthMuons(edm::Handle<std::vector<reco::GenParticle> >& genHandle)
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
  if ( (numberOfPromptMuons > 2) && (isZmumuSignal_) )
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
  unsigned int numberOfPromptMuons = promptMuonTruthCandidates_.size();
  //promptMuonTruthCandidates_ = {};

  // *** 1. Loop over MC particles
  for(size_t iGenParticle=0; iGenParticle<genHandle->size();iGenParticle++){

    // ** A. Prompt muons
    //getPromptTruthMuons( genHandle_ );

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

bool MuonIsolationAnalyzer::isGoodMuon(const reco::Muon& iMuon, edm::Handle<std::vector<reco::GenParticle> >& genHandle, edm::Handle<std::vector<reco::GenJet> >& genJetHandle, const SimVertex *genPV) const
{
  h_muon_cutflow_->Fill("All Muons", 1);
  
  // *** just check to make sure muon track available
  if (iMuon.pt() < 5.)
    return false;
  if (iMuon.track().isNull()) 
    return false;
  h_muon_cutflow_->Fill("track", 1);
  
  // *** 0. pT > 20 GeV ---> may need to change this for Bs->mumu studies, but keep for now to reproduce TDR results
  if (iMuon.pt() < 20.)
    return false;
  h_muon_cutflow_->Fill("pT > 20", 1);

  if (iMuon.pt() > 9999.)
    return false;
  h_muon_cutflow_->Fill("pT < 9999", 1);

  // *** 2. Loose ID
  //if ( iMuon.passed('CutBasedIdLoose'))//userFloat("CutBasedIdLoose") ) // CutBasedIdLoose, MvaLoose, others?
  if ( !(muon::isLooseMuon(iMuon)) )
    return false;
  h_muon_cutflow_->Fill("Loose ID", 1);
  
  // *** 1. |eta| < 2.4
  //if ( fabs(iMuon.eta()) > 2.4 )
  //  return false;
  //h_muon_cutflow_->Fill("|eta| < 2.4", 1);
  


  reco::TrackRef muonTrack = iMuon.track();
  float dz_sim  = std::abs( muonTrack->vz() - genPV->position().z() ); 
  float dxy_sim = sqrt ( pow(muonTrack->vx() - genPV->position().x(),2) + pow(muonTrack->vy() - genPV->position().y(),2) ); 

  float dz_4D = std::abs( muonTrack->dz(vertex4D.position()));
  float dxy_4D = std::abs( muonTrack->dxy(vertex4D.position()));
  // some plots
  if (fabs(iMuon.eta())<1.5){
    h_muon_dxy_BTL_->Fill( dxy_sim );
    h_muon_dz_BTL_->Fill( dz_sim );
  }

  // *** 4. z0 cut
  if ( dz_4D > dz_muonVertex )
    return false;
  h_muon_cutflow_->Fill("z0 < 0.5 cm", 1);
  
  // *** 3. d0 cut
  if ( dxy_4D > dxy_muonVertex )
    return false;
  h_muon_cutflow_->Fill("d0 < 0.2 cm", 1);
  
  // end

  if (vertex4D.isFake() || vertex3D.isFake() )
    return false;
  h_muon_cutflow_->Fill("nofake", 1);
  
  // calculate some booleans
  bool recoMuonMatchedToPromptTruth = isPromptMuon(iMuon, genHandle);
  bool recoMuonMatchedToGenJet      = isMatchedToGenJet(iMuon, genJetHandle);
  bool recoMuonFromTau              = isFromTau(iMuon, genHandle);

  // *** 5A. accept only prompt muons if Z->mumu signal
  if (isZmumuSignal_) {

    if (!recoMuonMatchedToPromptTruth)
      return false;
    h_muon_cutflow_->Fill("Signal Prompt Muon", 1);
  }
  // *** 5B. reject prompt muons if ttbar background
  else if (!isZmumuSignal_) {

    // ** . Muon is "good" non-prompt if !truthMatched && genJetMatched && !tauMatched
    if ( !recoMuonMatchedToPromptTruth && recoMuonMatchedToGenJet && !recoMuonFromTau) {
      h_muon_cutflow_->Fill("Non-prompt Bkg Muon", 1);
    }
    else 
      return false;
  }

  if (fabs(vertex4D.z()-genPV->position().z()) > 0.01 ) return false;
  h_muon_cutflow_->Fill("4D_gen", 1);
  if (fabs(vertex3D.z()-genPV->position().z()) > 0.01 ) return false;
  h_muon_cutflow_->Fill("3D_gen", 1);

  return true;

}

// ------------ function to one-line-ize deltaR logic ------------
bool MuonIsolationAnalyzer::passDeltaR( float coneSize, float eta1, float phi1, float eta2, float phi2) const
{
  
  float d_eta = eta1 - eta2;
  float d_phi = phi1 - phi2;
  float dR = sqrt( d_eta*d_eta + d_phi*d_phi);            

  if (dR > 0.0 && dR < coneSize)
    return true;
  else
    return false;

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
