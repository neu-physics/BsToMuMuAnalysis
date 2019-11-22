

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "FWCore/Framework/interface/IOVSyncValue.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/Common/interface/Handle.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "BsToMuMuAnalysis/MuonIsolationAnalyzer/interface/MuonIsolationAnalyzer.h"


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
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("BSTag"))),
  trackTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeValueMapTag" ) ) ),
  trackTimeErrToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackTimeErrValueMapTag" ) ) ),
  trackFastSimTimeToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackFastSimTimeValueMapTag" ) ) ),
  trackFastSimTimeErrToken_( consumes<ValueMap<float> >( iConfig.getParameter<InputTag>( "TrackFastSimTimeErrValueMapTag" ) ) ),
  triggerResultsToken_(consumes<edm::TriggerResults> (iConfig.getUntrackedParameter<InputTag>("TriggerTag"))),
  isZmumuSignal_(iConfig.getParameter<bool>("isZmumuSignal")),
  isBmumu_(iConfig.getParameter<bool>("isBmumu")),
  processName_(iConfig.getUntrackedParameter<string>("processName"))
  //lastFourFileID_(iConfig.getParameter<string>("lastFourFileID"))
{

   //now do what ever initialization is needed
  evInfo = new eventInfo;
  gRandom = new TRandom();
  gRandom2 = new TRandom();
  gRandom3 = new TRandom();
  //---outputs

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
  
   int n_evt = iEvent.id().event();
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

   // *** Load Trigger paths (not used for now)
   iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
   if (!triggerResultsHandle_.isValid()) {
     cout << "****Error in getting TriggerResults product from Event!" << endl;
     return;
   }
   std::vector<std::string> triggerNames = hltConfig_.triggerNames();
   for (unsigned int iPath=0; iPath<triggerNames.size();++iPath){
     std::string pathName = triggerNames[iPath];
     //std::cout << "HLT path name: " << pathName << std::endl;
   }

   //****load magnetic field
   ESHandle<MagneticField> bFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
   
   //****load BeamSpot
   reco::BeamSpot beamSpot;
   iEvent.getByToken(beamSpotToken_, beamSpotHandle_);
   if ( beamSpotHandle_.isValid() )
   {
     beamSpot = *beamSpotHandle_;
   } 
   else
   {
     std::cout << "No beam spot available from EventSetup " << std::endl;
   }


   initEventStructure();
   //---get truth PV
   const SimVertex *genPV = NULL;// = SimVertex();
   
   if(genVertexHandle_.isValid()){
     const vector<SimVertex>& genVertices = *genVertexHandle_;   
     genPV = &(genVertices.at(0));   
   }
   else {
     std::cout << "oh actually sometimes we don't find a gen pv" << std::endl;
     return;
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
   
   
   if(!isBmumu_){
     if (vtx_index_3D == -1) // no good PV found
       vtx_index_3D = 0;
     if (vtx_index_4D == -1)
       vtx_index_4D = 0;
     vertex3D = (*vertex3DHandle_)[vtx_index_3D];
     vertex4D = (*vertex4DHandle_)[vtx_index_4D];
             
     //vertex = (*genPV);
     //vertex = (genVertexHandle_.product()->at(0));
     

     if(vertex4D.isFake())
       return;
     h_event_cutflow_->Fill("4D vtx not fake", 1);
     if (vertex3D.isFake())
       return;
     h_event_cutflow_->Fill("3D vtx not fake", 1); 
     
     if (fabs(vertex4D.z()-genPV->position().z()) > 0.01 )
       return;
     h_event_cutflow_->Fill("dz(4D vtx,genPV)<0.01cm", 1);
     
     if (fabs(vertex3D.z()-genPV->position().z()) > 0.01 )
       return;
     h_event_cutflow_->Fill("dz(3D vtx,genPV)<0.01cm", 1);
   }


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
   //if ((!isZmumuSignal_)&&(!isBmumu_))
   //  FindtWbMuons(genHandle_);

   //for(size_t iGenParticle=0; iGenParticle<genHandle_->size();iGenParticle++){
   //  //if((fabs((*genHandle_)[iGenParticle].pdgId())>=500) && (fabs((*genHandle_)[iGenParticle].pdgId())<1000)){
   //  if(fabs((*genHandle_)[iGenParticle].pdgId())==13){
   //    const reco::Candidate * genCandidate = &(*genHandle_)[iGenParticle] ;
   //    std::cout << "PdgId: " << (*genHandle_)[iGenParticle].pdgId() << " status: " << (*genHandle_)[iGenParticle].status() << " # mothers: " << genCandidate->numberOfMothers() << std::endl;
   //    FindMother(genCandidate,1);
   //    //for (size_t iMother=0; iMother<genCandidate->numberOfMothers(); iMother++){
   //    //  std::cout << "   mother: " << genCandidate->mother(iMother)->pdgId() << std::endl;
   //    //}
   // 
   //  }
   //}

   //auto BDaughterVector = BDaughters(genHandle_);

   if (!isZmumuSignal_ && !isBmumu_) // background ttbar
     decayMode = ttbarDecayMode( genHandle_ );
   h_ttbarDecayMode_->Fill( decayMode );
   
   
   // *** 4. Load reco muons
   iEvent.getByToken(muonsToken_, muonsHandle_);
   iEvent.getByToken(pfCandToken_, pfCandHandle_);
   auto muons = *muonsHandle_.product();
   bool hasGoodMuon = false;
   //bool MatchBcand = false;
   //bool m_inv_match = false;
   double pT_max_pos = 0;
   double pT_max_neg = 0;
   reco::Muon mu1;
   reco::Muon mu2;
   
   nPromptMuons = 0;
   nNonPromptMuons = 0;
   //std::cout << "# muons: " <<  muons.size() << std::endl;
   for(unsigned iMuon = 0; iMuon < muons.size(); ++iMuon)   {
     auto muon = muons.at(iMuon);
     if ( !isGoodMuon( muon, genHandle_, genJetHandle_, genPV )) continue;
     //std::cout << "muon: " << muon.pt() << " charge: " << muon.charge() << std::endl;
     if (isBmumu_){
       if ((muon.charge()==1) && (muon.pt()>pT_max_pos)){
         pT_max_pos = muon.pt();
         mu1 = muon;
       }
       if ((muon.charge()==-1) && (muon.pt()>pT_max_neg)){
         pT_max_neg = muon.pt();
         mu2 = muon;
       }
     }

     h_event_nPromptMuons_->Fill( nPromptMuons ) ;
     h_event_nNonPromptMuons_->Fill( nNonPromptMuons ) ;
     h_muon_pT->Fill(muon.pt());     

     // ** A. "Analysis-like" isolation borrowed from ttH
     float pfRelIso03 = getMuonPFRelIso( muon );
     // ** B. "TDR-like" isolation using PFCandidates
     std::vector<float> pfCandIso03 = {0, 0, 0, 0}; // 0th: nominal isolation, 1st: no dxy cut, 2nd: nominal + 3sigma timing cut, 3rd: no dxy + 3sigma timing cut
     //int nPFCandidatesInCone = getMuonPFCandIso( muon, pfCandHandle_, genPV, pfCandIso03, trackFastSimTimeValueMap, trackFastSimTimeErrValueMap, 0.040, n_evt );
     getMuonPFCandIso( muon, pfCandHandle_, genPV, pfCandIso03, trackFastSimTimeValueMap, trackFastSimTimeErrValueMap, 0.040, n_evt );
     //if (nPFCandidatesInCone == 0) continue;      // FIXME: come up with way to skip filling if pfCand sum is 0 --> this may be what causes slightly higher inefficiency


     // Fill information about muon PF relIso (R=0.3)
     evInfo->nevt.push_back(n_evt);
     evInfo->muon_eta.push_back(muon.eta());
     evInfo->muon_pt.push_back(muon.pt());
     evInfo->muon_pfCand.push_back(pfCandIso03.at(0));
     evInfo->muon_pfCand_noDxy.push_back(pfCandIso03.at(1));
     evInfo->muon_pfCand_dt.push_back(pfCandIso03.at(2));
     evInfo->muon_pfCand_noDxy_dt.push_back(pfCandIso03.at(3));
     //std::cout << "nevt: " << n_evt << " eta: " << muon.eta() << " pt: " << muon.pt() << " iso: " << pfCandIso03.at(1) << " isodt: " << pfCandIso03.at(3) << std::endl;
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

   if(isBmumu_){
     if( (mu1.charge()==0) || (mu2.charge()==0) )  return;
     h_event_cutflow_->Fill("two opposite charged Muons",1);
     //calculate dca between leading and sub-leading opposite charged muons
     double dca = 0;
     GlobalPoint cp; //crossing point of two transient tracks
     std::pair<GlobalPoint,GlobalPoint> pts;
     reco::TrackRef tk1 = mu1.track();
     reco::TrackRef tk2 = mu2.track();
     //std::cout << "tracks: " << tk1->referencePoint() << " " << tk2->referencePoint() << std::endl;
     TransientTrack mu1TT(*tk1, &(*bFieldHandle));
     TransientTrack mu2TT(*tk2, &(*bFieldHandle));

     // ***************** vertex selection ***************************
     // refit the secondary vertex
     float mu_mass = 0.1057;
     float mu_mass_sigma = 4E-9f;
     KinematicParticleFactoryFromTransientTrack pFactory;
     vector<RefCountedKinematicParticle> kinParticles; //the container that contains the two muons(daughter of B meson)
     kinParticles.push_back(pFactory.particle(mu1TT,mu_mass,0.0f,0.0f,mu_mass_sigma));
     kinParticles.push_back(pFactory.particle(mu2TT,mu_mass,0.0f,0.0f,mu_mass_sigma));
     RefCountedKinematicTree kinTree;
     KinematicParticleVertexFitter kpvFitter;
     kinTree = kpvFitter.fit(kinParticles); // do the fit to find the secondary vertex
     RefCountedKinematicParticle B_kinPart = kinTree->currentParticle();

     AnalyticalImpactPointExtrapolator extrapolator(&(*bFieldHandle));
     double minDist = 999;
     // refit the PV using tracks but exclude the muons tracks because they should come from the secondary vertex
     for(size_t i = 0; i< vertex3DHandle_->size(); ++i){
       std::cout << "3D vtx: " << vertex3DHandle_->size() << std::endl;
       const auto& vtx = (*vertex3DHandle_)[i];
       std::pair<bool,Measurement1D> DecayL;
       AdaptiveVertexFitter avf;
       vector<TransientTrack> vrtxRefit;
       //std::cout << "Begin looping over all tracks: " << vtx.tracksSize() << std::endl;
       for (vector<TrackBaseRef>::const_iterator tIt=vtx.tracks_begin();tIt!=vtx.tracks_end();++tIt){
         TrackRef tref = tIt->castTo<TrackRef>();
         //std::cout << " given track: " << tref->charge() << std::endl;
         if( (tref==tk1) || (tref==tk2) ){
           //std::cout << "Found the same muon track :)" << std::endl;
           //std::cout << "Found the same muon track :)" << std::endl;
           continue;
         }
         //std::cout << "--not the muon track" << std::endl;
         TransientTrack trkTT(*tref, &(*bFieldHandle));
         vrtxRefit.push_back(trkTT);
       }
       if (vrtxRefit.size()<5) continue;
       //std::cout << "vtx has more than 5 tracks" << std::endl;
       TransientVertex newVtx;
       //newVtx = avf.vertex(vrtxRefit);
       newVtx = avf.vertex(vrtxRefit,beamSpot);
       if (!newVtx.isValid()) continue;
       //std::cout << "refit succeeded" << std::endl;
       Vertex rePV = reco::Vertex(newVtx); //refitted PV
       //std::cout << "transfer to vtx" << std::endl;
       TrajectoryStateOnSurface tsos = extrapolator.extrapolate(B_kinPart->currentState().freeTrajectoryState(), RecoVertex::convertPos(rePV.position()));
       //std::cout << "set up tsos" << std::endl;
       DecayL = IPTools::signedDecayLength3D(tsos, GlobalVector(0,0,1), rePV);
       //std::cout << "calc decay length" << std::endl;
       if (fabs(DecayL.second.value()) < minDist ){
         minDist = fabs(DecayL.second.value());
         vertex3D = rePV;
       }
     }

     minDist = 999;
     for(size_t i = 0; i< vertex4DHandle_->size(); ++i){
       std::cout << "4D vtx: " << vertex4DHandle_->size() << std::endl;
       const auto& vtx = (*vertex4DHandle_)[i];
       std::pair<bool,Measurement1D> DecayL;
       AdaptiveVertexFitter avf;
       vector<TransientTrack> vrtxRefit;
       for (vector<TrackBaseRef>::const_iterator tIt=vtx.tracks_begin();tIt!=vtx.tracks_end();++tIt){
         TrackRef tref = tIt->castTo<TrackRef>();
         if( (tref==tk1) || (tref==tk2) ){
           std::cout << "Found the same muon track :)" << std::endl;
           continue;
         }
         TransientTrack trkTT(*tref, &(*bFieldHandle));
         vrtxRefit.push_back(trkTT);
       }
       if (vrtxRefit.size()<5) continue;
       TransientVertex newVtx;
       newVtx = avf.vertex(vrtxRefit);
       if (!newVtx.isValid()) continue;
       Vertex rePV = reco::Vertex(newVtx); //refitted PV
       TrajectoryStateOnSurface tsos = extrapolator.extrapolate(B_kinPart->currentState().freeTrajectoryState(), RecoVertex::convertPos(rePV.position()));
       DecayL = IPTools::signedDecayLength3D(tsos, GlobalVector(0,0,1), rePV);
       if (fabs(DecayL.second.value()) < minDist ){
         minDist = fabs(DecayL.second.value());
         vertex4D = rePV;
       }
     }

     if(vertex4D.isFake())
       return;
     h_event_cutflow_->Fill("4D vtx not fake", 1);
     if (vertex3D.isFake())
       return;
     h_event_cutflow_->Fill("3D vtx not fake", 1); 
     
     if (fabs(vertex4D.z()-genPV->position().z()) > 0.01 )
       return;
     h_event_cutflow_->Fill("dz(4D vtx,genPV)<0.01cm", 1);
     
     if (fabs(vertex3D.z()-genPV->position().z()) > 0.01 )
       return;
     h_event_cutflow_->Fill("dz(3D vtx,genPV)<0.01cm", 1);

     TrajectoryStateClosestToPoint mu1TS = mu1TT.impactPointTSCP();
     TrajectoryStateClosestToPoint mu2TS = mu2TT.impactPointTSCP();
     if (mu1TS.isValid() && mu2TS.isValid()) {
       ClosestApproachInRPhi cApp;
       cApp.calculate(mu1TS.theState(), mu2TS.theState());
       if(cApp.status()){
         dca = cApp.distance();
         //pts = cApp.points();
         //std::cout << "cAPP: " << cApp.crossingPoint().x() << " " << cApp.crossingPoint().y() << " " << cApp.crossingPoint().z() << " " << cApp.distance() << std::endl;
       }
       else{
         dca = 999;
       }
     }
     else{
       dca = 999;
     }
     //std::cout << "dca: " << dca << std::endl;
     if(dca>0.08)  return;
     h_event_cutflow_->Fill("dca<0.08cm",1);
     TLorentzVector p1;
     TLorentzVector p2;
     TLorentzVector ptotal;
     p1.SetPtEtaPhiE(mu1.pt(),mu1.eta(),mu1.phi(),mu1.energy());
     p2.SetPtEtaPhiE(mu2.pt(),mu2.eta(),mu2.phi(),mu2.energy());
     ptotal = p1+p2;
     double m_inv = ptotal.M();
     if ((m_inv < 4.8) || (m_inv > 6.0))  return;
     h_event_cutflow_->Fill("4.8<m_inv<6.0",1);

     //B2MuMu isolation cutflows
     double IsoCut[3] = {0.1545,0.0875,0.0525};
     double IsoCut_200PU[3] = {0.2925,0.1855,0.1365};
     double IsoCut_MTD[3] = {0.1545,0.0875,0.0515};
     double IsoCut_200PU_MTD[3] = {0.2315,0.1395,0.0985};
     TString CutName[3] = {"95%WP(iso<0.1545)","90%WP(iso<0.0875)","85%WP(iso<0.0525)"};
     TString CutName_200PU[3] = {"95%WP(iso<0.2925)(PU)","90%WP(iso<0.1855)(PU)","85%WP(iso<0.1365)(PU)"};
     TString CutName_MTD[3] = {"95%WP(iso<0.1545)(MTD)","90%WP(iso<0.0875)(MTD)","85%WP(iso<0.0515)(MTD)"};
     TString CutName_200PU_MTD[3] = {"95%WP(iso<0.2315)(PU MTD)","90%WP(iso<0.1395)(PU MTD)","85%WP(iso<0.0985)(PU MTD)"};
     std::vector<float> mu1_Iso = {0,0,0,0};
     std::vector<float> mu2_Iso = {0,0,0,0};
     getMuonPFCandIso( mu1, pfCandHandle_, genPV, mu1_Iso, trackFastSimTimeValueMap, trackFastSimTimeErrValueMap, 0.040, n_evt );
     getMuonPFCandIso( mu2, pfCandHandle_, genPV, mu2_Iso, trackFastSimTimeValueMap, trackFastSimTimeErrValueMap, 0.040, n_evt );
     //std::cout << "anti-mu: pt: " << mu1.pt() << " charge: " << mu1.charge() << std::endl;
     //std::cout << "Isolation: ";
     //for(auto isola:mu1_Iso){
     //  std::cout << isola;
     //}
     //std::cout << std::endl;
     //std::cout << "mu: pt: " << mu2.pt() << " charge: " << mu2.charge() << std::endl;
     //std::cout << "Isolation: ";
     //for(auto isola:mu2_Iso){
     //  std::cout << isola;
     //}
     //std::cout << std::endl;


     for(int iCut=0;iCut<3;++iCut){
       if( (mu1_Iso.at(1)<IsoCut[iCut]) && (mu2_Iso.at(1)<IsoCut[iCut]) ){
         h_event_cutflow_->Fill(CutName[iCut],1);
       }
     }
     for(int iCut=0;iCut<3;++iCut){
       if( (mu1_Iso.at(3)<IsoCut_MTD[iCut]) && (mu2_Iso.at(3)<IsoCut_MTD[iCut]) ){
         h_event_cutflow_->Fill(CutName_MTD[iCut],1);
       }
     }
     for(int iCut=0;iCut<3;++iCut){
       if( (mu1_Iso.at(1)<IsoCut_200PU[iCut]) && (mu2_Iso.at(1)<IsoCut_200PU[iCut]) ){
         h_event_cutflow_->Fill(CutName_200PU[iCut],1);
       }
     }
     for(int iCut=0;iCut<3;++iCut){
       if( (mu1_Iso.at(3)<IsoCut_200PU_MTD[iCut]) && (mu2_Iso.at(3)<IsoCut_200PU_MTD[iCut]) ){
         h_event_cutflow_->Fill(CutName_200PU_MTD[iCut],1);
       }
     }

   }

  // if (MatchBcand)
  //   h_event_cutflow_->Fill("dca<0.08cm",1);
  // if (m_inv_match)
  //   h_event_cutflow_->Fill("4.8<m_inv<6.0",1);
   if ((hasGoodMuon) && (!isBmumu_))
     h_event_cutflow_->Fill(">= 1 Good Muon", 1);

   eventTree->Fill();
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
//void
//MuonIsolationAnalyzer::beginJob(const edm::Run& iRun, const edm::EventSetup& iSetup)
void MuonIsolationAnalyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");

  eventTree = fileService->make<TTree>( "tree_isolation", "tree_isolation" );
  eventTree->Branch( "nevt",                  &evInfo->nevt);
  eventTree->Branch( "muon_eta",              &evInfo->muon_eta);
  eventTree->Branch( "muon_pt",               &evInfo->muon_pt);  
  eventTree->Branch( "muon_pfCand",           &evInfo->muon_pfCand);
  eventTree->Branch( "muon_pfCand_dt",        &evInfo->muon_pfCand_dt);
  eventTree->Branch( "muon_pfCand_noDxy",     &evInfo->muon_pfCand_noDxy);
  eventTree->Branch( "muon_pfCand_noDxy_dt",  &evInfo->muon_pfCand_noDxy_dt);

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
  h_pfCandidate_cutflow_4_ = fileService->make<TH1D>("h_pfCandidate_cutflow_4", "h_pfCandidate_cutflow_4", 7, 0, 7);
  h_pfCandidate_cutflow_20_ = fileService->make<TH1D>("h_pfCandidate_cutflow_20", "h_pfCandidate_cutflow_20", 7, 0, 7);

}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonIsolationAnalyzer::endJob()
{
}

void
MuonIsolationAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (!(hltConfig_.init(iRun,iSetup,processName_,changed))){
    std::cout << "Warning, didn't find trigger process HLT,\t" << processName_ << std::endl;
    return;
  }
}

void
MuonIsolationAnalyzer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

void MuonIsolationAnalyzer::initEventStructure()
{
  evInfo->nevt.clear();
  evInfo->muon_eta.clear();
  evInfo->muon_pt.clear();
  evInfo->muon_pfCand.clear();
  evInfo->muon_pfCand_dt.clear();
  evInfo->muon_pfCand_noDxy.clear();
  evInfo->muon_pfCand_noDxy_dt.clear();
}

int MuonIsolationAnalyzer::getMuonPFCandIso(const reco::Muon& iMuon, edm::Handle<std::vector<reco::PFCandidate> >& pfCandHandle, const SimVertex *genPV, std::vector<float> &isolations, edm::Handle<ValueMap<float> > trackFastSimTimeValueMap, edm::Handle<ValueMap<float> > trackFastSimTimeErrValueMap, double timeResolution, int n_evt) const
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
     if(iMuon.pt()>4 && iMuon.pt()<20)
       h_pfCandidate_cutflow_4_->Fill("PF Candidate", 1);
     if(iMuon.pt()>20)
       h_pfCandidate_cutflow_20_->Fill("PF Candidate", 1);
     thisCandPassesDxy = false;

     // skip neutrals
     if (pfCandidate.charge()==0) 
       continue;
     h_pfCandidate_cutflow_->Fill("Charge != 0", 1);
     if(iMuon.pt()>4 && iMuon.pt()<20)
       h_pfCandidate_cutflow_4_->Fill("Charge !=0", 1);
     if(iMuon.pt()>20)
       h_pfCandidate_cutflow_20_->Fill("Charge !=0", 1);

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

     if ( std::abs(pfCandidate.eta()) < 1.48 && pfCandidate.pt() < 0.7 ) continue;
     if ( std::abs(pfCandidate.eta()) > 1.48 && pfCandidate.pt() < 0.4 ) continue;
     // calculate dxy/dz 
     float dz_sim  = std::abs( pfTrack->vz() - genPV->position().z() ); 
     float dxy_sim = sqrt ( pow(pfTrack->vx() - genPV->position().x(),2) + pow(pfTrack->vy() - genPV->position().y(),2) ); 
     if (!isBmumu_){
       float dz4D = pfTrack->dz(vertex4D.position());
       float dz3D = std::abs( pfTrack->dz(vertex3D.position()) );
       float dzmu = std::abs( pfTrack->dz(vertex4D.position()) - iMuon.track()->dz(vertex4D.position()) );
       //float dxy4D = std::abs( pfTrack->dxy( vertex4D.position() ));
       //float dz_sim = fabs(dz4D + vertex4D.z() - genPV->position().z());
       if ( ( dz_sim < 1.0 || (std::abs(dz4D) < 1.0) || dz3D < 1.0 || dzmu < 1.0) )
       {
         h_pfCandidate_cutflow_->Fill("loose dz", 1);
       }
       else
         continue;
       //h_pfCandidate_cutflow_->Fill("loose dz", 1); 
     }

     float dr = deltaR(iMuon.eta(),iMuon.phi(),pfCandidate.eta(),pfCandidate.phi());
     if (!(dr>0 && dr<isoCone))  continue;
     h_pfCandidate_cutflow_->Fill("dr", 1);

     if ( dxy_sim < dxy_pfCandVertex ){
       thisCandPassesDxy = true;
       //h_pfCandidate_cutflow_->Fill("dxy < 0.02", 1);
     }
     else
       thisCandPassesDxy = false;
     
     // kinematic cuts on PF candidate
     if ( fabs(pfCandidate.eta())<1.5){ // BTL acceptance
     //if ( fabs(pfCandidate.eta())<1.48){ // BTL acceptance
       //h_pfCandidate_cutflow_->Fill("In BTL Volume", 1);
     
       if (pfTrack->pt() < 0.7)
	       continue;
       //h_pfCandidate_cutflow_->Fill("pT > 0.7 GeV", 1);
       
       h_muon_pfCandIso03_dxy_BTL_->Fill( dxy_sim );
       h_muon_pfCandIso03_dz_BTL_->Fill( dz_sim );
       
     }
     h_pfCandidate_cutflow_->Fill("BTL pt", 1);
     if(iMuon.pt()>4 && iMuon.pt()<20)
       h_pfCandidate_cutflow_4_->Fill("BTL pt", 1);
     if(iMuon.pt()>20)
       h_pfCandidate_cutflow_20_->Fill("BTL pt", 1);
     //else if ( fabs(track->eta())>1.5 && fabs(track->eta())<2.8) { // ETL acceptance
     //else if ( fabs(pfCandidate.eta())>1.5 && fabs(pfCandidate.eta())<2.8) { // ETL acceptance
     //else if ( fabs(pfCandidate.eta())>1.48 && fabs(pfCandidate.eta())<2.8) { // ETL acceptance
     if ( fabs(pfCandidate.eta())>1.5 ) { // ETL acceptance
       if (pfTrack->pt() < 0.4)
	       continue;       
     }
     h_pfCandidate_cutflow_->Fill("ETL pt", 1);
     
     if ( dz_sim >= dz_pfCandVertex )
       continue;
     if ( fabs(iMuon.eta()) < 1.5 ){
       h_pfCandidate_cutflow_->Fill("dz", 1);
       if(iMuon.pt()>4 && iMuon.pt()<20)
         h_pfCandidate_cutflow_4_->Fill("dz", 1);
       if(iMuon.pt()>20)
         h_pfCandidate_cutflow_20_->Fill("dz", 1);
     }

     // sum candidates in cone
     //if ( passDeltaR( isoCone, iMuon.eta(), iMuon.phi(), pfCandidate.eta(), pfCandidate.phi()) ) {
	   numberOfAssociatedPFCandidates++;

     if ( timeResolution!=-1) { // timeResolution either passed by user (!=-1) or defaul (==-1, don't use)
       double pfcandtimeFastSim = (*trackFastSimTimeValueMap)[pfTrack];
	     double pfcandtimeErrFastSim = (*trackFastSimTimeErrValueMap)[pfTrack];
       double targetTimeResol = timeResolution;
	     //double defaultTimeResol  = 0.;
       double defaultTimeResolFastSim  = pfcandtimeErrFastSim;
	     //if ( pfCandidate.isTimeValid() ) 
	     //  defaultTimeResol = 0.035; // only using mtd5 samples I think..?

	     //double extra_resol = 0.;
	     double extra_resol_FastSim = 0.;
       //if ( targetTimeResol > defaultTimeResol) 
	     //  extra_resol = sqrt(targetTimeResol*targetTimeResol - defaultTimeResol*defaultTimeResol); 

       //if ( targetTimeResol > defaultTimeResolFastSim)
       //  extra_resol_FastSim = sqrt(targetTimeResol*targetTimeResol - defaultTimeResolFastSim*defaultTimeResolFastSim);

       if ( defaultTimeResolFastSim < 0.035 )
         extra_resol_FastSim = sqrt( 0.035*0.035 - defaultTimeResolFastSim*defaultTimeResolFastSim);

	     double dtsim = 0.;
	     double pfCandidateTime = -999.;
       //if ( pfCandidate.isTimeValid() && !isnan( pfCandidate.time() )) {
	     // -- emulate BTL and ETL efficiency
	     bool keepTrack = true;

	     // introduce inefficiency loss to mirror more realistic detector behaviour, [BBT 08-23-19: COMMENT OUT FOR NOW]
       double rndEff = gRandom2->Uniform(0.,1.);
	     if ( std::abs(pfCandidate.eta()) < 1.5 && rndEff > btlEfficiency ) keepTrack = false; 
	     if ( std::abs(pfCandidate.eta()) > 1.5 && rndEff > etlEfficiency ) keepTrack = false; 

       if ( pfcandtimeErrFastSim !=-1 ){
         double rndFastSim = gRandom->Gaus(0., extra_resol_FastSim);
         pfCandidateTime = pfcandtimeFastSim + rndFastSim;
       }
	     if(pfCandidateTime!=-999) {
         // -- extra smearing to emulate different time resolution
         double extra_smearing = sqrt((targetTimeResol*targetTimeResol - 0.035*0.035));
         pfCandidateTime = pfCandidateTime + gRandom3->Gaus(0,extra_smearing);
         //dtsim = std::abs(pfCandidateTime - genPV->position().t()*1000000000.);
	     }
       dtsim = std::abs(pfCandidateTime - genPV->position().t()*1000000000.);
	     //else
	       //dtsim = 0.;

	     // *** A. Fill all candidates
	     // ** 0. Nominal
	     if(thisCandPassesDxy)
	       isolations.at(0) += pfCandidate.pt();
	     // ** 1. No Dxy
	     isolations.at(1) += pfCandidate.pt(); 
	     // *** B. Fill candidate times if within 3 sigma
	     //if (dtsim < 3.*targetTimeResol){
       //if ((keepTrack && pfCandidateTime!=-999 && dtsim > 3.*targetTimeResol)){if (fabs(iMuon.eta()) < 1.5)  std::cout << "false: " << keepTrack << pfCandidateTime << dtsim << targetTimeResol << std::endl;}
       //else{
       if (!(keepTrack && pfCandidateTime!=-999 && dtsim > 3.*targetTimeResol)){
         if ( fabs(iMuon.eta()) < 1.5 ){
           h_pfCandidate_cutflow_->Fill("time", 1);
           if(iMuon.pt()>4 && iMuon.pt()<20)
             h_pfCandidate_cutflow_4_->Fill("time", 1);
           if(iMuon.pt()>20)
             h_pfCandidate_cutflow_20_->Fill("time", 1);
         }
	       // ** 2. Nominal + 3sigma cut on timing
	       if(thisCandPassesDxy)
	         isolations.at(2) += pfCandidate.pt(); // nominal + 3sigma timing
	       // ** 3. No Dxy + 3sigma cut on timing
	       isolations.at(3) += pfCandidate.pt();  // no Dxy + 3sigma timing
	     }
       else{
         if (fabs(iMuon.eta())<1.5){
           //std::cout << "false: " << keepTrack << pfCandidateTime << dtsim << targetTimeResol << std::endl;
         }
       }

           } // end timing
           else { // no check on timing
	     // ** 0. Nominal
	     if(thisCandPassesDxy)
	       isolations.at(0) += pfCandidate.pt();
	     // ** 1. No Dxy
	     isolations.at(1) += pfCandidate.pt(); 
	     
           }
         //}// end of if condition for within isoCone
     
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

bool MuonIsolationAnalyzer::isGoodMuon(const reco::Muon& iMuon, edm::Handle<std::vector<reco::GenParticle> >& genHandle, edm::Handle<std::vector<reco::GenJet> >& genJetHandle, const SimVertex *genPV, bool ForMatchUse) const
{
  if (!ForMatchUse)
    h_muon_cutflow_->Fill("All Muons", 1);
  
  // *** just check to make sure muon track available
  if (iMuon.track().isNull()) 
    return false;
  if (!ForMatchUse)
    h_muon_cutflow_->Fill("muon has track", 1);
  
  // *** 0. pT > 20 GeV ---> use 4GeV for Bs->mumu studies, but keep for now to reproduce TDR results
  if (iMuon.pt() < 4.)
    return false;
  h_muon_cutflow_->Fill("pT > 4", 1);
  //if (iMuon.pt() < 20.)
  //  return false;
  //h_muon_cutflow_->Fill("pT > 20", 1);

  // *** 2. Loose ID
  if ( !(muon::isLooseMuon(iMuon)) )
    return false;
  if (!ForMatchUse)
    h_muon_cutflow_->Fill("Loose ID", 1);
  
  // *** 1. |eta| < 2.4
  //if ( fabs(iMuon.eta()) > 2.4 )
  //  return false;
  //h_muon_cutflow_->Fill("|eta| < 2.4", 1);
  if (isBmumu_ && fabs(iMuon.eta()) >1.4 )
    return false;
  if (isBmumu_ && !ForMatchUse)
    h_muon_cutflow_->Fill("eta<1.4", 1);

  reco::TrackRef muonTrack = iMuon.track();
  float dz_sim  = std::abs( muonTrack->vz() - genPV->position().z() ); 
  float dxy_sim = sqrt ( pow(muonTrack->vx() - genPV->position().x(),2) + pow(muonTrack->vy() - genPV->position().y(),2) ); 

  float dz_4D = std::abs( muonTrack->dz(vertex4D.position()));
  float dxy_4D = std::abs( muonTrack->dxy(vertex4D.position()));
  // some plots
  if ((!ForMatchUse) && (fabs(iMuon.eta()))<1.5){
    h_muon_dxy_BTL_->Fill( dxy_sim );
    h_muon_dz_BTL_->Fill( dz_sim );
  }

  if (!isBmumu_){
    // *** 4. z0 cut
    if ( dz_4D > dz_muonVertex )
      return false;
    if (!ForMatchUse)
      h_muon_cutflow_->Fill("dz(4D vtx,muon) < 0.5 cm", 1);
    
    // *** 3. d0 cut
    if ( dxy_4D > dxy_muonVertex )
      return false;
    if (!ForMatchUse)
      h_muon_cutflow_->Fill("dxy(4D vtx,muon) < 0.2 cm", 1);
  }
  // end

  // calculate some booleans
  bool recoMuonMatchedToPromptTruth = isPromptMuon(iMuon, genHandle);
  bool recoMuonMatchedToGenJet      = isMatchedToGenJet(iMuon, genJetHandle);
  bool recoMuonFromTau              = isFromTau(iMuon, genHandle);

  // *** 5A. accept only prompt muons if Z->mumu signal
  if (isZmumuSignal_) {

    if (!recoMuonMatchedToPromptTruth)
      return false;
    if (!ForMatchUse)
      h_muon_cutflow_->Fill("Signal Prompt Muon", 1);
  }
  // *** 5B. reject prompt muons if ttbar background
  else if ((!isZmumuSignal_) && (!isBmumu_) ) {

    // *** 5B-1. Muon is "good" non-prompt if !truthMatched && genJetMatched && !tauMatched
    if ( !recoMuonMatchedToPromptTruth && recoMuonMatchedToGenJet && !recoMuonFromTau) {
      if (!ForMatchUse)
        h_muon_cutflow_->Fill("Non-prompt Bkg Muon", 1);
    }
    else 
      return false;
  }

  // *** 5C. accept only muons matched with generated muons from Bs/d decay
  if (isBmumu_){
    auto BDaughterVector = BDaughters(genHandle);
    bool isMtahcedToDaughter = isMatchedToBMMDaughter(iMuon, BDaughterVector);
    if (isMtahcedToDaughter){
      std::cout << "muon matched :)" << std::endl;
      if (!ForMatchUse)
        h_muon_cutflow_->Fill("Muon from B", 1);
    }
    else
      return false;
  }

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
