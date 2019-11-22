#ifndef _UTILS_H
#define _UTILS_H

// system include
#include <memory>
#include <cstdlib>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

//#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"

//#include "Geometry/CommonTopologies/interface/Topology.h"
//#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
//#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
//#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
//#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "DataFormats/Common/interface/View.h"

using namespace std;
using namespace edm;
using namespace reco;
/* Martina version, 08-15-19
bool isPromptMuon(const reco::Muon &muon, const edm::View<reco::GenParticle>& genParticles);
bool isMatchedToGenJet(const reco::Muon &muon, const edm::View<reco::GenJet>& genJet);
bool isFromTau(const reco::Muon &muon, const edm::View<reco::GenParticle>& genParticles);
*/
// Ben version, 08-16-19
bool isPromptMuon(const reco::Muon &muon, edm::Handle<std::vector<reco::GenParticle> >& genHandle);
bool isMatchedToGenJet(const reco::Muon &muon, edm::Handle<std::vector<reco::GenJet> >& genJetHandle);
bool isFromTau(const reco::Muon &muon, edm::Handle<std::vector<reco::GenParticle> >& genHandle);

bool isMatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles);
bool isUnmatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles);

bool isMatchedToBMMDaughter(const reco::Muon& muon, std::vector<const reco::Candidate*> TrueMVector);
std::vector<const reco::Candidate*> BDaughters(edm::Handle<std::vector<reco::GenParticle> >& genHandle);
void FindMother(const reco::Candidate* Part, int rank);
bool PrintLeptonDaughters(edm::Handle<std::vector<reco::GenParticle> >& genHandle, const reco::Candidate * W_new, unsigned int id);
void FindtWbMuons(edm::Handle<std::vector<reco::GenParticle> >& genHandle);
const reco::Candidate * GetObjectBeforeDecay( const reco::Candidate * particle );
double GetdR(double eta1, double phi1, double eta2, double phi2);

#endif
