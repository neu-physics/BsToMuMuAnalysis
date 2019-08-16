// copied from M. Malberti: https://github.com/martinamalberti/PrecisionTiming/blob/isolation_studies/PTAnalysis/plugins/Utils.cc
// user include files                                                                                                                                                                                                                  
#include "BsToMuMuAnalysis/MuonIsolationAnalyzer/interface/Utils.h"

// --- matching to gen muon ----------------------------------------------------------------
bool isPromptMuon(const reco::Muon& muon, edm::Handle<std::vector<reco::GenParticle> >& genHandle)
{
  bool isPrompt = false;
  float drmin = 9999.;

  auto genParticles = *genHandle.product();
  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    auto genParticle = genParticles.at(ip);
    //const reco::GenParticle& genp = genParticles[ip];
    if ( std::abs(genParticle.pdgId()) != 13) continue;
    if (genParticle.status() != 1 || !genParticle.isLastCopy() ) continue; // -- from Simone                                                                                                                                                          
    if ( !genParticle.isPromptFinalState() ) continue;
    double dr = deltaR(muon,genParticle);
    if (dr > 0.2){
      continue;
    }
    //else{
    //  isPrompt=true;
    //  break;
    //}
    if (dr < drmin && dr < 0.2){
      drmin = dr;
      isPrompt=true;   
    }
  }

  return isPrompt;
}
//-------------------------------------------------------------------------------------------



// --- matching to gen jet ------------------------------------------------------------------
bool isMatchedToGenJet(const reco::Muon& muon, edm::Handle<std::vector<reco::GenJet> >& genJetHandle)
{
  bool isMatched = false;

  auto genJets = *genJetHandle.product();
  for( unsigned int iGenJet = 0; iGenJet < genJetHandle->size(); ++iGenJet ) {
    auto genJet = genJets.at(iGenJet);

    if ( genJet.pt() < 15.0  || genJet.hadEnergy()/genJet.energy() < 0.3) continue;
    double dr = deltaR(muon,genJet);
    if (dr > 0.3){
      continue;
    }
    else{
      isMatched=true;
      break;
    }
  }

  return isMatched;
}
//-------------------------------------------------------------------------------------------



// --- matching to muons from tau decays ----------------------------------------------------
bool isFromTau(const reco::Muon& muon, edm::Handle<std::vector<reco::GenParticle> >& genHandle)
{

  bool fromTau = false;

  auto genParticles = *genHandle.product();
  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    auto genParticle = genParticles.at(ip);

    if ( std::abs(genParticle.pdgId()) != 13) continue;
    if ( !genParticle.isDirectPromptTauDecayProductFinalState() ) continue;
    double dr = deltaR(muon,genParticle);
    if (dr > 0.2){
      continue;
    }
    else{
      fromTau=true;
      break;
    }
  }

  return fromTau;
}
//-------------------------------------------------------------------------------------------




//--- matching pfcands to gen level particles (from PV) ------------------------------------------------------- 
bool isMatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles){

  bool isMatched = false;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    // -- stable particle
    if (genp.status() !=1 ) continue;
    // -- charged or neutral particle
    if (genp.charge()!=pfcand.charge()) continue;
    // -- pt matching 
    float dpt = std::abs((pfcand.pt()-genp.pt())/genp.pt());
    if ( dpt > 0.1) continue;
    // -- deltaR matching
    double dr = deltaR(pfcand,genp);
    if (dr > 0.05){
      continue;
    }
    else{
      isMatched=true;
      break;
    }
  }

  return isMatched;

}
//-------------------------------------------------------------------------------------------



// -- not matching to gen particles (looser deltaR) -------------------------------------------------------------
bool isUnmatchedToGenParticle(const reco::PFCandidate &pfcand, const edm::View<reco::GenParticle>& genParticles){

  bool isUnmatched = true;

  for(unsigned int ip=0; ip < genParticles.size(); ip++ ){
    const reco::GenParticle& genp = genParticles[ip];
    // -- stable particle
    if (genp.status() !=1 ) continue;
    // -- charged or neutral particle
    if (genp.charge()!=pfcand.charge()) continue;
    // -- not a muon
    if (genp.pdgId() == 13 ) continue;
    // -- pt matching
    //float dpt = std::abs((pfcand.pt()-genp.pt())/genp.pt());
    //if ( dpt > 0.2) continue; 
    // -- deltaR matching
    double dr = deltaR(pfcand,genp);
    if (dr > 0.15){
      continue;
    }
    else{
      isUnmatched=false;
      break;
    }
  }

  return isUnmatched;

}
//-------------------------------------------------------------------------------------------




