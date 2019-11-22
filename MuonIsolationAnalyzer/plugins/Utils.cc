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

bool isMatchedToBMMDaughter(const reco::Muon& muon, std::vector<const reco::Candidate*> TrueMVector)
{
  bool Matched = false;
  double dR_min = 999;
  for (auto* TureMu:TrueMVector){
    double dR = deltaR(muon, *TureMu);
    if (dR>0.2)   continue;
    if (dR<0.2 && dR<dR_min){
      dR_min = dR;
      Matched = true;
    }
  }
  return Matched;
}

std::vector<const reco::Candidate*> BDaughters(edm::Handle<std::vector<reco::GenParticle> >& genHandle)
{
  std::vector<const reco::Candidate*> TrueMuons = {};
  for(size_t iGenParticle=0; iGenParticle<genHandle->size();iGenParticle++){
     if(fabs((*genHandle)[iGenParticle].pdgId())>500 && fabs((*genHandle)[iGenParticle].pdgId())<600){
       const reco::Candidate * iCand = &(*genHandle)[iGenParticle];
       //reco::Candidate * Muon1;
       //reco::Candidate * Muon2;
       std::vector<const reco::Candidate*> muons;
       if(iCand->numberOfDaughters()!=2)  continue;
       bool TwoMDaughter = true;
       for(size_t iDaughter=0; iDaughter<iCand->numberOfDaughters(); ++iDaughter){
         if( fabs(iCand->daughter(iDaughter)->pdgId())!=13 ){
           TwoMDaughter = false;
           break;
         }
         else{
           //reco::Candidate muon = *(iCand->daughter(iDaughter));
           muons.push_back((iCand->daughter(iDaughter)));
           //TrueMuons.push_back(iCand->daughter(iDaughter));
         }
       }
       if ((TwoMDaughter)&&(muons.size()==2)&&(muons[0]->charge()==(-muons[1]->charge()))){
         for (auto* mu:muons){
           TrueMuons.push_back(mu);
         }
       }
     }
  }
  return TrueMuons;
}

void FindMother(const reco::Candidate* Part, int rank)
{
  for(size_t iMother=0;iMother<Part->numberOfMothers();++iMother){
    std::cout << "rank: " << rank << " particle: "<< Part->pdgId() <<" mother: " << Part->mother(iMother)->pdgId() << " status: " << Part->mother(iMother)->status() << "  # mothers: " << Part->numberOfMothers() << std::endl;
    if((Part->mother(iMother)->numberOfMothers())){
      FindMother(Part->mother(iMother),rank+1);
    }
  }

}

bool PrintLeptonDaughters(edm::Handle<std::vector<reco::GenParticle> >& genHandle, const reco::Candidate * W_new, unsigned int id)
{
  //std::cout << "Printing lepton Daughters " << id << std::endl;
  bool FoundDaughter = false;
  for(size_t iD=0; iD<W_new->numberOfDaughters(); ++iD){
    if(fabs(W_new->daughter(iD)->pdgId())==id){
      reco::GenParticle mu;
      double dr_min = 999;
      //std::cout << "----Muon from the W of t->Wb decay: " << std::endl;
      for(size_t iG=0; iG<genHandle->size();++iG){
        if (fabs((*genHandle)[iG].pdgId())!=id)  continue;
        double dr = GetdR((W_new->daughter(iD))->eta(),(W_new->daughter(iD))->phi(),(*genHandle)[iG].eta(),(*genHandle)[iG].phi());
        //std::cout << "iG: " << iG << " " << (W_new->daughter(iD))->eta() << " " << (W_new->daughter(iD))->phi() << "  the gen Mu: " << (*genHandle)[iG].eta() << "  " << (*genHandle)[iG].phi() <<  " dR: " << dr << std::endl;
        if(dr<0.1 && dr<dr_min){
          dr_min = dr;
          mu = (*genHandle)[iG];
        }
      }
      if(dr_min!=999){
        if(!FoundDaughter) FoundDaughter = true;
        std::cout <<"The cand lepton: " << "  Prompt: " << mu.isPromptFinalState() << " pT: " << mu.pt() << " eta: " << mu.eta() << " phi: " << mu.phi() << std::endl;
        for (size_t lD=0; lD<W_new->daughter(iD)->numberOfDaughters(); ++lD){
          std::cout << "  Lepton daughters: " << W_new->daughter(iD)->daughter(lD)->pdgId() << std::endl;
          if (fabs(W_new->daughter(iD)->daughter(lD)->pdgId())==id){
            const reco::Candidate * lep_daughter = GetObjectBeforeDecay(W_new->daughter(iD)->daughter(lD));
            for(size_t iL=0; iL<genHandle->size();++iL){
              if (fabs((*genHandle)[iL].pdgId())!=id)  continue;
              double dr = GetdR(lep_daughter->eta(),lep_daughter->phi(),(*genHandle)[iL].eta(),(*genHandle)[iL].phi());
              if(dr==0){
                std::cout <<"The cand lepton daughter: " << "  Prompt: " << mu.isPromptFinalState() << std::endl;
              }
            }
          }
        }
        //FindMother(W_new->daughter(iD),0);
      }
    }
  }
  return FoundDaughter;
}

void FindtWbMuons(edm::Handle<std::vector<reco::GenParticle> >& genHandle)
{
  const reco::Candidate * W;
  const unsigned int id = 13;
  //const reco::GenParticle * W; 
  for(size_t iGen=0; iGen<genHandle->size();++iGen){
    if(fabs((*genHandle)[iGen].pdgId())==6){
      if((*genHandle)[iGen].numberOfDaughters()!=2)  continue;
      if ( fabs((*genHandle)[iGen].daughter(0)->pdgId())==24 && fabs((*genHandle)[iGen].daughter(1)->pdgId())==5 ){
        W = (*genHandle)[iGen].daughter(0);
        //std::cout << "----Found t->Wb decay" << std::endl;
      } 
      else if (fabs((*genHandle)[iGen].daughter(0)->pdgId())==5 && fabs((*genHandle)[iGen].daughter(1)->pdgId())==24){
        W = (*genHandle)[iGen].daughter(1);
        //std::cout << "----Found t->Wb decay" << std::endl;
      }
      else
        continue;
      //if(!PrintLeptonDaughters(genHandle, W, id)){
      if(1){
        //std::cout << "W->W: " << std::endl;
        const reco::Candidate * W_new = GetObjectBeforeDecay(W);
        PrintLeptonDaughters(genHandle, W_new, id);
      }
    }
  }
}

const reco::Candidate * GetObjectBeforeDecay( const reco::Candidate * particle ){

  for ( unsigned int i = 0 ; i <  particle -> numberOfDaughters(); i++ ){
    if( particle -> daughter( i ) -> pdgId()  ==  particle -> pdgId() ){

      return GetObjectBeforeDecay( particle -> daughter (i) );

    } // end if
  } // end for

  return particle ;

}


double GetdR(double eta1, double phi1, double eta2, double phi2)
{
  double d_eta = eta1-eta2;
  double d_phi = std::abs(phi1-phi2);
  //std::cout << "calculate dR: " << eta1 << "  " << phi1 << "  " << eta2 << "  " << phi2 << " deta: " << d_eta << " dphi: " << d_phi << std::endl;
  if (d_phi>M_PI){
    d_phi -= 2.0*M_PI;
  }
  return std::sqrt(d_eta*d_eta + d_phi*d_phi);
}

//-------------------------------------------------------------------------------------------




