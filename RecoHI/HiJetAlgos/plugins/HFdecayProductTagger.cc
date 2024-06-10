// -*- C++ -*-
//
// Package:    HFdecayProductTagger
// Class:      HFdecayProductTagger
//

/**\class HFdecayProductTagger HFdecayProductTagger.cc
* @brief Given a list genParticles, produce a genParticle collection of final state particles,  
* where the decay products of heavy flavour hadrons have a hijacked status code
*
* The description of the run-time parameters can be found at fillDescriptions()
*
* The description of the products can be found at HFdecayProductTagger()
*/

// Original Author:  Matthew Nguyen, LLR


// system include files
#include <memory>
#include <utility>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CommonTools/CandUtils/interface/pdgIdUtils.h"
#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"
#include "TLorentzVector.h"



class HFdecayProductTagger : public edm::global::EDProducer<> {
public:
  explicit HFdecayProductTagger(const edm::ParameterSet& cfg);
  ~HFdecayProductTagger() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  // ----------member data ---------------------------
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

  bool isFinalB(const reco::Candidate &particle) const;  
  bool isFromB(const reco::Candidate &particle) const;  

  bool isFinalC(const reco::Candidate &particle) const;  
  bool isFromC(const reco::Candidate &particle) const;  

  reco::GenParticleCollection addDaughters(const reco::Candidate &particle, reco::GenParticleCollection daughterCollection, int hfCode) const;
  //void visible( reco::Candidate::PolarLorentzVector &v, const reco::Candidate &particle, int doCharge) const;

  bool tagBorC_; // true to tag B, false to tag C

};

HFdecayProductTagger::HFdecayProductTagger(const edm::ParameterSet& cfg)
  : genParticlesToken_(consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("genParticles"))),
	  tagBorC_(cfg.getParameter<bool>("tagBorC")){

  produces<std::vector<pat::PackedGenParticle>>();
  
	}


// ------------ method called to produce the data  ------------
void HFdecayProductTagger::produce(edm::StreamID, edm::Event &evt, const edm::EventSetup &setup) const {
  //std::cout << "Event being processed..." << std::endl;
  // std::cout << "HFdecayProductTagger produce" << std::endl;

  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByToken(genParticlesToken_, genParticles);
  // std::cout << "nb of input particles, HFdecayProductTagger: " << (*genParticles).size() << std::endl;

  // Create output collection
  auto outputCollection = std::make_unique<std::vector<pat::PackedGenParticle>>();
  //std::cout << "type of outputCollection: " << typeid(outputCollection).name() << std::endl;

  // Keep track of HF's
  int hfCode = 100;
  if (tagBorC_) {
    // Tag fromB's
    // std::cout << "Tagging fromB's as usual" << std::endl;
    for (const reco::GenParticle& genPart : *genParticles) {
      if(isFinalB(genPart)){
        // std::cout << "Found a B, adding its daughters" << std::endl;
        // Do NOT add the B
        //outputCollection->push_back(genPart);

        // Add the daughters to the output collection
        reco::GenParticleCollection daughterCollection = {};
        daughterCollection = addDaughters(genPart, daughterCollection, hfCode);
        // std::cout << "Found " << daughterCollection.size() << " daugthers"<< std::endl;
        // std::cout << "Trying to convert to packed " << std::endl;

        // auto testDaughter = pat::PackedGenParticle(daughterCollection[0], reco::GenParticleRef());
        // std::cout << "... Done " << std::endl;

        // testDaughter.unpack();
        // std::cout << "daughter pt " << daughterCollection[0].pt()
        //           << " packed daughter pt " << testDaughter.pt()
        //           << " daughter eta " << daughterCollection[0].eta()
        //           << " packed daughter eta " << testDaughter.eta()
        //           << " daughter phi " << daughterCollection[0].phi()
        //           << " packed daughter phi " << testDaughter.phi()
        //           << std::endl;


        for (const reco::GenParticle daughter : daughterCollection){
          // std::cout << "daughter charge " << daughter.charge() << " and status " << daughter.status() << std::endl;
          // auto packedDaughter = pat::PackedGenParticle(daughter, reco::GenParticleRef());
          pat::PackedGenParticle packedDaughter(daughter, reco::GenParticleRef());
          // std::cout << "packed daughter charge " << packedDaughter.charge() << " and status " << packedDaughter.status() << std::endl;
          outputCollection->push_back(packedDaughter);
        }
        // (*outputCollection).insert((*outputCollection).end(), daughterCollection.begin(), daughterCollection.end());
        hfCode += 100;
      }
      if (genPart.status() != 1) continue;
      // Skip b daughters
      if (isFromB(genPart)) {
        //std::cout << "Found from B!" << std::endl;
        continue;	
      }
      // Add the rest of the particles			 
      pat::PackedGenParticle packedGenParticle(genPart, reco::GenParticleRef());
      outputCollection->push_back(packedGenParticle);
    } // end loop over particles
  } else {
	  // Tag fromC's
	  //std::cout << "Tagging formC's" << std::endl;
    for (const reco::GenParticle& genPart : *genParticles) {
      if(isFinalC(genPart)){
        // Do NOT add the C
        //outputCollection->push_back(genPart);

        // Add the daughters to the output collection
        reco::GenParticleCollection daughterCollection = {};
        daughterCollection = addDaughters(genPart, daughterCollection, hfCode);
        for (const reco::GenParticle daughter : daughterCollection){
          pat::PackedGenParticle packedDaughter(daughter, reco::GenParticleRef());
          outputCollection->push_back(packedDaughter);
        }
        // (*outputCollection).insert((*outputCollection).end(), daughterCollection.begin(), daughterCollection.end());
        hfCode += 100;
      }
      // Skip non-stable particles
      if (genPart.status() != 1) continue;
      // Skip C daughters
      if (isFromC(genPart)) {
      //std::cout << "Found from C!" << std::endl;
      continue;	
      }
      // Add the rest of the particles			 
      pat::PackedGenParticle packedGenParticle(genPart, reco::GenParticleRef());
      outputCollection->push_back(packedGenParticle);
    } // end loop over particles
  } // end if tagBorC
  // std::cout << "nb of output particles in HFdecayProductTagger: " << (*outputCollection).size() << std::endl;
  // count charged particles in collection
  // int counts = 0;
  // for (auto particle : *outputCollection) {
  //  if (particle.charge() != 0) counts++;
  // }
  // std::cout << "out oh which: " << counts << " charged" << std::endl;
  evt.put(std::move(outputCollection));
}

bool HFdecayProductTagger::isFinalB(const reco::Candidate &particle) const 
{
  if (!CandMCTagUtils::hasBottom(particle)) return false;

  // If it is a b hadron, check if any of the daughters is also a b hadron
  unsigned int npart = particle.numberOfDaughters();
  for (size_t i = 0; i < npart; ++i) {
    // If any of them do, then it is not a final B
    if (CandMCTagUtils::hasBottom(*particle.daughter(i))) return false;
  }
  return true;
}


bool HFdecayProductTagger::isFromB( const reco::Candidate &particle) const {
  bool fromB = false;

  // Check if any of the mothers are final B
  unsigned int npart = particle.numberOfMothers();
  for (size_t i = 0; i < npart; ++i) {
    const reco::Candidate &mom = *particle.mother(i);
    if (CandMCTagUtils::hasBottom(mom)){
      if(isFinalB(mom)){
        fromB = true;
        break;
      }
      else{
        fromB = false;
        break;	
      }
    }
    else fromB = isFromB(mom);
  }
  return fromB;
}

bool HFdecayProductTagger::isFinalC(const reco::Candidate &particle) const {
  if (!CandMCTagUtils::hasCharm(particle) ) return false;

  // check if any of the daughters is also a c hadron
    unsigned int npart=particle.numberOfDaughters();
    
    for (size_t i = 0; i < npart; ++i) {
      if (CandMCTagUtils::hasCharm(*particle.daughter(i))) return false;
    }
    
    return true;
}


bool HFdecayProductTagger::isFromC( const reco::Candidate &particle) const{

  bool fromC = false;

  unsigned int npart=particle.numberOfMothers();
  for (size_t i = 0; i < npart; ++i) {
    const reco::Candidate &mom = *particle.mother(i);
    if (CandMCTagUtils::hasCharm(mom)){
      if(isFinalC(mom)){
		fromC = true;
		break;
      }
      else{
		fromC = false;
		break;	
      }
    }
    else fromC = isFromC(mom);
  }
  return fromC;
}


/*
void HFdecayProductTagger::visible
(
 reco::Candidate::PolarLorentzVector &v, const reco::Candidate &particle, int doCharge) const
{


  unsigned int npart=particle.numberOfDaughters();
  
  for (size_t i = 0; i < npart; ++i) {
    if(particle.daughter(i)->status()==1){
      int charge = particle.daughter(i)->charge();
      if(doCharge == 1 && charge ==0) continue;
      int pdgid = abs(particle.daughter(i)->pdgId());
      if(doCharge == 0 && charge !=0) continue;
      if(pdgid == 12 || pdgid == 14 || pdgid == 16) continue;     

      reco::Candidate::PolarLorentzVector vTemp(0.,0.,0.,0.);
      vTemp.SetPt(particle.daughter(i)->pt());
      vTemp.SetEta(particle.daughter(i)->eta());
      vTemp.SetPhi(particle.daughter(i)->phi());
      vTemp.SetM(particle.daughter(i)->mass());
      v+=vTemp;
      //std::cout<<" adding a particle:  pt= "<<particle.daughter(i)->pt()<<" pdg id "<<particle.daughter(i)->pdgId()<<std::endl;
    }
    else{
      visible(v,*particle.daughter(i),doCharge);
    }
  }
  
}
*/
reco::GenParticleCollection HFdecayProductTagger::addDaughters(const reco::Candidate &particle, 
                                                               reco::GenParticleCollection collection, 
                                                               int hfCode) const {

  reco::GenParticleCollection daughterCollection = collection;

  // Get daughters
  unsigned int ndaught = particle.numberOfDaughters();

  for (size_t i = 0; i < ndaught; i++) {
    const reco::Candidate& daughter = *particle.daughter(i);
    // Add only stable daughters
    if (daughter.status() == 1) {
      //std::cout << "Found a daughter!" << std::endl;
      daughterCollection.push_back(reco::GenParticle(daughter.charge(), daughter.p4(), daughter.vertex(), daughter.pdgId(), hfCode, true));
    }
	  // Get daughters of daughter
    unsigned int ndaughtdaught = daughter.numberOfDaughters();
    if (ndaughtdaught > 0) {
      daughterCollection = addDaughters(daughter, daughterCollection, hfCode);
    }
  }
  //std::cout << "Daughters found : " << daughterCollection.size() << std::endl;
  return daughterCollection;
}
			    

void HFdecayProductTagger::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
  desc.add<bool>("tagBorC", true);
  descriptions.add("HFdecayProductTagger", desc);
}


DEFINE_FWK_MODULE(HFdecayProductTagger);


