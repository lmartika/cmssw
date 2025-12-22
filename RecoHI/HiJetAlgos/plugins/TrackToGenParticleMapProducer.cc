// WIP 

// -*- C++ -*-
//
// Package:    TrackToGenParticleMapProducer
// Class:      TrackToGenParticleMapProducer
//

/**\class TrackToGenParticleMapProducer TrackToGenParticleMapProducer.cc
* @brief Produce a mapping between the tracks and the gen particles in the event
* given a jet collection and a list of final state generated particles
*
* The description of the run-time parameters can be found at fillDescriptions()
* The description of the products can be found at TrackToGenParticleMapProducer()
*/

// Original Author:  Lida Kalipoliti, LLR


// system include files
#include <memory>
#include <utility>
#include <vector>
#include <algorithm>
#include <map>

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

// Producer specific 
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "AnalysisDataFormats/TrackInfo/interface/TrackToGenParticleMap.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"


class TrackToGenParticleMapProducer : public edm::global::EDProducer<> {
public:
  explicit TrackToGenParticleMapProducer(const edm::ParameterSet& cfg);
  ~TrackToGenParticleMapProducer() override = default;

  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  // Definitions 
  typedef pat::Jet JetType; 
  typedef reco::GenJet GenJetType; 

  typedef reco::TrackToGenParticleMap MapType;

  typedef MapType::key_type TrackTypePtr; 
  typedef TrackTypePtr::value_type TrackType;
  
  typedef MapType::mapped_type GenTypePtr; 
  typedef GenTypePtr::value_type GenType;  
  
  // Configurables
  edm::EDGetTokenT<std::vector<GenType>> genParticlesToken_;
  edm::EDGetTokenT<std::vector<JetType>> inputJetsToken_;

  double maxDR2_;
  double minRelPt_;
  double maxRelPt_;

  // Private functions
  GenTypePtr findMatch(TrackTypePtr, edm::Handle<std::vector<GenType>>) const;
};

TrackToGenParticleMapProducer::TrackToGenParticleMapProducer(const edm::ParameterSet& cfg) {
  // Initialize configurables
  genParticlesToken_ = consumes<std::vector<GenType>>(cfg.getParameter<edm::InputTag>("genParticleSrc"));
  inputJetsToken_ = consumes<std::vector<JetType>>(cfg.getParameter<edm::InputTag>("jetSrc"));


  maxDR2_ = cfg.getUntrackedParameter<double>("maxDR2", 0.0004);
  minRelPt_ = cfg.getUntrackedParameter<double>("minRelPt", 0.8);
  maxRelPt_ = cfg.getUntrackedParameter<double>("maxRelPt", 1.2);
  
    /*  maxDR2_ = cfg.getUntrackedParameter<double>("maxDR2", 0.001);
  //  maxDR2_ = cfg.getUntrackedParameter<double>("maxDR2", 0.0004);
  minRelPt_ = cfg.getUntrackedParameter<double>("minRelPt", 0.6);
  maxRelPt_ = cfg.getUntrackedParameter<double>("maxRelPt", 1.4);  */

  produces<MapType>("trackToGenParticleMap");
  produces<MapType>("genConstitToGenParticleMap");
}


// ------------ method called to produce the data  ------------
void TrackToGenParticleMapProducer::produce(edm::StreamID, edm::Event &evt, const edm::EventSetup &setup) const {
  //  std::cout << "In TrackToGenParticleMapProducer produce" << std::endl;

  edm::Handle<std::vector<JetType>> inputJets;
  evt.getByToken(inputJetsToken_, inputJets);

  // Grab the gen particles
  edm::Handle<std::vector<GenType>> genParticles;
  evt.getByToken(genParticlesToken_, genParticles);

  // Create output collections
  std::unique_ptr<MapType> trackToGenParticleMap = std::make_unique<MapType>();
  std::unique_ptr<MapType> genConstitToGenParticleMap = std::make_unique<MapType>();

  /*  for (size_t igen = 0; igen < (*genParticles).size(); igen++) {
    edm::Ref<std::vector<GenType>> genParticleRef_(genParticles, int(igen));
    GenTypePtr genParticlePtr = edm::refToPtr(genParticleRef_);
   
    if (genParticlePtr->status() > 1) std::cout << "printing gen B dec products. charge, pt, status: " << genParticlePtr->charge() << " " << genParticlePtr->pt() << " " << genParticlePtr->status() << std::endl;
  } */

  
  for (const JetType &jet : *inputJets) {
    //    std::cout << " ----- New jet in map producer, with pt, eta, phi " << jet.pt() << " " << jet.eta() << " " << jet.phi() << std::endl;

    std::map<GenTypePtr,int> counter;
    std::map<GenTypePtr, std::pair<double, TrackTypePtr>> tracker;
      
    // Go over the charged jet constituents (aka tracks)
    for (const TrackTypePtr &constitPtr : jet.getJetConstituents()) {
      if (constitPtr->charge() == 0) continue;

      //      std::cout << "constitPtr in reco jet consts " << constitPtr->charge() << " " << constitPtr->pt() << " " <<  constitPtr->eta() << std::endl;
      // Look for match in charged gen particles
      GenTypePtr matchGenParticle = findMatch(constitPtr, genParticles);

      // To check double-matches; fill a map of gen particles so that gen ptcl points to the closest track
      if (matchGenParticle.isNonnull()) {
	double dr = reco::deltaR2(*constitPtr, *matchGenParticle);
	
	if (tracker.find(matchGenParticle) != tracker.end()) {
	  if ( dr < tracker[matchGenParticle].first ) tracker[matchGenParticle] = std::pair<double,TrackTypePtr>(dr,constitPtr);
	}
	else {
	  tracker[matchGenParticle] = 	  std::pair<double,TrackTypePtr>(dr,constitPtr);
	}

      }

      // If want unique matches (smallest dr), instead of this, go through the tracker-list
      ///////////
      if (matchGenParticle.isNonnull()) 
	      trackToGenParticleMap->insert(std::pair<TrackTypePtr, GenTypePtr>(constitPtr, matchGenParticle));
      
    } // end reco jet constituents loop

    // Then fill the map per-jet, following for the unique matches
    /*    for (auto it = tracker.begin(); it !=tracker.end(); it++) { 
      //      std::cout << "Testing: " << it->first->pt() << " " << it->second.second->pt() << std::endl;
      trackToGenParticleMap->insert(std::pair<TrackTypePtr, GenTypePtr>(it->second.second, it->first));
      } */
 
    
    // grab the gen jet
    const GenJetType *genJet = jet.genJet();
    if (!genJet) continue;
    //    std::cout << "Found a gen jet in map producer" << std::endl;

    // Go over the charged gen jet constituents
    //    std::cout << "Gen jet: " << std::endl;
    //    std::cout << " ----- Corresponding genjet in map producer, with pt, eta, phi " << genJet->pt() << " " << genJet->eta() << " " << genJet->phi() << std::endl;
    for (const TrackTypePtr &constitPtr : genJet->getJetConstituents()) {
      if (constitPtr->charge() == 0) continue;

      //      std::cout << "constitPtr in gen consts " << constitPtr->charge() << " " << constitPtr->pt() << " " <<  constitPtr->eta() << std::endl;
            GenTypePtr matchGenParticle = findMatch(constitPtr, genParticles);
            if (matchGenParticle.isNonnull()) 
	      genConstitToGenParticleMap->insert(std::pair<TrackTypePtr, GenTypePtr>(constitPtr, matchGenParticle));
    } // end gen jet constituents loop



    // [TEST]: Go over the SV tracks
    /*    std::string svTagInfoLabel_ = "pfInclusiveSecondaryVertexFinder";
    const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());
    int nsv = svTagInfo->nVertices(); 
    for (int isv = 0; isv < nsv; isv++) {
      const std::vector<reco::CandidatePtr> svTracks = svTagInfo->vertexTracks(isv);
      for (auto svTrkPtr : svTracks) {
        if (svTrkPtr->charge() == 0) continue;
        // Look if the particle is already in the map
        if (trackToGenParticleMap->find(svTrkPtr) != trackToGenParticleMap->end()) {
          continue;
        }

        // If not, look for match in charged gen particles 
        GenTypePtr matchGenParticle = findMatch(svTrkPtr, genParticles);
        if (matchGenParticle.isNonnull()) { 
	  //          trackToGenParticleMap->insert(std::pair<TrackTypePtr, GenTypePtr>(svTrkPtr, matchGenParticle));
	  //	  std::cout << "Added a SV track" << std::endl;
	}
      } // end sv trk loop
      } // end sv loop  */
  } // end jet loop

  evt.put(std::move(trackToGenParticleMap), "trackToGenParticleMap");
  evt.put(std::move(genConstitToGenParticleMap), "genConstitToGenParticleMap");
}

  TrackToGenParticleMapProducer::GenTypePtr 
  TrackToGenParticleMapProducer::findMatch(TrackTypePtr cand, 
                                           edm::Handle<std::vector<GenType>> genParticles) const {
  double minDR2 = std::numeric_limits<double>::max();
  GenTypePtr matchGenParticle;
  std::map<int,int> counter;
  //  int matchidx = 0;
  
  for (size_t igen = 0; igen < (*genParticles).size(); igen++) {
    edm::Ref<std::vector<GenType>> genParticleRef_(genParticles, int(igen));
    GenTypePtr genParticlePtr = edm::refToPtr(genParticleRef_);

    if (genParticlePtr->charge() == 0) continue;

    double DR2  = reco::deltaR2(*cand, *genParticlePtr);
    if (DR2 > maxDR2_) continue;

    //    std::cout << "matching  to gen ptcl with charge, pt, eta, phi, massm, status, DR: " << genParticlePtr->charge() << " " << genParticlePtr->pt() << " " << genParticlePtr->eta() << " " << genParticlePtr->phi() << " " << genParticlePtr->mass() << " " << genParticlePtr->status() << " " << DR2 << " igen " << igen << std::endl;

  double relPt = cand->pt() / genParticlePtr->pt();
  if (relPt < minRelPt_ || relPt > maxRelPt_) continue;

    if (DR2 < minDR2) {
      minDR2 = DR2;
      matchGenParticle = genParticlePtr;
      //      matchidx = igen;
    }
  } // end gen particle loop

  /*  if (matchGenParticle.isNonnull())   {
    
    std::cout << "matched to gen ptcl with charge, pt, eta, phi, massm, status, DR: " << matchGenParticle->charge() << " " << matchGenParticle->pt() << " " << matchGenParticle->eta() << " " << matchGenParticle->phi() << " " << matchGenParticle->mass() << " " << matchGenParticle->status() << " " << minDR2 << " matchidx " << std::endl;
    std::cout << "cand was: charge, pt, eta, phi massm: " << cand->charge() << " " << cand->pt() << " " << cand->eta() << " " << cand->phi() << " " << cand->mass() << std::endl;

  }
  else {
    std::cout << "No match gen ptcl" << std::endl;
    //    std::cout << "cand was: charge, pt, eta, phi massm: " << cand->charge() << " " << cand->pt() << " " << cand->eta() << " " << cand->phi() << " " << cand->mass() << std::endl;
    } */
  
  return matchGenParticle;
}

void TrackToGenParticleMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("jet constituent to gen particle map");
  desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("genParticleSrc", edm::InputTag("packedGenParticles"));

  /*   desc.addUntracked<double>("maxDR2", 0.001);
  desc.addUntracked<double>("minRelPt", 0.6);
  desc.addUntracked<double>("maxRelPt", 1.4);
  */

  desc.addUntracked<double>("maxDR2", 0.0004);
  desc.addUntracked<double>("minRelPt", 0.8);
  desc.addUntracked<double>("maxRelPt", 1.2);

  descriptions.add("TrackToGenParticleMapProducer", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(TrackToGenParticleMapProducer);
