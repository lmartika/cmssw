
/* 
 
 Aggregate B hadrons in PF Collection.
  
*/

#include <memory>
#include <vector>
#include <cmath>
#include <map>
#include <random>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackToGenParticleMap.h"
#include "CommonTools/MVAUtils/interface/TMVAEvaluator.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"


class aggregatedPFCollection : public edm::global::EDProducer<> {
public:
  explicit aggregatedPFCollection(const edm::ParameterSet&);
  ~aggregatedPFCollection() override = default;

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  
  
  // ------------- member data ----------------------------
  edm::EDGetTokenT<pat::JetCollection> jetSrc_;
  edm::EDGetTokenT<pat::JetCollection> matchTag_;
  
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> constitSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>>  packedConstitSrc_;
  edm::EDGetTokenT<reco::TrackToGenParticleMap> candToGenParticleMapToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<int> CentralityBinTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> primaryVerticesToken_;
  edm::Handle<std::vector<reco::Vertex>> primaryVertices;  

  std::unique_ptr<TMVAEvaluator> tmvaTagger;

  bool isMC_;

  bool writeConstits_;
  bool doGenJets_;
  bool chargedOnly_;
  bool domatch_;
  
  double rParam_;
  double ptCut_; // For tracks in aggregation
  double trkInefRate_; // 0 by default
  double jetPtCut_;// skip low pT jets
  double jetEtaCut_ = 3.0;// skip forward jets

  bool aggregateHF_;
  bool withTruthInfo_;
  bool withCuts_;
  bool withTMVA_;
  edm::FileInPath tmva_path_;
  std::vector<std::string> tmva_variable_names_;
  std::vector<std::string> tmva_spectator_names_;
  std::string ipTagInfoLabel_;
  std::string svTagInfoLabel_;
};

aggregatedPFCollection::aggregatedPFCollection(const edm::ParameterSet& iConfig) {
  // Get configuration parameters
    isMC_ = iConfig.getParameter<bool>("isMC");
    chargedOnly_ = iConfig.getParameter<bool>("chargedOnly");
    aggregateHF_ = iConfig.getParameter<bool>("aggregateHF"); 
    doGenJets_ = iConfig.getParameter<bool>("doGenJets");
    domatch_ = iConfig.getParameter<bool>("domatch");
  
    ptCut_ = iConfig.getParameter<double>("ptCut");
    jetPtCut_ = iConfig.getParameter<double>("jetPtCut");
    trkInefRate_ = iConfig.getParameter<double>("trkInefRate");

    if (aggregateHF_) {
        withTruthInfo_ = iConfig.getParameter<bool>("aggregateWithTruthInfo");
        withCuts_ = iConfig.getParameter<bool>("aggregateWithCuts");
        withTMVA_ = iConfig.getParameter<bool>("aggregateWithTMVA");

        if (withTMVA_) {
            tmva_path_ = iConfig.getParameter<edm::FileInPath>("tmva_path");
            tmva_variable_names_ = iConfig.getParameter<std::vector<std::string>>("tmva_variables");
            tmva_spectator_names_ = iConfig.getParameter<std::vector<std::string>>("tmva_spectators");
        }
    } 


    // Get labels
    ipTagInfoLabel_ = iConfig.getParameter<std::string>("ipTagInfoLabel");
    svTagInfoLabel_ = iConfig.getParameter<std::string>("svTagInfoLabel");
  
    // Get tokens
    jetSrc_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"));
    matchTag_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("matchTag"));
    constitSrc_ = consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"));
    packedConstitSrc_ = consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"));
    CentralityBinTag_ = consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc"));    
    primaryVerticesToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices")));

    if (aggregateHF_ && isMC_) candToGenParticleMapToken_ = consumes<reco::TrackToGenParticleMap>(iConfig.getParameter<edm::InputTag>("candToGenParticleMap"));

    // Initialize objects 
    if (aggregateHF_ && withTMVA_) {
        tmvaTagger = std::make_unique<TMVAEvaluator>();
        tmvaTagger->initialize("Color:Silent:Error",
                               "BDTG",
                                tmva_path_.fullPath(),
                                tmva_variable_names_,
                                tmva_spectator_names_,
                                false,
                                false);
    }

    produces<reco::PFCandidateCollection>();
//    produces<reco::PFCandidateCollection>("newPFCollectionHFgen");
//    produces<reco::PFCandidateCollection>("newPFCollectionHFreco");

}

void aggregatedPFCollection::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  //     std::cout << "In aggregatedPFCollection::produce" << std::endl;

    auto newPFCandCollection = std::make_unique<reco::PFCandidateCollection>();

    //    auto newPFCandCollectionHF = std::make_unique<reco::PFCandidateCollection>();

    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetSrc_, jets);    

    // Matched jets for taginfos and TMVA aggregation 
    edm::Handle<pat::JetCollection> matchedjets;
    iEvent.getByToken(matchTag_, matchedjets);
    // -- For aggregation -- //
    edm::Handle<reco::TrackToGenParticleMap> candToGenParticleMap;

    if (aggregateHF_ && isMC_) {
        iEvent.getByToken(candToGenParticleMapToken_, candToGenParticleMap);
    } 

    //centrality
    edm::Handle<int> cbin_;
    iEvent.getByToken(CentralityBinTag_, cbin_);
    int hiBin = *cbin_;    
  
    edm::Handle<std::vector<reco::Vertex>> primaryVertices;
    iEvent.getByToken(primaryVerticesToken_, primaryVertices);

    std::vector<fastjet::PseudoJet> jetConstituents = {};

    // std::cout << "Matched jets: " << matchedjets->size() << std::endl; // debug

    for(unsigned int j = 0; j < jets->size(); ++j){

        const pat::Jet& jet = (*jets)[j];

	if ( jet.pt() < jetPtCut_ ) continue;
	if ( jet.eta() > jetEtaCut_ ) continue;

        if (aggregateHF_) {
	  
            if (doGenJets_ && isMC_) {

                //  std::cout << "------->Aggregating HF for gen jet" << std::endl;               
                reco::PFCandidate outputPseudoHF;
                std::vector<reco::PFCandidate> constituentsNoHF;

                const reco::GenJet *genJet = jet.genJet();

                // Particle collection to aggregate into pseudo-Bs
                //std::vector<edm::Ptr<reco::Candidate>> inputJetConstituents = genJet->getJetConstituents();
                std::map<int, std::vector<edm::Ptr<reco::Candidate>>> hfConstituentsMap;
                reco::Candidate::PolarLorentzVector totalPseudoHF(0., 0., 0., 0.);

                //int HFjet = 0;
                // Go over gen particles

                if (genJet) {
		  //		  std::cout << "Gen pt: " << genJet->pt() << std::endl;
                    for (const edm::Ptr<reco::Candidate> &constit : (*genJet).getJetConstituents()) {
   
		        if(chargedOnly_ && constit->charge() == 0) continue;

		        bool isNeutrino = (constit->pdgId() == 12); // nue
                        isNeutrino &= (constit->pdgId() == 14); // numu
                        isNeutrino &= (constit->pdgId() == 16); // nutau
                        isNeutrino &= (constit->pdgId() == 18); // nutau'
                        if (isNeutrino) {
                            std::cout << "found a neutrino" << std::endl;
                            continue;
                        }

			if (constit->charge() == 0 or constit->pt() < ptCut_) {   // Neutrals and low pT tracks are not used in aggregation, let's fill them directly in the collection
			  reco::PFCandidate constituentsPF;
			  constituentsPF.setCharge(constit->charge());
			  constituentsPF.setP4(constit->p4());
			  constituentsPF.setPdgId(constit->pdgId());
			  constituentsNoHF.push_back(constituentsPF);
			  newPFCandCollection->push_back(constituentsPF);
			  //			  std::cout << "Filling a gen const with ch " << constit->charge() << " " << constit->pt() <<" id: " << constit->pdgId() << std::endl;
			  continue;
			}

                        //std::cout << "---- Got gen level constit with charge " <<  constit->charge() << std::endl;
 
                        // Get status of matched gen particle
                        int status = 1;
                        if ((*candToGenParticleMap).find(constit) != (*candToGenParticleMap).end()) {
                            edm::Ptr<pat::PackedGenParticle> matchGenParticle = (*candToGenParticleMap).at(constit);
                            status = matchGenParticle->status();
                            //std::cout << "---- constit status in the gen map " << status << std::endl;               
                        }
                        // Add particle to output collection or from HF map
                        if (status >= 100) {
                            hfConstituentsMap[status].push_back(constit);
                        }

                        else if (status == 1) {
                            //add for PF candidate collection output
                            reco::PFCandidate constituentsPF;
                            constituentsPF.setCharge(constit->charge());
                            constituentsPF.setP4(constit->p4());
                            constituentsPF.setPdgId(constit->pdgId());
                            constituentsNoHF.push_back(constituentsPF);
                            newPFCandCollection->push_back(constituentsPF);
                        }

                    } // end constit loop 

                    // Aggregate particles coming from HF decays into pseudo-B/D's and add them to the collection
                    for (auto it = hfConstituentsMap.begin(); it != hfConstituentsMap.end(); ++it) {
		      //                        reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
                        for (edm::Ptr<reco::Candidate> hfConstituent : it->second) {
                            reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
                            productLorentzVector.SetPt(hfConstituent->pt());
                            productLorentzVector.SetEta(hfConstituent->eta());
                            productLorentzVector.SetPhi(hfConstituent->phi());
                            productLorentzVector.SetM(hfConstituent->mass());
			    //     pseudoHF += productLorentzVector;
                            totalPseudoHF += productLorentzVector;                           
                            reco::PFCandidate daughter;
                            daughter.setP4(productLorentzVector);
                            outputPseudoHF.addDaughter(daughter);
                        }
                    } // end map loop
                    

                    if (hfConstituentsMap.size() > 0) {
                        outputPseudoHF.setP4(totalPseudoHF);
                        outputPseudoHF.setPdgId(211);
                        outputPseudoHF.setCharge(-5);
                        outputPseudoHF.setMass(totalPseudoHF.mass()*-1);
     
                        newPFCandCollection->push_back(outputPseudoHF);
			//			std::cout << "Filling a gen aggregatedB with pt " << outputPseudoHF.pt() << std::endl;
			//                        newPFCandCollectionHF->insert(newPFCandCollectionHF->end(), constituentsNoHF.begin(), constituentsNoHF.end());
			//                        newPFCandCollectionHF->push_back(outputPseudoHF);
                    }
                }
            }
                
            else {

	      // std::cout << "------->Aggregating HF for reco jet" << std::endl;
	      
	      int matchIndex = -1;
	      if (domatch_) {   // 
		double drMin = 100;
		for (unsigned int imatch = 0; imatch < matchedjets->size(); ++imatch) {
		  const pat::Jet& mjet = (*matchedjets)[imatch];
		  double dr = deltaR(jet, mjet);
		  if (dr < drMin) {
		    drMin = dr;
		    matchIndex = imatch;
		  }
		}
	      }
	      
	      //    std::cout << "MatchIndex: " << matchIndex << std::endl;

                reco::TrackToGenParticleMap recoMap = isMC_ ? *candToGenParticleMap : reco::TrackToGenParticleMap();	        
		//		const pat::Jet& mjet = (*matchedjets)[matchIndex];
		const pat::Jet& mjet = (domatch_ ? (*matchedjets)[matchIndex] : jet);
		//		std::cout << "Test jet pT:s " << mjet.pt() << " " << jet.pt() << std::endl;
                std::vector<edm::Ptr<reco::Candidate>> inputJetConstituents = jet.getJetConstituents();
		//                std::vector<reco::PFCandidate> droppedTracks = {};
                reco::PFCandidate outputPseudoHF;
                std::vector<reco::PFCandidate> constituentsNoHF;

                // Particle collection to aggregate into pseudo-Bs
                std::map<int, std::vector<edm::Ptr<reco::Candidate>>> hfConstituentsMap;
                reco::Candidate::PolarLorentzVector totalPseudoHF(0., 0., 0., 0.);

                // Grab the IP and SV tag info from the jet
		/*	std::cout << "Grab tag infos" << std::endl;
                std::cout << "HAS IPtag " << ipTagInfoLabel_.c_str() << " " << jet.hasTagInfo(ipTagInfoLabel_.c_str()) << std::endl;
                std::cout << "HAS svtag " << svTagInfoLabel_.c_str() << " " << jet.hasTagInfo(svTagInfoLabel_.c_str()) << std::endl;

		std::cout << "Grab tag infos for mjet" << std::endl;
                std::cout << "HAS IPtag " << ipTagInfoLabel_.c_str() << " " << mjet.hasTagInfo(ipTagInfoLabel_.c_str()) << std::endl;
                std::cout << "HAS svtag " << svTagInfoLabel_.c_str() << " " << mjet.hasTagInfo(svTagInfoLabel_.c_str()) << std::endl; */

                const reco::CandIPTagInfo *ipTagInfo = mjet.tagInfoCandIP(ipTagInfoLabel_.c_str());
                const std::vector<reco::btag::TrackIPData> ipData = ipTagInfo->impactParameterData();
                const std::vector<edm::Ptr<reco::Candidate>> ipTracks = ipTagInfo->selectedTracks();

                const reco::CandSecondaryVertexTagInfo *svTagInfo = mjet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());

                for (const edm::Ptr<reco::Candidate> &constit : mjet.getJetConstituents()) {
                    //std::cout << "constit loop" << std::endl;
		  if (chargedOnly_ && constit->charge() == 0) continue;
		  if (constit->charge() == 0 or constit->pt() < ptCut_) {   // Neutrals and low pT tracks are not used in aggregation, let's fill them directly in the collection
		    reco::PFCandidate constituentsPF;
		    constituentsPF.setCharge(constit->charge());
		    constituentsPF.setP4(constit->p4());
		    constituentsPF.setPdgId(constit->pdgId());
		    constituentsNoHF.push_back(constituentsPF);
		    newPFCandCollection->push_back(constituentsPF);
		    // std::cout << "Filling a const with ch " << constit->charge() << " " << constit->pt() << std::endl;
		    continue;
		  }
		  
                    // Look for particle in ipTracks
                    auto itIPTrack = std::find(ipTracks.begin(), ipTracks.end(), constit);
                    if (itIPTrack == ipTracks.end()) {
		      // std::cout << "Didn't find a track in ipTrack list: " << constit->charge() << " " << constit->pt() << std::endl;
		      reco::PFCandidate constituentsPF;
		      constituentsPF.setCharge(constit->charge());
		      constituentsPF.setP4(constit->p4());
		      constituentsPF.setPdgId(constit->pdgId());
		      constituentsNoHF.push_back(constituentsPF);
		      newPFCandCollection->push_back(constituentsPF);
		      continue;
		    }
		    
                    // For track inefficiency uncertainty 
                    /* const double range_from = 0;
                    const double range_to = 1;
                    std::random_device rand_dev;
                    std::mt19937 generator(rand_dev());
                    std::uniform_real_distribution<double> distr(range_from, range_to);
                    double rand_uniform = distr(generator);
                    if (rand_uniform < trkInefRate_) {
                        reco::PFCandidate droppedTrack;
                        droppedTrack.setCharge(constit->charge());
                        droppedTrack.setP4(constit->p4());
                        droppedTrack.setPdgId(constit->pdgId());
                        droppedTracks.push_back(droppedTrack);
                        continue;
                    } */

                    int status = 1;

                    if (isMC_ && withTruthInfo_) {
                        if (recoMap.find(constit) != recoMap.end()) {
                            edm::Ptr<pat::PackedGenParticle> matchGenParticle = recoMap.at(constit);
                            status = matchGenParticle->status();
                            //std::cout << " matchgenptcl status = " << status << " pdgid " << matchGenParticle->pdgId() << std::endl;
                        }
                    }
                    else {
                        // Initialize values 
                        const double missing_value = -9999;
                        float ip3dSig = missing_value;
                        float ip2dSig = missing_value;
                        float trkdz = missing_value;
                        float distanceToJetAxis = missing_value;
                        //bool isLepton = false;
                        bool inSV = false;

                        float svtxdls = missing_value;
                        float svtxdls2d = missing_value;
                        float svtxm = missing_value;
                        float svtxmcorr = missing_value;
                        float svtxNtrk = missing_value;
                        float svtxchi2 = missing_value;
                        float svtxTrkPt = missing_value;
                        //float svtxTrkPtOverSv = missing_value;

                        float jtpt = jet.pt();
                        
                        // Get IP info 
                        int trkIPIndex = itIPTrack - ipTracks.begin();
                        const reco::btag::TrackIPData trkIPdata = ipData[trkIPIndex];
                        ip3dSig = trkIPdata.ip3d.significance();
                        ip2dSig = trkIPdata.ip2d.significance();
                        distanceToJetAxis = trkIPdata.distanceToJetAxis.value();
                        //int pdg = constit->pdgId();
                        //isLepton = (std::abs(pdg) == 11) || (std::abs(pdg) == 13);

                        const reco::Track *constitTrack = constit->bestTrack();
                        if (constitTrack) {
			    //                            std::cout << "testTrack dz " << testTrack->dz() << std::endl;
                            trkdz = constitTrack->dz(primaryVertices->at(0).position());
                        } else {
                            trkdz = missing_value;
                        }
  
                        // if nan go back to missing_value
                        if (ip3dSig != ip3dSig) ip3dSig = missing_value;
                        if (ip2dSig != ip2dSig) ip2dSig = missing_value;
                        if (distanceToJetAxis != distanceToJetAxis) distanceToJetAxis = missing_value;
                        if (trkdz != trkdz) trkdz = missing_value;

                        // Get SV info           
                        for (uint ivtx = 0; ivtx < svTagInfo->nVertices(); ivtx++) {
                            std::vector<edm::Ptr<reco::Candidate>> isvTracks = svTagInfo->vertexTracks(ivtx);
                            auto itSVTrack = std::find(isvTracks.begin(), isvTracks.end(), constit);
                            if (itSVTrack == isvTracks.end()) continue;
                            
                            inSV = true;

                            svtxNtrk = (float) svTagInfo->nVertexTracks(ivtx);

                            Measurement1D m1D = svTagInfo->flightDistance(ivtx, 0);
                            svtxdls = m1D.significance();

                            Measurement1D m2D = svTagInfo->flightDistance(ivtx, 2);
                            svtxdls2d = m2D.significance();

                            const reco::VertexCompositePtrCandidate svtx = svTagInfo->secondaryVertex(ivtx);
                            svtxm = svtx.p4().mass();

                            double svtxpt = svtx.p4().pt();
                            svtxTrkPt = svtxpt;
                            //svtxTrkPtOverSv = constit->pt() / svtxpt;
        
                            double sinth = svtx.p4().Vect().Unit().Cross((svTagInfo->flightDirection(ivtx)).unit()).Mag2();
                            sinth = sqrt(sinth);
                            double underRoot = std::pow(svtxm, 2) + (std::pow(svtxpt, 2) * std::pow(sinth, 2));
                            svtxmcorr = std::sqrt(underRoot) + (svtxpt * sinth);

                            //svtxnormchi2 = svtx.vertexNormalizedChi2();
                            svtxchi2 = svtx.vertexChi2();  
                            break;
                        } // end vtx loop

                        if (withCuts_) {
                            if (inSV || (ip3dSig > 2.5)) {
                                status = 100;
                            }
                        }

                        else if (withTMVA_) {
                            std::map<std::string, float> inputs;
                            inputs["trkIp3dSig"] = ip3dSig;
                            inputs["trkIp2dSig"] = ip2dSig;
                            inputs["trkDistToAxis"] = distanceToJetAxis;
                            inputs["trkDz"] = trkdz;
                            inputs["svtxdls"] = svtxdls;
                            inputs["svtxdls2d"] = svtxdls2d;
                            inputs["svtxm"] = svtxm;
                            inputs["svtxmcorr"] = svtxmcorr;
                            inputs["svtxchi2"] = svtxchi2;
                            inputs["svtxNtrk"] = svtxNtrk;
                            inputs["svtxpt"] = svtxTrkPt;
                            inputs["jtpt"] = jtpt;
                            inputs["hiBin"] = hiBin;

                            float prediction = -99.;

                            prediction = tmvaTagger->evaluate(inputs);

                            if (prediction > -0.3) {
                                status = 100;
                            }

                        } // endif 
                    } // endif *not* with truth info


                    // Add particle to output collection or from HF map
                    if (status == 1) {
                        fastjet::PseudoJet outConstit(constit->px(), constit->py(), constit->pz(), constit->energy());
                        //add for PF candidate collection output
                        reco::PFCandidate constituentsPF;
                        constituentsPF.setCharge(constit->charge());
                        constituentsPF.setP4(constit->p4());
                        constituentsPF.setPdgId(constit->pdgId());
			//                        constituentsNoHF.push_back(constituentsPF);
                        newPFCandCollection->push_back(constituentsPF);
                    }
                    else if (status >= 100) {
                        hfConstituentsMap[status].push_back(constit);
			//			std::cout << "HF Constit pt: " << constit->pt() << " id: " << constit->pdgId() << " mass " << constit->mass() << std::endl;
                    }
    
                } // end jet constituents loop
   
                // Aggregate particles coming from HF decays into pseudo-B/C's  and add them to the collection
                for (auto itTrackFromHF = hfConstituentsMap.begin(); itTrackFromHF != hfConstituentsMap.end(); itTrackFromHF++) {
		  //                    reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
                    for (edm::Ptr<reco::Candidate> hfConstituent : itTrackFromHF->second) {
                        reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
                        productLorentzVector.SetPt(hfConstituent->pt());
                        productLorentzVector.SetEta(hfConstituent->eta());
                        productLorentzVector.SetPhi(hfConstituent->phi());
                        productLorentzVector.SetM(hfConstituent->mass());
			//      pseudoHF += productLorentzVector;
                        totalPseudoHF += productLorentzVector;

                        reco::PFCandidate daughter;
                        daughter.setP4(productLorentzVector);
                        outputPseudoHF.addDaughter(daughter);
                    }
                } // end tracks from B loop

                if (hfConstituentsMap.size() > 0) {
                    outputPseudoHF.setP4(totalPseudoHF);
                    outputPseudoHF.setPdgId(211); // Charged hadron
                    //std::cout << "Setting id: " << outputPseudoHF.pdgId() << " mass " << outputPseudoHF.mass() << std::endl;
                    outputPseudoHF.setCharge(-5); // Set charge to avoid problems with chargedOnly selections
                    outputPseudoHF.setMass(totalPseudoHF.mass()*(-1));
 
                    newPFCandCollection->push_back(outputPseudoHF);

		}
		//  newPFCandCollectionHF->insert(newPFCandCollectionHF->end(), constituentsNoHF.begin(), constituentsNoHF.end());
                //    newPFCandCollectionHF->push_back(outputPseudoHF);
		    // }
            } 
        }

        // This part is directly from the pp analysis
        else {
            // std::cout << "\tNot aggregating" << std::endl;
            std::vector<edm::Ptr<reco::Candidate>> constituents = {}; 
            if (doGenJets_ && isMC_) {
                const reco::GenJet *genJet = jet.genJet();
                if (genJet) constituents = genJet->getJetConstituents();
            } 
            else {
                constituents = jet.getJetConstituents();
            }

            for (edm::Ptr<reco::Candidate> constituent : constituents) {
                if ((chargedOnly_) && (constituent->charge() == 0)) continue;
                if (constituent->pt() < ptCut_) continue;
                jetConstituents.push_back(fastjet::PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
            }
        }

    } // end jet loop

    //iEvent.put(std::move(newPFCandCollectionHF));
    iEvent.put(std::move(newPFCandCollection));
  
}


void aggregatedPFCollection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Aggregated PF Collection");

  // Input collections
  desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("matchTag", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("candToGenParticleMap", edm::InputTag("TrackToGenParticleMapProducer"));

  desc.add<edm::InputTag>("CentralityBinSrc", edm::InputTag("centralityBin","HFtowers"));

  // Configuration parameters
  desc.add<bool>("isMC", true);
  desc.add<bool>("chargedOnly", false);
  
  desc.add<double>("ptCut", 1.);
  desc.add<double>("jetPtCut", 50.);
  desc.add<double>("trkInefRate", 0.);

  desc.add<bool>("doGenJets", false);
  desc.add<bool>("domatch", true);

  desc.add<bool>("aggregateHF", true);
  desc.add<bool>("aggregateWithTruthInfo", true);
  desc.add<bool>("aggregateWithCuts", false);
  desc.add<bool>("aggregateWithTMVA", false);
  desc.add<edm::FileInPath>("tmva_path", edm::FileInPath("RecoHI/HiJetAlgos/data/TMVAClassification_BDTG.weights.xml"));
  //  desc.add<edm::FileInPath>("tmva_path", edm::FileInPath(""));

  desc.add<std::vector<std::string>>("tmva_variables", {});
  desc.add<std::vector<std::string>>("tmva_spectators", {});
  // Tag info labels
  desc.add<std::string>("ipTagInfoLabel", "pfImpactParameter");
  desc.add<std::string>("svTagInfoLabel", "pfInclusiveSecondaryVertexFinder");

  descriptions.add("aggregatedPFCands", desc);
}

using aggregatedPFCands = aggregatedPFCollection;

// define this as a plug-in
DEFINE_FWK_MODULE(aggregatedPFCands);
