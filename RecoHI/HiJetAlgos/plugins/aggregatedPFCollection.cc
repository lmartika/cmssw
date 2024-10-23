/*

  Perform dynamic grooming on jets
  Cannot be done in FastjetJetProducer, as this algorithm is not part of fastjet
  I borrowed liberally from CompoundJetProducer
  If I was a better programmer, I would just call that class

  -Matt Nguyen, 02-02-2022 (Groundhog's day) 

*/

#include <memory>
#include <vector>
#include <cmath>
#include <map>
#include <random>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"


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

#include "HeavyIonsAnalysis/JetAnalysis/interface/HiInclusiveJetAnalyzer.h"

//template <class T>

class aggregatedPFCollection : public edm::EDProducer {
public:
    explicit aggregatedPFCollection(const edm::ParameterSet&);
    ~aggregatedPFCollection();

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;


    // ------------- member data ----------------------------
    //edm::EDGetTokenT<std::vector<>> jetSrc_;
    edm::EDGetTokenT<pat::JetCollection> jetSrc_;
    edm::EDGetTokenT<std::vector<reco::PFCandidate>> constitSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>>  packedConstitSrc_;
    edm::EDGetTokenT<reco::TrackToGenParticleMap> candToGenParticleMapToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

    std::unique_ptr<TMVAEvaluator> tmvaTagger;

    bool isMC_;

    bool chargedOnly_;
    bool doGenJets_;
  
    double ptCut_;
    double trkInefRate_; // 0 by default

    bool aggregateHF_;
    bool withTruthInfo_;
    bool withCuts_;
    bool withXGB_;
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
  
    ptCut_ = iConfig.getParameter<double>("ptCut");
    trkInefRate_ = iConfig.getParameter<double>("trkInefRate");

    if (aggregateHF_) {
        withTruthInfo_ = iConfig.getParameter<bool>("aggregateWithTruthInfo");
        withCuts_ = iConfig.getParameter<bool>("aggregateWithCuts");
        withXGB_ = iConfig.getParameter<bool>("aggregateWithXGB");
        withTMVA_ = iConfig.getParameter<bool>("aggregateWithTMVA");
        if (withXGB_) {
            //xgb_path_ = iConfig.getParameter<edm::FileInPath>("xgb_path");
        }

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
    constitSrc_ = consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"));
    packedConstitSrc_ = consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("constitSrc"));

    if (aggregateHF_ && isMC_) candToGenParticleMapToken_ = consumes<reco::TrackToGenParticleMap>(iConfig.getParameter<edm::InputTag>("candToGenParticleMap"));

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


}

aggregatedPFCollection::~aggregatedPFCollection()
{

   // do anything here that needs to be done at desctruction time
   //    // (e.g. close files, deallocate resources etc.)
   //
   //    }
   //
}

void aggregatedPFCollection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

    auto newPFCandCollection = std::make_unique<reco::PFCandidateCollection>();

    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetSrc_, jets);    


    //    edm::Handle<std::vector<reco::PFCandidate>> pfcands;
    // bool isPF = iEvent.getByToken(constitSrc_, pfcands);

    // edm::Handle<edm::View<pat::PackedCandidate>> pfcandsPacked;
    // bool isPackedPF = iEvent.getByToken(packedConstitSrc_, pfcandsPacked);

    // -- For aggregation -- //
    edm::Handle<reco::TrackToGenParticleMap> candToGenParticleMap;

    if (aggregateHF_ && isMC_) {
        iEvent.getByToken(candToGenParticleMapToken_, candToGenParticleMap);
    } 

  
    std::vector<fastjet::PseudoJet> jetConstituents = {};
    reco::PFCandidate pseudoHF;


    for(unsigned int j = 0; j < jets->size(); ++j){

        const pat::Jet& jet = (*jets)[j];

        if (aggregateHF_) {
      
            if (doGenJets_ && isMC_) {

                // std::cout << "------->Aggregating HF for gen jet" << std::endl;               
                reco::PFCandidate outputPseudoHF;

                const reco::GenJet *genJet = jet.genJet();

                // Particle collection to aggregate into pseudo-Bs
                std::vector<edm::Ptr<reco::Candidate>> inputJetConstituents = genJet->getJetConstituents();
                std::map<int, std::vector<edm::Ptr<reco::Candidate>>> hfConstituentsMap;
                reco::Candidate::PolarLorentzVector totalPseudoHF(0., 0., 0., 0.);

                // Go over gen particles
                for (edm::Ptr<reco::Candidate> constit : inputJetConstituents) {            
                    if(chargedOnly_ && constit->charge() == 0) continue;
                    if(constit->pt() < ptCut_) continue;

                    bool isNeutrino = (constit->pdgId() == 12); // nue
                    isNeutrino &= (constit->pdgId() == 14); // numu
                    isNeutrino &= (constit->pdgId() == 16); // nutau
                    isNeutrino &= (constit->pdgId() == 18); // nutau'
                    if (isNeutrino) {
                        std::cout << "found a neutrino" << std::endl;
                        continue;
                    }

                    // Get status of matched gen particle
                    edm::Ptr<pat::PackedGenParticle> matchGenParticle = (*candToGenParticleMap).at(constit);
                    int status = matchGenParticle->status();

                    // Add particle to output collection or from HF map
                    if (status == 1) {
                        //add for PF candidate collection output
                        reco::PFCandidate constituentsPF;
                        constituentsPF.setCharge(constit->charge());
                        constituentsPF.setP4(constit->p4());
                        constituentsPF.setPdgId(constit->pdgId());
                        newPFCandCollection->push_back(constituentsPF);
                    } 
                    else if (status >= 100) {
                        hfConstituentsMap[status].push_back(constit);
                    }

                } // end constit loop 

                // Aggregate particles coming from HF decays into pseudo-B/D's and add them to the collection
                for (auto it = hfConstituentsMap.begin(); it != hfConstituentsMap.end(); ++it) {
                    reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
                    for (edm::Ptr<reco::Candidate> hfConstituent : it->second) {
                        reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
                        productLorentzVector.SetPt(hfConstituent->pt());
                        productLorentzVector.SetEta(hfConstituent->eta());
                        productLorentzVector.SetPhi(hfConstituent->phi());
                        productLorentzVector.SetM(hfConstituent->mass());
                        pseudoHF += productLorentzVector;
                        totalPseudoHF += productLorentzVector;                           
                        reco::PFCandidate daughter;
                        daughter.setP4(productLorentzVector);
                        outputPseudoHF.addDaughter(daughter);
                    }
                } // end map loop
                    
                outputPseudoHF.setP4(totalPseudoHF);
                newPFCandCollection->push_back(outputPseudoHF);
            }
                
            else {

		std::cout << "------->Aggregating HF for reco jet" << std::endl; 

                reco::TrackToGenParticleMap recoMap = isMC_ ? *candToGenParticleMap : reco::TrackToGenParticleMap();	        

                std::vector<edm::Ptr<reco::Candidate>> inputJetConstituents = jet.getJetConstituents();
                std::vector<reco::PFCandidate> droppedTracks = {};
                reco::PFCandidate outputPseudoHF;


                // Particle collection to aggregate into pseudo-Bs
                std::map<int, std::vector<edm::Ptr<reco::Candidate>>> hfConstituentsMap;
                reco::Candidate::PolarLorentzVector totalPseudoHF(0., 0., 0., 0.);

                // Grab the IP and SV tag info from the jet
                const reco::CandIPTagInfo *ipTagInfo = jet.tagInfoCandIP(ipTagInfoLabel_.c_str());
                const std::vector<reco::btag::TrackIPData> ipData = ipTagInfo->impactParameterData();
                const std::vector<edm::Ptr<reco::Candidate>> ipTracks = ipTagInfo->selectedTracks();

                const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());

                for (const edm::Ptr<reco::Candidate> constit : jet.getJetConstituents()) {
                    if (chargedOnly_ && constit->charge() == 0) continue;
                    if (constit->pt() < ptCut_) continue;

                    // Look for particle in ipTracks
                    //auto itIPTrack = std::find(ipTracks.begin(), ipTracks.end(), constit);
                    //if (itIPTrack == ipTracks.end()) continue;
		    
		    reco::CandidatePtr itIPTrack;
		    int trkIPIndex = -1;
		    for (auto iterIPTrack : ipTracks) {
		      float eps = 1e-5;                                                                                                                                
		      trkIPIndex++;
		        
		      if (std::abs(constit->eta()-iterIPTrack->eta())>eps) continue;
		      else if (std::abs(constit->phi()-iterIPTrack->phi())>eps) continue;
		      else if (std::abs(constit->pt()-iterIPTrack->pt())>eps) continue;
		      itIPTrack = iterIPTrack;
		    }

		    if (!(itIPTrack)) continue; 

                    // For track inefficiency uncertainty 
                    const double range_from = 0;
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
                    }

                    int status = 1;

                    if (isMC_ && withTruthInfo_) {
                        if (recoMap.find(constit) != recoMap.end()) {
                            edm::Ptr<pat::PackedGenParticle> matchGenParticle = recoMap.at(constit);
                            status = matchGenParticle->status();
                        }
                    }
                    else {
                        // Initialize values 
                        const double missing_value = -1000000.;

                        float ip3dSig = missing_value;
                        float ip2dSig = missing_value;
                        float distanceToJetAxis = missing_value;
			//                        bool isLepton = false;
                        bool inSV = false;

                        float svtxdls = missing_value;
                        float svtxdls2d = missing_value;
                        float svtxm = missing_value;
                        float svtxmcorr = missing_value;
                        float svtxNtrk = missing_value;
                        float svtxnormchi2 = missing_value;
                        float svtxTrkPtOverSv = missing_value;

                        float jtpt = jet.pt();
                        
                        // Get IP info 
                        const reco::btag::TrackIPData trkIPdata = ipData[trkIPIndex];
                        ip3dSig = trkIPdata.ip3d.significance();
                        ip2dSig = trkIPdata.ip2d.significance();
                        distanceToJetAxis = trkIPdata.distanceToJetAxis.value();
			//                         int pdg = constit->pdgId();
			//                        isLepton = (std::abs(pdg) == 11) || (std::abs(pdg) == 13);

                        // if nan go back to missing_value
                        if (ip3dSig != ip3dSig) ip3dSig = missing_value;
                        if (ip2dSig != ip2dSig) ip2dSig = missing_value;
                        if (distanceToJetAxis != distanceToJetAxis) distanceToJetAxis = missing_value;

                        // Get SV info
            
                        for (uint ivtx = 0; ivtx < svTagInfo->nVertices(); ivtx++) {
                            std::vector<edm::Ptr<reco::Candidate>> isvTracks = svTagInfo->vertexTracks(ivtx);
                            // auto itSVTrack = std::find(isvTracks.begin(), isvTracks.end(), constit);
                            // if (itSVTrack == isvTracks.end()) continue;
			    
			    reco::CandidatePtr itSVTrack;
			    for (auto iterSVTrack : isvTracks) {
			      float eps = 1e-5;

			      if (std::abs(constit->eta()-iterSVTrack->eta())>eps) continue;
			      else if (std::abs(constit->phi()-iterSVTrack->phi())>eps) continue;
                              else if (std::abs(constit->pt()-iterSVTrack->pt())>eps) continue;

			      itSVTrack = iterSVTrack;
			    }

			    if (!(itSVTrack)) continue; 

                            inSV = true;

                            svtxNtrk = (float) svTagInfo->nVertexTracks(ivtx);

                            Measurement1D m1D = svTagInfo->flightDistance(ivtx, 0);
                            svtxdls = m1D.significance();

                            Measurement1D m2D = svTagInfo->flightDistance(ivtx, 2);
                            svtxdls2d = m2D.significance();

                            const reco::VertexCompositePtrCandidate svtx = svTagInfo->secondaryVertex(ivtx);
                            svtxm = svtx.p4().mass();

                            double svtxpt = svtx.p4().pt();
                            svtxTrkPtOverSv = constit->pt() / svtxpt;
        
                            double sinth = svtx.p4().Vect().Unit().Cross((svTagInfo->flightDirection(ivtx)).unit()).Mag2();
                            sinth = sqrt(sinth);
                            double underRoot = std::pow(svtxm, 2) + (std::pow(svtxpt, 2) * std::pow(sinth, 2));
                            svtxmcorr = std::sqrt(underRoot) + (svtxpt * sinth);

                            svtxnormchi2 = svtx.vertexNormalizedChi2();
                            svtxTrkPtOverSv = constit->pt() / svtxpt;

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
                            inputs["svtxdls"] = svtxdls;
                            inputs["svtxdls2d"] = svtxdls2d;
                            inputs["svtxm"] = svtxm;
                            inputs["svtxmcorr"] = svtxmcorr;
                            inputs["svtxnormchi2"] = svtxnormchi2;
                            inputs["svtxNtrk"] = svtxNtrk;
                            inputs["svtxTrkPtOverSv"] = svtxTrkPtOverSv;
                            inputs["jtpt"] = jtpt;

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
                        newPFCandCollection->push_back(constituentsPF);
                    }
                    else if (status >= 100) {
                        hfConstituentsMap[status].push_back(constit);
                    }
    
                } // end jet constituents loop
    
                // Aggregate particles coming from HF decays into pseudo-B/C's  and add them to the collection
                for (auto itTrackFromHF = hfConstituentsMap.begin(); itTrackFromHF != hfConstituentsMap.end(); itTrackFromHF++) {
		  std::cout << "In aggregation loop" << std::endl;
                    reco::Candidate::PolarLorentzVector pseudoHF(0., 0., 0., 0.);
                    for (edm::Ptr<reco::Candidate> hfConstituent : itTrackFromHF->second) {
                        reco::Candidate::PolarLorentzVector productLorentzVector(0., 0., 0., 0.);
                        productLorentzVector.SetPt(hfConstituent->pt());
                        productLorentzVector.SetEta(hfConstituent->eta());
                        productLorentzVector.SetPhi(hfConstituent->phi());
                        productLorentzVector.SetM(hfConstituent->mass());
                        pseudoHF += productLorentzVector;
                        totalPseudoHF += productLorentzVector;

                        reco::PFCandidate daughter;
                        daughter.setP4(productLorentzVector);
                        outputPseudoHF.addDaughter(daughter);
                    }

                } // end tracks from B loop

                outputPseudoHF.setP4(totalPseudoHF);
                newPFCandCollection->push_back(outputPseudoHF);

            } 
        }

        else {
            std::cout << "\tNot aggregating" << std::endl;
            std::vector<edm::Ptr<reco::Candidate>> constituents = {}; 
            if (doGenJets_ && isMC_) {
                const reco::GenJet *genJet = jet.genJet();
                if (genJet) constituents = genJet->getJetConstituents();
            } 
            else {
                constituents = jet.getJetConstituents();
		std::cout << "try jet pt " << jet.pt() << " " << std::endl;
            }

	    for (edm::Ptr<reco::Candidate> constituent : constituents) {
	      //std::cout << "Test const " << constituent->pt() << std::endl;
	      if ((chargedOnly_) && (constituent->charge() == 0)) continue;
	      if (constituent->pt() < ptCut_) continue;
	      jetConstituents.push_back(fastjet::PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
		}  
        }

    } // end jet loop
    
    iEvent.put(std::move(newPFCandCollection));
  
}


//template <class T>
void aggregatedPFCollection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("New PF Collection");

  // Input collections
  //   desc.add<edm::InputTag>("jetSrc", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("jetSrc", edm::InputTag("updatedPatJets"));
  desc.add<edm::InputTag>("constitSrc", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("candToGenParticleMap", edm::InputTag("TrackToGenParticleMapProducer"));

  // Configuration parameters
  desc.add<bool>("isMC", true);
  desc.add<bool>("writeConstits", false);
  desc.add<bool>("doLateKt", false);
  desc.add<bool>("chargedOnly", false);
  
  desc.add<double>("zcut", 0.1);
  desc.add<double>("beta", 0.0);
  desc.add<double>("dynktcut", 1.0);
  desc.add<double>("ktcut", 1.0);
  desc.add<double>("rParam", 0.4);
  desc.add<double>("ptCut", 1.);
  desc.add<double>("trkInefRate", 0.); // between 0 and 1

  desc.add<bool>("doGenJets", false);

  desc.add<bool>("aggregateHF", false);
  desc.add<bool>("aggregateWithTruthInfo", true);
  desc.add<bool>("aggregateWithCuts", false);
  desc.add<bool>("aggregateWithXGB", false);
  desc.add<bool>("aggregateWithTMVA", false);
  // desc.add<edm::FileInPath>("xgb_path", edm::FileInPath("RecoHI/HiJetAlgos/data/dummy.model"));
  // desc.add<edm::FileInPath>("tmva_path", edm::FileInPath("RecoHI/HiJetAlgos/data/dummy.weights.xml"));
  desc.add<std::vector<std::string>>("tmva_variables", {});
  desc.add<std::vector<std::string>>("tmva_spectators", {});
  // Tag info labels
  desc.add<std::string>("ipTagInfoLabel", "pfImpactParameter");
  desc.add<std::string>("svTagInfoLabel", "pfInclusiveSecondaryVertexFinder");

  descriptions.add("aggregatedPFCands", desc);
}

// ------------ method called once each job just before starting event loop  ------------
void
aggregatedPFCollection::beginJob() {
}

// // ------------ method called once each job just after ending the event loop  ------------
void
aggregatedPFCollection::endJob() {
}
 
using aggregatedPFCands = aggregatedPFCollection;

// define this as a plug-in
DEFINE_FWK_MODULE(aggregatedPFCands);
