// -*- C++ -*-
//
// Package:    RhoAnalyser/RandomConeAnalysis
// Class:      RandomConeAnalysis
//
/**\class RandomConeAnalysis RandomConeAnalysis.cc RhoAnalyser/RandomConeAnalysis/plugins/RandomConeAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]

*/
//
// Original Author:  Vangelis Vladimirov (copied from FJ analyzer)
//         Created:  Thu, 27 Jul 2023 12:13:58 GMT
//
//

#include "DataFormats/Common/interface/Handle.h"
 #include "DataFormats/Common/interface/View.h"
 #include "DataFormats/JetReco/interface/Jet.h"
 #include "DataFormats/Math/interface/deltaR.h"
 #include "DataFormats/PatCandidates/interface/PackedCandidate.h"
 #include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
 #include "DataFormats/Provenance/interface/EventID.h"
 #include "FWCore/Framework/interface/Event.h"
 #include "FWCore/Framework/interface/MakerMacros.h"
 #include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
 #include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
 #include "FWCore/ParameterSet/interface/ParameterSetDescriptionFiller.h"
 #include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
 #include "DataFormats/JetReco/interface/JetCollection.h"
 #include "FWCore/Framework/interface/stream/EDProducer.h"
 #include "FWCore/Utilities/interface/StreamID.h"
 #include "FWCore/Utilities/interface/EDGetToken.h"
 #include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "DataFormats/Candidate/interface/Candidate.h"
 #include "DataFormats/Candidate/interface/CandidateFwd.h"
 #include "DataFormats/Common/interface/Association.h"
 #include "DataFormats/Common/interface/RefVector.h"
 #include "DataFormats/Math/interface/deltaPhi.h"
 #include "DataFormats/PatCandidates/interface/CovarianceParameterization.h"
 #include "DataFormats/VertexReco/interface/Vertex.h"
 #include "DataFormats/VertexReco/interface/VertexFwd.h"
 #include "FWCore/Utilities/interface/thread_safety_macros.h"

#include "TF1.h"
 #include "TF2.h"
 #include "TH1.h"
#include "TMath.h"
 #include "Math/ProbFuncMathCore.h"
 #include "TMinuitMinimizer.h"
 
 #include <algorithm>
 #include <cmath>
 #include <memory>
 #include <string>
 #include <utility>
 #include <vector>


#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

//ROOT include files
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "Math/ProbFunc.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

// user include files

// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace edm;
using namespace reco;
using namespace pat;

float delta_phi(float phi1, float phi2){
    float result = phi1 - phi2;
    while ( result < -M_PI )
      {
        result += 2.*M_PI;
      }
    while ( result > M_PI )
      {
        result -= 2.*M_PI;
      }
    return result;
}

class RandomConeAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RandomConeAnalysis(const edm::ParameterSet& iConfig);
  ~RandomConeAnalysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  // void fillDescriptionsFromRandomConeAnalysis(edm::ParameterSetDescription& desc);

private:
  void beginJob() override;
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;
  void DefineRandomCone(const std::vector<pat::PackedCandidate> pfCandidateColl);
  void RhoPhiModulation(const double object_phi, const double minFlowChi2Prob_, const double maxFlowChi2Prob_, const edm::Handle<std::vector<double>>  rhoFlowFitParams, double &rhoModulationFactor);
  void ConeAxisModulation(const edm::Handle<std::vector<double>> rhoFlowFitParams, const bool useModulatedRho_ );
  double getModulatedRhoFactor(const double phi, const edm::Handle<std::vector<double>> flowComponents);
  void coneSegmentCorrection(edm::Handle<std::vector<double>> rhoFlowFitParams);
  // ----------member data ---------------------------
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidateToken_;
  // edm::EDGetTokenT<reco::PFCandidateCollection>          pfCandidateLabel_;
  edm::EDGetTokenT<std::vector<double>>                  etaToken_;
  edm::EDGetTokenT<std::vector<double>>                  rhoToken_;
  // edm::EDGetTokenT<std::vector<double>>                  rhomToken_;
  edm::EDGetTokenT<std::vector<double>>                  rhoPhiModParamsToken_;

  bool    useModulatedRho_;
  double minFlowChi2Prob_;  /// flowFit chi2/ndof minimum compatability requirement
  double maxFlowChi2Prob_;  /// flowFit chi2/ndof maximum compatability requirement
  double rngConeRadius_;

  TTree *tree_;
  edm::Service<TFileService> fs_;

  struct RHO {
        std::vector<double> etaMin;
        std::vector<double> etaMax;
        std::vector<double> rho;
        std::vector<double> rhoPhiMod;
        double cone_eta;
        double cone_phi;
        double cone_pt_raw;
        double cone_pt_axisCorr;
        double cone_pt_bandsCorr;
  };
  RHO rhoObj_;
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
RandomConeAnalysis::RandomConeAnalysis(const edm::ParameterSet& iConfig)
  :pfCandidateToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSource"))),
  etaToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("etaMap"))),
  rhoToken_(consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("rho"))),

  useModulatedRho_(iConfig.getParameter<bool>("useModulatedRho")),
  minFlowChi2Prob_(iConfig.getParameter<double>("minFlowChi2Prob")),
  maxFlowChi2Prob_(iConfig.getParameter<double>("maxFlowChi2Prob")),
  rngConeRadius_(iConfig.getParameter<double>("rngConeRadius"))
{
  if(useModulatedRho_)
    rhoPhiModParamsToken_ = consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("rhoFlowFitParams"));
    
  //now do what ever initialization is needed
}

RandomConeAnalysis::~RandomConeAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void RandomConeAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  rhoObj_ = {};

  edm::Handle<std::vector<double>> etaRanges;
  edm::Handle<std::vector<double>> rho;
  edm::Handle<std::vector<double>> rhoFlowFitParams;
  auto const& pfCandidates = iEvent.get(pfCandidateToken_);

  iEvent.getByToken(etaToken_, etaRanges);
  iEvent.getByToken(rhoToken_, rho);
  if(useModulatedRho_)
    iEvent.getByToken(rhoPhiModParamsToken_, rhoFlowFitParams);
  // iEvent.getByToken(pfCandidateLabel_, pfCandidates);

  int neta = (int)etaRanges->size();
  for(int ieta = 0; ieta<(neta-1); ieta++) {
    rhoObj_.etaMin.push_back(etaRanges->at(ieta));
    rhoObj_.etaMax.push_back(etaRanges->at(ieta+1));
    rhoObj_.rho.push_back(rho->at(ieta));
  }

  //random cone generation, raw 
  DefineRandomCone(pfCandidates);
  //axis correction to cone pt
  ConeAxisModulation(rhoFlowFitParams, useModulatedRho_);
  //segment correction to cone pt
  coneSegmentCorrection(rhoFlowFitParams);
  tree_->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void RandomConeAnalysis::beginJob() {
  // please remove this method if not needed
  TString jetTagTitle = "HiPuRho Jet background analysis tree";
  tree_ = fs_->make<TTree>("t",jetTagTitle.Data());
  // tree_->Branch("etaMin",&(rhoObj_.etaMin));
  // tree_->Branch("etaMax",&(rhoObj_.etaMax));
  // tree_->Branch("rho",&(rhoObj_.rho));
  tree_->Branch("rcone_eta",&(rhoObj_.cone_eta));
  tree_->Branch("rcone_phi",&(rhoObj_.cone_phi));
  tree_->Branch("rcone_pt_raw",&(rhoObj_.cone_pt_raw));
  tree_->Branch("rcone_pt_axisCorr",&(rhoObj_.cone_pt_axisCorr));
  tree_->Branch("rcone_pt_bandsCorr",&(rhoObj_.cone_pt_bandsCorr));

}

// ------------ method called once each job just after ending the event loop  ------------
void RandomConeAnalysis::endJob() {
  // please remove this method if not needed
}


void RandomConeAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  edm::ParameterSetDescription desc;
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  desc.add<edm::InputTag>("etaMap", {"hiPuRhoProducer", "mapEtaEdges"});
  desc.add<edm::InputTag>("rho", {"hiPuRhoProducer", "mapToRho"});
  desc.add<double>("minFlowChi2Prob", 0.05);
  desc.add<double>("maxFlowChi2Prob", 0.95);
  desc.add<double>("rngConeRadius", 0.4);
  desc.add<bool>("useModulatedRho", true);
  desc.add<edm::InputTag>("rhoFlowFitParams", {"hiFJRhoFlowModulationProducer", "rhoFlowFitParams"});
  desc.add<edm::InputTag>("pfCandSource", edm::InputTag("packedPFCandidates"));
  // desc.add<edm::InputTag>("particleFlowInput", edm::InputTag("cleanedParticleFlow"));
  // descriptions.add("RandomConeAnalysis", desc);
  descriptions.addWithDefaultLabel(desc);
}


void RandomConeAnalysis::DefineRandomCone(const std::vector<pat::PackedCandidate> pfCandidateColl){
//define random cone from a flat distribution in eta [-2+x;2-x] (where x is the cone radius) and phi [0;2pi], extract raw p_t of PF candidates
    TRandom *r1 = new TRandom3(0);
    TRandom *r2 = new TRandom3(0);
    //jets don't follow uniform distribution in eta, but for now keep it this way
    rhoObj_.cone_phi = r1->Uniform(0,2*M_PI);
    rhoObj_.cone_eta = r2->Uniform(-2 + rngConeRadius_, 2 - rngConeRadius_);
    //pfCandidate(LorentzVector) does not have Delta_R function...
    for(unsigned icand=0;icand<pfCandidateColl.size(); icand++) {
        const pat::PackedCandidate& pfCandidate = pfCandidateColl.at(icand);
        float Delta_R = sqrt(pow(delta_phi(pfCandidate.phi(), rhoObj_.cone_phi),2) + pow(pfCandidate.eta()-rhoObj_.cone_eta,2));
        if( Delta_R < rngConeRadius_ ){
            rhoObj_.cone_pt_raw  += pfCandidate.pt();
        }
    }
}


void RandomConeAnalysis::RhoPhiModulation(double object_phi, double minFlowChi2Prob_, double maxFlowChi2Prob_, edm::Handle<std::vector<double> >  rhoFlowFitParams, double &rhoModulationFactor){
//take an object phi as input and calculate the modulation factor for that specific phi
  bool minProb = false;
  bool maxProb = false;
  if(rhoFlowFitParams->empty()) return;
  int chi2index = rhoFlowFitParams->size() - 3;
  double val = ROOT::Math::chisquared_cdf_c(rhoFlowFitParams->at(chi2index), rhoFlowFitParams->at(chi2index+1));
  minProb = val > minFlowChi2Prob_;
  maxProb = val < maxFlowChi2Prob_;
  if (minProb && maxProb)
    //gMRF defined in CSJetProducer
    rhoModulationFactor = getModulatedRhoFactor(object_phi, rhoFlowFitParams);
}


void RandomConeAnalysis::ConeAxisModulation(const edm::Handle<std::vector<double> > rhoFlowFitParams, const bool useModulatedRho_ ){
//take the map of rho=rho(eta) and find the rho of the random cone axis, also modulate by phi if set
  double rhoModulationFactor = 1.;

  if(useModulatedRho_)
    RhoPhiModulation(rhoObj_.cone_phi, minFlowChi2Prob_, maxFlowChi2Prob_, rhoFlowFitParams, rhoModulationFactor);
  // //find out which eta bin the random cone axis is located in
  // std::cout << "New modulation factor " << rhoModulationFactor << "for " << useModulatedRho_ << std::endl;
  size_t cone_eta_bin = 0;
  for( size_t i{0}; i < rhoObj_.etaMin.size(); ++i ){
    if(rhoObj_.etaMin.at(i) < rhoObj_.cone_eta && rhoObj_.etaMax.at(i) > rhoObj_.cone_eta){
      cone_eta_bin = i;
      break;
    }
  }
  rhoObj_.cone_pt_axisCorr = rhoObj_.cone_pt_raw - rngConeRadius_*rngConeRadius_*M_PI*rhoObj_.rho.at(cone_eta_bin)*rhoModulationFactor;
}


double RandomConeAnalysis::getModulatedRhoFactor(double phi, edm::Handle<std::vector<double>> flowComponents) {
  //the function is a protected member of CSJetProducer, not sure if it's better to make RandomConeAnalysis inherit from that class, or this is OK?
  // get the rho modulation as function of phi
  // flow modulation fit is done in RecoHI/HiJetAlgos/plugins/HiFJRhoFlowModulationProducer
  int nFlow = (flowComponents->size()-4) / 2;
  double modulationValue = 1;
  double firstFlowComponent = flowComponents->at(nFlow*2+3);
  for(int iFlow = 0; iFlow < nFlow; iFlow++){
    modulationValue += 2.0 * flowComponents->at(2*iFlow+1) * cos((iFlow+firstFlowComponent) * (phi - flowComponents->at(2*iFlow+2)) );
  }
  return modulationValue;
}

void RandomConeAnalysis::coneSegmentCorrection(edm::Handle<std::vector<double>> rhoFlowFitParams){
  //divide the random cone into different strips of eta and calculate the separate contributions from each strip to the total cone correction factor
  //find the two bins which contain the first/last segments of the cone, and its axis
  double rhoModulationFactor = 1.;
  if(useModulatedRho_)
    RhoPhiModulation(rhoObj_.cone_phi, minFlowChi2Prob_, maxFlowChi2Prob_, rhoFlowFitParams, rhoModulationFactor);
  size_t cone_min_idx = 0;
  size_t cone_max_idx = 0;
  size_t cone_axis_idx = 0;
  for(size_t i{0}; i<rhoObj_.etaMin.size(); ++i){
    //finding etamin and centre should continue the loop, finding the max should break out as we ascend in eta
    if(rhoObj_.etaMin.at(i) < rhoObj_.cone_eta-rngConeRadius_ && rhoObj_.etaMax.at(i) > rhoObj_.cone_eta-rngConeRadius_){
      cone_min_idx = i;
      continue;
    }
    if(rhoObj_.etaMin.at(i) < rhoObj_.cone_eta+rngConeRadius_ && rhoObj_.etaMax.at(i) > rhoObj_.cone_eta+rngConeRadius_){
      cone_max_idx = i;
      break;
    }
    if(rhoObj_.etaMin.at(i) < rhoObj_.cone_eta && rhoObj_.etaMax.at(i) > rhoObj_.cone_eta){
      cone_axis_idx = i;
      continue;
    }
  }
  //now loop over the relevant bands containing the cone and calculate the areas
  //area of first and last eta band of the cone are equal to the segment (sector - triangle) defined by the etaMax and etaMin lines respectively
  //eta bin containing cone axis divided into two about the centre
  //rest of bins defined as difference between segments defined by etaMax and etaMin
  std::vector<double> areas(rhoObj_.etaMin.size(),0);
  std::vector<double> area_times_rho(rhoObj_.etaMin.size(),0);
  for(size_t i{cone_min_idx}; i <= cone_max_idx; ++i){
    double height_min = abs(rhoObj_.cone_eta - rhoObj_.etaMin.at(i));
    double height_max = abs(rhoObj_.cone_eta - rhoObj_.etaMax.at(i));
    double theta_min = 2*acos(height_min/rngConeRadius_);
    double theta_max = 2*acos(height_max/rngConeRadius_);
    // std::cout << height_min << " " << height_max << " " << theta_min << " " << theta_max << " variables" << std::endl;
    if(i==cone_min_idx){
      areas.at(i) = 0.5*rngConeRadius_*rngConeRadius_*(theta_max-sin(theta_max));
    }
    else if(i==cone_max_idx){
      areas.at(i) = 0.5*rngConeRadius_*rngConeRadius_*(theta_min-sin(theta_min));
    }
    else if(i==cone_axis_idx){
      areas.at(i) = 0.5*rngConeRadius_*rngConeRadius_*(2*M_PI-theta_min+sin(theta_min)-theta_max+sin(theta_max));
    }
    else{
      areas.at(i) = abs(0.5*rngConeRadius_*rngConeRadius_*(theta_max-theta_min-sin(theta_max)+sin(theta_min)));
    }
    // std::cout << areas.at(i) << " area at position " << i << std::endl;
    area_times_rho.at(i) = areas.at(i)*rhoObj_.rho.at(i)*rhoModulationFactor;
  }
  double total_area = 0;
  double total_area_times_rho = 0;
  for(size_t i{0}; i<areas.size(); ++i){
    total_area += areas.at(i);
    total_area_times_rho += area_times_rho.at(i);
  }
  // std::cout << "deviation from the total area is " << 0.4*0.4*M_PI-total_area << std::endl;
  rhoObj_.cone_pt_bandsCorr = rhoObj_.cone_pt_raw - total_area_times_rho;

}

//define this as a plug-in
DEFINE_FWK_MODULE(RandomConeAnalysis);
