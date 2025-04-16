/*
  Based on the jet response analyzer
  Modified by Matt Nguyen, November 2010
*/

#include "HeavyIonsAnalysis/JetAnalysis/interface/HiInclusiveJetAnalyzer.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "AnalysisDataFormats/TrackInfo/interface/TrackToGenParticleMap.h"

using namespace std;
using namespace edm;
using namespace reco;

HiInclusiveJetAnalyzer::HiInclusiveJetAnalyzer(const edm::ParameterSet& iConfig) {
  doMatch_ = iConfig.getUntrackedParameter<bool>("matchJets", false);
  jetTag_ = consumes<pat::JetCollection>(iConfig.getParameter<InputTag>("jetTag"));
  originalCSTag_ = consumes<pat::JetCollection>(iConfig.getParameter<InputTag>("originalCSTag"));
  caloJetTag_ = consumes<reco::CaloJetCollection>(iConfig.getParameter<InputTag>("caloJetTag"));
  matchTag_ = consumes<pat::JetCollection>(iConfig.getUntrackedParameter<InputTag>("matchTag"));

  runSubstructure =  iConfig.getUntrackedParameter<bool>("runSubstructure", false);
  doChargedConstOnly_ =  iConfig.getUntrackedParameter<bool>("doChargedConstOnly", true);
  doPFjetID =  iConfig.getUntrackedParameter<bool>("doPFjetID", false);

  useQuality_ = iConfig.getUntrackedParameter<bool>("useQuality", true);
  trackQuality_ = iConfig.getUntrackedParameter<string>("trackQuality", "highPurity");

  jetName_ = iConfig.getUntrackedParameter<string>("jetName");
  doGenTaus_ = iConfig.getUntrackedParameter<bool>("doGenTaus", false);
  doGenSym_ = iConfig.getUntrackedParameter<bool>("doGenSym", false);
  doSubJets_ = iConfig.getUntrackedParameter<bool>("doSubJets", false);
  doJetConstituents_ = iConfig.getUntrackedParameter<bool>("doJetConstituents", false);
  doCaloJets_ = iConfig.getUntrackedParameter<bool>("doCaloJets", true);

  matchOriginal_ = iConfig.getUntrackedParameter<bool>("matchOriginal", true);
   
  doGenSubJets_ = iConfig.getUntrackedParameter<bool>("doGenSubJets", false);
  if (doGenSubJets_)
    subjetGenTag_ = consumes<reco::JetView>(iConfig.getUntrackedParameter<InputTag>("subjetGenTag"));

  //reWTA reclustering
  doWTARecluster_ = iConfig.getUntrackedParameter<bool>("doWTARecluster", false);

  if (doGenTaus_) {
    tokenGenTau1_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("genTau1"));
    tokenGenTau2_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("genTau2"));
    tokenGenTau3_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("genTau3"));
  }

  if (doGenSym_) {
    tokenGenSym_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("genSym"));
    tokenGenDroppedBranches_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("genDroppedBranches"));
  }

  isMC_ = iConfig.getUntrackedParameter<bool>("isMC", false);
  useHepMC_ = iConfig.getUntrackedParameter<bool>("useHepMC", false);
  fillGenJets_ = iConfig.getUntrackedParameter<bool>("fillGenJets", false);
  
  doHiJetID_ = iConfig.getUntrackedParameter<bool>("doHiJetID", false);
  doStandardJetID_ = iConfig.getUntrackedParameter<bool>("doStandardJetID", false);

  rParam = iConfig.getParameter<double>("rParam");
  hardPtMin_ = iConfig.getUntrackedParameter<double>("hardPtMin", 4);
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
  jetAbsEtaMax_ = iConfig.getUntrackedParameter<double>("jetAbsEtaMax", 5.1);

  if (isMC_) {
    genjetTag_ = consumes<edm::View<reco::GenJet>>(iConfig.getParameter<InputTag>("genjetTag"));
    if (useHepMC_)
      eventInfoTag_ = consumes<HepMCProduct>(iConfig.getParameter<InputTag>("eventInfoTag"));
    eventGenInfoTag_ = consumes<GenEventInfoProduct>(iConfig.getParameter<InputTag>("eventInfoTag"));
    jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( iConfig.getParameter<edm::InputTag>("jetFlavourInfos") );
  }
  useRawPt_ = iConfig.getUntrackedParameter<bool>("useRawPt", true);
  doLegacyBtagging_ = iConfig.getUntrackedParameter<bool>("doLegacyBtagging", true);
  doCandidateBtagging_ = iConfig.getUntrackedParameter<bool>("doCandidateBtagging", true);
  useNewBtaggers_ = iConfig.getUntrackedParameter<bool>("useNewBtaggers", false);

  pfCandidateLabel_ =
      consumes<edm::View<pat::PackedCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateLabel"));

  if (isMC_)
    genParticleSrc_ =
        consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"));

  if (doLegacyBtagging_) {
    trackCHEBJetTags_ = "trackCountingHighEffBJetTags";
    trackCHPBJetTags_ = "trackCountingHighPurBJetTags";
    jetPBJetTags_ = "jetProbabilityBJetTags";
    jetBPBJetTags_ = "jetBProbabilityBJetTags";
    simpleSVHighEffBJetTags_ = "simpleSecondaryVertexHighEffBJetTags";
    simpleSVHighPurBJetTags_ = "simpleSecondaryVertexHighPurBJetTags";
    combinedSVV2BJetTags_ = "combinedSecondaryVertexV2BJetTags";
  }
  if (doCandidateBtagging_ && useNewBtaggers_) {
    for (const std::string label : {"pfJetProbabilityBJetTag", "pfDeepCSVJetTags", "pfDeepFlavourJetTags", "pfParticleTransformerAK4JetTags", "pfUnifiedParticleTransformerAK4JetTags"}) {
      const auto& tag = iConfig.getUntrackedParameter<string>(label, "");
      if (tag.empty())
        continue;
      else if (label == "pfJetProbabilityBJetTag")
        jetTaggers_["pfJP"].emplace("probb", consumes<JetTagCollection>(tag));
      else if (label == "pfDeepCSVJetTags")
        for (const auto& cat : {"probb", "probbb"})
           jetTaggers_["deepCSV"].emplace(cat, consumes<JetTagCollection>(tag+":"+cat));
      else if (label == "pfDeepFlavourJetTags")
        for (const auto& cat : {"probb", "probbb", "problepb"})
          jetTaggers_["deepFlavour"].emplace(cat, consumes<JetTagCollection>(tag+":"+cat));
      else if (label == "pfParticleTransformerAK4JetTags")
        for (const auto& cat : {"probb", "probbb", "problepb"})
          jetTaggers_["particleTransformer"].emplace(cat, consumes<JetTagCollection>(tag+":"+cat));
      else if (label == "pfUnifiedParticleTransformerAK4JetTags")
        for (const auto& cat : {"probb", "probbb", "problepb"})
				/*	      for (const auto& cat : {"probb", "probbb", "problepb", "probc", "probg", "probu", "probd", "probs",
                                "probtaup1h0p", "probtaup1h1p", "probtaup1h2p", "probtaup3h0p", "probtaup3h1p", "probtaum1h0p", "probtaum1h1p", "probtaum1h2p", "probtaum3h0p", "probtaum3h1p",
                                "probele", "probmu", "ptcorr", "ptnu"}) */
          jetTaggers_["unifiedParticleTransformer"].emplace(cat, consumes<JetTagCollection>(tag+":"+cat));
    }
  }
  else if (doCandidateBtagging_) {
    deepCSVJetTags_ = jetName_ + "pfDeepCSVJetTags:probb";
    pfJPJetTags_ = jetName_ + "pfJetProbabilityBJetTags";
  }
  doSubEvent_ = false;

  if (isMC_) {
    genPtMin_ = iConfig.getUntrackedParameter<double>("genPtMin", 10);
    doSubEvent_ = iConfig.getUntrackedParameter<bool>("doSubEvent", false);
  }

    ptCut = iConfig.getUntrackedParameter<double>("ptCut",1.);
  trkInefRate_ = iConfig.getUntrackedParameter<double>("trkInefRate",0.);

  doSvtx_ = iConfig.getUntrackedParameter<bool>("doSvtx",false);
  if (doSvtx_) {
    svTagInfoLabel_ = iConfig.getUntrackedParameter<std::string>("svTagInfoLabel");
  }

  doTracks_ = iConfig.getUntrackedParameter<bool>("doTracks", false);

  if (doTracks_) {
    trkPtCut_ = iConfig.getUntrackedParameter<double>("trkPtCut", 1.);
    ipTagInfoLabel_ = iConfig.getUntrackedParameter<std::string>("ipTagInfoLabel");
    if (isMC_) {
      trackToGenParticleMapToken_ = consumes<reco::TrackToGenParticleMap>(iConfig.getUntrackedParameter<edm::InputTag>("trackToGenParticleMap", edm::InputTag("TrackToGenParticleMapProducer", "trackToGenParticleMap")));
    }
  }
  primaryVerticesToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices", edm::InputTag("offlineSlimmedPrimaryVertices")));
  
}

HiInclusiveJetAnalyzer::~HiInclusiveJetAnalyzer() {}

void HiInclusiveJetAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& es) {}
void HiInclusiveJetAnalyzer::endRun(const edm::Run& run, const edm::EventSetup& es) {}

void HiInclusiveJetAnalyzer::beginJob() {
  string jetTagTitle = jetTagLabel_.label() + " Jet Analysis Tree";
  t = fs1->make<TTree>("t", jetTagTitle.c_str());

  t->Branch("run", &jets_.run, "run/I");
  t->Branch("evt", &jets_.evt, "evt/I");
  t->Branch("lumi", &jets_.lumi, "lumi/I");
  t->Branch("nref", &jets_.nref, "nref/I");
  t->Branch("rawpt", jets_.rawpt, "rawpt[nref]/F");
  t->Branch("jtpt", jets_.jtpt, "jtpt[nref]/F");
  t->Branch("jtptCS", jets_.jtptCS, "jtptCS[nref]/F");
  t->Branch("jtCSdr", jets_.jtCSdr, "jtCSdr[nref]/F");
  t->Branch("jteta", jets_.jteta, "jteta[nref]/F");
  //  t->Branch("jty", jets_.jty, "jty[nref]/F");
  t->Branch("jtphi", jets_.jtphi, "jtphi[nref]/F");
  /*  t->Branch("jtpu", jets_.jtpu, "jtpu[nref]/F");
  t->Branch("jtm", jets_.jtm, "jtm[nref]/F");
  t->Branch("jtarea", jets_.jtarea, "jtarea[nref]/F"); */

  //  t->Branch("nCSjets", &jets_.nCSjets, "nCSjets/I");

  t->Branch("massHF", jets_.massHF, "massHF[nref]/F");
  if (isMC_)  t->Branch("massHFgen", jets_.massHFgen, "massHFgen[nref]/F");
  
  if(doCaloJets_){
    t->Branch("ncalo", &jets_.ncalo, "ncalo/I");
    t->Branch("calopt", jets_.calopt, "calopt[ncalo]/F");
    t->Branch("caloeta", jets_.caloeta, "caloeta[ncalo]/F");
    t->Branch("calophi", jets_.calophi, "calophi[ncalo]/F");
  }


  if (doPFjetID) {
    t->Branch("jtPfCHF", jets_.jtPfCHF, "jtPfCHF[nref]/F");
    t->Branch("jtPfNHF", jets_.jtPfNHF, "jtPfNHF[nref]/F");
    t->Branch("jtPfCEF", jets_.jtPfCEF, "jtPfCEF[nref]/F");
    t->Branch("jtPfNEF", jets_.jtPfNEF, "jtPfNEF[nref]/F");
    t->Branch("jtPfMUF", jets_.jtPfMUF, "jtPfMUF[nref]/F");
    
    t->Branch("jtPfCHM", jets_.jtPfCHM, "jtPfCHM[nref]/I");
    t->Branch("jtPfNHM", jets_.jtPfNHM, "jtPfNHM[nref]/I");
    t->Branch("jtPfCEM", jets_.jtPfCEM, "jtPfCEM[nref]/I");
    t->Branch("jtPfNEM", jets_.jtPfNEM, "jtPfNEM[nref]/I");
    t->Branch("jtPfMUM", jets_.jtPfMUM, "jtPfMUM[nref]/I");
  }
  
  if (doSubJets_) {
    t->Branch("jtSubJetPt", &jets_.jtSubJetPt);
    t->Branch("jtSubJetEta", &jets_.jtSubJetEta);
    t->Branch("jtSubJetPhi", &jets_.jtSubJetPhi);
    t->Branch("jtSubJetM", &jets_.jtSubJetM);
    t->Branch("jtsym", jets_.jtsym, "jtsym[nref]/F");
    t->Branch("jtdroppedBranches", jets_.jtdroppedBranches, "jtdroppedBranches[nref]/I");
  }

  if (doJetConstituents_) {
    t->Branch("jtConstituentsId", &jets_.jtConstituentsId);
    t->Branch("jtConstituentsE", &jets_.jtConstituentsE);
    t->Branch("jtConstituentsPt", &jets_.jtConstituentsPt);
    t->Branch("jtConstituentsEta", &jets_.jtConstituentsEta);
    t->Branch("jtConstituentsPhi", &jets_.jtConstituentsPhi);
    t->Branch("jtConstituentsM", &jets_.jtConstituentsM);
    t->Branch("jtSDConstituentsId", &jets_.jtSDConstituentsId);
    t->Branch("jtSDConstituentsE", &jets_.jtSDConstituentsE);
    t->Branch("jtSDConstituentsPt", &jets_.jtSDConstituentsPt);
    t->Branch("jtSDConstituentsEta", &jets_.jtSDConstituentsEta);
    t->Branch("jtSDConstituentsPhi", &jets_.jtSDConstituentsPhi);
    t->Branch("jtSDConstituentsM", &jets_.jtSDConstituentsM);
  }
  // jet ID information, jet composition
  if (doHiJetID_) {
    t->Branch("trackMax", jets_.trackMax, "trackMax[nref]/F");
    t->Branch("trackSum", jets_.trackSum, "trackSum[nref]/F");
    t->Branch("trackN", jets_.trackN, "trackN[nref]/I");
    t->Branch("trackHardSum", jets_.trackHardSum, "trackHardSum[nref]/F");
    t->Branch("trackHardN", jets_.trackHardN, "trackHardN[nref]/I");

    t->Branch("chargedMax", jets_.chargedMax, "chargedMax[nref]/F");
    t->Branch("chargedSum", jets_.chargedSum, "chargedSum[nref]/F");
    t->Branch("chargedN", jets_.chargedN, "chargedN[nref]/I");
    t->Branch("chargedHardSum", jets_.chargedHardSum, "chargedHardSum[nref]/F");
    t->Branch("chargedHardN", jets_.chargedHardN, "chargedHardN[nref]/I");

    t->Branch("photonMax", jets_.photonMax, "photonMax[nref]/F");
    t->Branch("photonSum", jets_.photonSum, "photonSum[nref]/F");
    t->Branch("photonN", jets_.photonN, "photonN[nref]/I");
    t->Branch("photonHardSum", jets_.photonHardSum, "photonHardSum[nref]/F");
    t->Branch("photonHardN", jets_.photonHardN, "photonHardN[nref]/I");

    t->Branch("neutralMax", jets_.neutralMax, "neutralMax[nref]/F");
    t->Branch("neutralSum", jets_.neutralSum, "neutralSum[nref]/F");
    t->Branch("neutralN", jets_.neutralN, "neutralN[nref]/I");

    t->Branch("eMax", jets_.eMax, "eMax[nref]/F");
    t->Branch("eSum", jets_.eSum, "eSum[nref]/F");
    t->Branch("eN", jets_.eN, "eN[nref]/I");

    t->Branch("muMax", jets_.muMax, "muMax[nref]/F");
    t->Branch("muSum", jets_.muSum, "muSum[nref]/F");
    t->Branch("muN", jets_.muN, "muN[nref]/I");
  }

  if (doStandardJetID_) {
    t->Branch("fHPD", jets_.fHPD, "fHPD[nref]/F");
    t->Branch("fRBX", jets_.fRBX, "fRBX[nref]/F");
    t->Branch("n90", jets_.n90, "n90[nref]/I");
    t->Branch("fSubDet1", jets_.fSubDet1, "fSubDet1[nref]/F");
    t->Branch("fSubDet2", jets_.fSubDet2, "fSubDet2[nref]/F");
    t->Branch("fSubDet3", jets_.fSubDet3, "fSubDet3[nref]/F");
    t->Branch("fSubDet4", jets_.fSubDet4, "fSubDet4[nref]/F");
    t->Branch("restrictedEMF", jets_.restrictedEMF, "restrictedEMF[nref]/F");
    t->Branch("nHCAL", jets_.nHCAL, "nHCAL[nref]/I");
    t->Branch("nECAL", jets_.nECAL, "nECAL[nref]/I");
    t->Branch("apprHPD", jets_.apprHPD, "apprHPD[nref]/F");
    t->Branch("apprRBX", jets_.apprRBX, "apprRBX[nref]/F");
    t->Branch("n2RPC", jets_.n2RPC, "n2RPC[nref]/I");
    t->Branch("n3RPC", jets_.n3RPC, "n3RPC[nref]/I");
    t->Branch("nRPC", jets_.nRPC, "nRPC[nref]/I");

    t->Branch("fEB", jets_.fEB, "fEB[nref]/F");
    t->Branch("fEE", jets_.fEE, "fEE[nref]/F");
    t->Branch("fHB", jets_.fHB, "fHB[nref]/F");
    t->Branch("fHE", jets_.fHE, "fHE[nref]/F");
    t->Branch("fHO", jets_.fHO, "fHO[nref]/F");
    t->Branch("fLong", jets_.fLong, "fLong[nref]/F");
    t->Branch("fShort", jets_.fShort, "fShort[nref]/F");
    t->Branch("fLS", jets_.fLS, "fLS[nref]/F");
    t->Branch("fHFOOT", jets_.fHFOOT, "fHFOOT[nref]/F");
  }

  // Jet ID
  if (doMatch_) {
    t->Branch("mjtPt", jets_.mjtPt, "mjtPt[nref]/F");
    t->Branch("mjtRawPt", jets_.mjtRawPt, "mjtRawPt[nref]/F");
    //    t->Branch("mjtPu", jets_.mjtPu, "mjtPu[nref]/F");
    t->Branch("mjtR", jets_.mjtR, "mjtR[nref]/F");
    if (isMC_) {
      t->Branch("mjtHadronFlavor", jets_.mjtHadronFlavor, "mjtHadronFlavor[nref]/I");
      t->Branch("mjtPartonFlavor", jets_.mjtPartonFlavor, "mjtPartonFlavor[nref]/I");
      t->Branch("mjtNbHad", jets_.mjtNbHad, "mjtNbHad[nref]/I");
      t->Branch("mjtNcHad", jets_.mjtNcHad, "mjtNcHad[nref]/I");
      t->Branch("mjtNbPar", jets_.mjtNbPar, "mjtNbPar[nref]/I");
      t->Branch("mjtNcPar", jets_.mjtNcPar, "mjtNcPar[nref]/I");
}
  }

  // b-jet discriminators
  if (doLegacyBtagging_) {
    t->Branch("discr_ssvHighEff", jets_.discr_ssvHighEff, "discr_ssvHighEff[nref]/F");
    t->Branch("discr_ssvHighPur", jets_.discr_ssvHighPur, "discr_ssvHighPur[nref]/F");
    t->Branch("discr_csvV2", jets_.discr_csvV2, "discr_csvV2[nref]/F");
    t->Branch("discr_muByIp3", jets_.discr_muByIp3, "discr_muByIp3[nref]/F");
    t->Branch("discr_muByPt", jets_.discr_muByPt, "discr_muByPt[nref]/F");
    t->Branch("discr_prob", jets_.discr_prob, "discr_prob[nref]/F");
    t->Branch("discr_probb", jets_.discr_probb, "discr_probb[nref]/F");
    t->Branch("discr_tcHighEff", jets_.discr_tcHighEff, "discr_tcHighEff[nref]/F");
    t->Branch("discr_tcHighPur", jets_.discr_tcHighPur, "discr_tcHighPur[nref]/F");

    t->Branch("mue", jets_.mue, "mue[nref]/F");
    t->Branch("mupt", jets_.mupt, "mupt[nref]/F");
    t->Branch("mueta", jets_.mueta, "mueta[nref]/F");
    t->Branch("muphi", jets_.muphi, "muphi[nref]/F");
    t->Branch("mudr", jets_.mudr, "mudr[nref]/F");
    t->Branch("muptrel", jets_.muptrel, "muptrel[nref]/F");
    t->Branch("muchg", jets_.muchg, "muchg[nref]/I");
  }
  if (doCandidateBtagging_ && useNewBtaggers_) {
    for (const auto& tg : jetTaggers_) {
      auto& discr = jets_discr_[tg.first];
      t->Branch(("discr_"+tg.first).c_str(), discr["b"].data(), ("discr_"+tg.first+"[nref]/F").c_str());
      if (tg.first == "unifiedParticleTransformer") {
	//        t->Branch(("discr_"+tg.first+"_probtau").c_str(), discr["probtau"].data(), ("discr_"+tg.first+"_probtau[nref]/F").c_str());
        for (const auto& c : tg.second)
          if (c.first.rfind("probtau",0)!=0)
            t->Branch(("discr_"+tg.first+"_"+c.first).c_str(), discr[c.first].data(), ("discr_"+tg.first+"_"+c.first+"[nref]/F").c_str());
      }
    }
  }
  else if (doCandidateBtagging_) {
    t->Branch("discr_deepCSV", jets_.discr_deepCSV, "discr_deepCSV[nref]/F");
    t->Branch("discr_pfJP", jets_.discr_pfJP, "discr_pfJP[nref]/F");
  }

  if (runSubstructure) {
    t->Branch("jt_z_SD",jets_.jt_z_SD,"jt_z_SD[nref]/F");
    t->Branch("jt_rg_SD",jets_.jt_rg_SD,"jt_rg_SD[nref]/F");
    t->Branch("jt_ktg_SD",jets_.jt_ktg_SD,"jt_ktg_SD[nref]/F");
    t->Branch("jt_split_SD",jets_.jt_split_SD,"jt_split_SD[nref]/I");
    t->Branch("jt_hasHF_SD",jets_.jt_hasHF_SD,"jt_hasHF_SD[nref]/O");

    t->Branch("jt_z_latekt",jets_.jt_z_latekt,"jt_z_latekt[nref]/F");
    t->Branch("jt_rg_latekt",jets_.jt_rg_latekt,"jt_rg_latekt[nref]/F");
    t->Branch("jt_ktg_latekt",jets_.jt_ktg_latekt,"jt_ktg_latekt[nref]/F");
    t->Branch("jt_split_latekt",jets_.jt_split_latekt,"jt_split_latekt[nref]/I");
    t->Branch("jt_hasHF_latekt",jets_.jt_hasHF_latekt,"jt_hasHF_latekt[nref]/O");

    if (isMC_) {
      t->Branch("ref_z_SD",jets_.ref_z_SD,"ref_z_SD[nref]/F");
      t->Branch("ref_rg_SD",jets_.ref_rg_SD,"ref_rg_SD[nref]/F");
      t->Branch("ref_ktg_SD",jets_.ref_ktg_SD,"ref_ktg_SD[nref]/F");
      t->Branch("ref_split_SD",jets_.ref_split_SD,"ref_split_SD[nref]/I");
      t->Branch("ref_hasHF_SD",jets_.jt_hasHF_SD,"jt_hasHF_SD[nref]/O");
      
      t->Branch("ref_z_latekt",jets_.ref_z_latekt,"ref_z_latekt[nref]/F");
      t->Branch("ref_rg_latekt",jets_.ref_rg_latekt,"ref_rg_latekt[nref]/F");
      t->Branch("ref_ktg_latekt",jets_.ref_ktg_latekt,"ref_ktg_latekt[nref]/F");
      t->Branch("ref_split_latekt",jets_.ref_split_latekt,"ref_split_latekt[nref]/I");
      t->Branch("ref_hasHF_latekt",jets_.jt_hasHF_latekt,"jt_hasHF_latekt[nref]/O");

      t->Branch("jt_isClosestToTruth_latekt", jets_.jt_isClosestToTruth_latekt,"jt_isClosestToTruth_latekt[nref]/O");
      t->Branch("ref_isClosestToReco_latekt", jets_.ref_isClosestToReco_latekt,"ref_isClosestToReco_latekt[nref]/O");
      t->Branch("jt_ref_dR_latekt", jets_.jt_ref_dR_latekt,"jt_ref_dR_latekt[nref]/F");
    
      t->Branch("jt_isClosestToTruth_SD", jets_.jt_isClosestToTruth_SD,"jt_isClosestToTruth_SD[nref]/O");
      t->Branch("ref_isClosestToReco_SD", jets_.ref_isClosestToReco_SD,"ref_isClosestToReco_SD[nref]/O");
      t->Branch("jt_ref_dR_SD", jets_.jt_ref_dR_SD,"jt_ref_dR_SD[nref]/F");

    }
  }
  
  if (isMC_) {
    if (useHepMC_) {
      t->Branch("beamId1", &jets_.beamId1, "beamId1/I");
      t->Branch("beamId2", &jets_.beamId2, "beamId2/I");
    }

    //    t->Branch("pthat", &jets_.pthat, "pthat/F");

    // Only matched gen jets
    t->Branch("refpt", jets_.refpt, "refpt[nref]/F");
    t->Branch("refeta", jets_.refeta, "refeta[nref]/F");
    //    t->Branch("refy", jets_.refy, "refy[nref]/F");
    t->Branch("refphi", jets_.refphi, "refphi[nref]/F");
    //    t->Branch("refm", jets_.refm, "refm[nref]/F");
    //    t->Branch("refarea", jets_.refarea, "refarea[nref]/F");

    if (doGenTaus_) {
      t->Branch("reftau1", jets_.reftau1, "reftau1[nref]/F");
      t->Branch("reftau2", jets_.reftau2, "reftau2[nref]/F");
      t->Branch("reftau3", jets_.reftau3, "reftau3[nref]/F");
    }
    t->Branch("refdphijt", jets_.refdphijt, "refdphijt[nref]/F");
    t->Branch("refdrjt", jets_.refdrjt, "refdrjt[nref]/F");
    // matched parton
    //    t->Branch("refparton_pt", jets_.refparton_pt, "refparton_pt[nref]/F");
    //    t->Branch("refparton_flavor", jets_.refparton_flavor, "refparton_flavor[nref]/I");
    t->Branch("refparton_flavorForB", jets_.refparton_flavorForB, "refparton_flavorForB[nref]/I");

    if (doGenSubJets_) {
      t->Branch("refptG", jets_.refptG, "refptG[nref]/F");
      t->Branch("refetaG", jets_.refetaG, "refetaG[nref]/F");
      t->Branch("refphiG", jets_.refphiG, "refphiG[nref]/F");
      t->Branch("refmG", jets_.refmG, "refmG[nref]/F");
      t->Branch("refSubJetPt", &jets_.refSubJetPt);
      t->Branch("refSubJetEta", &jets_.refSubJetEta);
      t->Branch("refSubJetPhi", &jets_.refSubJetPhi);
      t->Branch("refSubJetM", &jets_.refSubJetM);
      t->Branch("refsym", jets_.refsym, "refsym[nref]/F");
      t->Branch("refdroppedBranches", jets_.refdroppedBranches, "refdroppedBranches[nref]/I");
    }

    if (doJetConstituents_) {
      t->Branch("refConstituentsId", &jets_.refConstituentsId);
      t->Branch("refConstituentsE", &jets_.refConstituentsE);
      t->Branch("refConstituentsPt", &jets_.refConstituentsPt);
      t->Branch("refConstituentsEta", &jets_.refConstituentsEta);
      t->Branch("refConstituentsPhi", &jets_.refConstituentsPhi);
      t->Branch("refConstituentsM", &jets_.refConstituentsM);
      t->Branch("refSDConstituentsId", &jets_.refSDConstituentsId);
      t->Branch("refSDConstituentsE", &jets_.refSDConstituentsE);
      t->Branch("refSDConstituentsPt", &jets_.refSDConstituentsPt);
      t->Branch("refSDConstituentsEta", &jets_.refSDConstituentsEta);
      t->Branch("refSDConstituentsPhi", &jets_.refSDConstituentsPhi);
      t->Branch("refSDConstituentsM", &jets_.refSDConstituentsM);
    }

    /*    t->Branch("genChargedSum", jets_.genChargedSum, "genChargedSum[nref]/F");
    t->Branch("genHardSum", jets_.genHardSum, "genHardSum[nref]/F");
    t->Branch("signalChargedSum", jets_.signalChargedSum, "signalChargedSum[nref]/F");
    t->Branch("signalHardSum", jets_.signalHardSum, "signalHardSum[nref]/F"); */

    if (doSubEvent_) {
      t->Branch("subid", jets_.subid, "subid[nref]/I");
    }

    if (fillGenJets_) {
      // For all gen jets, matched or unmatched
      t->Branch("ngen", &jets_.ngen, "ngen/I");
      t->Branch("genmatchindex", jets_.genmatchindex, "genmatchindex[ngen]/I");
      t->Branch("genpt", jets_.genpt, "genpt[ngen]/F");
      t->Branch("geneta", jets_.geneta, "geneta[ngen]/F");
      //      t->Branch("geny", jets_.geny, "geny[ngen]/F");
      if (doGenTaus_) {
        t->Branch("gentau1", jets_.gentau1, "gentau1[ngen]/F");
        t->Branch("gentau2", jets_.gentau2, "gentau2[ngen]/F");
        t->Branch("gentau3", jets_.gentau3, "gentau3[ngen]/F");
      }
      t->Branch("genphi", jets_.genphi, "genphi[ngen]/F");
      // t->Branch("genm", jets_.genm, "genm[ngen]/F");
      t->Branch("gendphijt", jets_.gendphijt, "gendphijt[ngen]/F");
      t->Branch("gendrjt", jets_.gendrjt, "gendrjt[ngen]/F");


      if (doGenSubJets_) {
        t->Branch("genptG", jets_.genptG, "genptG[ngen]/F");
        t->Branch("genetaG", jets_.genetaG, "genetaG[ngen]/F");
        t->Branch("genphiG", jets_.genphiG, "genphiG[ngen]/F");
        t->Branch("genmG", jets_.genmG, "genmG[ngen]/F");
        t->Branch("genSubJetPt", &jets_.genSubJetPt);
        t->Branch("genSubJetEta", &jets_.genSubJetEta);
        t->Branch("genSubJetPhi", &jets_.genSubJetPhi);
        t->Branch("genSubJetM", &jets_.genSubJetM);
        t->Branch("gensym", jets_.gensym, "gensym[ngen]/F");
        t->Branch("gendroppedBranches", jets_.gendroppedBranches, "gendroppedBranches[ngen]/I");
      }

      if (doJetConstituents_) {
        t->Branch("genConstituentsId", &jets_.genConstituentsId);
        t->Branch("genConstituentsE", &jets_.genConstituentsE);
        t->Branch("genConstituentsPt", &jets_.genConstituentsPt);
        t->Branch("genConstituentsEta", &jets_.genConstituentsEta);
        t->Branch("genConstituentsPhi", &jets_.genConstituentsPhi);
        t->Branch("genConstituentsM", &jets_.genConstituentsM);
        t->Branch("genSDConstituentsId", &jets_.genSDConstituentsId);
        t->Branch("genSDConstituentsE", &jets_.genSDConstituentsE);
        t->Branch("genSDConstituentsPt", &jets_.genSDConstituentsPt);
        t->Branch("genSDConstituentsEta", &jets_.genSDConstituentsEta);
        t->Branch("genSDConstituentsPhi", &jets_.genSDConstituentsPhi);
        t->Branch("genSDConstituentsM", &jets_.genSDConstituentsM);
      }

      if (doSubEvent_) {
        t->Branch("gensubid", jets_.gensubid, "gensubid[ngen]/I");
      }
    }
  }

  if (doTracks_) {
    t->Branch("jtNtrk", jets_.jtNtrk, "jtNtrk[nref]/I");
    t->Branch("ntrk", &jets_.ntrk, "ntrk/I");
    t->Branch("trkJetId", jets_.trkJetId, "trkJetId[ntrk]/I");
    t->Branch("trkSvtxId", jets_.trkSvtxId, "trkSvtxId[ntrk]/I");
    t->Branch("trkPt", jets_.trkPt, "trkPt[ntrk]/F");
    t->Branch("trkEta", jets_.trkEta, "trkEta[ntrk]/F");
    t->Branch("trkPhi", jets_.trkPhi, "trkPhi[ntrk]/F");
    t->Branch("trkIp3d", jets_.trkIp3d, "trkIp3d[ntrk]/F");
    t->Branch("trkIp3dSig", jets_.trkIp3dSig, "trkIp3dSig[ntrk]/F");
    t->Branch("trkIp2d", jets_.trkIp2d, "trkIp2d[ntrk]/F");
    t->Branch("trkIp2dSig", jets_.trkIp2dSig, "trkIp2dSig[ntrk]/F");
    t->Branch("trkDistToAxisSig", jets_.trkDistToAxisSig, "trkDistToAxisSig[ntrk]/F");
    t->Branch("trkDistToAxis", jets_.trkDistToAxis, "trkDistToAxis[ntrk]/F");
    t->Branch("trkIpProb3d", jets_.trkIpProb3d, "trkIpProb3d[ntrk]/F");
    t->Branch("trkIpProb2d", jets_.trkIpProb2d, "trkIpProb2d[ntrk]/F");
    t->Branch("trkDz", jets_.trkDz, "trkDz[ntrk]/F");
    t->Branch("trkPdgId", jets_.trkPdgId, "trkPdgId[ntrk]/I");
    t->Branch("trkMatchSta", jets_.trkMatchSta, "trkMatchSta[ntrk]/I");

    t->Branch("jtptCh", jets_.jtptCh, "jtptCh[nref]/F");
    if (isMC_) {
      t->Branch("refptCh", jets_.refptCh, "refptCh[nref]/F");
      t->Branch("refNtrk", jets_.refNtrk, "refNtrk[nref]/I");
    }
  }

  if (doSvtx_) {
    t->Branch("jtNsvtx", jets_.jtNsvtx, "jtNsvtx[nref]/I");
    t->Branch("nsvtx", &jets_.nsvtx, "nsvtx/I");
    t->Branch("svtxJetId", jets_.svtxJetId, "svtxJetId[nsvtx]/I");
    t->Branch("svtxNtrk", jets_.svtxNtrk, "svtxNtrk[nsvtx]/I");
    t->Branch("svtxdl", jets_.svtxdl, "svtxdl[nsvtx]/F");
    t->Branch("svtxdls", jets_.svtxdls, "svtxdls[nsvtx]/F");
    t->Branch("svtxdl2d", jets_.svtxdl2d, "svtxdl2d[nsvtx]/F");
    t->Branch("svtxdls2d", jets_.svtxdls2d, "svtxdls2d[nsvtx]/F");
    t->Branch("svtxm", jets_.svtxm, "svtxm[nsvtx]/F");
    t->Branch("svtxmcorr", jets_.svtxmcorr, "svtxmcorr[nsvtx]/F");
    t->Branch("svtxpt", jets_.svtxpt, "svtxpt[nsvtx]/F");
    t->Branch("svtxnormchi2", jets_.svtxnormchi2, "svtxnormchi2[nsvtx]/F");
    t->Branch("svtxchi2", jets_.svtxchi2, "svtxchi2[nsvtx]/F");
    /*    
    t->Branch("ntrkInSvtxNotInJet", &jets_.ntrkInSvtxNotInJet, "ntrkInSvtxNotInJet/I");
    t->Branch("trkInSvtxNotInJetSvId", jets_.trkInSvtxNotInJetSvId, "trkInSvtxNotInJetSvId[ntrkInSvtxNotInJet]/I");
    t->Branch("trkInSvtxNotInJetOtherJetId", jets_.trkInSvtxNotInJetOtherJetId, "trkInSvtxNotInJetOtherJetId[ntrkInSvtxNotInJet]/I");
    t->Branch("trkInSvtxNotInJetMatchSta", jets_.trkInSvtxNotInJetMatchSta, "trkInSvtxNotInJetMatchSta[ntrkInSvtxNotInJet]/I");
    t->Branch("trkInSvtxNotInJetPt", jets_.trkInSvtxNotInJetPt, "trkInSvtxNotInJetPt[ntrkInSvtxNotInJet]/F");
    t->Branch("trkInSvtxNotInJetEta", jets_.trkInSvtxNotInJetEta, "trkInSvtxNotInJetEta[ntrkInSvtxNotInJet]/F");
    t->Branch("trkInSvtxNotInJetPhi", jets_.trkInSvtxNotInJetPhi, "trkInSvtxNotInJetPhi[ntrkInSvtxNotInJet]/F");
    */
  }

  if (doLegacyBtagging_) {
    /* clear arrays */
    memset(jets_.discr_csvV2, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_muByIp3, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_muByPt, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_prob, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_probb, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_tcHighEff, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_tcHighPur, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_ssvHighEff, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_ssvHighPur, 0, MAXJETS * sizeof(float));
  }
  if (doCandidateBtagging_ && useNewBtaggers_) {
    for (auto& t : jets_discr_)
      for (auto& c : t.second)
        memset(c.second.data(), 0, MAXJETS * sizeof(float));
  }
  else if (doCandidateBtagging_) {
    memset(jets_.discr_deepCSV, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pfJP, 0, MAXJETS * sizeof(float));
  }
}

void HiInclusiveJetAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup) {
  int event = iEvent.id().event();
  int run = iEvent.id().run();
  int lumi = iEvent.id().luminosityBlock();

  jets_.run = run;
  jets_.evt = event;
  jets_.lumi = lumi;

  LogDebug("HiInclusiveJetAnalyzer") << "START event: " << event << " in run " << run << endl;

  // loop the events
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetTag_, jets);

  edm::Handle<reco::CaloJetCollection> calojets;
  if(doCaloJets_)iEvent.getByToken(caloJetTag_, calojets);

  edm::Handle<pat::JetCollection> matchedjets;
  iEvent.getByToken(matchTag_, matchedjets);

  edm::Handle<pat::JetCollection> originalCSjets;
  if (matchOriginal_) iEvent.getByToken(originalCSTag_, originalCSjets);

  if (doGenSubJets_)
    iEvent.getByToken(subjetGenTag_, gensubjets_);
  if (doGenSym_) {
    iEvent.getByToken(tokenGenSym_, genSymVM_);
    iEvent.getByToken(tokenGenDroppedBranches_, genDroppedBranchesVM_);
  }

  edm::Handle<edm::View<pat::PackedCandidate>> pfCandidates;
  iEvent.getByToken(pfCandidateLabel_, pfCandidates);
  edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfos;

  edm::Handle<reco::TrackToGenParticleMap> trackToGenParticleMap;
  //  edm::Handle<reco::TrackToGenParticleMap> genConstitToGenParticleMap;
    
  if (isMC_) {
    edm::Handle<reco::GenParticleCollection> genparts;
    iEvent.getByToken(genParticleSrc_, genparts);
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos );
    if (doTracks_) iEvent.getByToken(trackToGenParticleMapToken_, trackToGenParticleMap);
  }

  iEvent.getByToken(primaryVerticesToken_, primaryVertices);

  std::map<std::string, std::map<std::string, edm::Handle<JetTagCollection>>> jetTaggers;
  if (doCandidateBtagging_ && useNewBtaggers_) {
    for (const auto& t : jetTaggers_)
      for (const auto& c : t.second)
        jetTaggers[t.first].emplace(c.first, iEvent.getHandle(c.second));
  }

  // FILL JRA TREE
  jets_.nref = 0;
  jets_.ncalo = 0;
  
  jets_.nvtx = primaryVertices->size();
  if (doTracks_) jets_.ntrk = 0;
  if (doSvtx_) {
    jets_.nsvtx = 0;
    jets_.ntrkInSvtxNotInJet = 0;
  }
  int nsvtxCounterForTracks = 0;
  
  if (doJetConstituents_) {
    jets_.jtConstituentsId.clear();
    jets_.jtConstituentsE.clear();
    jets_.jtConstituentsPt.clear();
    jets_.jtConstituentsEta.clear();
    jets_.jtConstituentsPhi.clear();
    jets_.jtConstituentsM.clear();
    jets_.jtSDConstituentsE.clear();
    jets_.jtSDConstituentsPt.clear();
    jets_.jtSDConstituentsEta.clear();
    jets_.jtSDConstituentsPhi.clear();
    jets_.jtSDConstituentsM.clear();

    jets_.refConstituentsId.clear();
    jets_.refConstituentsE.clear();
    jets_.refConstituentsPt.clear();
    jets_.refConstituentsEta.clear();
    jets_.refConstituentsPhi.clear();
    jets_.refConstituentsM.clear();
    jets_.refSDConstituentsE.clear();
    jets_.refSDConstituentsPt.clear();
    jets_.refSDConstituentsEta.clear();
    jets_.refSDConstituentsPhi.clear();
    jets_.refSDConstituentsM.clear();

    jets_.genConstituentsId.clear();
    jets_.genConstituentsE.clear();
    jets_.genConstituentsPt.clear();
    jets_.genConstituentsEta.clear();
    jets_.genConstituentsPhi.clear();
    jets_.genConstituentsM.clear();
    jets_.genSDConstituentsE.clear();
    jets_.genSDConstituentsPt.clear();
    jets_.genSDConstituentsEta.clear();
    jets_.genSDConstituentsPhi.clear();
    jets_.genSDConstituentsM.clear();
  }

  auto getTag = [](const edm::Handle<reco::JetTagCollection> &bTags,const pat::Jet &jet) {
    float tagValue(-999),maxDR(3.1415);
    for (const auto &t : *bTags) {
      auto const dR = deltaR(jet, *(t.first));
      if (dR>maxDR) continue;
      maxDR=dR;
      tagValue=t.second;
    }
    if(maxDR>0.4) tagValue=-999;
    return tagValue;
  };

  // Count original jets with pt > cut, eta < cut
  /*   int countoriginal = 0;
  //  std::cout << "Check original jets " << std::endl;
  for (unsigned int j = 0; j < originalCSjets->size(); ++j) {
    const pat::Jet& jet = (*originalCSjets)[j];
    if (jet.pt() < jetPtMin_) continue;
    if (std::abs(jet.eta()) > jetAbsEtaMax_) continue;
    //    std::cout << "original jet with pt: " << jet.pt() << endl;
    countoriginal++; 
  }
  jets_.nCSjets = countoriginal; */
  
  for (unsigned int j = 0; j < jets->size(); ++j) {
    const pat::Jet& jet = (*jets)[j];

    auto pt = useRawPt_ ? jet.correctedJet("Uncorrected").pt() : jet.pt();
    if (pt < jetPtMin_)
      continue;
    if (std::abs(jet.eta()) > jetAbsEtaMax_)
      continue;

    int matchCSIndex = -1;
    jets_.jtptCS[jets_.nref] = -999;
    jets_.jtCSdr[jets_.nref] = -999;
    
    if (matchOriginal_) {   // 
      double drMin = 100;
      for (unsigned int imatch = 0; imatch < originalCSjets->size(); ++imatch) {
	const pat::Jet& mjet = (*originalCSjets)[imatch];
	double dr = deltaR(jet, mjet);
	if (dr < drMin) {
	  drMin = dr;
	  matchCSIndex = imatch;
	}
      }
      const pat::Jet& mjet = (*originalCSjets)[matchCSIndex];
      jets_.jtptCS[jets_.nref] = mjet.pt();
      jets_.jtCSdr[jets_.nref] = drMin;
    }
    

    if (doCandidateBtagging_ && useNewBtaggers_) {
      for (const auto& t : jetTaggers) {
        auto& discr = jets_discr_.at(t.first);
        if (t.first == "pfJP")
          discr["b"][jets_.nref] = getTag(t.second.at("probb"),jet);
        else if (t.first == "deepCSV")
          discr["b"][jets_.nref]  = getTag(t.second.at("probb"),jet)+getTag(t.second.at("probbb"),jet);
        else if (t.first == "deepFlavour" || t.first == "particleTransformer" || t.first == "unifiedParticleTransformer")
          discr["b"][jets_.nref] = getTag(t.second.at("probb"),jet)+getTag(t.second.at("probbb"),jet)+getTag(t.second.at("problepb"),jet);
        if (t.first == "unifiedParticleTransformer") {
          /* float tag(0.0);
	            for (const auto& n : {"probtaup1h0p", "probtaup1h1p", "probtaup1h2p", "probtaup3h0p", "probtaup3h1p", "probtaum1h0p", "probtaum1h1p", "probtaum1h2p", "probtaum3h0p", "probtaum3h1p"})
            tag += getTag(t.second.at(n), jet);
	    discr["probtau"][jets_.nref] = tag; */
          for (const auto& c : t.second)
            if (c.first.rfind("probtau",0)!=0)
              discr[c.first][jets_.nref] = getTag(c.second,jet);
        }
      }
    }
    else if (doCandidateBtagging_) {
      jets_.discr_deepCSV[jets_.nref] = jet.bDiscriminator(deepCSVJetTags_);
      jets_.discr_pfJP[jets_.nref] = jet.bDiscriminator(pfJPJetTags_);
    }
    if (doLegacyBtagging_) {
      jets_.discr_ssvHighEff[jets_.nref] = jet.bDiscriminator(simpleSVHighEffBJetTags_);
      jets_.discr_ssvHighPur[jets_.nref] = jet.bDiscriminator(simpleSVHighPurBJetTags_);
      jets_.discr_csvV2[jets_.nref] = jet.bDiscriminator(combinedSVV2BJetTags_);
      jets_.discr_prob[jets_.nref] = jet.bDiscriminator(jetPBJetTags_);
      jets_.discr_probb[jets_.nref] = jet.bDiscriminator(jetBPBJetTags_);
      jets_.discr_tcHighEff[jets_.nref] = jet.bDiscriminator(trackCHEBJetTags_);
      jets_.discr_tcHighPur[jets_.nref] = jet.bDiscriminator(trackCHPBJetTags_);

      const edm::View<pat::PackedCandidate>* pfCandidateColl = &(*pfCandidates);
      int pfMuonIndex = getPFJetMuon(jet, pfCandidateColl);

      if (pfMuonIndex >= 0) {
        const pat::PackedCandidate muon = pfCandidates->at(pfMuonIndex);
        jets_.mupt[jets_.nref] = muon.pt();
        jets_.mueta[jets_.nref] = muon.eta();
        jets_.muphi[jets_.nref] = muon.phi();
        jets_.mue[jets_.nref] = muon.energy();
        jets_.mudr[jets_.nref] = reco::deltaR(jet, muon);
        jets_.muptrel[jets_.nref] = getPtRel(muon, jet);
        jets_.muchg[jets_.nref] = muon.charge();
      } else {
        jets_.mupt[jets_.nref] = 0.0;
        jets_.mueta[jets_.nref] = 0.0;
        jets_.muphi[jets_.nref] = 0.0;
        jets_.mue[jets_.nref] = 0.0;
        jets_.mudr[jets_.nref] = 9.9;
        jets_.muptrel[jets_.nref] = 0.0;
        jets_.muchg[jets_.nref] = 0;
      }
    }

    if (doHiJetID_) {
      // Jet ID variables

      jets_.muMax[jets_.nref] = 0;
      jets_.muSum[jets_.nref] = 0;
      jets_.muN[jets_.nref] = 0;

      jets_.eMax[jets_.nref] = 0;
      jets_.eSum[jets_.nref] = 0;
      jets_.eN[jets_.nref] = 0;

      jets_.neutralMax[jets_.nref] = 0;
      jets_.neutralSum[jets_.nref] = 0;
      jets_.neutralN[jets_.nref] = 0;

      jets_.photonMax[jets_.nref] = 0;
      jets_.photonSum[jets_.nref] = 0;
      jets_.photonN[jets_.nref] = 0;
      jets_.photonHardSum[jets_.nref] = 0;
      jets_.photonHardN[jets_.nref] = 0;

      jets_.chargedMax[jets_.nref] = 0;
      jets_.chargedSum[jets_.nref] = 0;
      jets_.chargedN[jets_.nref] = 0;
      jets_.chargedHardSum[jets_.nref] = 0;
      jets_.chargedHardN[jets_.nref] = 0;

      jets_.trackMax[jets_.nref] = 0;
      jets_.trackSum[jets_.nref] = 0;
      jets_.trackN[jets_.nref] = 0;
      jets_.trackHardSum[jets_.nref] = 0;
      jets_.trackHardN[jets_.nref] = 0;

      jets_.genChargedSum[jets_.nref] = 0;
      jets_.genHardSum[jets_.nref] = 0;

      jets_.signalChargedSum[jets_.nref] = 0;
      jets_.signalHardSum[jets_.nref] = 0;

      jets_.subid[jets_.nref] = -1;

      for (unsigned int icand = 0; icand < pfCandidates->size(); ++icand) {
        const pat::PackedCandidate& t = (*pfCandidates)[icand];

        if (!t.hasTrackDetails())
          continue;

        reco::Track const& track = t.pseudoTrack();

        if (useQuality_) {
          bool goodtrack = track.quality(reco::TrackBase::qualityByName(trackQuality_));
          if (!goodtrack)
            continue;
        }

        double dr = deltaR(jet, track);
        if (dr < rParam) {
          double ptcand = track.pt();
          jets_.trackSum[jets_.nref] += ptcand;
          jets_.trackN[jets_.nref] += 1;

          if (ptcand > hardPtMin_) {
            jets_.trackHardSum[jets_.nref] += ptcand;
            jets_.trackHardN[jets_.nref] += 1;
          }
          if (ptcand > jets_.trackMax[jets_.nref])
            jets_.trackMax[jets_.nref] = ptcand;
        }
      }

      reco::PFCandidate converter = reco::PFCandidate();
      for (unsigned int icand = 0; icand < pfCandidates->size(); ++icand) {
        const pat::PackedCandidate& track = (*pfCandidates)[icand];
        double dr = deltaR(jet, track);
        if (dr < rParam) {
          double ptcand = track.pt();
          int pfid = converter.translatePdgIdToType(track.pdgId());

          switch (pfid) {
            case 1:
              jets_.chargedSum[jets_.nref] += ptcand;
              jets_.chargedN[jets_.nref] += 1;
              if (ptcand > hardPtMin_) {
                jets_.chargedHardSum[jets_.nref] += ptcand;
                jets_.chargedHardN[jets_.nref] += 1;
              }
              if (ptcand > jets_.chargedMax[jets_.nref])
                jets_.chargedMax[jets_.nref] = ptcand;
              break;

            case 2:
              jets_.eSum[jets_.nref] += ptcand;
              jets_.eN[jets_.nref] += 1;
              if (ptcand > jets_.eMax[jets_.nref])
                jets_.eMax[jets_.nref] = ptcand;
              break;

            case 3:
              jets_.muSum[jets_.nref] += ptcand;
              jets_.muN[jets_.nref] += 1;
              if (ptcand > jets_.muMax[jets_.nref])
                jets_.muMax[jets_.nref] = ptcand;
              break;

            case 4:
              jets_.photonSum[jets_.nref] += ptcand;
              jets_.photonN[jets_.nref] += 1;
              if (ptcand > hardPtMin_) {
                jets_.photonHardSum[jets_.nref] += ptcand;
                jets_.photonHardN[jets_.nref] += 1;
              }
              if (ptcand > jets_.photonMax[jets_.nref])
                jets_.photonMax[jets_.nref] = ptcand;
              break;

            case 5:
              jets_.neutralSum[jets_.nref] += ptcand;
              jets_.neutralN[jets_.nref] += 1;
              if (ptcand > jets_.neutralMax[jets_.nref])
                jets_.neutralMax[jets_.nref] = ptcand;
              break;

            default:
              break;
          }
        }
      }
    }

    int matchIndex = -1;
    if (doMatch_) {
      // Alternative reconstruction matching (PF for calo, calo for PF)

      double drMin = 100;
      for (unsigned int imatch = 0; imatch < matchedjets->size(); ++imatch) {
        const pat::Jet& mjet = (*matchedjets)[imatch];

        double dr = deltaR(jet, mjet);
        if (dr < drMin) {
          jets_.mjtPt[jets_.nref] = mjet.pt();

          jets_.mjtRawPt[jets_.nref] = mjet.correctedJet("Uncorrected").pt();
	  //	  jets_.mjtPu[jets_.nref] = mjet.pileup();
          if (isMC_) {
            jets_.mjtHadronFlavor[jets_.nref] = mjet.hadronFlavour();
            jets_.mjtPartonFlavor[jets_.nref] = mjet.partonFlavour();

	    for (const JetFlavourInfoMatching& jetFlavourInfoMatching : *jetFlavourInfos) {
	      if (deltaR(mjet.p4(), jetFlavourInfoMatching.first->p4()) < 1e-6) {
		JetFlavourInfo jetInfo = jetFlavourInfoMatching.second;
		const GenParticleRefVector &bHadronsInJet = jetInfo.getbHadrons();
		const GenParticleRefVector &cHadronsInJet = jetInfo.getcHadrons();

		jets_.mjtNbHad[jets_.nref] = bHadronsInJet.size();
		jets_.mjtNcHad[jets_.nref] = cHadronsInJet.size();

		const GenParticleRefVector &partonsInJet = jetInfo.getPartons();

		int nb = 0, nc = 0;
		for (GenParticleRefVector::const_iterator it = partonsInJet.begin(); it != partonsInJet.end(); ++it) {
		  int parFlav = (*it)->pdgId();
		  if (abs(parFlav)==5) nb++;
		  if (abs(parFlav)==4) nc++;
		}

		jets_.mjtNbPar[jets_.nref] = nb;
		jets_.mjtNcPar[jets_.nref] = nc;
		
		break;
	      }
	    } // end loop over flavour info
          }

          jets_.mjtR[jets_.nref] = dr;
          drMin = dr;
	  matchIndex = imatch;
        }
      }
    }


   if (doSvtx_ && matchIndex>=0 ) {

      const pat::Jet& mjet = (*matchedjets)[matchIndex];
      if(mjet.hasTagInfo(svTagInfoLabel_.c_str())){
	const reco::CandSecondaryVertexTagInfo *svTagInfo = mjet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());
	
	int nsv = svTagInfo->nVertices();
	jets_.jtNsvtx[jets_.nref] = 0;
	for (int isv = 0; isv < nsv; isv++) {
	  int ijetSvtx = jets_.nsvtx + isv;
	  jets_.svtxNtrk[ijetSvtx] = svTagInfo->nVertexTracks(isv);
	  
	  Measurement1D dl3d = svTagInfo->flightDistance(isv);
	  jets_.svtxdl[ijetSvtx] = dl3d.value();
	  jets_.svtxdls[ijetSvtx] = dl3d.significance();
	  
	  Measurement1D dl2d = svTagInfo->flightDistance(isv, 2);
	  jets_.svtxdl2d[ijetSvtx] = dl2d.value();
	  jets_.svtxdls2d[ijetSvtx] = dl2d.significance();
	  
	  const VertexCompositePtrCandidate svtx = svTagInfo->secondaryVertex(isv);
	  double svtxM = svtx.p4().mass();
	  double svtxPt = svtx.p4().pt();
	  double normalizedChi2 = svtx.vertexNormalizedChi2();
	  double Chi2 = svtx.vertexChi2();
	  
	  //mCorr=srqt(m^2+p^2sin^2(th)) + p*sin(th)
	  double sinth = svtx.p4().Vect().Unit().Cross((svTagInfo->flightDirection(isv)).unit()).Mag2();
	  sinth = sqrt(sinth);
	  double underRoot = std::pow(svtxM, 2) + (std::pow(svtxPt, 2) * std::pow(sinth, 2));
	  double svtxMcorr = std::sqrt(underRoot) + (svtxPt * sinth);
	  
	  jets_.svtxnormchi2[ijetSvtx] = normalizedChi2;
	  jets_.svtxchi2[ijetSvtx] = Chi2;
	  jets_.svtxm[ijetSvtx] = svtxM;
	  jets_.svtxmcorr[ijetSvtx] = svtxMcorr;
	  jets_.svtxpt[ijetSvtx] = svtxPt;
	  
	  jets_.svtxJetId[ijetSvtx] = jets_.nref;
	  
	  const std::vector<reco::CandidatePtr> svTracks = svTagInfo->vertexTracks(isv);
	  
	} // sv loop
	jets_.jtNsvtx[jets_.nref] = nsv;
	jets_.nsvtx += nsv;
      } // endif doSvtx_
    }
    //for (auto label : jet.tagInfoLabels())  std::cout << label << std::endl;
    

    if (doTracks_ && matchIndex>=0 ) {

      const pat::Jet& mjet = (*matchedjets)[matchIndex];
      if (mjet.hasTagInfo(ipTagInfoLabel_.c_str())) {

	jets_.jtNtrk[jets_.nref] = 0;
	jets_.jtptCh[jets_.nref] = 0.;
	const reco::CandIPTagInfo *ipTagInfo = mjet.tagInfoCandIP(ipTagInfoLabel_.c_str());
	const std::vector<reco::btag::TrackIPData> ipData = ipTagInfo->impactParameterData();
	const std::vector<reco::CandidatePtr> ipTracks = ipTagInfo->selectedTracks();
	
	reco::Candidate::PolarLorentzVector chJet(0., 0., 0., 0.);
	
	// For debugging
        // for (auto itIPTrack : ipTracks)  std::cout << " ip track: " << itIPTrack->pt() << " " <<  itIPTrack->eta() << " " <<  itIPTrack->phi() << " " <<  std::endl;
		  
	for (const reco::CandidatePtr &constit : mjet.getJetConstituents()) {
	  // std::cout << "new jet constit with pt, eta, phi " << constit->pt() << " "  << constit->eta() << " "  << constit->phi() << " " << constit->charge() <<  std::endl;
	  if (constit->charge() == 0) continue;
	  if (constit->pt() < trkPtCut_) continue;
	  
	  auto itIPTrack = std::find(ipTracks.begin(), ipTracks.end(), constit);
	  if (itIPTrack == ipTracks.end()) continue;
	  
	  // Check if the track was dropped from the aggregation
	  /* if (isMC_) {
	     bool drop = false;
	     double eps = 1e-4;
	     // std::cout << "before droppedTracks" << std::endl;
	     for (size_t idropped=0; idropped<droppedTracks->size(); idropped++) {
	     reco::PFCandidate droppedTrack = (*droppedTracks)[idropped];
	     if (std::abs(constit->eta()-droppedTrack.eta())>eps) continue;
	     if (std::abs(constit->phi()-droppedTrack.phi())>eps) continue;
	     if (std::abs(constit->pt()-droppedTrack.pt())>eps) continue;
	     drop = true;
	     }
	     if (drop) continue;
	     // std::cout << "after droppedTracks" << std::endl;
	     } */
	  
	  int ijetTrack = jets_.ntrk + jets_.jtNtrk[jets_.nref];
	  int itrk = itIPTrack - ipTracks.begin();
	  
	  reco::Candidate::PolarLorentzVector constitV(0., 0., 0., 0.);
	  constitV.SetPt(constit->pt());
	  constitV.SetEta(constit->eta());
	  constitV.SetPhi(constit->phi());
	  constitV.SetM(constit->mass());
	  chJet += constitV;
	  
	  const reco::btag::TrackIPData trkIPData = ipData[itrk];
	  
	  jets_.trkJetId[ijetTrack] = jets_.nref;  
	  
	  jets_.trkPt[ijetTrack] = constit->pt();
	  jets_.trkEta[ijetTrack] = constit->eta();
	  jets_.trkPhi[ijetTrack] = constit->phi();
	  
	  jets_.trkIp3d[ijetTrack] = trkIPData.ip3d.value();
	  jets_.trkIp3dSig[ijetTrack] = trkIPData.ip3d.significance();
	  
	  jets_.trkIp2d[ijetTrack] = trkIPData.ip2d.value();
	  jets_.trkIp2dSig[ijetTrack] = trkIPData.ip2d.significance();
	  
	  jets_.trkIpProb3d[ijetTrack] = ipTagInfo->probabilities(0)[itrk];
	  jets_.trkIpProb2d[ijetTrack] = ipTagInfo->probabilities(1)[itrk];
	  
	  jets_.trkDistToAxis[ijetTrack] = trkIPData.distanceToJetAxis.value();
	  jets_.trkDistToAxisSig[ijetTrack] = trkIPData.distanceToJetAxis.significance();
	  
	  jets_.trkSvtxId[ijetTrack] = -1;
	
	  if (doSvtx_ && matchIndex >=0 && mjet.hasTagInfo(svTagInfoLabel_.c_str())) {
	    if(mjet.hasTagInfo(svTagInfoLabel_.c_str())){
	      const reco::CandSecondaryVertexTagInfo *svTagInfo = mjet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());
	      int nsv = svTagInfo->nVertices();
	      for (int isv = 0; isv < nsv; isv++) {
		const std::vector<reco::CandidatePtr> svTracks = svTagInfo->vertexTracks(isv);
		auto itSVTrack = std::find(svTracks.begin(), svTracks.end(), constit); // TODO: replace
		if (itSVTrack == svTracks.end()) continue;
		jets_.trkSvtxId[ijetTrack] = nsvtxCounterForTracks + isv;
	      } // end sv loop for tracks
	    } // end doSvtx_
	  }
	  Int_t status = -1; // default, no match
	  if (isMC_ && trackToGenParticleMap->find(constit) != trackToGenParticleMap->end()) { 
	    edm::Ptr<pat::PackedGenParticle> matchGenParticle = trackToGenParticleMap->at(constit);
	    status = matchGenParticle->status();
	  }
	  jets_.trkMatchSta[ijetTrack] = status;
	  jets_.trkPdgId[ijetTrack] = constit->pdgId();
	  
	  const reco::Track *constitTrack = constit->bestTrack();
	  if (constitTrack) {
	    // std::cout << "track exists " << std::endl;
	    // std::cout << "testTrack dz " << testTrack->dz() << std::endl;
	    jets_.trkDz[ijetTrack] = constitTrack->dz(primaryVertices->at(0).position());
	  } else {
	    jets_.trkDz[ijetTrack] = -100000.;
	  }
	  
	  jets_.jtNtrk[jets_.nref]++;
	} // jet constituent loop
	//      std::cout << "pt from constituents" << ptcounter << std::endl;
	
	jets_.jtptCh[jets_.nref] = chJet.pt();
	nsvtxCounterForTracks += jets_.jtNsvtx[jets_.nref];
	jets_.ntrk += jets_.jtNtrk[jets_.nref];
      } // endif doTracks_
      
    }
 

    jets_.rawpt[jets_.nref] = jet.correctedJet("Uncorrected").pt();
    jets_.jtpt[jets_.nref] = jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();
    //    jets_.jty[jets_.nref] = jet.eta();
    //    jets_.jtpu[jets_.nref] = jet.pileup();
    //    jets_.jtm[jets_.nref] = jet.mass();
    //    jets_.jtarea[jets_.nref] = jet.jetArea();


    //------------------------------------------------------------------

    jets_.jtsym[jets_.nref] = -999.;
    jets_.jtdroppedBranches[jets_.nref] = -999;

    if (doSubJets_)
      analyzeSubjets(jet);

    if (jet.hasUserFloat(jetName_ + "Jets:sym"))
      jets_.jtsym[jets_.nref] = jet.userFloat(jetName_ + "Jets:sym");
    if (jet.hasUserInt(jetName_ + "Jets:droppedBranches"))
      jets_.jtdroppedBranches[jets_.nref] = jet.userInt(jetName_ + "Jets:droppedBranches");

    if (doPFjetID) {
      if (jet.isPFJet()) {
	jets_.jtPfCHF[jets_.nref] = jet.chargedHadronEnergyFraction();
	jets_.jtPfNHF[jets_.nref] = jet.neutralHadronEnergyFraction();
	jets_.jtPfCEF[jets_.nref] = jet.chargedEmEnergyFraction();
	jets_.jtPfNEF[jets_.nref] = jet.neutralEmEnergyFraction();
	jets_.jtPfMUF[jets_.nref] = jet.muonEnergyFraction();

	jets_.jtPfCHM[jets_.nref] = jet.chargedHadronMultiplicity();
	jets_.jtPfNHM[jets_.nref] = jet.neutralHadronMultiplicity();
	jets_.jtPfCEM[jets_.nref] = jet.electronMultiplicity();
	jets_.jtPfNEM[jets_.nref] = jet.photonMultiplicity();
	jets_.jtPfMUM[jets_.nref] = jet.muonMultiplicity();
      } else {
	jets_.jtPfCHF[jets_.nref] = 0;
	jets_.jtPfNHF[jets_.nref] = 0;
	jets_.jtPfCEF[jets_.nref] = 0;
	jets_.jtPfNEF[jets_.nref] = 0;
	jets_.jtPfMUF[jets_.nref] = 0;
	
	jets_.jtPfCHM[jets_.nref] = 0;
	jets_.jtPfNHM[jets_.nref] = 0;
	jets_.jtPfCEM[jets_.nref] = 0;
	jets_.jtPfNEM[jets_.nref] = 0;
	jets_.jtPfMUM[jets_.nref] = 0;
      }
    }
    
    //    if(isMC_){

    //      for(UInt_t i = 0; i < genparts->size(); ++i){
    // const reco::GenParticle& p = (*genparts)[i];
    // if ( p.status()!=1 || p.charge()==0) continue;
    // double dr = deltaR(jet,p);
    // if(dr < rParam){
    //   double ppt = p.pt();
    //   jets_.genChargedSum[jets_.nref] += ppt;
    //   if(ppt > hardPtMin_) jets_.genHardSum[jets_.nref] += ppt;
    //   if(p.collisionId() == 0){
    //     jets_.signalChargedSum[jets_.nref] += ppt;
    //     if(ppt > hardPtMin_) jets_.signalHardSum[jets_.nref] += ppt;
    //   }
    // }
    //      }
    //    }

    //    IterativeDeclusteringRec(groom_type, groom_combine, jet, sub1Hyb, sub2Hyb);
    if (runSubstructure) IterativeDeclusteringRec(0, 1, jet);
    
    if (isMC_) {
      const reco::GenJet* genjet = jet.genJet();

      if (genjet) {
        jets_.refpt[jets_.nref] = genjet->pt();
        jets_.refeta[jets_.nref] = genjet->eta();
        jets_.refphi[jets_.nref] = genjet->phi();
	//        jets_.refm[jets_.nref] = genjet->mass();
	//        jets_.refarea[jets_.nref] = genjet->jetArea();
	//        jets_.refy[jets_.nref] = genjet->eta();
        jets_.refdphijt[jets_.nref] = reco::deltaPhi(jet.phi(), genjet->phi());
        jets_.refdrjt[jets_.nref] = reco::deltaR(jet.eta(), jet.phi(), genjet->eta(), genjet->phi());

	if (runSubstructure) {
	  IterativeDeclusteringGen(0, 1, *genjet);

	  if (doSplitMatching_) {
	    TruthRecoRecoTruthMatching_SD();
	    TruthRecoRecoTruthMatching_latekt();
	    //vangi's
	    jets_.jtJetSplits = {};
	    jets_.refJetSplits = {};	    
	  }
	}

        if (doSubEvent_) {
          const GenParticle* gencon = genjet->getGenConstituent(0);
          jets_.subid[jets_.nref] = gencon->collisionId();
        }

        if (doGenSubJets_)
          analyzeRefSubjets(*genjet);

      } else {
        jets_.refpt[jets_.nref] = -999.;
        jets_.refeta[jets_.nref] = -999.;
        jets_.refphi[jets_.nref] = -999.;
        jets_.refm[jets_.nref] = -999.;
        jets_.refarea[jets_.nref] = -999.;
        jets_.refy[jets_.nref] = -999.;
        jets_.refdphijt[jets_.nref] = -999.;
        jets_.refdrjt[jets_.nref] = -999.;

        if (doJetConstituents_) {
          jets_.refConstituentsId.emplace_back(1, -999);
          jets_.refConstituentsE.emplace_back(1, -999);
          jets_.refConstituentsPt.emplace_back(1, -999);
          jets_.refConstituentsEta.emplace_back(1, -999);
          jets_.refConstituentsPhi.emplace_back(1, -999);
          jets_.refConstituentsM.emplace_back(1, -999);

          jets_.refSDConstituentsId.emplace_back(1, -999);
          jets_.refSDConstituentsE.emplace_back(1, -999);
          jets_.refSDConstituentsPt.emplace_back(1, -999);
          jets_.refSDConstituentsEta.emplace_back(1, -999);
          jets_.refSDConstituentsPhi.emplace_back(1, -999);
          jets_.refSDConstituentsM.emplace_back(1, -999);
        }

        if (doGenSubJets_) {
          jets_.refptG[jets_.nref] = -999.;
          jets_.refetaG[jets_.nref] = -999.;
          jets_.refphiG[jets_.nref] = -999.;
          jets_.refmG[jets_.nref] = -999.;
          jets_.refsym[jets_.nref] = -999.;
          jets_.refdroppedBranches[jets_.nref] = -999;

          jets_.refSubJetPt.emplace_back(1, -999);
          jets_.refSubJetEta.emplace_back(1, -999);
          jets_.refSubJetPhi.emplace_back(1, -999);
          jets_.refSubJetM.emplace_back(1, -999);
        }
      }
      jets_.reftau1[jets_.nref] = -999.;
      jets_.reftau2[jets_.nref] = -999.;
      jets_.reftau3[jets_.nref] = -999.;

      jets_.refparton_flavorForB[jets_.nref] = jet.partonFlavour();

      //      if(jet.genParton()){
      // // matched partons
      // const reco::GenParticle & parton = *jet.genParton();

      // jets_.refparton_pt[jets_.nref] = parton.pt();
      // jets_.refparton_flavor[jets_.nref] = parton.pdgId();

      //      } else {
      //      jets_.refparton_pt[jets_.nref] = -999;
      //  jets_.refparton_flavor[jets_.nref] = -999;
      //      }
    }

    jets_.nref++;
  }

  if (isMC_ and fillGenJets_) {
    if (useHepMC_) {
      edm::Handle<HepMCProduct> hepMCProduct;
      iEvent.getByToken(eventInfoTag_, hepMCProduct);
      const HepMC::GenEvent* MCEvt = hepMCProduct->GetEvent();

      std::pair<HepMC::GenParticle*, HepMC::GenParticle*> beamParticles = MCEvt->beam_particles();
      jets_.beamId1 = (beamParticles.first != 0) ? beamParticles.first->pdg_id() : 0;
      jets_.beamId2 = (beamParticles.second != 0) ? beamParticles.second->pdg_id() : 0;
    }

    /*    edm::Handle<GenEventInfoProduct> hEventInfo;
    iEvent.getByToken(eventGenInfoTag_, hEventInfo);

    // binning values and qscale appear to be equivalent, but binning values not always present
    jets_.pthat = hEventInfo->qScale(); */

    edm::Handle<edm::View<reco::GenJet>> genjets;
    iEvent.getByToken(genjetTag_, genjets);

    //get gen-level n-jettiness
    edm::Handle<edm::ValueMap<float>> genTau1s;
    edm::Handle<edm::ValueMap<float>> genTau2s;
    edm::Handle<edm::ValueMap<float>> genTau3s;
    if (doGenTaus_) {
      iEvent.getByToken(tokenGenTau1_, genTau1s);
      iEvent.getByToken(tokenGenTau2_, genTau2s);
      iEvent.getByToken(tokenGenTau3_, genTau3s);
    }

    jets_.ngen = 0;

    for (unsigned int igen = 0; igen < genjets->size(); ++igen) {
      const reco::GenJet& genjet = (*genjets)[igen];
      float genjet_pt = genjet.pt();

      float tau1 = -999.;
      float tau2 = -999.;
      float tau3 = -999.;
      Ptr<reco::GenJet> genJetPtr = genjets->ptrAt(igen);
      if (doGenTaus_) {
        tau1 = (*genTau1s)[genJetPtr];
        tau2 = (*genTau2s)[genJetPtr];
        tau3 = (*genTau3s)[genJetPtr];
      }

      // find matching patJet if there is one
      jets_.gendrjt[jets_.ngen] = -1.0;
      jets_.genmatchindex[jets_.ngen] = -1;

      for (int ijet = 0; ijet < jets_.nref; ++ijet) {
        // poor man's matching, someone fix please

        double deltaPt = fabs(genjet.pt() - jets_.refpt[ijet]);  //Note: precision of this ~ .0001, so cut .01
        double deltaEta = fabs(
            genjet.eta() -
            jets_.refeta
                [ijet]);  //Note: precision of this is  ~.0000001, but keep it low, .0001 is well below cone size and typical pointing resolution
        double deltaPhi = fabs(reco::deltaPhi(
            genjet.phi(),
            jets_.refphi
                [ijet]));  //Note: precision of this is  ~.0000001, but keep it low, .0001 is well below cone size and typical pointing resolution

        if (deltaPt < 0.01 && deltaEta < .0001 && deltaPhi < .0001) {
          if (genjet_pt > genPtMin_) {
            jets_.genmatchindex[jets_.ngen] = (int)ijet;
            jets_.gendphijt[jets_.ngen] = reco::deltaPhi(jets_.refphi[ijet], genjet.phi());
            jets_.gendrjt[jets_.ngen] =
                sqrt(pow(jets_.gendphijt[jets_.ngen], 2) + pow(fabs(genjet.eta() - jets_.refeta[ijet]), 2));
          }
          if (doGenTaus_) {
            jets_.reftau1[ijet] = tau1;
            jets_.reftau2[ijet] = tau2;
            jets_.reftau3[ijet] = tau3;
          }
          break;
        }
      }

      //-------------------------------------------------

      // threshold to reduce size of output in minbias PbPb
      if (genjet_pt > genPtMin_) {
        jets_.genpt[jets_.ngen] = genjet_pt;
        jets_.geneta[jets_.ngen] = genjet.eta();
        jets_.genphi[jets_.ngen] = genjet.phi();
        //        jets_.genm[jets_.ngen] = genjet.mass();
	//        jets_.geny[jets_.ngen] = genjet.eta();

        if (doGenTaus_) {
          jets_.gentau1[jets_.ngen] = tau1;
          jets_.gentau2[jets_.ngen] = tau2;
          jets_.gentau3[jets_.ngen] = tau3;
        }

        if (doGenSubJets_)
          analyzeGenSubjets(genjet);

        if (doSubEvent_) {
          const GenParticle* gencon = genjet.getGenConstituent(0);
          jets_.gensubid[jets_.ngen] = gencon->collisionId();
        }
        jets_.ngen++;
      }
    }
  }
  
  if(doCaloJets_){
    for (unsigned int j = 0; j < calojets->size(); ++j) {
      const reco::Jet& jet = (*calojets)[j];
      jets_.calopt[jets_.ncalo] = jet.pt();
      jets_.caloeta[jets_.ncalo] = jet.eta();
      jets_.calophi[jets_.ncalo] = jet.phi();
      jets_.ncalo++;
    }
  }

  t->Fill();

  //memset(&jets_,0,sizeof jets_);
  jets_ = {0};
}



void HiInclusiveJetAnalyzer::IterativeDeclusteringRec(double groom_type, double groom_combine, const reco::Jet& jet)
{
  
  Int_t nsplit = 0;

  double z = 0;
  
  double zg_SD = 0;
  double zg_latekt = 0;

  double ktg_SD = 0;
  double ktg_latekt = 0;

  double rg_SD = 0;
  double rg_latekt = 0;

  Int_t SD_split = -1; 
  Int_t latekt_split = -1;

  double jet_radius_ca = 1.0;

  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
  // Reclustering jet constituents with new algorithm
  
  try{
    std::vector<fastjet::PseudoJet> particles = {};                         
    auto daughters = jet.getJetConstituents();

    for (auto it = daughters.begin(); it!=daughters.end(); ++it){

      // std::cout << "scan jet consts in declustering, mass: " << (**it).mass() << " charge: " << (**it).charge() << " id: "<< (**it).pdgId() << " pt: " << (**it).pt() << std::endl;
      if (doChargedConstOnly_ && (**it).charge()==0) continue;
      
      //      std::cout << "scan jet consts in declustering, mass: " << (**it).mass() << " charge: " << (**it).charge() << " id: "<< (**it).pdgId() << " pt: " << (**it).pt() << std::endl;
      if ((**it).mass() < 0 and (**it).charge() < -4) {      // Assumes charge is set in aggregator
	jets_.massHF[jets_.nref] = -((**it).mass());
	//	std::cout << "HF mass is " << -((**it).mass()) << std::endl;
      }
      //if we want only charged constituents and the daughter charge is 0, skip it
      
      if ((**it).pt()<1) continue; //Particle pt cut

      double PFE_scale = 1.;

      if (isMC_){ //if it is MC, rescale the 4-momentum of the particles by pfCCES(+-1%)
        if ((**it).charge()!=0){ //if Charged candidate

          if (doPFChargedEnergyScaleVar_ == 1.)       PFE_scale = 1. + 0.01;
          else if (doPFChargedEnergyScaleVar_ == -1.) PFE_scale = 1. - 0.01;
          else if (doPFChargedEnergyScaleVar_ == 0.)  PFE_scale = 1.;
          else cout << "you should not be here (Charged)" << endl;
        }

        else if ((**it).pdgId()==130){//If Neutral candidate

          if(doPFNeutralEnergyScaleVar_ == 1.)       PFE_scale = 1. + 0.05;
          else if(doPFNeutralEnergyScaleVar_ == -1.) PFE_scale = 1. - 0.05;
          else if(doPFNeutralEnergyScaleVar_ == 0.)  PFE_scale = 1.;
          else cout << "you should not be here (Neutral)" << endl;
        }

        else if ((**it).pdgId()==22){  //If Gamma candidate

          if (doPFGammaEnergyScaleVar_ == 1.)         PFE_scale = 1. + 0.03;
          else if (doPFGammaEnergyScaleVar_ == -1.)   PFE_scale = 1. - 0.03;
          else if (doPFGammaEnergyScaleVar_ == 0.)    PFE_scale = 1.;
          else cout << "you should not be here (Gamma)" << endl;
        }

        else cout << "Found no charged, charged or photon candidaties; pdgID: " << (**it).pdgId() << std::endl;
      }
      particles.push_back(fastjet::PseudoJet((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale));
    }

    if (particles.size() == 0){ 
      jets_.jt_split_SD[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_split_latekt[jets_.nref] = std::numeric_limits<int>::min();
    }

    if (particles.size()!=0 ) {
      fastjet::ClusterSequence csiter(particles, jet_def);
      std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
      output_jets = sorted_by_pt(output_jets);

      fastjet::PseudoJet jj = output_jets[0];
      fastjet::PseudoJet j1;
      fastjet::PseudoJet j2;

      if(!jj.has_parents(j1,j2)) {
        jets_.jt_split_SD[jets_.nref] = std::numeric_limits<int>::min();
        jets_.jt_split_latekt[jets_.nref] = std::numeric_limits<int>::min();
      }
      
      int stopSD = 0;
      bool flagHF = false, flagHFSD = false, flagHFkt = false;
      
      while (jj.has_parents(j1,j2)) {
        if(j1.perp() < j2.perp()) std::swap(j1,j2);
	//        vector <fastjet::PseudoJet> constitj1 = sorted_by_pt(j1.constituents());

	std::vector<fastjet::PseudoJet> j1constits = j1.constituents();
	for (size_t icon = 0; icon < j1constits.size(); icon++) {
	  fastjet::PseudoJet constit = j1constits[icon];
	  if (abs(constit.m() + jets_.massHF[jets_.nref]) < 1e-3) {
	    flagHF = true; // the leading prong has the HF; massHF and constit.m have opposite signs
	    //std::cout << "pt  " << constit.pt() << " mass " << constit.m() << " " << flagHF << std::endl;
	    break;
	  }
	}
	//	std::cout << " ---------- " << std::endl;
	
        double delta_R = j1.delta_R(j2);
        if (doSplitMatching_ && isMC_) {
          jets_.jtJetSplits.push_back(j2);
        }
        double k_t = j2.perp()*delta_R;
        z = j2.perp()/(j1.perp()+j2.perp());

      //  std::cout << "Reco split " << nsplit << " with k_T=" << k_t << " z=" << z << " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl; 

        if (((groom_combine == 0) and (groom_type == 1) and (z > SDcut) and (stopSD == 0)) or ((groom_combine == 1) and (z > SDcut) and (stopSD == 0))) { 
          stopSD = 1;
          zg_SD = z;
          rg_SD  = delta_R;
          ktg_SD = k_t;
          SD_split = nsplit;
	  flagHFSD = flagHF;
        }
        
        if (((groom_combine == 0) and (groom_type == 0) and (k_t > latektcut)) or ((groom_combine == 1) and (k_t > latektcut))) {
          zg_latekt = z;
          rg_latekt  = delta_R;
          ktg_latekt = k_t;
          latekt_split = nsplit;
	  flagHFkt = flagHF;
        }
        jj = j1;
        nsplit = nsplit+1;
      }
    //  std::cout << "end results: " << zg_SD << " " << zg_latekt << std::endl;
    jets_.jt_z_SD[jets_.nref] = zg_SD;
    jets_.jt_rg_SD[jets_.nref] = rg_SD;
    jets_.jt_ktg_SD[jets_.nref] = ktg_SD;
    jets_.jt_split_SD[jets_.nref] = SD_split;
    jets_.jt_hasHF_SD[jets_.nref] = flagHFSD;

    jets_.jt_z_latekt[jets_.nref] = zg_latekt;
    jets_.jt_rg_latekt[jets_.nref] = rg_latekt;
    jets_.jt_ktg_latekt[jets_.nref] = ktg_latekt;
    jets_.jt_split_latekt[jets_.nref] = latekt_split;
    jets_.jt_hasHF_latekt[jets_.nref] = flagHFkt;
    }
  } 

  catch (fastjet::Error const&){ 
    cout << "Fastjet error" << endl;
  }
  catch (Int_t MyNum) {
    cout<<"catch neutralN = "<<jets_.neutralN[jets_.nref]<<endl;
    jets_.jt_z_SD[jets_.nref] = 0;
    jets_.jt_rg_SD[jets_.nref] = 0;
    jets_.jt_ktg_SD[jets_.nref] = 0;
    jets_.jt_z_latekt[jets_.nref] = 0;
    jets_.jt_rg_latekt[jets_.nref] = 0;
    jets_.jt_ktg_latekt[jets_.nref] = 0;
  }
  
}

void HiInclusiveJetAnalyzer::IterativeDeclusteringGen(double groom_type, double groom_combine,const reco::GenJet& jet)
{
  double nsplit = 0;

  double z = 0;
  double zg_SD = 0;
  double zg_latekt = 0;

  double ktg_SD = 0;
  double ktg_latekt = 0;

  double rg_SD = 0;
  double rg_latekt = 0;

  Int_t SD_split = -1; 
  Int_t latekt_split = -1; 

  double jet_radius_ca = 1.0;
  
	
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
    // Reclustering jet constituents with new algorithm
  try{
    std::vector<fastjet::PseudoJet> particles = {};                         
    auto daughters = jet.getJetConstituents();
    // fastjet::PseudoJet tmp_jet;


    for(auto it = daughters.begin(); it!=daughters.end(); ++it){
      // tmp_jet += fastjet::PseudoJet((**it).px(), (**it).py(), (**it).pz(), (**it).energy());
      //      std::cout << "scan gen jet consts, mass: " << (**it).mass() << " charge: " << (**it).charge() << " id: "<< (**it).pdgId() << std::endl;

      //if we want only charged constituents and the daughter charge is 0, skip it
      if (doChargedConstOnly_ && (**it).charge()==0) continue;
    
      if ((**it).mass() < 0 and (**it).charge() < -4) {      // Assumes charge is set in aggregator
	jets_.massHFgen[jets_.nref] = -((**it).mass());
	//	std::cout << "HF mass is " << -((**it).mass()) << std::endl;
      }
	    
    
      if ((**it).pt()<1) continue; //Particle Pt cut
      // std::cout << "scan gen jet consts, mass: " << (**it).mass() << " charge: " << (**it).charge() << std::endl;
      // cout<<"pdg Id = "<< (**it).pdgId()<< ", pt = "<< (**it).pt() << ", eta = "<< (**it).eta() << endl;
      particles.push_back(fastjet::PseudoJet((**it).px(), (**it).py(), (**it).pz(), (**it).energy()));
    }
      //  cout<< "tmp pt = " << tmp_jet.perp() << endl;

    if (particles.size() == 0){ 
      jets_.ref_split_SD[jets_.nref] = std::numeric_limits<int>::min();
      jets_.ref_split_latekt[jets_.nref] = std::numeric_limits<int>::min();
    }

    if (particles.size() != 0 ) {
      fastjet::ClusterSequence csiter(particles, jet_def);
      std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
      output_jets = sorted_by_pt(output_jets);

      fastjet::PseudoJet jj = output_jets[0];
      fastjet::PseudoJet j1;
      fastjet::PseudoJet j2;

      if(!jj.has_parents(j1,j2)) {
        jets_.ref_split_SD[jets_.nref] = std::numeric_limits<int>::min();
        jets_.ref_split_latekt[jets_.nref] = std::numeric_limits<int>::min();
      }

      int stopSD = 0;
      bool flagHF = false, flagHFSD = false, flagHFkt = false;
        
      while (jj.has_parents(j1,j2)) {
        if (j1.perp() < j2.perp()) std::swap(j1,j2); // j1 is the hardest prong
	//        vector <fastjet::PseudoJet> constitj1 = sorted_by_pt(j1.constituents()); // Vector containing j1 costituents

	std::vector<fastjet::PseudoJet> j1constits = j1.constituents();
	for (size_t icon = 0; icon < j1constits.size(); icon++) {
	  fastjet::PseudoJet constit = j1constits[icon];
	  if (abs(constit.m() + jets_.massHF[jets_.nref]) < 1e-3) {
	    flagHF = true;
	    //std::cout << "pt  " << constit.pt() << " mass " << constit.m() << " " << flagHF << std::endl;
	    break;
	  }
	}
	//	std::cout << " ---------- " << std::endl;
	
	double delta_R = j1.delta_R(j2);
        if (doSplitMatching_ && isMC_) {
	  //vangi's
          jets_.refJetSplits.push_back(j2);
        }
        double k_t = j2.perp()*delta_R;
        z = j2.perp()/(j1.perp()+j2.perp());   

      //  std::cout << "Truth split " << nsplit << " with k_T=" << k_t << " z=" << z << " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl; 

        if (((groom_combine == 0) and (groom_type == 1) and (z > SDcut) and (stopSD == 0)) or ((groom_combine == 1) and (z > SDcut) and (stopSD == 0))) { 
          stopSD = 1;
          zg_SD = z;
          rg_SD  = delta_R;
          ktg_SD = k_t;
          SD_split = nsplit;
	  flagHFSD = flagHF;
        }
        
        if (((groom_combine == 0) and (groom_type == 0) and (k_t > latektcut)) or ((groom_combine == 1) and (k_t > latektcut))) {
	  zg_latekt = z;
          rg_latekt  = delta_R;
          ktg_latekt = k_t;
          latekt_split = nsplit;
	  flagHFkt = flagHF;
        }
        jj = j1;
        nsplit = nsplit+1;
      }
      jets_.ref_z_SD[jets_.nref] = zg_SD;
      jets_.ref_rg_SD[jets_.nref] = rg_SD;
      jets_.ref_ktg_SD[jets_.nref] = ktg_SD;
      jets_.ref_split_SD[jets_.nref] = SD_split;
      jets_.ref_hasHF_SD[jets_.nref] = flagHFSD;
    
      jets_.ref_z_latekt[jets_.nref] = zg_latekt;
      jets_.ref_rg_latekt[jets_.nref] = rg_latekt;
      jets_.ref_ktg_latekt[jets_.nref] = ktg_latekt;
      jets_.ref_split_latekt[jets_.nref] = latekt_split;
      jets_.ref_hasHF_latekt[jets_.nref] = flagHFkt;
    }

  }
  catch (fastjet::Error const&){
    cout << "Fastjet error in gen level" << endl;
  }
  catch (Int_t MyNum){
    cout<<"MyNum catch"<<endl;
    jets_.ref_z_SD[jets_.nref] = 0;
    jets_.ref_rg_SD[jets_.nref] = 0;
    jets_.ref_ktg_SD[jets_.nref] = 0;
    jets_.ref_z_latekt[jets_.nref] = 0;
    jets_.ref_rg_latekt[jets_.nref] = 0;
    jets_.ref_ktg_latekt[jets_.nref] = 0;
  }
}

//vangi's
void HiInclusiveJetAnalyzer::RecoTruthSplitMatching(std::vector<fastjet::PseudoJet> &allSplitsLevel1, fastjet::PseudoJet &toMatchLevel2, bool *bool_array, int *splitLevel1){
  float mindR = std::numeric_limits<float>::max();
  size_t closestInLevel1 = 0;
   
  for (size_t i{0}; i < allSplitsLevel1.size(); ++i) {
    float dR = allSplitsLevel1.at(i).delta_R(toMatchLevel2);
    if (mindR > dR) {
      closestInLevel1 = i;
      mindR = dR;
    }
  }
  //  std::cout << "Looped over " <<  allSplitsLevel1.size() << " splittings, mindR was " << mindR << std::endl;
  
  if (static_cast<int>(closestInLevel1) == splitLevel1[jets_.nref] ) bool_array[jets_.nref] = true;
  else bool_array[jets_.nref] = false;
  // std::cout << "bool array at jet: " << bool_array[jets_.nref] << std::endl;
}

void HiInclusiveJetAnalyzer::TruthRecoRecoTruthMatching_latekt(){
  // std::cout << "# reco latekt split: "<<jets_.jt_split_latekt[jets_.nref] << ", # truth latekt split " << jets_.ref_split_latekt[jets_.nref] << std::endl;
  if( (jets_.jt_split_latekt[jets_.nref] == std::numeric_limits<int>::min())
      or (jets_.ref_split_latekt[jets_.nref] == std::numeric_limits<int>::min())
      or (jets_.jtJetSplits.size() == 0)
      or (jets_.refJetSplits.size() == 0)) {
    
    jets_.ref_isClosestToReco_latekt[jets_.nref] = false;
    jets_.jt_isClosestToTruth_latekt[jets_.nref] = false;
    jets_.jt_ref_dR_latekt[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }

  if ( (jets_.jt_split_latekt[jets_.nref] == -1) and (jets_.ref_split_latekt[jets_.nref] == -1) ){
    // std::cout<<"\nuntgged?\n"<<std::endl;
    jets_.ref_isClosestToReco_latekt[jets_.nref] = true;
    jets_.jt_isClosestToTruth_latekt[jets_.nref] = true;
    jets_.jt_ref_dR_latekt[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }
  else if ( ((jets_.jt_split_latekt[jets_.nref] == -1) and (jets_.ref_split_latekt[jets_.nref] != -1)) or ((jets_.jt_split_latekt[jets_.nref] != -1) and (jets_.ref_split_latekt[jets_.nref] == -1)) ){
    // std::cout<<"\nuntgged latekt?\n"<<std::endl;
    jets_.ref_isClosestToReco_latekt[jets_.nref] = false;
    jets_.jt_isClosestToTruth_latekt[jets_.nref] = false;
    jets_.jt_ref_dR_latekt[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }
  else {
    fastjet::PseudoJet latekt_R_split = jets_.jtJetSplits.at(jets_.jt_split_latekt[jets_.nref]);
    fastjet::PseudoJet latekt_T_split = jets_.refJetSplits.at(jets_.ref_split_latekt[jets_.nref]);
    // std::cout <<"Reco: eta = "<< latekt_R_split.eta() << ", phi = " << latekt_R_split.phi() << ", DeltaR = "<< latekt_R_split.delta_R(latekt_T_split) <<" latekt reco  splitting in matching" << std::endl;
 
    // std::cout << "Angle between latekt splits is dR = " << latekt_R_split.delta_R(latekt_T_split) << std::endl;
    // for(size_t i{0};i<jets_.jtJetConstituent.size();++i) std::cout<< " reco "<< i <<" : DeltaR = "<< jets_.jtJetConstituent.at(i).delta_R(latekt_T_split)<<", k_T = "<< jets_.jtJetConstituent.at(i).pt()<<" GeV"<< std::endl;
    jets_.jt_ref_dR_latekt[jets_.nref] = latekt_R_split.delta_R(latekt_T_split);
    RecoTruthSplitMatching(jets_.refJetSplits, latekt_R_split, jets_.ref_isClosestToReco_latekt, jets_.ref_split_latekt);
    RecoTruthSplitMatching(jets_.jtJetSplits,  latekt_T_split, jets_.jt_isClosestToTruth_latekt, jets_.jt_split_latekt);
  }
}

void HiInclusiveJetAnalyzer::TruthRecoRecoTruthMatching_SD(){
  //  std::cout << "# reco SD split : "<<jets_.jt_split_SD[jets_.nref] << ", # truth SD split " << jets_.ref_split_SD[jets_.nref] << " consts " <<  jets_.jtJetSplits.size() << " " << jets_.refJetSplits.size() << std::endl;
  if ((jets_.jt_split_SD[jets_.nref] == std::numeric_limits<int>::min()) or (jets_.ref_split_SD[jets_.nref] == std::numeric_limits<int>::min()) or (jets_.jtJetSplits.size() == 0) or  (jets_.refJetSplits.size() == 0) ) {
    jets_.ref_isClosestToReco_SD[jets_.nref] = false;
    jets_.jt_isClosestToTruth_SD[jets_.nref] = false;
    jets_.jt_ref_dR_SD[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }
  
  if( (jets_.jt_split_SD[jets_.nref] == -1) and (jets_.ref_split_SD[jets_.nref] == -1)) {
    /// std::cout<<"\nuntgged?\n"<<std::endl;
    jets_.ref_isClosestToReco_SD[jets_.nref] = true;
    jets_.jt_isClosestToTruth_SD[jets_.nref] = true;
    jets_.jt_ref_dR_SD[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }
  else if (((jets_.jt_split_SD[jets_.nref] == -1) and (jets_.ref_split_SD[jets_.nref] != -1)) or ((jets_.jt_split_SD[jets_.nref] != -1) and (jets_.ref_split_SD[jets_.nref] == -1))) {
    //   std::cout<<"\nuntgged 2 SD?\n"<<std::endl;
    jets_.ref_isClosestToReco_SD[jets_.nref] = false;
    jets_.jt_isClosestToTruth_SD[jets_.nref] = false;
    jets_.jt_ref_dR_SD[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }
  else {
    //  std::cout<<"Run SD split matching"<<std::endl;
    fastjet::PseudoJet SD_R_split = jets_.jtJetSplits.at(jets_.jt_split_SD[jets_.nref]);
    fastjet::PseudoJet SD_T_split = jets_.refJetSplits.at(jets_.ref_split_SD[jets_.nref]);
    jets_.jt_ref_dR_SD[jets_.nref] = SD_R_split.delta_R(SD_T_split);
    RecoTruthSplitMatching(jets_.refJetSplits, SD_R_split, jets_.ref_isClosestToReco_SD, jets_.ref_split_SD);
    RecoTruthSplitMatching(jets_.jtJetSplits,  SD_T_split, jets_.jt_isClosestToTruth_SD,  jets_.jt_split_SD);
  }
}


int HiInclusiveJetAnalyzer::getPFJetMuon(const pat::Jet& pfJet,
                                         const edm::View<pat::PackedCandidate>* pfCandidateColl) {
  int pfMuonIndex = -1;
  float ptMax = 0.;

  for (unsigned icand = 0; icand < pfCandidateColl->size(); icand++) {
    const pat::PackedCandidate& pfCandidate = pfCandidateColl->at(icand);
    int id = pfCandidate.pdgId();
    if (abs(id) != 3)
      continue;

    if (reco::deltaR(pfJet, pfCandidate) > 0.5)
      continue;

    double pt = pfCandidate.pt();
    if (pt > ptMax) {
      ptMax = pt;
      pfMuonIndex = (int)icand;
    }
  }

  return pfMuonIndex;
}

double HiInclusiveJetAnalyzer::getPtRel(const pat::PackedCandidate& lep, const pat::Jet& jet)

{
  float lj_x = jet.p4().px();
  float lj_y = jet.p4().py();
  float lj_z = jet.p4().pz();

  // absolute values squared
  float lj2 = lj_x * lj_x + lj_y * lj_y + lj_z * lj_z;
  float lep2 = lep.px() * lep.px() + lep.py() * lep.py() + lep.pz() * lep.pz();

  // projection vec(mu) to lepjet axis
  float lepXlj = lep.px() * lj_x + lep.py() * lj_y + lep.pz() * lj_z;

  // absolute value squared and normalized
  float pLrel2 = lepXlj * lepXlj / lj2;

  // lep2 = pTrel2 + pLrel2
  float pTrel2 = lep2 - pLrel2;

  return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetAnalyzer::analyzeSubjets(const reco::Jet& jet) {
  std::vector<float> sjpt;
  std::vector<float> sjeta;
  std::vector<float> sjphi;
  std::vector<float> sjm;
  if (jet.numberOfDaughters() > 0) {
    for (unsigned k = 0; k < jet.numberOfDaughters(); ++k) {
      const reco::Candidate& dp = *jet.daughter(k);
      sjpt.push_back(dp.pt());
      sjeta.push_back(dp.eta());
      sjphi.push_back(dp.phi());
      sjm.push_back(dp.mass());
    }
  } else {
    sjpt.push_back(-999.);
    sjeta.push_back(-999.);
    sjphi.push_back(-999.);
    sjm.push_back(-999.);
  }
  jets_.jtSubJetPt.push_back(sjpt);
  jets_.jtSubJetEta.push_back(sjeta);
  jets_.jtSubJetPhi.push_back(sjphi);
  jets_.jtSubJetM.push_back(sjm);
}

//--------------------------------------------------------------------------------------------------
int HiInclusiveJetAnalyzer::getGroomedGenJetIndex(const reco::GenJet& jet) const {
  //Find closest soft-dropped gen jet
  double drMin = 100;
  int imatch = -1;
  for (unsigned int i = 0; i < gensubjets_->size(); ++i) {
    const reco::Jet& mjet = (*gensubjets_)[i];

    double dr = deltaR(jet, mjet);
    if (dr < drMin) {
      imatch = i;
      drMin = dr;
    }
  }
  return imatch;
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetAnalyzer::analyzeRefSubjets(const reco::GenJet& jet) {
  //Find closest soft-dropped gen jet
  int imatch = getGroomedGenJetIndex(jet);
  double dr = 999.;
  if (imatch > -1) {
    const reco::Jet& mjet = (*gensubjets_)[imatch];
    dr = deltaR(jet, mjet);
  }

  jets_.refptG[jets_.nref] = -999.;
  jets_.refetaG[jets_.nref] = -999.;
  jets_.refphiG[jets_.nref] = -999.;
  jets_.refmG[jets_.nref] = -999.;
  jets_.refsym[jets_.nref] = -999.;
  jets_.refdroppedBranches[jets_.nref] = -999;

  std::vector<float> sjpt;
  std::vector<float> sjeta;
  std::vector<float> sjphi;
  std::vector<float> sjm;
  if (imatch > -1 && dr < 0.4) {
    const reco::Jet& mjet = (*gensubjets_)[imatch];
    jets_.refptG[jets_.nref] = mjet.pt();
    jets_.refetaG[jets_.nref] = mjet.eta();
    jets_.refphiG[jets_.nref] = mjet.phi();
    jets_.refmG[jets_.nref] = mjet.mass();

    if (mjet.numberOfDaughters() > 0) {
      for (unsigned k = 0; k < mjet.numberOfDaughters(); ++k) {
        const reco::Candidate& dp = *mjet.daughter(k);
        sjpt.push_back(dp.pt());
        sjeta.push_back(dp.eta());
        sjphi.push_back(dp.phi());
        sjm.push_back(dp.mass());
      }
    }
    if (doGenSym_) {
      Ptr<reco::Jet> genJetPtr = gensubjets_->ptrAt(imatch);
      float gensym = (*genSymVM_)[genJetPtr];
      jets_.refsym[jets_.nref] = gensym;
      int db = (*genDroppedBranchesVM_)[genJetPtr];
      jets_.refdroppedBranches[jets_.nref] = db;
    }
  } else {
    jets_.refptG[jets_.nref] = -999.;
    jets_.refetaG[jets_.nref] = -999.;
    jets_.refphiG[jets_.nref] = -999.;
    jets_.refmG[jets_.nref] = -999.;

    sjpt.push_back(-999.);
    sjeta.push_back(-999.);
    sjphi.push_back(-999.);
    sjm.push_back(-999.);
  }

  jets_.refSubJetPt.push_back(sjpt);
  jets_.refSubJetEta.push_back(sjeta);
  jets_.refSubJetPhi.push_back(sjphi);
  jets_.refSubJetM.push_back(sjm);
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetAnalyzer::analyzeGenSubjets(const reco::GenJet& jet) {
  //Find closest soft-dropped gen jet
  int imatch = getGroomedGenJetIndex(jet);
  double dr = 999.;
  if (imatch > -1) {
    const reco::Jet& mjet = (*gensubjets_)[imatch];
    dr = deltaR(jet, mjet);
  }

  jets_.genptG[jets_.ngen] = -999.;
  jets_.genetaG[jets_.ngen] = -999.;
  jets_.genphiG[jets_.ngen] = -999.;
  jets_.genmG[jets_.ngen] = -999.;
  jets_.gensym[jets_.ngen] = -999.;
  jets_.gendroppedBranches[jets_.ngen] = -999;

  std::vector<float> sjpt;
  std::vector<float> sjeta;
  std::vector<float> sjphi;
  std::vector<float> sjm;
  std::vector<float> sjarea;
  if (imatch > -1 && dr < 0.4) {
    const reco::Jet& mjet = (*gensubjets_)[imatch];
    jets_.genptG[jets_.ngen] = mjet.pt();
    jets_.genetaG[jets_.ngen] = mjet.eta();
    jets_.genphiG[jets_.ngen] = mjet.phi();
    jets_.genmG[jets_.ngen] = mjet.mass();

    if (mjet.numberOfDaughters() > 0) {
      for (unsigned k = 0; k < mjet.numberOfDaughters(); ++k) {
        const reco::Candidate& dp = *mjet.daughter(k);
        sjpt.push_back(dp.pt());
        sjeta.push_back(dp.eta());
        sjphi.push_back(dp.phi());
        sjm.push_back(dp.mass());
        //sjarea.push_back(dp.castTo<reco::JetRef>()->jetArea());
      }
    }
    if (doGenSym_) {
      Ptr<reco::Jet> genJetPtr = gensubjets_->ptrAt(imatch);
      float gensym = (*genSymVM_)[genJetPtr];
      jets_.gensym[jets_.ngen] = gensym;
      int db = (*genDroppedBranchesVM_)[genJetPtr];
      jets_.gendroppedBranches[jets_.ngen] = db;
    }
  } else {
    jets_.genptG[jets_.ngen] = -999.;
    jets_.genetaG[jets_.ngen] = -999.;
    jets_.genphiG[jets_.ngen] = -999.;
    jets_.genmG[jets_.ngen] = -999.;

    sjpt.push_back(-999.);
    sjeta.push_back(-999.);
    sjphi.push_back(-999.);
    sjm.push_back(-999.);
    sjarea.push_back(-999.);
  }

  jets_.genSubJetPt.push_back(sjpt);
  jets_.genSubJetEta.push_back(sjeta);
  jets_.genSubJetPhi.push_back(sjphi);
  jets_.genSubJetM.push_back(sjm);
  jets_.genSubJetArea.push_back(sjarea);
}

DEFINE_FWK_MODULE(HiInclusiveJetAnalyzer);
