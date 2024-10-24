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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
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
  matchTag_ = consumes<pat::JetCollection>(iConfig.getUntrackedParameter<InputTag>("matchTag"));

  useQuality_ = iConfig.getUntrackedParameter<bool>("useQuality", true);
  trackQuality_ = iConfig.getUntrackedParameter<string>("trackQuality", "highPurity");

  jetName_ = iConfig.getUntrackedParameter<string>("jetName");
  doGenTaus_ = iConfig.getUntrackedParameter<bool>("doGenTaus", false);
  doGenSym_ = iConfig.getUntrackedParameter<bool>("doGenSym", false);
  doSubJets_ = iConfig.getUntrackedParameter<bool>("doSubJets", false);
  doJetConstituents_ = iConfig.getUntrackedParameter<bool>("doJetConstituents", false);
  doGenSubJets_ = iConfig.getUntrackedParameter<bool>("doGenSubJets", false);
  if (doGenSubJets_)
    subjetGenTag_ = consumes<reco::JetView>(iConfig.getUntrackedParameter<InputTag>("subjetGenTag"));

  //reWTA reclustering
  doWTARecluster_ = iConfig.getUntrackedParameter<bool>("doWTARecluster", false);

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
  }
  useRawPt_ = iConfig.getUntrackedParameter<bool>("useRawPt", true);

  doLegacyBtagging_ = iConfig.getUntrackedParameter<bool>("doLegacyBtagging", true);
  doCandidateBtagging_ = iConfig.getUntrackedParameter<bool>("doCandidateBtagging", false);

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
  if (doCandidateBtagging_) {
    pnetBvsAllJetTags_ = "pfParticleNetAK4DiscriminatorsJetTags:BvsAll";
    pnetProbBJetTags_ = "pfParticleNetAK4JetTags:probb";
    pnetProbBBJetTags_ = "pfParticleNetAK4JetTags:probbb";
    pnetProbCJetTags_ = "pfParticleNetAK4JetTags:probc";
    pnetProbCCJetTags_ = "pfParticleNetAK4JetTags:probcc";
    pnetProbGJetTags_ = "pfParticleNetAK4JetTags:probg";
    pnetProbUDSJetTags_ = "pfParticleNetAK4JetTags:probuds";
    pnetProbPUJetTags_ = "pfParticleNetAK4JetTags:probpu";
    pnetProbUNDEFJetTags_ = "pfParticleNetAK4JetTags:probundef";
    pfJPJetTags_ = jetName_ + "pfJetProbabilityBJetTags";
  }
  doSubEvent_ = false;

  if (isMC_) {
    genPtMin_ = iConfig.getUntrackedParameter<double>("genPtMin", 10);
    doSubEvent_ = iConfig.getUntrackedParameter<bool>("doSubEvent", false);
  }

  //  ptCut = iConfig.getParameter<double>("ptCut",1.);
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

void HiInclusiveJetAnalyzer::beginJob() {
  string jetTagTitle = jetTagLabel_.label() + " Jet Analysis Tree";
  t = fs1->make<TTree>("t", jetTagTitle.c_str());

  t->Branch("run", &jets_.run, "run/I");
  t->Branch("evt", &jets_.evt, "evt/I");
  t->Branch("lumi", &jets_.lumi, "lumi/I");
  t->Branch("nref", &jets_.nref, "nref/I");
  t->Branch("rawpt", jets_.rawpt, "rawpt[nref]/F");
  t->Branch("jtpt", jets_.jtpt, "jtpt[nref]/F");
  t->Branch("jteta", jets_.jteta, "jteta[nref]/F");
  t->Branch("jty", jets_.jty, "jty[nref]/F");
  t->Branch("jtphi", jets_.jtphi, "jtphi[nref]/F");
  t->Branch("jtpu", jets_.jtpu, "jtpu[nref]/F");
  t->Branch("jtm", jets_.jtm, "jtm[nref]/F");
  t->Branch("jtarea", jets_.jtarea, "jtarea[nref]/F");

  t->Branch("massHF", jets_.massHF, "massHF[nref]/F");

  t->Branch("jtdyn_split",jets_.jtdyn_split,"jtdyn_split[nref]/I");
  // t->Branch("jtdyn_theta",jets_.jtdyn_theta,"jtdyn_theta[nref]/F");
  t->Branch("jtdyn_deltaR",jets_.jtdyn_deltaR,"jtdyn_deltaR[nref]/F");
  t->Branch("jtdyn_kt",jets_.jtdyn_kt,"jtdyn_kt[nref]/F");
  t->Branch("jtdyn_z",jets_.jtdyn_z,"jtdyn_z[nref]/F");
  t->Branch("jt_intjet_multi", jets_.jt_intjet_multi,"jt_intjet_multi[nref]/I");
  t->Branch("jt_girth", jets_.jt_girth, "jt_girth[nref]/F");
  
  
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
    t->Branch("matchedPt", jets_.matchedPt, "matchedPt[nref]/F");
    t->Branch("matchedRawPt", jets_.matchedRawPt, "matchedRawPt[nref]/F");
    t->Branch("matchedPu", jets_.matchedPu, "matchedPu[nref]/F");
    t->Branch("matchedR", jets_.matchedR, "matchedR[nref]/F");
    if (isMC_) {
      t->Branch("matchedHadronFlavor", jets_.matchedHadronFlavor, "matchedHadronFlavor[nref]/I");
      t->Branch("matchedPartonFlavor", jets_.matchedPartonFlavor, "matchedPartonFlavor[nref]/I");
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
  if(doCandidateBtagging_){
    t->Branch("discr_pnetBvsAll", jets_.discr_pnetBvsAll, "discr_pnetBvsAll[nref]/F");
    t->Branch("discr_pnetProbB", jets_.discr_pnetProbB, "discr_pnetProbB[nref]/F");
    t->Branch("discr_pnetProbBB", jets_.discr_pnetProbBB, "discr_pnetProbBB[nref]/F");
    t->Branch("discr_pnetProbC", jets_.discr_pnetProbC, "discr_pnetProbC[nref]/F");
    t->Branch("discr_pnetProbCC", jets_.discr_pnetProbCC, "discr_pnetProbCC[nref]/F");
    t->Branch("discr_pnetProbG", jets_.discr_pnetProbG, "discr_pnetProbG[nref]/F");
    t->Branch("discr_pnetProbUDS", jets_.discr_pnetProbUDS, "discr_pnetProbUDS[nref]/F");
    t->Branch("discr_pnetProbPU", jets_.discr_pnetProbPU, "discr_pnetProbPU[nref]/F");
    t->Branch("discr_pnetProbUNDEF", jets_.discr_pnetProbUNDEF, "discr_pnetProbUNDEF[nref]/F");
    t->Branch("discr_pfJP", jets_.discr_pfJP, "discr_pfJP[nref]/F");
    }
  if (isMC_) {
    if (useHepMC_) {
      t->Branch("beamId1", &jets_.beamId1, "beamId1/I");
      t->Branch("beamId2", &jets_.beamId2, "beamId2/I");
    }

    t->Branch("pthat", &jets_.pthat, "pthat/F");

    // Only matched gen jets
    t->Branch("refpt", jets_.refpt, "refpt[nref]/F");
    t->Branch("refeta", jets_.refeta, "refeta[nref]/F");
    t->Branch("refy", jets_.refy, "refy[nref]/F");
    t->Branch("refphi", jets_.refphi, "refphi[nref]/F");
    t->Branch("refm", jets_.refm, "refm[nref]/F");
    t->Branch("refarea", jets_.refarea, "refarea[nref]/F");

    // Simple subjet stuff
    t->Branch("refdyn_split",jets_.refdyn_split,"refdyn_split[nref]/I");
    // t->Branch("refdyn_theta",jets_.refdyn_theta,"refdyn_theta[nref]/F");
    t->Branch("refdyn_deltaR",jets_.refdyn_deltaR,"refdyn_deltaR[nref]/F");

    t->Branch("refdyn_kt",jets_.refdyn_kt,"refdyn_kt[nref]/F");
    t->Branch("refdyn_z",jets_.refdyn_z,"refdyn_z[nref]/F");
    t->Branch("ref_intjet_multi", jets_.ref_intjet_multi,"ref_intjet_multi[nref]/I");
    t->Branch("ref_girth", jets_.ref_girth, "ref_girth[nref]/F");
    t->Branch("jtdyn_isClosestToTruth", jets_.jtdyn_isClosestToTruth,"jtdyn_isClosestToTruth[nref]/O");
    t->Branch("refdyn_isClosestToReco", jets_.refdyn_isClosestToReco,"refdyn_isClosestToReco[nref]/O");
    t->Branch("jtdyn_refdyn_dR", jets_.jtdyn_refdyn_dR,"jtdyn_refdyn_dR[nref]/F");
    /////////////////

    
    t->Branch("refdphijt", jets_.refdphijt, "refdphijt[nref]/F");
    t->Branch("refdrjt", jets_.refdrjt, "refdrjt[nref]/F");
    // matched parton
    t->Branch("refparton_pt", jets_.refparton_pt, "refparton_pt[nref]/F");
    t->Branch("refparton_flavor", jets_.refparton_flavor, "refparton_flavor[nref]/I");
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

    t->Branch("genChargedSum", jets_.genChargedSum, "genChargedSum[nref]/F");
    t->Branch("genHardSum", jets_.genHardSum, "genHardSum[nref]/F");
    t->Branch("signalChargedSum", jets_.signalChargedSum, "signalChargedSum[nref]/F");
    t->Branch("signalHardSum", jets_.signalHardSum, "signalHardSum[nref]/F");

    if (doSubEvent_) {
      t->Branch("subid", jets_.subid, "subid[nref]/I");
    }

    if (fillGenJets_) {
      // For all gen jets, matched or unmatched
      t->Branch("ngen", &jets_.ngen, "ngen/I");
      t->Branch("genmatchindex", jets_.genmatchindex, "genmatchindex[ngen]/I");
      t->Branch("genpt", jets_.genpt, "genpt[ngen]/F");
      t->Branch("geneta", jets_.geneta, "geneta[ngen]/F");
      t->Branch("geny", jets_.geny, "geny[ngen]/F");
      t->Branch("genphi", jets_.genphi, "genphi[ngen]/F");
      t->Branch("genm", jets_.genm, "genm[ngen]/F");
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
    
    t->Branch("ntrkInSvtxNotInJet", &jets_.ntrkInSvtxNotInJet, "ntrkInSvtxNotInJet/I");
    t->Branch("trkInSvtxNotInJetSvId", jets_.trkInSvtxNotInJetSvId, "trkInSvtxNotInJetSvId[ntrkInSvtxNotInJet]/I");
    t->Branch("trkInSvtxNotInJetOtherJetId", jets_.trkInSvtxNotInJetOtherJetId, "trkInSvtxNotInJetOtherJetId[ntrkInSvtxNotInJet]/I");
    t->Branch("trkInSvtxNotInJetMatchSta", jets_.trkInSvtxNotInJetMatchSta, "trkInSvtxNotInJetMatchSta[ntrkInSvtxNotInJet]/I");
    t->Branch("trkInSvtxNotInJetPt", jets_.trkInSvtxNotInJetPt, "trkInSvtxNotInJetPt[ntrkInSvtxNotInJet]/F");
    t->Branch("trkInSvtxNotInJetEta", jets_.trkInSvtxNotInJetEta, "trkInSvtxNotInJetEta[ntrkInSvtxNotInJet]/F");
    t->Branch("trkInSvtxNotInJetPhi", jets_.trkInSvtxNotInJetPhi, "trkInSvtxNotInJetPhi[ntrkInSvtxNotInJet]/F");
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
  if (doCandidateBtagging_) {
    memset(jets_.discr_pnetBvsAll, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbB, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbBB, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbC, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbCC, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbG, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbUDS, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbPU, 0, MAXJETS * sizeof(float));
    memset(jets_.discr_pnetProbUNDEF, 0, MAXJETS * sizeof(float));
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

  edm::Handle<pat::JetCollection> matchedjets;
  iEvent.getByToken(matchTag_, matchedjets);


  if (doGenSubJets_)
    iEvent.getByToken(subjetGenTag_, gensubjets_);
  if (doGenSym_) {
    iEvent.getByToken(tokenGenSym_, genSymVM_);
    iEvent.getByToken(tokenGenDroppedBranches_, genDroppedBranchesVM_);
  }

  edm::Handle<edm::View<pat::PackedCandidate>> pfCandidates;
  iEvent.getByToken(pfCandidateLabel_, pfCandidates);

  
  edm::Handle<reco::TrackToGenParticleMap> trackToGenParticleMap;
  //  edm::Handle<reco::TrackToGenParticleMap> genConstitToGenParticleMap;
  if (isMC_ and doTracks_) {
    edm::Handle<reco::GenParticleCollection> genparts;
    iEvent.getByToken(genParticleSrc_, genparts);

    // Track-gen ptcl tagging
    iEvent.getByToken(trackToGenParticleMapToken_, trackToGenParticleMap);
    //    std::cout << "TracktoGen Map size: " << trackToGenParticleMap->size() << std::endl;
  }

  iEvent.getByToken(primaryVerticesToken_, primaryVertices);
  // FILL JRA TREE
  jets_.nref = 0;
  jets_.nvtx = primaryVertices->size();
  /*  std::cout << "PV:" <<  primaryVertices->size() << std::endl;
      std::cout << "PV:" <<  primaryVertices->at(0).position() << std::endl; */

  //  std::cout << "Number of tracks: " << jets_.ntrk << std::endl; 
  if (doTracks_) jets_.ntrk = 0;
  if (doSvtx_) {
    jets_.nsvtx = 0;
    jets_.ntrkInSvtxNotInJet = 0;
  }
  int nsvtxCounterForTracks = 0;

  // For aggregation
  auto pseudoHFCollection = std::make_unique<std::vector<reco::PFCandidate>>();
  auto droppedTrackCollection = std::make_unique<std::vector<reco::PFCandidate>>();
  
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

  for (unsigned int j = 0; j < jets->size(); ++j) {
    const pat::Jet& jet = (*jets)[j];
    
    auto pt = useRawPt_ ? jet.correctedJet("Uncorrected").pt() : jet.pt();
    if (pt < jetPtMin_)
      continue;
    if (std::abs(jet.eta()) > jetAbsEtaMax_)
      continue;


    ///////////////
    // DEBUG SVTX:
    if (doSvtx_) {
      //      std::cout << "Try getting taginfo with label " << svTagInfoLabel_.c_str() << endl;
      const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());
      // const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex();
      if (jet.hasTagInfo(svTagInfoLabel_.c_str())) {
	int nsv = svTagInfo->nVertices();
	std::cout << "got secondary vtx info, nsv = " << svTagInfo->nVertices() << " " << svTagInfo->nVertexTracks() << std::endl;
      
	jets_.nsvtx += nsv;
	jets_.jtNsvtx[jets_.nref] = nsv;

	if (nsv > 0) std::cout << "jet has nsv = " << nsv << std::endl;
      }
    }
    /////////////

    if (doCandidateBtagging_){
      jets_.discr_pnetBvsAll[jets_.nref]=jet.bDiscriminator(pnetBvsAllJetTags_);
      jets_.discr_pnetProbB[jets_.nref]=jet.bDiscriminator(pnetProbBJetTags_);
      jets_.discr_pnetProbBB[jets_.nref]=jet.bDiscriminator(pnetProbBBJetTags_);
      jets_.discr_pnetProbC[jets_.nref]=jet.bDiscriminator(pnetProbCJetTags_);
      jets_.discr_pnetProbCC[jets_.nref]=jet.bDiscriminator(pnetProbCCJetTags_);
      jets_.discr_pnetProbG[jets_.nref]=jet.bDiscriminator(pnetProbGJetTags_);
      jets_.discr_pnetProbUDS[jets_.nref]=jet.bDiscriminator(pnetProbUDSJetTags_);
      jets_.discr_pnetProbPU[jets_.nref]=jet.bDiscriminator(pnetProbPUJetTags_);
      jets_.discr_pnetProbUNDEF[jets_.nref]=jet.bDiscriminator(pnetProbUNDEFJetTags_);
      jets_.discr_pfJP[jets_.nref]=jet.bDiscriminator(pfJPJetTags_);
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

    // std::cout << "New jet with nref " << jets_.nref << std::endl;
    // std::cout << svTagInfos_ << std::endl;
    // std::cout << "start of svtx" << std::endl;
    if (doSvtx_ && jet.hasTagInfo(svTagInfoLabel_.c_str())) {
      const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());
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

        //mCorr=srqt(m^2+p^2sin^2(th)) + p*sin(th)
        double sinth = svtx.p4().Vect().Unit().Cross((svTagInfo->flightDirection(isv)).unit()).Mag2();
        sinth = sqrt(sinth);
        double underRoot = std::pow(svtxM, 2) + (std::pow(svtxPt, 2) * std::pow(sinth, 2));
        double svtxMcorr = std::sqrt(underRoot) + (svtxPt * sinth);

        jets_.svtxnormchi2[ijetSvtx] = normalizedChi2;
        jets_.svtxm[ijetSvtx] = svtxM;
        jets_.svtxmcorr[ijetSvtx] = svtxMcorr;
        jets_.svtxpt[ijetSvtx] = svtxPt;

        jets_.svtxJetId[ijetSvtx] = jets_.nref;

        const std::vector<reco::CandidatePtr> svTracks = svTagInfo->vertexTracks(isv);
        // std::cout << "nsvtracks = " << svTracks.size() << std::endl;
        for (auto svTrk : svTracks) {
          // look for track in this jet's tracks 
          const std::vector<reco::CandidatePtr> jetTracks = jet.getJetConstituents();
          auto itJetTrack = std::find(jetTracks.begin(), jetTracks.end(), svTrk);
          if (itJetTrack == jetTracks.end()) {
            jets_.trkInSvtxNotInJetSvId[jets_.ntrkInSvtxNotInJet] = ijetSvtx;
            jets_.trkInSvtxNotInJetPt[jets_.ntrkInSvtxNotInJet] = svTrk->pt();
            jets_.trkInSvtxNotInJetEta[jets_.ntrkInSvtxNotInJet] = svTrk->eta();
            jets_.trkInSvtxNotInJetPhi[jets_.ntrkInSvtxNotInJet] = svTrk->phi();

            // look for track in the other jets in the event
            jets_.trkInSvtxNotInJetOtherJetId[jets_.ntrkInSvtxNotInJet] = -1;
            for (unsigned int ji = 0; ji < jets->size(); ++ji) {
              if (ji == j) continue; // skip the same jet
              const pat::Jet& tempJet = (*jets)[ji];
              const std::vector<reco::CandidatePtr> tempJetTracks = tempJet.getJetConstituents();
              auto itTempJetTrack = std::find(tempJetTracks.begin(), tempJetTracks.end(), svTrk);
              if (itTempJetTrack == tempJetTracks.end()) {
                continue;
              }
              jets_.trkInSvtxNotInJetOtherJetId[jets_.ntrkInSvtxNotInJet] = ji;
              break;
            } // other jet loop

            // get status of track 
            int status = -1;
            if (isMC_ && trackToGenParticleMap->find(svTrk) != trackToGenParticleMap->end()) { 
	      edm::Ptr<pat::PackedGenParticle> matchGenParticle = trackToGenParticleMap->at(svTrk);
              status = matchGenParticle->status();
            }
            jets_.trkInSvtxNotInJetMatchSta[jets_.ntrkInSvtxNotInJet] = status;

            // std::cout << "jets_.ntrkInSvtxNotInJet = " << jets_.ntrkInSvtxNotInJet << std::endl;
            jets_.ntrkInSvtxNotInJet++;
          }
        } // sv track loop
      } // sv loop
      jets_.jtNsvtx[jets_.nref] = nsv;
      jets_.nsvtx += nsv;
    } // endif doSvtx_
    std::cout << "end of svtx" << std::endl;
    //std::cout << "Jet has tag infos: " << std::endl;
    for (auto label : jet.tagInfoLabels())  std::cout << label << std::endl;
    

    //cout << "New jet and its tracks, pt jet and uncorrected " << jet.pt() << " " << jet.correctedJet("Uncorrected").pt()  << endl;
    if (doTracks_ && jet.hasTagInfo(ipTagInfoLabel_.c_str())) {
      jets_.jtNtrk[jets_.nref] = 0;
      jets_.jtptCh[jets_.nref] = 0.;
      //      std::cout << "Has IP tag info" << std::endl;
      const reco::CandIPTagInfo *ipTagInfo = jet.tagInfoCandIP(ipTagInfoLabel_.c_str());
      const std::vector<reco::btag::TrackIPData> ipData = ipTagInfo->impactParameterData();
      const std::vector<reco::CandidatePtr> ipTracks = ipTagInfo->selectedTracks();

      reco::Candidate::PolarLorentzVector chJet(0., 0., 0., 0.);

      // For debugging
      // for (auto itIPTrack : ipTracks) {
      //	  std::cout << " ip track: " << itIPTrack->pt() << " " <<  itIPTrack->eta() << " " <<  itIPTrack->phi() << " " <<  std::endl;
      // 	}
	
      //  float ptcounter = 0;
      for (const reco::CandidatePtr constit : jet.getJetConstituents()) {
        // std::cout << "new jet constit with pt, eta, phi " << constit->pt() << " "  << constit->eta() << " "  << constit->phi() << " " << constit->charge() <<  std::endl;
	//	ptcounter += constit->pt();
        if (constit->charge() == 0) continue;
        if (constit->pt() < trkPtCut_) continue;

	// Find IPTrack that matches to a jet constitute track
        reco::CandidatePtr itIPTrack;

	int itrk = -1; // Counter for track index, probably there is something smarter to do. Used for fetching track IP data.
	for (auto iterIPTrack : ipTracks) {
	  float eps = 1e-5;                                                                                                                                                                                          itrk++;
	  
          if (std::abs(constit->eta()-iterIPTrack->eta())>eps) continue;                                                                                                                                   
          else if (std::abs(constit->phi()-iterIPTrack->phi())>eps) continue;                                                                                                                                        else if (std::abs(constit->pt()-iterIPTrack->pt())>eps) continue;                                                                                                                                       
          itIPTrack = iterIPTrack;
	}

	if (!(itIPTrack)) continue; 
	
        // TODO
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
        if (doSvtx_ && jet.hasTagInfo(svTagInfoLabel_.c_str())) {
          const reco::CandSecondaryVertexTagInfo *svTagInfo = jet.tagInfoCandSecondaryVertex(svTagInfoLabel_.c_str());
          int nsv = svTagInfo->nVertices();
          for (int isv = 0; isv < nsv; isv++) {
            const std::vector<reco::CandidatePtr> svTracks = svTagInfo->vertexTracks(isv);
            auto itSVTrack = std::find(svTracks.begin(), svTracks.end(), constit); // TODO: replace
            if (itSVTrack == svTracks.end()) continue;
            jets_.trkSvtxId[ijetTrack] = nsvtxCounterForTracks + isv;
          } // end sv loop for tracks
        } // end doSvtx_

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

      for (unsigned int icand = 0; icand < pfCandidates->size(); ++icand) {
        const pat::PackedCandidate& track = (*pfCandidates)[icand];
        double dr = deltaR(jet, track);
        if (dr < rParam) {
          double ptcand = track.pt();
          int pfid = track.pdgId();

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

    if (doMatch_) {
      // Alternative reconstruction matching (PF for calo, calo for PF)

      double drMin = 100;
      for (unsigned int imatch = 0; imatch < matchedjets->size(); ++imatch) {
        const pat::Jet& mjet = (*matchedjets)[imatch];

        double dr = deltaR(jet, mjet);
        if (dr < drMin) {
          jets_.matchedPt[jets_.nref] = mjet.pt();

          jets_.matchedRawPt[jets_.nref] = mjet.correctedJet("Uncorrected").pt();
          jets_.matchedPu[jets_.nref] = mjet.pileup();
          if (isMC_) {
            jets_.matchedHadronFlavor[jets_.nref] = mjet.hadronFlavour();
            jets_.matchedPartonFlavor[jets_.nref] = mjet.partonFlavour();
          }

          jets_.matchedR[jets_.nref] = dr;
          drMin = dr;
        }
      }
    }

    jets_.rawpt[jets_.nref] = jet.correctedJet("Uncorrected").pt();
    jets_.jtpt[jets_.nref] = jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();
    jets_.jty[jets_.nref] = jet.rapidity();
    jets_.jtpu[jets_.nref] = jet.pileup();
    jets_.jtm[jets_.nref] = jet.mass();
    jets_.jtarea[jets_.nref] = jet.jetArea();

    jets_.jtsym[jets_.nref] = -999.;
    jets_.jtdroppedBranches[jets_.nref] = -999;

   ///// substructure stuff

    jets_.jtdyn_split[jets_.nref] = 0;
    // jets_.jtdyn_theta[jets_.nref] = 0;
    jets_.jtdyn_kt[jets_.nref] = 0;
    jets_.jtdyn_z[jets_.nref] = 0;

    fastjet::PseudoJet *sub1Gen = new fastjet::PseudoJet();
    fastjet::PseudoJet *sub2Gen = new fastjet::PseudoJet();
    fastjet::PseudoJet *sub1Hyb = new fastjet::PseudoJet();
    fastjet::PseudoJet *sub2Hyb = new fastjet::PseudoJet();
    //IterativeDeclusteringRec(groom_type, groom_combine, jet, sub1Hyb, sub2Hyb);
        IterativeDeclusteringRec(0.0, 0.0, jet, sub1Hyb, sub2Hyb);
    // IterativeDeclustering(jetConstituents, pseudoHF);
   
    jets_.refpt[jets_.nref] = 0;
    jets_.refeta[jets_.nref] = 0;
    jets_.refphi[jets_.nref] = 0;
    jets_.refsym[jets_.nref] = 0.;
        // jets_.refrg[jets_.nref] = 0;
        // jets_.refdynkt[jets_.nref] = 0;
    // jets_.refangu[jets_.nref] = 0; 
        // jets_.refdyn_pt1[jets_.nref] = 0;
    // jets_.refdyn_var[jets_.nref] = 0;
    jets_.refdyn_split[jets_.nref] = 0;
    // jets_.refdyn_theta[jets_.nref] = 0;
    jets_.refdyn_kt[jets_.nref] = 0;
    jets_.refdyn_z[jets_.nref] = 0;
    jets_.jtdyn_isClosestToTruth[jets_.nref] = 0;
    jets_.refdyn_isClosestToReco[jets_.nref] = 0;
    jets_.jtdyn_refdyn_dR[jets_.nref] = 0;
    jets_.refsub11[jets_.nref] = 0;
    jets_.refsub12[jets_.nref] = 0;
    jets_.refsub21[jets_.nref] = 0;
    jets_.refsub22[jets_.nref] = 0;
    // std::cout << jets_.jtJetConstituent.size() << " " << jets_.refJetConstituent.size() << " sizes of consts" << std::endl; 

    
    if (isMC_){
      const reco::GenJet * genjet = jet.genJet();
      if(!genjet) continue;

      //TODO: correct map has to be debugged
      
      if(jet.genParton()){
        const reco::GenParticle & parton = *jet.genParton();
        jets_.refparton_pt[jets_.nref] = parton.pt();
        jets_.refparton_flavor[jets_.nref] = parton.pdgId();
      }
      else {
        jets_.refparton_pt[jets_.nref] = -999;
        jets_.refparton_flavor[jets_.nref] = -999;
      }
      jets_.refpt[jets_.nref] = genjet->pt();
      jets_.refeta[jets_.nref] = genjet->eta();
      jets_.refphi[jets_.nref] = genjet->phi();

      //cout<<"jet daughters gen"<<genjet->numberOfDaughters()<<endl;
      //      if(dopthatcut) if(pthat<0.35*genjet->pt()) continue;

      //      IterativeDeclusteringGen(groom_type, groom_combine, *genjet, sub1Gen, sub2Gen);
      // Todo: aggregation outputs into declustering
      IterativeDeclusteringGen(0.0, 0.0, *genjet, sub1Gen, sub2Gen);
      
    }
    delete sub1Gen;
    delete sub2Gen;
     delete sub1Hyb;
    delete sub2Hyb;
    /////////////////////
    
    if (doSubJets_)
      analyzeSubjets(jet);

    if (jet.hasUserFloat(jetName_ + "Jets:sym"))
      jets_.jtsym[jets_.nref] = jet.userFloat(jetName_ + "Jets:sym");
    if (jet.hasUserInt(jetName_ + "Jets:droppedBranches"))
      jets_.jtdroppedBranches[jets_.nref] = jet.userInt(jetName_ + "Jets:droppedBranches");

    /*    if (jet.isPFJet()) {
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
      } */

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

    if (isMC_) {
      const reco::GenJet* genjet = jet.genJet();

      if (genjet) {
        jets_.refpt[jets_.nref] = genjet->pt();
        jets_.refeta[jets_.nref] = genjet->eta();
        jets_.refphi[jets_.nref] = genjet->phi();
        jets_.refm[jets_.nref] = genjet->mass();
        jets_.refarea[jets_.nref] = genjet->jetArea();
        jets_.refy[jets_.nref] = genjet->eta();
        jets_.refdphijt[jets_.nref] = reco::deltaPhi(jet.phi(), genjet->phi());
        jets_.refdrjt[jets_.nref] = reco::deltaR(jet.eta(), jet.phi(), genjet->eta(), genjet->phi());

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

      jets_.refparton_flavorForB[jets_.nref] = jet.partonFlavour();

      //      if(jet.genParton()){
      // // matched partons
      // const reco::GenParticle & parton = *jet.genParton();

      // jets_.refparton_pt[jets_.nref] = parton.pt();
      // jets_.refparton_flavor[jets_.nref] = parton.pdgId();

      //      } else {
      jets_.refparton_pt[jets_.nref] = -999;
      jets_.refparton_flavor[jets_.nref] = -999;
      //      }
    }
    jets_.nref++;
  }


  if (isMC_) {
    if (useHepMC_) {
      edm::Handle<HepMCProduct> hepMCProduct;
      iEvent.getByToken(eventInfoTag_, hepMCProduct);
      const HepMC::GenEvent* MCEvt = hepMCProduct->GetEvent();

      std::pair<HepMC::GenParticle*, HepMC::GenParticle*> beamParticles = MCEvt->beam_particles();
      jets_.beamId1 = (beamParticles.first != 0) ? beamParticles.first->pdg_id() : 0;
      jets_.beamId2 = (beamParticles.second != 0) ? beamParticles.second->pdg_id() : 0;
    }

    edm::Handle<GenEventInfoProduct> hEventInfo;
    iEvent.getByToken(eventGenInfoTag_, hEventInfo);

    // binning values and qscale appear to be equivalent, but binning values not always present
    jets_.pthat = hEventInfo->qScale();

    edm::Handle<edm::View<reco::GenJet>> genjets;
    iEvent.getByToken(genjetTag_, genjets);

    jets_.ngen = 0;

    for (unsigned int igen = 0; igen < genjets->size(); ++igen) {
      const reco::GenJet& genjet = (*genjets)[igen];
      float genjet_pt = genjet.pt();

      Ptr<reco::GenJet> genJetPtr = genjets->ptrAt(igen);

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
          
          break;
        }
      }

      // threshold to reduce size of output in minbias PbPb
      if (genjet_pt > genPtMin_) {
        jets_.genpt[jets_.ngen] = genjet_pt;
        jets_.geneta[jets_.ngen] = genjet.eta();
        jets_.genphi[jets_.ngen] = genjet.phi();
        jets_.genm[jets_.ngen] = genjet.mass();
        jets_.geny[jets_.ngen] = genjet.eta();

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

  t->Fill();

  //memset(&jets_,0,sizeof jets_);
  jets_ = {0};
}



void HiInclusiveJetAnalyzer::IterativeDeclusteringRec(double groom_type, double groom_combine, const reco::Jet& jet, fastjet::PseudoJet *sub1, fastjet::PseudoJet *sub2)
{
  //  TRandom *r1 = new TRandom3(0);
  std::cout << "running declustering" << std::endl;
  int intjet_multi = 0;
  float jet_girth = 0;
  Int_t nsplit = 0;
  double dyn_kt = std::numeric_limits<double>::min();
  int dyn_split = 0;
  double z = 0;
  double dyn_deltaR = 0;
  // double dyn_var = std::numeric_limits<double>::min();
  double dyn_z = 0;
  double jet_radius_ca = 1.0;
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);

  fastjet::PseudoJet myjet;
  fastjet::PseudoJet mypart;
  myjet.reset(jet.p4().px(),jet.p4().py(),jet.p4().pz(),jet.p4().e());
  // Reclustering jet constituents with new algorithm
  try{
    std::vector<fastjet::PseudoJet> particles = {};                         
    auto daughters = jet.getJetConstituents();
    //        std::cout << "Number of pfCand " << pfCandidates.size() << std::endl;
        // Geometrical PF Candidate x Jet Constituent Matching - Added by Bharadwaj - Apr 2023
        // poor man's matching, someone fix please
    // std::vector<int> vec_jet_consituent_charge;
    for (auto it = daughters.begin(); it!=daughters.end(); ++it){
      //if we want only charged constituents and the daughter charge is 0, skip it
      if (doChargedConstOnly_ && (**it).charge()==0) continue;
      double PFE_scale = 1.;
      //if it is MC, rescale the 4-momentum of the charged particles (we accept only them above) by pfCCES(+-1%)
      //      if(isMC_) PFE_scale = pfChargedCandidateEnergyScale_;
      //vary tracking efficiency - drop ~4% of particles within the jet
      //      if(isMC_ && doTrackVariation_){
      //  // std::cout << "doing the track variation" << std::endl;
      //  if(r1->Uniform(0,1)<0.05) continue;
      //  }
      // std::cout << "Rescaling charged pfCand energy by " << PFE_scale << std::endl;
      particles.push_back(fastjet::PseudoJet((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale));
      mypart.reset((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale);
      intjet_multi++;
      jet_girth += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();
    }
    std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
    if (particles.empty()){
      // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(123);
    }
    // std::cout << "Clustering " << particles.size() << " number of reco particles" << std::endl;
    fastjet::ClusterSequence csiter(particles, jet_def);
    std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
    output_jets = sorted_by_pt(output_jets);
    fastjet::PseudoJet jj = output_jets[0];
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet highest_splitting;
    if(!jj.has_parents(j1,j2)){
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(124);
    }
    while(jj.has_parents(j1,j2)){
      /*if(j1.perp() < j2.perp()) std::swap(j1,j2);
      double delta_R = j1.delta_R(j2);
      if(doHardestSplitMatching_ && isMC_) jets_.jtJetConstituent.push_back(j2);
      double k_t = j2.perp()*delta_R;
      z = j2.perp()/(j1.perp()+j2.perp());
      // double dyn = z*(1-z)*j2.perp()*pow(delta_R/rParam,mydynktcut);
      // double dyn = 1./output_jets[0].perp()*z*(1-z)*jj.perp()*pow(delta_R/rParam,mydynktcut);
      // std::cout << "Reco split " << nsplit << " with k_T=" << k_t << " z=" << z << " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl;
      if(k_t > dyn_kt){
        // dyn_var = dyn;
        highest_splitting = j2;

        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
      }*/
      if (j1.perp() < j2.perp()) std::swap(j1,j2); // j1 is the harder subjet
      z = j2.perp()/(j1.perp()+j2.perp());
      if (z > 0.1) { // TODO: do not hard code these things..
	double delta_R = j1.delta_R(j2);	
	double k_t = j2.perp()*delta_R;
	//	highest_splitting = j2; // is this needed?
        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
	break;
      }
      jj = j1;
      nsplit = nsplit+1;
    }
    // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest reco splitting eta phi at " << dyn_split << std::endl;
    // jets_.jtdyn_var[jets_.nref] = dyn_var;
    jets_.jtdyn_split[jets_.nref] = dyn_split;
    jets_.jtdyn_deltaR[jets_.nref] = dyn_deltaR;
    jets_.jtdyn_kt[jets_.nref] = dyn_kt;
    jets_.jtdyn_z[jets_.nref] = dyn_z;
    jets_.jt_intjet_multi[jets_.nref] = intjet_multi;
    jets_.jt_girth[jets_.nref] = jet_girth;
  } 
  catch (fastjet::Error const&) { /*return -1;*/ }
  catch (Int_t MyNum){
    if(MyNum == 123)
      std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all reco jet split variables to numeric min." << std::endl;
    if(MyNum == 124)
      std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
      }
}
void HiInclusiveJetAnalyzer::IterativeDeclusteringGen(double groom_type, double groom_combine,const reco::GenJet& jet,fastjet::PseudoJet *sub1,fastjet::PseudoJet *sub2)
{
  int intjet_multi = 0;
  float jet_girth = 0;
  double nsplit = 0;
  // double dyn_theta = 0;
  double dyn_kt = std::numeric_limits<double>::min();
  int dyn_split = 0;
  double dyn_deltaR = 0;
  double z = 0;
  double dyn_z = 0;
  double jet_radius_ca = 1.0;
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
  fastjet::PseudoJet myjet;
  fastjet::PseudoJet mypart;
  myjet.reset(jet.p4().px(),jet.p4().py(),jet.p4().pz(),jet.p4().e());  
    // Reclustering jet constituents with new algorithm
  try{
    std::vector<fastjet::PseudoJet> particles = {};                         
    auto daughters = jet.getJetConstituents();
    for(auto it = daughters.begin(); it!=daughters.end(); ++it){
      //if we want only charged constituents and the daughter charge is 0, skip it
      if(doChargedConstOnly_ && (**it).charge()==0) continue;
      particles.push_back(fastjet::PseudoJet((**it).px(), (**it).py(), (**it).pz(), (**it).energy()));
      mypart.reset((**it).px(), (**it).py(), (**it).pz(), (**it).energy());
      intjet_multi++;
      jet_girth += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();
    }

    // std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
    if(particles.empty()){
      // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.refdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.ref_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(123);
    }
    // std::cout << "Clustering " << particles.size() << " number of truth particles" << std::endl;
    fastjet::ClusterSequence csiter(particles, jet_def);
    std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
    output_jets = sorted_by_pt(output_jets);
    fastjet::PseudoJet jj = output_jets[0];
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet highest_splitting;
    if(!jj.has_parents(j1,j2)){
      jets_.refdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.refdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.ref_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(124);
    }
    // Edit this to SD -> 
    while(jj.has_parents(j1,j2)){
      if(j1.perp() < j2.perp()) std::swap(j1,j2); // j1 is the harder subjet
      z = j2.perp()/(j1.perp()+j2.perp());
      if (z > 0.1) {
	double delta_R = j1.delta_R(j2);	
	double k_t = j2.perp()*delta_R;
	//	highest_splitting = j2; // is this needed?
        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
	break;
      }
     
      //if(doHardestSplitMatching_ && isMC_) jets_.refJetConstituent.push_back(j2); // is this needed?
      //       // std::cout << "Truth split " << nsplit << " with k_T=" << k_t << " z=" << z <<  " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl;
      
      jj = j1; // We follow the hard prong if SD condition not satisfied
      nsplit = nsplit+1;
    }
    // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest truth splitting eta phi at " << dyn_split << std::endl;
    jets_.refdyn_split[jets_.nref] = dyn_split;
    jets_.refdyn_deltaR[jets_.nref] = dyn_deltaR;
    jets_.refdyn_kt[jets_.nref] = dyn_kt;
    jets_.refdyn_z[jets_.nref] = dyn_z;
    jets_.ref_intjet_multi[jets_.nref] = intjet_multi;
    jets_.ref_girth[jets_.nref] = jet_girth;
  } 
  catch (fastjet::Error const&) { /*return -1;*/ }
  catch (Int_t MyNum){
    if(MyNum == 123)
           std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all gen jet split variables to numeric min." << std::endl;
    if(MyNum == 124)
            std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
      } 
}

// Iteratice declustering for aggregated jets, inputs are the aggregated constituent collection and aggregated HF candidate. TODO: check how to handle genJets
void HiInclusiveJetAnalyzer::IterativeDeclustering(std::vector<fastjet::PseudoJet> jetConstituents, reco::PFCandidate pseudoHF)
{
  //  TRandom *r1 = new TRandom3(0);
  int intjet_multi = 0;
  float jet_girth = 0;
  Int_t nsplit = 0;
  double dyn_kt = std::numeric_limits<double>::min();
  int dyn_split = 0;
  double z = 0;
  double dyn_deltaR = 0;
  // double dyn_var = std::numeric_limits<double>::min();
  double dyn_z = 0;
  double jet_radius_ca = 1.0;
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
  fastjet::PseudoJet myjet;
  fastjet::PseudoJet mypart;
  //  myjet.reset(jet.p4().px(),jet.p4().py(),jet.p4().pz(),jet.p4().e());
  // TODO: Follow the HF
  // Reclustering jet constituents with new algorithm
  try{
    
    //std::vector<fastjet::PseudoJet> particles = {};                         // Different thing to be done with GEN jets
    //auto daughters = jetConstituents; //jet.getJetConstituents();
        // std::cout << "Number of pfCand " << pfCandidates.size() << std::endl;
        // Geometrical PF Candidate x Jet Constituent Matching - Added by Bharadwaj - Apr 2023
        // poor man's matching, someone fix please
    // std::vector<int> vec_jet_consituent_charge;
    /*   for(auto it = daughters.begin(); it!=daughters.end(); ++it){  // Should this in general be skipped for genjets (right now configured PFEscale = 1)? -> could add isgenjet boolean?
      //if we want only charged constituents and the daughter charge is 0, skip it
      if(doChargedConstOnly_ && (**it).charge()==0) continue;
      double PFE_scale = 1.;
      //if it is MC, rescale the 4-momentum of the charged particles (we accept only them above) by pfCCES(+-1%)
      //      if(isMC_) PFE_scale = pfChargedCandidateEnergyScale_;
      //vary tracking efficiency - drop ~4% of particles within the jet
      //  if(isMC_ && doTrackVariation_){ 
      //  // std::cout << "doing the track variation" << std::endl;
      //  if(r1->Uniform(0,1)<0.05) continue;
      // }
      // std::cout << "Rescaling charged pfCand energy by " << PFE_scale << std::endl;
        particles.push_back(fastjet::PseudoJet((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale));
      mypart.reset((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale);
      intjet_multi++;
      jet_girth += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();
    } */
    // std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
    
    if (jetConstituents.size() == 0) {
    //if(particles.empty()){
      // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(123);
    }
    // std::cout << "Clustering " << particles.size() << " number of reco particles" << std::endl;
    fastjet::ClusterSequence csiter(jetConstituents, jet_def);
    std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
    output_jets = sorted_by_pt(output_jets);
    fastjet::PseudoJet jj = output_jets[0];
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet highest_splitting;
    
    if(!jj.has_parents(j1,j2)){
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();

      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(124);
    }
    while(jj.has_parents(j1,j2)){
      if (j1.perp() < j2.perp()) std::swap(j1,j2); // j1 is the harder subjet
      z = j2.perp()/(j1.perp()+j2.perp());
      if (z > 0.1) { // TODO: SD vs latekt, define threshold and method in the python configuration
	double delta_R = j1.delta_R(j2);	
	double k_t = j2.perp()*delta_R;
	//	highest_splitting = j2; // is this needed?
        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
	break;
      }
      jj = j1; // Follow the hard prong if SD not satisfied
      nsplit = nsplit+1;
    }
    // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest reco splitting eta phi at " << dyn_split << std::endl;
    // jets_.jtdyn_var[jets_.nref] = dyn_var;
    jets_.jtdyn_split[jets_.nref] = dyn_split;
    jets_.jtdyn_deltaR[jets_.nref] = dyn_deltaR;
    jets_.jtdyn_kt[jets_.nref] = dyn_kt;
    jets_.jtdyn_z[jets_.nref] = dyn_z;
    jets_.jt_intjet_multi[jets_.nref] = intjet_multi;
    jets_.jt_girth[jets_.nref] = jet_girth;
  } 
  catch (fastjet::Error const&) { /*return -1;*/ }
  catch (Int_t MyNum){
    if(MyNum == 123)
         std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all reco jet split variables to numeric min." << std::endl;
    if(MyNum == 124)
         std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
      }
}

//maybe there is a more elegant way than the one below for matching...
void HiInclusiveJetAnalyzer::RecoTruthSplitMatching(std::vector<fastjet::PseudoJet> &constituents_level1, fastjet::PseudoJet &hardest_level2, bool *bool_array, int *hardest_level1_split){
    //for now only include geometric matching, maybe consider pt/z, Lund plane location, etc...
  float min_dR = std::numeric_limits<float>::max();
  size_t closest_level1 = 0;
  // std::cout << "Starting loop over " << constituents_level1.size() << " particles" << std::endl;
  for(size_t i{0};i<constituents_level1.size();++i){
    float dR = constituents_level1.at(i).delta_R(hardest_level2);
    if(min_dR > dR){
      closest_level1 = i;
      min_dR = dR;
    }
  }
  // std::cout << "Compare particle " << static_cast<int>(closest_level1) << " with hardest split " << hardest_level1_split[jets_.nref] << std::endl;
  if(static_cast<int>(closest_level1) == hardest_level1_split[jets_.nref] ){
    bool_array[jets_.nref] = true;
  }
  else{
    // std::cout << "Sorry, closest pair is " << min_dR << " away with index " << static_cast<int>(closest_level1) << " as opposed to " << hardest_level1_split[jets_.nref] << std::endl;
    bool_array[jets_.nref] = false;
  }
}
void HiInclusiveJetAnalyzer::TruthRecoRecoTruthMatching(){
  // std::cout << jets_.jtdyn_split[jets_.nref] << " " << jets_.refdyn_split[jets_.nref] << " numbers of highest splits" << std::endl;
  if( jets_.jtdyn_split[jets_.nref] == std::numeric_limits<int>::min() || jets_.refdyn_split[jets_.nref] == std::numeric_limits<int>::min() || jets_.jtJetConstituent.size() == 0 || jets_.refJetConstituent.size() == 0 ){
    jets_.refdyn_isClosestToReco[jets_.nref] = false;

    jets_.jtdyn_isClosestToTruth[jets_.nref] = false;
    jets_.jtdyn_refdyn_dR[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }
  //mind how the split number is defined in the reclustering
  fastjet::PseudoJet hardest_R_split = jets_.jtJetConstituent.at(jets_.jtdyn_split[jets_.nref]);
  fastjet::PseudoJet hardest_T_split = jets_.refJetConstituent.at(jets_.refdyn_split[jets_.nref]);
  // std::cout << hardest_R_split.eta() << " " << hardest_R_split.phi() << " hardest reco  splitting in matching" << std::endl;
  // std::cout << hardest_T_split.eta() << " " << hardest_T_split.phi() << " hardest truth splitting in matching" << std::endl;
  // std::cout << "Angle between hardest splits is dR = " << hardest_R_split.delta_R(hardest_T_split) << std::endl;
  jets_.jtdyn_refdyn_dR[jets_.nref] = hardest_R_split.delta_R(hardest_T_split);
  // std::cout << "truth loop" << std::endl;
  RecoTruthSplitMatching(jets_.refJetConstituent, hardest_R_split, jets_.refdyn_isClosestToReco, jets_.refdyn_split);
  // std::cout << "reco loop" << std::endl;
  RecoTruthSplitMatching(jets_.jtJetConstituent,  hardest_T_split, jets_.jtdyn_isClosestToTruth, jets_.jtdyn_split);
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
