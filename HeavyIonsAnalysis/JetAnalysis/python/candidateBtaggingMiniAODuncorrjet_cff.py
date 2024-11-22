import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
ak2PFJets = ak4PFJets.clone(jetPtMin = 1.0, rParam = 0.2, src = 'packedPFCandidates')
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJetsNoNu
ak2GenJetsNoNuNonAgg = ak4GenJetsNoNu.clone(src = 'packedGenParticlesSignal', rParam = 0.2)
from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetGenJetMatch
ak2PatJetGenJetMatch = patJetGenJetMatch.clone(src = "ak2PFJets", matched = "ak2GenJetsNoNuNonAgg", maxDeltaR = 0.2)
from RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi import pfImpactParameterTagInfos
pfImpactParameterTagInfos.jets = "ak2PFJets"
pfImpactParameterTagInfos.candidates = "packedPFCandidates"
pfImpactParameterTagInfos.primaryVertex = "offlineSlimmedPrimaryVertices"
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

from RecoBTau.JetTagComputer.jetTagRecord_cfi import *
from RecoBTag.ImpactParameter.candidateJetProbabilityComputer_cfi import *
from RecoBTag.ImpactParameter.pfJetProbabilityBJetTags_cfi import *

from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi import pfInclusiveSecondaryVertexFinderTagInfos
inclusiveCandidateVertexFinder.primaryVertices  = "offlineSlimmedPrimaryVertices"
inclusiveCandidateVertexFinder.tracks= "packedPFCandidates"
inclusiveCandidateVertexFinder.minHits = 0
inclusiveCandidateVertexFinder.minPt = 0.8
candidateVertexArbitrator.tracks = "packedPFCandidates"
candidateVertexArbitrator.primaryVertices = "offlineSlimmedPrimaryVertices"
pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = "inclusiveCandidateSecondaryVertices"
from PhysicsTools.PatAlgos.producersHeavyIons.heavyIonJets_cff import hiSignalGenParticles
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
selectedHadronsAndPartons.particles = "hiSignalGenParticles"
from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
ak2JetFlavourInfos = ak4JetFlavourInfos.clone(jets = "ak2PFJets", rParam=0.2)


from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
ak2PFpatJets = patJets.clone(
    jetSource = "ak2PFJets",
    tagInfoSources = cms.VInputTag('pfImpactParameterTagInfos','pfInclusiveSecondaryVertexFinderTagInfos'),
    discriminatorSources = cms.VInputTag("pfJetProbabilityBJetTags"),
    genJetMatch = cms.InputTag("ak2PatJetGenJetMatch"),
    JetFlavourInfoSource = "ak2JetFlavourInfos", 
    addAssociatedTracks = False,
    addJetCorrFactors = False,
    addGenPartonMatch = False,
    addGenJetMatch = True,
    addDiscriminators = True,
    getJetMCFlavour = True,
    useLegacyJetMCFlavour = False,
    addBTagInfo = True,
    addTagInfos = True,
    addJetFlavourInfo = True
)


unsubCandidateBtagging = cms.Sequence(
    ak2PFJets +
    hiSignalGenParticles +
    ak2GenJetsNoNuNonAgg +
    ak2PatJetGenJetMatch +
    selectedHadronsAndPartons +
    ak2JetFlavourInfos +
    pfImpactParameterTagInfos +
    pfJetProbabilityBJetTags + 
    inclusiveCandidateVertexFinder +
    candidateVertexMerger +
    candidateVertexArbitrator +
    inclusiveCandidateSecondaryVertices +
    pfInclusiveSecondaryVertexFinderTagInfos +
    ak2PFpatJets 
)
