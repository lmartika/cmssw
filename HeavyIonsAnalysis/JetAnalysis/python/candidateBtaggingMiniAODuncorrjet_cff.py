import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
ak2PFJets = ak4PFJets.clone(jetPtMin = 1.0, rParam = 0.2, src = 'packedPFCandidates')
from RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi import pfImpactParameterTagInfos
pfImpactParameterTagInfos.jets = "ak2PFJets"
pfImpactParameterTagInfos.candidates = "packedPFCandidates"
pfImpactParameterTagInfos.primaryVertex = "offlineSlimmedPrimaryVerticesRecovery"
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from RecoBTau.JetTagComputer.jetTagRecord_cfi import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi import pfInclusiveSecondaryVertexFinderTagInfos
inclusiveCandidateVertexFinder.primaryVertices  = "offlineSlimmedPrimaryVerticesRecovery"
inclusiveCandidateVertexFinder.tracks= "packedPFCandidates"
inclusiveCandidateVertexFinder.minHits = 0
inclusiveCandidateVertexFinder.minPt = 0.8
candidateVertexArbitrator.tracks = "packedPFCandidates"
candidateVertexArbitrator.primaryVertices = "offlineSlimmedPrimaryVerticesRecovery"
pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = "inclusiveCandidateSecondaryVertices"

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
ak2PFpatJets = patJets.clone(
    jetSource = cms.InputTag("ak2PFJets"),
    tagInfoSources = cms.VInputTag('pfImpactParameterTagInfos','pfInclusiveSecondaryVertexFinderTagInfos'),
    addAssociatedTracks = False,
    addJetCorrFactors = False,
    addGenPartonMatch = False,
    addGenJetMatch = False,
    addDiscriminators = False,
    getJetMCFlavour = False,
    useLegacyJetMCFlavour = False,
    addBTagInfo = True,
    addTagInfos = True
)


unsubCandidateBtagging = cms.Sequence(
    ak2PFJets +
    pfImpactParameterTagInfos +
    inclusiveCandidateVertexFinder +
    candidateVertexMerger +
    candidateVertexArbitrator +
    inclusiveCandidateSecondaryVertices +
    pfInclusiveSecondaryVertexFinderTagInfos +
    ak2PFpatJets 
)
