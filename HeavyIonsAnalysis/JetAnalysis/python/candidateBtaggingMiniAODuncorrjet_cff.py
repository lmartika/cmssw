import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
ak2PFJets = ak4PFJets.clone(jetPtMin = 1.0, rParam = 0.2, src = 'packedPFCandidates')
from RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi import pfImpactParameterTagInfos
pfImpactParameterTagInfos.jets = "ak2PFJets"
pfImpactParameterTagInfos.candidates = "packedPFCandidates"
pfImpactParameterTagInfos.primaryVertex = "offlineSlimmedPrimaryVerticesRecovery"
#from RecoBTag.SecondaryVertex.pfSecondaryVertexTagInfos_cfi import pfSecondaryVertexTagInfos
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import inclusiveCandidateVertexFinder
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import candidateVertexMerger
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import candidateVertexArbitrator
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import inclusiveCandidateSecondaryVertices
from RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi import pfInclusiveSecondaryVertexFinderTagInfos
inclusiveCandidateVertexFinder.primaryVertices  = "offlineSlimmedPrimaryVerticesRecovery"
inclusiveCandidateVertexFinder.tracks= "packedPFCandidates"
candidateVertexArbitrator.tracks = "packedPFCandidates"
candidateVertexArbitrator.primaryVertices = "offlineSlimmedPrimaryVerticesRecovery"
pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = "slimmedSecondaryVertices"  

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
    addTagInfos = True
)


from RecoBTag.ONNXRuntime.pfParticleNetAK4_cff import pfParticleNetAK4TagInfos, pfParticleNetAK4JetTags
from RecoBTag.ONNXRuntime.pfParticleNetAK4DiscriminatorsJetTags_cfi import pfParticleNetAK4DiscriminatorsJetTags
# temporarily run PNET on unsubtracted jets.  FIXME!  -Matt
pfParticleNetAK4TagInfos.jets = "ak2PFpatJets"
#pfParticleNetAK4TagInfos.unsubjet_map = "unsubJets"
pfParticleNetAK4TagInfos.use_puppiP4 = False
pfParticleNetAK4TagInfos.pf_candidates = "packedPFCandidates"
pfParticleNetAK4TagInfos.puppi_value_map = ''
pfParticleNetAK4TagInfos.vertex_associator = ""
pfParticleNetAK4TagInfos.vertices = "offlineSlimmedPrimaryVerticesRecovery"

pfParticleNetAK4JetTags.src = "pfParticleNetAK4TagInfos"



unsubCandidateBtagging = cms.Sequence(
    ak2PFJets +
    pfImpactParameterTagInfos +
    #pfSecondaryVertexTagInfos +
    inclusiveCandidateVertexFinder +
    candidateVertexMerger +
    candidateVertexArbitrator +
    inclusiveCandidateSecondaryVertices +
    pfInclusiveSecondaryVertexFinderTagInfos +
    ak2PFpatJets +
    pfParticleNetAK4TagInfos +
    pfParticleNetAK4JetTags +  
    pfParticleNetAK4DiscriminatorsJetTags 
)
