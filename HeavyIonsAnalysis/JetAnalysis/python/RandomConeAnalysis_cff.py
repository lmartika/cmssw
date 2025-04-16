import FWCore.ParameterSet.Config as cms
from RecoHI.HiJetAlgos.PFTowers_cfi import PFTowers
# hiPuRho and hiFJRhoFlowModulation producers are both called by the CSJetProducer, so call rhoAnalysis after that and load the inputs from there
randomConeAnalysisR4 = cms.EDAnalyzer('RandomConeAnalysis',
                                etaMap             = cms.InputTag('hiPuRho','mapEtaEdges'),
                                rho                = cms.InputTag('hiPuRho','mapToRho'),
                                useModulatedRho = cms.bool(True),
                                minFlowChi2Prob = cms.double(0.05),
                                maxFlowChi2Prob = cms.double(0.95),
                                rngConeRadius   = cms.double(0.4),
                                pfCandSource = cms.InputTag('packedPFCandidates'),
                                rhoFlowFitParams   = cms.InputTag('hiFJRhoFlowModulationProducer', 'rhoFlowFitParams'),
)

randomConeAnalysisR2 = randomConeAnalysisR4.clone(
                                rngConeRadius   = cms.double(0.2),
)

randomConeAnalysisR3 = randomConeAnalysisR4.clone(
                                rngConeRadius   = cms.double(0.3),
)