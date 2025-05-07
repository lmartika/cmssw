import FWCore.ParameterSet.Config as cms
# from RecoHI.HiJetAlgos.PFTowers_cfi import PFTowers
# hiPuRho and hiFJRhoFlowModulation producers are both called by the CSJetProducer, so call rhoAnalysis after that and load the inputs from there
rhoAnalysis = cms.EDAnalyzer('RhoAnalysis',
                                etaMap             = cms.InputTag('hiPuRho','mapEtaEdges'),
                                rho                = cms.InputTag('hiPuRho','mapToRho'),
                                # useModulatedRho = cms.bool(True),
                                # minFlowChi2Prob = cms.double(0.05),
                                # maxFlowChi2Prob = cms.double(0.95),
                                # rngConeRadius   = cms.double(0.4),
                                # pfCandSource = cms.InputTag('packedPFCandidates'),
                                # rhoFlowFitParams   = cms.InputTag('hiFJRhoFlowModulationProducer', 'rhoFlowFitParams'),
)

# rhoAnalysisR2 = rhoAnalysisR4.clone(
#                                 rngConeRadius   = cms.double(0.2),
# )