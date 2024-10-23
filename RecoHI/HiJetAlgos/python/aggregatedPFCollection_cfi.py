import FWCore.ParameterSet.Config as cms

aggregatedPFCands = cms.EDProducer('aggregatedPFCands',
#                                    jetSrc = cms.InputTag('slimmedJets'),
                                    constitSrc = cms.InputTag('packedPFCandidates'),
                                    candToGenParticleMap = cms.InputTag("TrackToGenParticleMapProducer", "trackToGenParticleMap"),
                                    isMC = cms.bool(True),
                                    chargedOnly = cms.bool(False),
                                    ptCut = cms.double(1),
                                    trkInefRate = cms.double(0),
                                    doGenJets = cms.bool(False),
                                    aggregateHF = cms.bool(True),
                                    aggregateWithTruthInfo = cms.bool(True),
                                    aggregateWithCuts = cms.bool(False),
                                    aggregateWithXGB = cms.bool(False),
                                    aggregateWithTMVA = cms.bool(False),
                                    ipTagInfoLabel = cms.string("pfImpactParameter"),
                                    svTagInfoLabel = cms.string("pfInclusiveSecondaryVertexFinder"),
#                                    tmva_path = cms.FileInPath("RecoHI/HiJetAlgos/data/TMVAClassification_BDTG.weights.xml"),
                                    tmva_variables = cms.vstring("trkIp3dSig", "trkIp2dSig", "trkDistToAxis", "svtxdls", "svtxdls2d", "svtxm", "svtxmcorr", "svtxnormchi2", "svtxNtrk", "svtxTrkPtOverSv", "jtpt"),
                                    tmva_spectators = cms.vstring(),

                                    )             
