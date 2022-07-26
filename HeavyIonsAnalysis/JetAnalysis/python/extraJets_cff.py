import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.producersHeavyIons.heavyIonJets_cff import PackedPFTowers, hiPuRho, hiSignalGenParticles, allPartons
hiSignalGenParticles.src = "prunedGenParticles"
hiPuRho.src = 'PackedPFTowers' 
extraJetsData = cms.Sequence(PackedPFTowers + hiPuRho)
extraJetsMC = cms.Sequence(PackedPFTowers + hiPuRho + hiSignalGenParticles + allPartons)

