## HiForest Configuration
# Input: miniAOD
# Type: mc

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_pp_on_AA_cff import Run2_2018_pp_on_AA
from Configuration.ProcessModifiers.run2_miniAOD_pp_on_AA_103X_cff import run2_miniAOD_pp_on_AA_103X
process = cms.Process('HiForest', Run2_2018_pp_on_AA,run2_miniAOD_pp_on_AA_103X)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 112X, mc")

# import subprocess, os
# version = subprocess.check_output(
#     ['git', '-C', os.path.expandvars('$CMSSW_BASE/src'), 'describe', '--tags'])
# if version == '':
#     version = 'no git info'
# process.HiForestInfo.HiForestVersion = cms.string(version)

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
#        '/store/himc/HINPbPbSpring21MiniAOD/DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/MINIAODSIM/FixL1CaloGT_New_Release_112X_upgrade2018_realistic_HI_v9-v1/2520000/00147e49-765a-424c-a7e0-29860f11847d.root'
        'file:/eos/user/l/lamartik/testsamples/pbpbdijet2018/00147e49-765a-424c-a7e0-29860f11847d.root'
    ),
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(10)
#    input = cms.untracked.int32(2000)
    )

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic_hi', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("JPcalib_MC103X_2018PbPb_v4"),
             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
         )
])


###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD.root"))

# # edm output for debugging purposes
# process.output = cms.OutputModule(
#     "PoolOutputModule",
#     fileName = cms.untracked.string('HiForestEDM.root'),
#     outputCommands = cms.untracked.vstring(
#         'keep *',
#         )
#     )

# process.output_path = cms.EndPath(process.output)

###############################################################################

#############################
# Gen Analyzer
#############################
process.load('HeavyIonsAnalysis.EventAnalysis.HiGenAnalyzer_cfi')
# making cuts looser so that we can actually check dNdEta
process.HiGenParticleAna.ptMin = cms.untracked.double(0.4) # default is 5
process.HiGenParticleAna.etaMax = cms.untracked.double(5.) # default is 2.5

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi')
#process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_mc
process.hltobject.triggerNames = trigger_list_mc

################################
# electrons, photons, muons
SS2018PbPbMC = "HeavyIonsAnalysis/EGMAnalysis/data/SS2018PbPbMC.dat"
process.load('HeavyIonsAnalysis.EGMAnalysis.correctedElectronProducer_cfi')
process.correctedElectrons.correctionFile = SS2018PbPbMC

process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doGenParticles = cms.bool(True)
process.ggHiNtuplizer.doMuons = cms.bool(False)
process.ggHiNtuplizer.electronSrc = "correctedElectrons"
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_mc_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
#muons
process.load("HeavyIonsAnalysis.MuonAnalysis.unpackedMuons_cfi")
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.muonAnalyzer.doGen = cms.bool(True)


process.genJetSequence = cms.Sequence()

###############################################################################



###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.hltanalysis +
#    process.hltobject +
#    process.l1object +
#    process.trackSequencePbPb +
#    process.particleFlowAnalyser +
    process.hiEvtAnalyzer +
    process.HiGenParticleAna + 
    process.genJetSequence
#    process.ggHiNtuplizer +
    )

#customisation

addR2Jets = True
addR4Jets = False

addCandidateTagging = True

doTracks = False
doSvtx = True
doGenAnalysis = False

if addR2Jets or addR4Jets :
    process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
    process.load("RecoHI.HiJetAlgos.EventConstSub_cfi")
    process.forest += process.extraJetsMC

    from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupHeavyIonJets

    if addR2Jets :
        process.jetsR2 = cms.Sequence()
        setupHeavyIonJets('akCs2PF', process.jetsR2, process, isMC = 1, radius = 0.20, JECTag = 'AK2PF')
        process.akCs2PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.akCs2PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "akCs2PFpatJets", jetName = 'akCs2PF', genjetTag = "ak2GenJetsNoNu")      
        process.forest += process.jetsR2 * process.akCs2PFJetAnalyzer

        
    if addR4Jets :
        # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections
        process.jetsR4 = cms.Sequence()
        setupHeavyIonJets('akCs0PF', process.jetsR4, process, isMC = 1, radius = 0.40, JECTag = 'AK4PF')
        process.akCs0PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.akCs4PFJetAnalyzer.jetTag = 'akCs0PFpatJets'
        process.akCs4PFJetAnalyzer.jetName = 'akCs0PF'
        process.forest += process.jetsR4 * process.akCs4PFJetAnalyzer

# To apply PNET we have to remove the JEC first.  We could have just done this before the JEC is applied, but I copied this from the pp workflow
# Maybe something to clean up for the future
# Also note that it's only setup for R=0.4 for the moment. -> own edit to AK2


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

ipTagInfoLabel_ = "pfImpactParameter"
svTagInfoLabel_ = "pfInclusiveSecondaryVertexFinder"

if addCandidateTagging:
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")

    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK2PF', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators = ['pfCombinedSecondaryVertexV2BJetTags', 'pfDeepCSVDiscriminatorsJetTags:BvsAll', 'pfDeepCSVDiscriminatorsJetTags:CvsB', 'pfDeepCSVDiscriminatorsJetTags:CvsL'], ## to add discriminators,
        btagPrefix = 'TEST',
    )
    if addR2Jets : process.updatedPatJets.jetSource = 'akCs2PFpatJets'
    process.updatedPatJets.addJetCorrFactors = True
    process.undoPatJetCorrFactors = process.akCs2PFpatJetCorrFactors.clone(
        levels = cms.vstring(),
        src = cms.InputTag("akCs2PFpatJets"),
    )
    process.updatedPatJets.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("undoPatJetCorrFactors"))

    process.redoPatJetCorrFactors = process.akCs2PFpatJetCorrFactors.clone(
        src = cms.InputTag("updatedPatJets"),
    )

#SV needed for aggregation
    process.load("RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi")
    process.pfImpactParameterTagInfos.candidates  = "packedPFCandidates"
    process.pfImpactParameterTagInfos.primaryVertex = "offlineSlimmedPrimaryVertices"
#    process.pfImpactParameterTagInfos.primaryVertex = "hiSelectedVertex"
    process.pfImpactParameterTagInfos.jets = "updatedPatJets"

    process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi")
    process.pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = "slimmedSecondaryVertices"

    process.updatedCorrectedPatJets = process.updatedPatJets.clone(
        jetSource = "updatedPatJets",
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("redoPatJetCorrFactors")),
        tagInfoSources = cms.VInputTag(
            cms.InputTag(ipTagInfoLabel_ + "TagInfos"),
            cms.InputTag(svTagInfoLabel_ + "TagInfos")),
        discriminatorSources = cms.VInputTag(
            cms.InputTag('pfParticleNetAK4JetTags:probb'),
            cms.InputTag('pfParticleNetAK4JetTags:probbb'),
            cms.InputTag('pfParticleNetAK4JetTags:probc'),
            cms.InputTag('pfParticleNetAK4JetTags:probcc'),
            cms.InputTag('pfParticleNetAK4JetTags:probpu'),
            cms.InputTag('pfParticleNetAK4JetTags:probg'),
            cms.InputTag('pfParticleNetAK4JetTags:probuds'),
            cms.InputTag('pfParticleNetAK4JetTags:probundef'),
            cms.InputTag('pfParticleNetAK4DiscriminatorsJetTags:BvsAll'),
            cms.InputTag('pfParticleNetAK4DiscriminatorsJetTags:CvsL'),
            cms.InputTag('pfParticleNetAK4DiscriminatorsJetTags:QvsG'),
            cms.InputTag('pfParticleNetAK4DiscriminatorsJetTags:CvsB'),
       )
    )
    process.updatedCorrectedPatJets.addTagInfos = True

#    process.updatedCorrectedPatJets.tagInfoSources = cms.VInputTag(
#        cms.InputTag(ipTagInfoLabel_ + "TagInfos"),
#        cms.InputTag(svTagInfoLabel_ + "TagInfos"),
#    )

    process.forest.insert(-1,
                          process.undoPatJetCorrFactors*
                          process.updatedPatJets*
                          process.pfImpactParameterTagInfos *
                          process.pfInclusiveSecondaryVertexFinderTagInfos *
                          process.candidateBtagging*
                          process.redoPatJetCorrFactors*
                          process.updatedCorrectedPatJets
                      )
    process.akCs2PFJetAnalyzer.jetTag = "updatedCorrectedPatJets"
    process.akCs2PFJetAnalyzer.doCandidateBtagging = True

if doTracks:
    process.akCs2PFJetAnalyzer.doTracks = cms.untracked.bool(True)
    process.akCs2PFJetAnalyzer.ipTagInfoLabel = cms.untracked.string(ipTagInfoLabel_)
if doSvtx:
    process.akCs2PFJetAnalyzer.doSvtx = cms.untracked.bool(True)
    process.akCs2PFJetAnalyzer.svTagInfoLabel = cms.untracked.string(svTagInfoLabel_)


if doGenAnalysis:
## Track-Gen-matches, try to use AK2    
    process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
    process.genJetSequence += process.mergedGenParticles
    ## Produces a reco::GenParticleCollection named mergedGenParticles

    process.load("RecoHI.HiJetAlgos.HFdecayProductTagger_cfi")
    process.HFdecayProductTagger.genParticles = cms.InputTag("mergedGenParticles")
    process.HFdecayProductTagger.tagBorC = cms.bool(True) # tag B
    process.genJetSequence += process.HFdecayProductTagger

    taggedGenParticlesName_ = "HFdecayProductTagger"
    ## Produces a std::vector<pat::PackedGenParticle> named HFdecayProductTagger
    process.akCs2PFJetAnalyzer.genParticles = cms.untracked.InputTag(taggedGenParticlesName_)

    process.bDecayAna = process.HiGenParticleAna.clone(
        genParticleSrc = cms.InputTag(taggedGenParticlesName_),
        useRefVector = cms.untracked.bool(False),
        partonMEOnly = cms.untracked.bool(False),
        chargedOnly = True, 
        doHI = False,
        etaMax = cms.untracked.double(10),
        ptMin = cms.untracked.double(0),
        stableOnly = False
    )
    process.genJetSequence += process.bDecayAna

    process.load("RecoHI.HiJetAlgos.TrackToGenParticleMapProducer_cfi")

    process.TrackToGenParticleMapProducer.jetSrc = cms.InputTag("updatedCorrectedPatJets")  
    process.TrackToGenParticleMapProducer.genParticleSrc = cms.InputTag(taggedGenParticlesName_)
    process.forest.insert(-1,process.TrackToGenParticleMapProducer)

    
#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.pAna = cms.EndPath(process.skimanalysis)
