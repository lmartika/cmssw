### HiForest Configuration
# Input: miniAOD
# Type: mc

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
process = cms.Process('HiForest', Run3_pp_on_PbPb_2023)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, mc")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #'/store/group/phys_heavyions/jviinika/PythiaHydjetRun3_5p36TeV_dijet_ptHat15_100kEvents_miniAOD_2023_08_30/PythiaHydjetDijetRun3/PythiaHydjetRun3_dijet_ptHat15_5p36TeV_miniAOD/230830_165931/0000/pythiaHydjet_miniAOD_11.root'
        #'/store/mc/HINPbPbSpring23MiniAOD/TT-2Jets_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/MINIAODSIM/132X_mcRun3_2023_realistic_HI_v9-v3/120000/0004a160-9198-41fb-a74b-3b4c009dd6ac.root'
        '/store/group/phys_heavyions/lamartik/bjet/043213d2-944a-4e18-b1b5-ef71e93ef850.root'  # 2018 bjet sample
    ),
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
    )

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_mcRun3_2023_realistic_HI_v10', '')
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
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
#process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
#process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

#from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_mc
#process.hltobject.triggerNames = trigger_list_mc

################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_mc_cff')
# gen jets
process.genJetSequence = cms.Sequence()

################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")

###############################################################################

process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")

process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAODuncorrjet_cff")

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
#    process.hltanalysis +
#    process.hltobject +
#    process.l1object +
#    process.trackSequencePbPb +
#    process.particleFlowAnalyser +
    process.hiEvtAnalyzer +
    process.genJetSequence +
    process.extraJetsMC+
    process.unsubCandidateBtagging#+

    #    process.HiGenParticleAna 
    )


doTracks = True
doSvtx = True

doGenAnalysis = True
runAggregation = True
#customisation

matchJets = False             # Enables q/g and heavy flavor jet identification in MC
addCandidateTagging = False
#doHIJetID = True             # Fill jet ID and composition information branches
#doWTARecluster = False        # Add jet phi and eta for WTA axis

if doGenAnalysis:
    process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
    process.genJetSequence += process.mergedGenParticles
    ## Produces a reco::GenParticleCollection named mergedGenParticles

    process.load("RecoHI.HiJetAlgos.HFdecayProductTagger_cfi")
    process.HFdecayProductTagger.genParticles = cms.InputTag("mergedGenParticles")
    process.HFdecayProductTagger.tagBorC = cms.bool(True) # tag B
    process.genJetSequence += process.HFdecayProductTagger

    taggedGenParticlesName_ = "HFdecayProductTagger"
    ## Produces a std::vector<pat::PackedGenParticle> named HFdecayProductTagger

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

    process.TrackToGenParticleMapProducer.jetSrc = cms.InputTag("ak2PFpatJets")
    process.TrackToGenParticleMapProducer.genParticleSrc = cms.InputTag(taggedGenParticlesName_)
    process.forest += process.TrackToGenParticleMapProducer

    if runAggregation:
            process.load("RecoHI.HiJetAlgos.aggregatedPFCollection_cfi")
            process.aggregatedPFCands.aggregateHF = True
            process.aggregatedPFCands.jetSrc = "ak2PFpatJets"
            process.aggregatedPFCands.constitSrc = "packedPFCandidates"
            process.aggregatedPFCands.doGenJets = False
            process.aggregatedPFCands.aggregateWithTruthInfo = True
            process.aggregatedPFCands.candToGenParticleMap = ["TrackToGenParticleMapProducer", "trackToGenParticleMap"]

            process.aggregatedGenLevel  = process.aggregatedPFCands.clone(
                chargedOnly = cms.bool(True),
                aggregateHF = cms.bool(True),
                jetSrc = cms.InputTag("ak2PFpatJets"),
                constitSrc = cms.InputTag("packedGenParticles"),
                doGenJets = cms.bool(True),
                candToGenParticleMap = cms.InputTag("TrackToGenParticleMapProducer", "genConstitToGenParticleMap"),
            )

            process.forest +=process.aggregatedPFCands + process.aggregatedGenLevel 

            process.aggregatedJets = cms.Sequence()
            from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupHeavyIonJets
            setupHeavyIonJets('akCs2PF', process.aggregatedJets, process, isMC = 1, radius = 0.20, JECTag = 'AK2PF')
#            process.akCs2PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute'] 
            process.akCs2PFpatJetCorrFactors.levels = [] 
            process.akCs2PFpatJets.addJetCorrFactors = True  #check!  Matt
            process.akCs2PFpatJets.addBTagInfo = False
            process.akCs2PFpatJets.addDiscriminators = False
            process.akCs2PFJets.src = 'aggregatedPFCands'
            process.ak2GenJetsNoNu.src = 'aggregatedGenLevel'


            process.unsubJets = cms.EDProducer("JetMatcherDR",
                                               matched = cms.InputTag("ak2PFpatJets"),
                                               source = cms.InputTag("akCs2PFpatJets")
                                           )



#            updateJetCollection(
#                process,
#                jetSource = cms.InputTag('akCs2PFpatJets'),
#                jetCorrections = ('AK2PF', cms.vstring(), 'None'),
#                btagDiscriminators = ['pfParticleNetAK4JetTags:probb', 'pfParticleNetAK4JetTags:probbb'], ## to add discriminators,
#                btagPrefix = 'TEST',
#            )
#            process.updatedPatJetCorrFactors = process.akCs2PFpatJetCorrFactors.clone(src = "akCs2PFpatJets", levels = ['L2Relative', 'L3Absolute'] )
#            process.updatedPatJets.jetCorrFactorsSource = cms.VInputTag("updatedPatJetCorrFactors")
#            process.updatedPatJets.addBTagInfo = True
#            process.updatedPatJets.addTagInfos = False

#            process.akCs2PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "updatedPatJets", jetName = 'akCs2PF', genjetTag = "ak2GenJetsNoNu")
            process.akCs2PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "akCs2PFpatJets", jetName = 'akCs2PF', genjetTag = "ak2GenJetsNoNu")
#            process.akCs2PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "ak2PFpatJets", jetName = 'ak2PF', genjetTag = "ak2GenJetsNoNu")      
#            #process.akCs2PFJetAnalyzer.jetPtMin = cms.double(70.0)
#            process.akCs2PFJetAnalyzer.doCandidateBtagging = True  
            process.akCs2PFJetAnalyzer.unsubjet_map = cms.InputTag("unsubJets")
            process.akCs2PFJetAnalyzer.matchJets = cms.untracked.bool(True)
            process.akCs2PFJetAnalyzer.matchTag = cms.untracked.InputTag("ak2PFpatJets")
#            process.forest += process.aggregatedJets * process.unsubJets *process.pfParticleNetAK4TagInfos * process.pfParticleNetAK4JetTags * process.pfParticleNetAK4DiscriminatorsJetTags *process.updatedPatJetCorrFactors * process.updatedPatJets * process.akCs2PFJetAnalyzer
            process.forest += process.aggregatedJets * process.unsubJets * process.akCs2PFJetAnalyzer
            process.akCs2PFpatJets.addTagInfos = True
            process.akCs2PFpatJets.addBTagInfo = True
            
            #process.ak2PFJetAnalyzer = process.akCs2PFJetAnalyzer.clone(jetTag = "ak2PFpatJets", jetName = 'ak2PF', isMC = False, saveRawPt = False)
            process.ak2PFpatJets.tagInfoSources = ["pfImpactParameterTagInfos", "pfInclusiveSecondaryVertexFinderTagInfos"]
            if doTracks:
                process.akCs2PFJetAnalyzer.doTracks = cms.untracked.bool(True)
                process.akCs2PFJetAnalyzer.ipTagInfoLabel = cms.untracked.string("pfImpactParameter")
            if doSvtx:
                process.akCs2PFJetAnalyzer.doSvtx = cms.untracked.bool(True)
                process.akCs2PFJetAnalyzer.svTagInfoLabel = cms.untracked.string("pfInclusiveSecondaryVertexFinder")
            
#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)
process.pAna = cms.EndPath(process.skimanalysis)
