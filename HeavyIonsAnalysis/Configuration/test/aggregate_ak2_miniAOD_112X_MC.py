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

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #'file:/eos/user/l/lamartik/testsamples/pbpbbjet2018/043213d2-944a-4e18-b1b5-ef71e93ef850.root'
        "/store/group/phys_heavyions/lamartik/bjet/043213d2-944a-4e18-b1b5-ef71e93ef850.root"
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
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_mc_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")

process.genJetSequence = cms.Sequence()

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
    process.hiEvtAnalyzer+
#    process.HiGenParticleAna + 
    process.genJetSequence+
    process.extraJetsMC+
    process.unsubCandidateBtagging#+
    #process.updatedPatJets
)



doTracks = True
doSvtx = True

doGenAnalysis = True
runAggregation = True



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
    # Following probably not needed?
    #process.ak2PFJetAnalyzer.genParticles = cms.untracked.InputTag(taggedGenParticlesName_)

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
            #process.akCs2PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute'] 
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


            process.load("RecoBTag.ONNXRuntime.pfParticleNetAK4_cff")
            process.load("RecoBTag.ONNXRuntime.pfParticleNetAK4DiscriminatorsJetTags_cfi")
            process.pfParticleNetAK4TagInfos.jets = "akCs2PFpatJets" 
            process.pfParticleNetAK4TagInfos.unsubjet_map = "unsubJets"
            process.pfParticleNetAK4TagInfos.use_puppiP4 = False
            process.pfParticleNetAK4TagInfos.pf_candidates = "packedPFCandidates"
            process.pfParticleNetAK4TagInfos.puppi_value_map = ''
            process.pfParticleNetAK4TagInfos.vertex_associator = ""
            process.pfParticleNetAK4TagInfos.vertices = "offlineSlimmedPrimaryVerticesRecovery"

            updateJetCollection(
                process,
                jetSource = cms.InputTag('akCs2PFpatJets'),
                jetCorrections = ('AK2PF', cms.vstring(), 'None'),
                btagDiscriminators = ['pfParticleNetAK4JetTags:probb', 'pfParticleNetAK4JetTags:probbb'], ## to add discriminators,
                btagPrefix = 'TEST',
            )
            process.updatedPatJetCorrFactors = process.akCs2PFpatJetCorrFactors.clone(src = "akCs2PFpatJets", levels = ['L2Relative', 'L3Absolute'] )
            process.updatedPatJets.jetCorrFactorsSource = cms.VInputTag("updatedPatJetCorrFactors")
            process.updatedPatJets.addBTagInfo = True
            process.updatedPatJets.addTagInfos = False
            process.updatedPatJets.discriminatorSources =  ['pfParticleNetAK4JetTags:probb', 'pfParticleNetAK4JetTags:probb', 'pfParticleNetAK4JetTags:probbb', 'pfParticleNetAK4JetTags:probc', 'pfParticleNetAK4JetTags:probcc','pfParticleNetAK4JetTags:probpu','pfParticleNetAK4JetTags:probg','pfParticleNetAK4JetTags:probuds','pfParticleNetAK4JetTags:probundef','pfParticleNetAK4DiscriminatorsJetTags:BvsAll','pfParticleNetAK4DiscriminatorsJetTags:CvsL','pfParticleNetAK4DiscriminatorsJetTags:QvsG','pfParticleNetAK4DiscriminatorsJetTags:CvsB']

            process.akCs2PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "updatedPatJets", jetName = 'akCs2PF', genjetTag = "ak2GenJetsNoNu")      
            #process.akCs2PFJetAnalyzer.jetPtMin = cms.double(70.0)
            process.akCs2PFJetAnalyzer.doCandidateBtagging = True  
            process.akCs2PFJetAnalyzer.unsubjet_map = cms.InputTag("unsubJets")
            process.akCs2PFJetAnalyzer.matchJets = cms.untracked.bool(True)
            process.akCs2PFJetAnalyzer.matchTag = cms.untracked.InputTag("ak2PFpatJets")
            process.forest += process.aggregatedJets * process.unsubJets *process.pfParticleNetAK4TagInfos * process.pfParticleNetAK4JetTags * process.pfParticleNetAK4DiscriminatorsJetTags *process.updatedPatJetCorrFactors * process.updatedPatJets * process.akCs2PFJetAnalyzer
            process.akCs2PFpatJets.addTagInfos = False
            process.akCs2PFpatJets.addBTagInfo = False
            
            #process.ak2PFJetAnalyzer = process.akCs2PFJetAnalyzer.clone(jetTag = "ak2PFpatJets", jetName = 'ak2PF', isMC = False, saveRawPt = False)
            process.ak2PFpatJets.tagInfoSources = ["pfImpactParameterTagInfos", "pfInclusiveSecondaryVertexFinderTagInfos"]
            if doTracks:
                process.akCs2PFJetAnalyzer.doTracks = cms.untracked.bool(True)
                process.akCs2PFJetAnalyzer.ipTagInfoLabel = cms.untracked.string("pfImpactParameter")
            if doSvtx:
                process.akCs2PFJetAnalyzer.doSvtx = cms.untracked.bool(True)
                process.akCs2PFJetAnalyzer.svTagInfoLabel = cms.untracked.string("pfInclusiveSecondaryVertexFinder")

            #process.forest += process.ak2PFJetAnalyzer

    
#########################
# Event Selection
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.pAna = cms.EndPath(process.skimanalysis)

