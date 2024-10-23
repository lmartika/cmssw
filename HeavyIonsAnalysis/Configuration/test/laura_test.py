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
        'file:/eos/user/l/lamartik/testsamples/pbpbbjet2018/043213d2-944a-4e18-b1b5-ef71e93ef850.root'
    ),
)

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
#    input = cms.untracked.int32(-1)
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



###############################################################################
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
    process.genJetSequence
    )

#customisation

addR2Jets = True

addCandidateTagging = True
doTracks = True
doSvtx = True

doGenAnalysis = True
runAggregation = True

if addR2Jets :
    process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
    process.load("RecoHI.HiJetAlgos.EventConstSub_cfi")
    process.forest += process.extraJetsMC

    from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupHeavyIonJets

    if addR2Jets :
        process.jetsR2 = cms.Sequence()
        setupHeavyIonJets('akCs2PF', process.jetsR2, process, isMC = 1, radius = 0.20, JECTag = 'AK2PF')
#        process.akCs2PFpatJetCorrFactors.levels = ['L2Relative', 'L3Absolute']   # if JECs wanted
        process.akCs2PFpatJetCorrFactors.levels =  cms.vstring()
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAODuncorrjet_cff")
        process.akCs2PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "akCs2PFpatJets", jetName = 'akCs2PF', genjetTag = "ak2GenJetsNoNu")      
        process.forest += process.jetsR2 * process.akCs2PFJetAnalyzer
        process.akCs2PFJetAnalyzer.jetPtMin = cms.double(70.0)

# To apply PNET we have to remove the JEC first.  We could have just done this before the JEC is applied, but I copied this from the pp workflow
# Maybe something to clean up for the future
# Right now the problem is that B tagging is hard-coded to take updatedJets so maybe we will need the new collection after all?

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

ipTagInfoLabel_ = "pfImpactParameter"
svTagInfoLabel_ = "pfInclusiveSecondaryVertexFinder"

if addCandidateTagging:
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAODuncorrjet_cff")

    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = ('AK2PF', cms.vstring(), 'None'),
        btagDiscriminators = ['pfCombinedSecondaryVertexV2BJetTags', 'pfDeepCSVDiscriminatorsJetTags:BvsAll', 'pfDeepCSVDiscriminatorsJetTags:CvsB', 'pfDeepCSVDiscriminatorsJetTags:CvsL'], ## to add discriminators,
        btagPrefix = 'TEST',
    )
    if addR2Jets : process.updatedPatJets.jetSource = 'akCs2PFpatJets'
    process.updatedPatJets.addJetCorrFactors = False

#SV needed for aggregation
    process.load("RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi")
    process.pfImpactParameterTagInfos.candidates  = "packedPFCandidates"
    process.pfImpactParameterTagInfos.primaryVertex = "offlineSlimmedPrimaryVertices"
    process.pfImpactParameterTagInfos.jets = "akCs2PFpatJets"

    process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi")
    process.pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = "slimmedSecondaryVertices"

    process.updatedPatJets.addTagInfos = True
    process.updatedPatJets.tagInfoSources = [ipTagInfoLabel_+"TagInfos", svTagInfoLabel_+"TagInfos"]
    process.updatedPatJets.discriminatorSources =  ['pfParticleNetAK4JetTags:probb', 'pfParticleNetAK4JetTags:probb', 'pfParticleNetAK4JetTags:probbb', 'pfParticleNetAK4JetTags:probc', 'pfParticleNetAK4JetTags:probcc','pfParticleNetAK4JetTags:probpu','pfParticleNetAK4JetTags:probg','pfParticleNetAK4JetTags:probuds','pfParticleNetAK4JetTags:probundef','pfParticleNetAK4DiscriminatorsJetTags:BvsAll','pfParticleNetAK4DiscriminatorsJetTags:CvsL','pfParticleNetAK4DiscriminatorsJetTags:QvsG','pfParticleNetAK4DiscriminatorsJetTags:CvsB']

    
    process.forest.insert(-1,
                          process.pfImpactParameterTagInfos *
                          process.pfInclusiveSecondaryVertexFinderTagInfos *
                          process.candidateBtagging*
                          process.updatedPatJets
                      )
    process.akCs2PFJetAnalyzer.jetTag = "updatedPatJets"
    process.akCs2PFJetAnalyzer.doCandidateBtagging = True
    
    process.akCs2PFJetAnalyzer.doSubJets = False   # saves subjet kinematics per-jet


if doTracks:
    process.akCs2PFJetAnalyzer.doTracks = cms.untracked.bool(True)
    process.akCs2PFJetAnalyzer.ipTagInfoLabel = cms.untracked.string(ipTagInfoLabel_)
if doSvtx:
    process.akCs2PFJetAnalyzer.doSvtx = cms.untracked.bool(True)
    process.akCs2PFJetAnalyzer.svTagInfoLabel = cms.untracked.string(svTagInfoLabel_)

if doGenAnalysis:      ## Track-Gen-matches, try to use AK2    
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

    process.TrackToGenParticleMapProducer.jetSrc = cms.InputTag("updatedPatJets")  
    process.TrackToGenParticleMapProducer.genParticleSrc = cms.InputTag(taggedGenParticlesName_)
    process.forest.insert(-1,process.TrackToGenParticleMapProducer)

# At the moment together with maps; debug how to deliver candidates5B
    if runAggregation:
            process.load("RecoHI.HiJetAlgos.aggregatedPFCollection_cfi")
        ##    process.aggregatedPFCollection.tmva_path = ""
            process.aggregatedPFCands.aggregateHF = True
            process.aggregatedPFCands.candToGenParticleMap = ["TrackToGenParticleMapProducer", "trackToGenParticleMap"]
            process.forest.insert(-1,process.aggregatedPFCands)

            # TODO: change source collection
            process.aggregatedJets = cms.Sequence()
            setupHeavyIonJets('akCs2PF2', process.aggregatedJets, process, isMC = 1, radius = 0.20, JECTag = 'AK2PF')
           
            process.akCs2PF2Jets.src = 'aggregatedPFCands'      # To be debugged
#            process.akCs2PF2Jets.src = 'packedPFCandidates'   # This is the default thing

            process.akCs2PF2patJetCorrFactors.levels =  cms.vstring()
            process.akCs2PF2JetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "akCs2PF2patJets", jetName = 'akCs2PF2', genjetTag = "ak2GenJetsNoNu")      
            process.akCs2PF2JetAnalyzer.jetPtMin = cms.double(70.0)
#           process.forest += process.aggregatedJets * process.akCs2PF2JetAnalyzer
            process.forest.insert(-1,process.aggregatedJets * process.akCs2PF2JetAnalyzer)

    
#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.pAna = cms.EndPath(process.skimanalysis)
