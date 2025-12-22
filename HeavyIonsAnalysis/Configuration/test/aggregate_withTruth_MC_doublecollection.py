### HiForest Configuration
# Input: miniAOD
# Type: mc

doRun2 = True
import FWCore.ParameterSet.Config as cms

if doRun2:
    from Configuration.Eras.Era_Run2_2018_pp_on_AA_cff import Run2_2018_pp_on_AA
    from Configuration.ProcessModifiers.run2_miniAOD_pp_on_AA_103X_cff import run2_miniAOD_pp_on_AA_103X
    process = cms.Process('HiForest', Run2_2018_pp_on_AA,run2_miniAOD_pp_on_AA_103X)
else:
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
        '/store/group/phys_heavyions/jviinika/PythiaHydjetRun3_5p36TeV_dijet_ptHat15_100kEvents_miniAOD_2023_08_30/PythiaHydjetDijetRun3/PythiaHydjetRun3_dijet_ptHat15_5p36TeV_miniAOD/230830_165931/0000/pythiaHydjet_miniAOD_11.root'
    ),                            
)

if doRun2:
#    process.source.fileNames = cms.untracked.vstring('/store/himc/HINPbPbSpring21MiniAOD/Bjet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/MINIAODSIM/FixL1CaloGT_New_Release_112X_upgrade2018_realistic_HI_v9-v1/260000/0af4190a-74ea-4a92-a88e-96ee753d8ccb.root')
    process.source.fileNames = cms.untracked.vstring('/store/group/phys_heavyions/lamartik/bjet/043213d2-944a-4e18-b1b5-ef71e93ef850.root')
#    process.source.eventsToProcess = cms.untracked.VEventRange("1:1251612-1:1251612"); # example of B splitting in reco level, centr.161
    #    process.source.fileNames = cms.untracked.vstring('/store/group/phys_heavyions/lamartik/bjet/043213d2-944a-4e18-b1b5-ef71e93ef850.root')
#    process.source.eventsToProcess = cms.untracked.VEventRange("1:334986-1:334986"); # example of B splitting in reco level, centr.161
#    process.source.eventsToProcess = cms.untracked.VEventRange("1:1758572-1:1758572");  # cent 107, split + lost by TMVA
#    process.source.eventsToProcess = cms.untracked.VEventRange("1:1769063-1:1769063");  # cent 46, seems like b split over several jets, 2 b jet
#    process.source.eventsToProcess = cms.untracked.VEventRange("1:1758509-1:1758509");   # cent 2; no all bs in jets, tmva dropping 3/4 daughters
#    process.source.eventsToProcess = cms.untracked.VEventRange("1:1761370-1:1761370");   # cent 80


   
# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
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
if doRun2:   # auto:phase1_2018_realistic_hi
#    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_hi', '')
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic_hi', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("JPcalib_MC103X_2018PbPb_v4"),
             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
         )
])
if doRun2:
    process.GlobalTag.toGet.extend([
        cms.PSet(record = cms.string("GEMRecoGeometryRcd"),
                 tag = cms.string("GEMRECO_Geometry_131DD4hepV1_mc_v2"),
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

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_mc_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

#from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_mc
#process.hltobject.triggerNames = trigger_list_mc

################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_mc_cff')

###############################################################################

###############################
# rho and random cone stuff
process.load("HeavyIonsAnalysis.JetAnalysis.RhoAnalysis_cff")
process.load("HeavyIonsAnalysis.JetAnalysis.RandomConeAnalysis_cff")
process.load("RecoHI.HiJetAlgos.hiFJRhoFlowModulationProducer_cfi")
###############################

# ZDC analyzer
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018Producer_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018RecHit_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.zdcanalyzer_cfi')

process.zdcanalyzer.doZDCRecHit = False
process.zdcanalyzer.doZDCDigi = True
process.zdcanalyzer.zdcRecHitSrc = cms.InputTag("QWzdcreco")
process.zdcanalyzer.zdcDigiSrc = cms.InputTag("hcalDigis", "ZDC")
process.zdcanalyzer.calZDCDigi = False
process.zdcanalyzer.verbose = False

process.genJetSequence = cms.Sequence()
###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.hltanalysis +
#    process.hltobject +
#    process.l1object +
#    process.trackSequencePbPb +
#    process.particleFlowAnalyser +
    process.genJetSequence + 
    process.hiEvtAnalyzer #+
#    process.HiGenParticleAna +
#    process.ggHiNtuplizer +
#    process.zdcdigi +
#    process.QWzdcreco +
    #process.akCs4PFJetAnalyzer +
#    process.zdcanalyzer #+
#    process.unpackedMuons +
#    process.muonAnalyzer
    )

#customisation
process.particleFlowAnalyser.ptMin = 0.0

# Gen
process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.genJetSequence += process.mergedGenParticles
## Produces a reco::GenParticleCollection named mergedGenParticles

process.load("RecoHI.HiJetAlgos.HFdecayProductTagger_cfi")
process.HFdecayProductTagger.genParticles = cms.InputTag("mergedGenParticles")
process.HFdecayProductTagger.tagBorC = cms.bool(True) # tag B
process.genJetSequence += process.HFdecayProductTagger

taggedGenParticlesName_ = "HFdecayProductTagger"

process.bDecayAna = process.HiGenParticleAna.clone(
    genParticleSrc = cms.InputTag(taggedGenParticlesName_),
    useRefVector = cms.untracked.bool(False),
    partonMEOnly = cms.untracked.bool(False),
    chargedOnly = True, 
    doHI = False,
    etaMax = cms.untracked.double(30),
    ptMin = cms.untracked.double(0),
    stableOnly = False
)
# process.genJetSequence += process.bDecayAna

# Select the types of jets filled
matchJets = True             # Enables q/g and heavy flavor jet identification in MC
jetPtMin = 15
jetAbsEtaMax = 2.5
doCaloJets = False

# Toggle on/off for saving track info
doTracks = True
doSvtx = False

runAggregation = True

# Choose which additional information is added to jet trees
doHIJetID = False             # Fill jet ID and composition information branches
doWTARecluster = False        # Add jet phi and eta for WTA axis
doBtagging  =  True         # Note that setting to True increases computing time a lot

# 0 means use original mini-AOD jets, otherwise use R value, e.g., 3,4,8
jetLabel = "2"

# add candidate tagging, copy/paste to add other jet radii
from HeavyIonsAnalysis.JetAnalysis.deepNtupleSettingsFullAggregation_cff import candidateBtaggingMiniAOD
candidateBtaggingMiniAOD(process, isMC = True, jetPtMin = jetPtMin, jetCorrLevels = ['L2Relative', 'L3Absolute'], doBtagging = doBtagging, labelR = jetLabel, runAggregation = runAggregation)

# setup jet analyzer for inclusive jets
setattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive",process.akCs4PFJetAnalyzer.clone())
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetTag = 'patJetsAKCs'+jetLabel+'PFNotAggrJets'
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetName = 'akCs'+jetLabel+'PF'
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").genjetTag = "ak"+jetLabel+"GenJetsRecluster"
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").fillGenJets = cms.untracked.bool(False)   # save only gen matched to reco jets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").matchJets = matchJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").matchTag = 'patJetsAK'+jetLabel+'PFUnsubJets'

## Matching between the two CS jet collections
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").matchOriginal = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").originalCSTag = 'patJetsAKCs'+jetLabel+'PFNotAggrJets'

getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doHiJetID = doHIJetID
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").runSubstructure = cms.untracked.bool(True)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doWTARecluster = doWTARecluster
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doCaloJets = doCaloJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetPtMin = 50
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").genPtMin = 80

getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").latektcut = cms.untracked.double(1)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").SDktcut = cms.untracked.double(1)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").SDcut = cms.untracked.double(0.2)

getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").useRawPt = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doPFjetID = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").rParam = int(jetLabel)*0.1
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetFlavourInfos = "ak"+jetLabel+"PFUnsubJetFlavourInfos"
if doBtagging:
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").useNewBtaggers = True
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").pfJetProbabilityBJetTag = cms.untracked.string("pfJetProbabilityBJetTagsDeepFlavour") 
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").pfUnifiedParticleTransformerAK4JetTags = cms.untracked.string("pfUnifiedParticleTransformerAK4JetTagsDeepFlavour")

if doTracks:
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doTracks = cms.untracked.bool(True)
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").ipTagInfoLabel = cms.untracked.string("pfImpactParameter")
if doSvtx:
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doSvtx = cms.untracked.bool(True)
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").svTagInfoLabel = cms.untracked.string("pfInclusiveSecondaryVertexFinder")

process.forest += getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive")

# setup jet analyzer for aggregated jets
setattr(process,"akCs"+jetLabel+"PFJetAnalyzer",process.akCs4PFJetAnalyzer.clone())
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetTag = 'selectedUpdatedPatJetsDeepFlavour'
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetName = 'akCs'+jetLabel+'PF'
#getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").genjetTag = "ak"+jetLabel+"GenJetsWithNu"
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").genjetTag = "ak"+jetLabel+"aggregatedGenJetsNoNu"
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").fillGenJets = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").matchJets = matchJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").matchTag = 'patJetsAK'+jetLabel+'PFUnsubJets'

getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").matchOriginal = cms.untracked.bool(True)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").originalCSTag = 'patJetsAKCs'+jetLabel+'PFNotAggrJets'

getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doHiJetID = doHIJetID
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").runSubstructure = cms.untracked.bool(True)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doWTARecluster = doWTARecluster
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doCaloJets = doCaloJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetPtMin = 50
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").genPtMin = 80

getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").latektcut = cms.untracked.double(1)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").SDktcut = cms.untracked.double(1)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").SDcut = cms.untracked.double(0.2)

getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").useRawPt = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doPFjetID = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").rParam = int(jetLabel)*0.1
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetFlavourInfos = "ak"+jetLabel+"PFUnsubJetFlavourInfos"
if doBtagging:
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").useNewBtaggers = True
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").pfJetProbabilityBJetTag = cms.untracked.string("pfJetProbabilityBJetTagsDeepFlavour") 
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").pfUnifiedParticleTransformerAK4JetTags = cms.untracked.string("pfUnifiedParticleTransformerAK4JetTagsDeepFlavour")

if doTracks:
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doTracks = cms.untracked.bool(True)
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").ipTagInfoLabel = cms.untracked.string("pfImpactParameter")
if doSvtx:
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doSvtx = cms.untracked.bool(True)
    getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").svTagInfoLabel = cms.untracked.string("pfInclusiveSecondaryVertexFinder")


process.ak2PFXpatJets = cms.EDFilter("PatJetXSelector",
                                     src = cms.InputTag("selectedUpdatedPatJetsDeepFlavour"),
                                     cut = cms.string("pt > 0.0 && abs(rapidity()) < 3.")
                                     )
process.forest += getattr(process,"ak"+jetLabel+"PFXpatJets")

#process.akCs2PFJets.writeJetsWithConst = False

process.akCs2PFJetAnalyzer.jetTag = "ak2PFXpatJets"
process.forest += getattr(process,"akCs"+jetLabel+"PFJetAnalyzer")


#########################
# rho stuff, maybe not all needed
# process.forest += process.hiFJRhoFlowModulationProducer * process.rhoAnalysis * process.randomConeAnalysisR4 * process.randomConeAnalysisR2
process.forest += process.rhoAnalysis
#########################


#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)
process.pAna = cms.EndPath(process.skimanalysis)


#### Tag infos
process.patJetsAK2PFUnsubJets.addBTagInfo = True
process.patJetsAK2PFUnsubJets.addTagInfos = True
process.patJetsAK2PFUnsubJets.tagInfoSources = cms.VInputTag(["pfInclusiveSecondaryVertexFinderTagInfos","pfImpactParameterTagInfos"])

process.patJetsAK2PFUnsubJets.addDiscriminators = False

### true aggregation?
process.aggregatedPFCands.aggregateWithTruthInfo = True
process.aggregatedPFCands.aggregateWithCuts = False
process.aggregatedPFCands.aggregateWithTMVA = False

