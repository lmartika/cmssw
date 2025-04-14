### HiForest Configuration
# Input: miniAOD
# Type: data

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
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")

###############################################################################

# input files
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        '/store/group/phys_heavyions/jviinika/PythiaHydjetRun3_5p36TeV_dijet_ptHat15_100kEvents_miniAOD_2023_08_30/PythiaHydjetDijetRun3/PythiaHydjetRun3_dijet_ptHat15_5p36TeV_miniAOD/230830_165931/0000/pythiaHydjet_miniAOD_11.root'
    ),                            
)

if doRun2:
    process.source.fileNames = cms.untracked.vstring('/store/hidata/HIRun2018A/HIHardProbes/MINIAOD/PbPb18_MiniAODv1-v1/00000/034807c1-0ae5-4540-bb81-80ab2f3bc01d.root')


# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(4000)
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
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data_promptlike_hi', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
             tag = cms.string("JPcalib_Data103X_2018PbPb_v1"),
             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
         )
])

###############################################################################

# root output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("HiForestMiniAOD_data.root"))

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
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')

###############################################################################

###############################
# rho and random cone stuff
process.load("HeavyIonsAnalysis.JetAnalysis.RhoAnalysis_cff")
process.load("HeavyIonsAnalysis.JetAnalysis.RandomConeAnalysis_cff")
###############################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.hltanalysis +
#    process.hltobject +
#    process.l1object +
#    process.trackSequencePbPb +
#    process.particleFlowAnalyser +
    process.hiEvtAnalyzer #+

    )

# Select the types of jets filled
matchJets = True             # Enables q/g and heavy flavor jet identification in MC
jetPtMin = 15
jetAbsEtaMax = 2.5
doCaloJets = False

# Toggle on/off for saving track info
doTracks = False
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
candidateBtaggingMiniAOD(process, isMC = False, jetPtMin = jetPtMin, jetCorrLevels = ['L2Relative', 'L3Absolute'], doBtagging = doBtagging, labelR = jetLabel, runAggregation = runAggregation)

# setup jet analyzers
#### Inclusive jets
setattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive",process.akCs4PFJetAnalyzer.clone())
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetTag = 'patJetsAKCs'+jetLabel+'PFNotAggrJets'
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetName = 'akCs'+jetLabel+'PF'
#getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").genjetTag = "ak"+jetLabel+"GenJetsWithNu"
#getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").genjetTag = "ak"+jetLabel+"aggregatedGenJetsRecluster" 
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").matchJets = matchJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").matchTag = 'patJetsAK'+jetLabel+'PFUnsubJets'

## Matching between the two CS jet collections
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").matchOriginal = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").originalCSTag = 'patJetsAKCs'+jetLabel+'PFNotAggrJets'

getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doHiJetID = doHIJetID
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").runSubstructure = cms.untracked.bool(True)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doWTARecluster = doWTARecluster
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doCaloJets = doCaloJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetPtMin = 70
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").useRawPt = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").doPFjetID = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzerInclusive").rParam = int(jetLabel)*0.1

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


#### These are b-jets
setattr(process,"akCs"+jetLabel+"PFJetAnalyzer",process.akCs4PFJetAnalyzer.clone())
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetTag = 'selectedUpdatedPatJetsDeepFlavour'
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetName = 'akCs'+jetLabel+'PF'
#getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").genjetTag = "ak"+jetLabel+"GenJetsWithNu"
#getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").genjetTag = "ak"+jetLabel+"aggregatedGenJetsRecluster" 
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").matchJets = matchJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").matchTag = 'patJetsAK'+jetLabel+'PFUnsubJets'

## Matching between the two CS jet collections
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").matchOriginal = cms.untracked.bool(True)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").originalCSTag = 'patJetsAKCs'+jetLabel+'PFNotAggrJets'

getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doHiJetID = doHIJetID
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").runSubstructure = cms.untracked.bool(True)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doWTARecluster = doWTARecluster
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doCaloJets = doCaloJets
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetPtMin = 70
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").useRawPt = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").doPFjetID = cms.untracked.bool(False)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").jetAbsEtaMax = cms.untracked.double(jetAbsEtaMax)
getattr(process,"akCs"+jetLabel+"PFJetAnalyzer").rParam = int(jetLabel)*0.1

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
