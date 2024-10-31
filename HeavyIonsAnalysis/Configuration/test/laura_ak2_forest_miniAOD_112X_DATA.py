### HiForest Configuration
# Input: miniAOD
# Type: data

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_pp_on_AA_cff import Run2_2018_pp_on_AA
from Configuration.ProcessModifiers.run2_miniAOD_pp_on_AA_103X_cff import run2_miniAOD_pp_on_AA_103X
process = cms.Process('HiForest', Run2_2018_pp_on_AA,run2_miniAOD_pp_on_AA_103X)

###############################################################################

# HiForest info
process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 112X, data")

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
#        "/store/hidata/HIRun2018A/HISingleMuon/MINIAOD/PbPb18_MiniAODv1-v1/00000/00345f79-641f-4002-baf1-19ae8e83c48b.root"
#        "/store/hidata/HIRun2018A/HIHardProbes/MINIAOD/PbPb18_MiniAODv1-v1/00000/034807c1-0ae5-4540-bb81-80ab2f3bc01d.root"
        "file:/eos/user/l/lamartik/testsamples/pbpbdata2018/034807c1-0ae5-4540-bb81-80ab2f3bc01d.root"
    ),
)
#input file produced from:
#"file:/afs/cern.ch/work/r/rbi/public/forest/HIHardProbes_HIRun2018A-PromptReco-v2_AOD.root"

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
    )

###############################################################################

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data_promptlike_hi', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag

centralityTag = "CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"
process.HiForestInfo.info.append(centralityTag)

print('\n')
print('\033[31m~*~ CENTRALITY TABLE FOR 2018 PBPB DATA ~*~\033[0m')
print('\033[36m~*~ TAG: ' + centralityTag + ' ~*~\033[0m')
print('\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(
        record = cms.string("HeavyIonRcd"),
        tag = cms.string(centralityTag),
        label = cms.untracked.string("HFtowers"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        ),
    ])

process.GlobalTag.toGet.extend([
    cms.PSet(
        record = cms.string("BTagTrackProbability3DRcd"),
        tag = cms.string("JPcalib_Data103X_2018PbPb_v1"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
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

# event analysis
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.l1object_cfi')

from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data
process.hltobject.triggerNames = trigger_list_data

process.load('HeavyIonsAnalysis.EventAnalysis.particleFlowAnalyser_cfi')
################################
# electrons, photons, muons
SSHIRun2018A = "HeavyIonsAnalysis/EGMAnalysis/data/SSHIRun2018A.dat"
process.load('HeavyIonsAnalysis.EGMAnalysis.correctedElectronProducer_cfi')
process.correctedElectrons.correctionFile = SSHIRun2018A

process.load('HeavyIonsAnalysis.MuonAnalysis.unpackedMuons_cfi')
process.load("HeavyIonsAnalysis.MuonAnalysis.muonAnalyzer_cfi")
process.load('HeavyIonsAnalysis.EGMAnalysis.ggHiNtuplizer_cfi')
process.ggHiNtuplizer.doMuons = cms.bool(False)
process.ggHiNtuplizer.electronSrc = "correctedElectrons"
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
################################
# jet reco sequence
process.load('HeavyIonsAnalysis.JetAnalysis.akCs4PFJetSequence_pponPbPb_data_cff')
################################
# tracks
process.load("HeavyIonsAnalysis.TrackAnalysis.TrackAnalyzers_cff")
###############################################################################

# ZDC RecHit Producer
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018Producer_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018RecHit_cfi')

process.load('HeavyIonsAnalysis.ZDCAnalysis.zdcanalyzer_cfi')
process.zdcanalyzer.doZDCRecHit = True
process.zdcanalyzer.doZDCDigi = False
process.zdcanalyzer.zdcRecHitSrc = cms.InputTag("QWzdcreco")
process.zdcanalyzer.calZDCDigi = True
################################

###############################################################################
# main forest sequence
process.forest = cms.Path(
    process.HiForestInfo +
    process.hltanalysis +
#    process.hltobject +
#    process.l1object +
#    process.trackSequencePbPb +
#    process.particleFlowAnalyser +
    process.hiEvtAnalyzer
#    process.unpackedMuons +
#    process.correctedElectrons +
#    process.ggHiNtuplizer +
#    process.zdcdigi +
#    process.QWzdcreco +
#    process.zdcanalyzer +
#    process.muonAnalyzer 
    )

#customisation

addR2Jets = True
addR4Jets = False

addCandidateTagging = True

doTracks = False
doSvtx = True

if addR2Jets or addR4Jets :
    process.load("HeavyIonsAnalysis.JetAnalysis.extraJets_cff")
    process.forest += process.extraJetsData
    
    from HeavyIonsAnalysis.JetAnalysis.clusterJetsFromMiniAOD_cff import setupHeavyIonJets

    if addR2Jets :
        process.jetsR2 = cms.Sequence()
        setupHeavyIonJets('akCs2PF', process.jetsR2, process, isMC = 0, radius = 0.20, JECTag = 'AK2PF')
        process.akCs2PFpatJetCorrFactors.levels = ['L2Relative', 'L2L3Residual']
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.akCs2PFJetAnalyzer = process.akCs4PFJetAnalyzer.clone(jetTag = "akCs2PFpatJets", jetName = 'akCs2PF')
        process.forest += process.jetsR2 * process.akCs2PFJetAnalyzer

    if addR4Jets :
        # Recluster using an alias "0" in order not to get mixed up with the default AK4 collections
        process.jetsR4 = cms.Sequence()
        setupHeavyIonJets('akCs0PF', process.jetsR4, process, isMC = 0, radius = 0.40, JECTag = 'AK4PF')
        process.akCs0PFpatJetCorrFactors.levels = ['L2Relative', 'L2L3Residual']
        process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
        process.akCs4PFJetAnalyzer.jetTag = 'akCs0PFpatJets'
        process.akCs4PFJetAnalyzer.jetName = 'akCs0PF'
        process.forest += process.jetsR4 * process.akCs4PFJetAnalyzer

# To apply PNET we have to remove the JEC first.  We could have just done this before the JEC is applied, but I copied this from the pp workflow                                 
# Maybe something to clean up for the future                                                                                                                                     
# Also note that it's only setup for R=0.4 for the moment.  

ipTagInfoLabel_ = "pfImpactParameter"
svTagInfoLabel_ = "pfInclusiveSecondaryVertexFinder"

if addCandidateTagging:
    process.load("HeavyIonsAnalysis.JetAnalysis.candidateBtaggingMiniAOD_cff")
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

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

    process.akCs2PFJetAnalyzer.jetPtMin = cms.double(70.0)
    process.akCs2PFJetAnalyzer.doSubJets = False   # saves subjet kinematics per-jet

if doTracks:
    process.akCs2PFJetAnalyzer.doTracks = cms.untracked.bool(True)
    process.akCs2PFJetAnalyzer.ipTagInfoLabel = cms.untracked.string(ipTagInfoLabel_)
if doSvtx:
    process.akCs2PFJetAnalyzer.doSvtx = cms.untracked.bool(True)
    process.akCs2PFJetAnalyzer.svTagInfoLabel = cms.untracked.string(svTagInfoLabel_)


#########################
# Event Selection -> add the needed filters here
#########################

process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.pclusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
process.pprimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.pphfCoincFilter2Th4 = cms.Path(process.phfCoincFilter2Th4)
process.pAna = cms.EndPath(process.skimanalysis)


