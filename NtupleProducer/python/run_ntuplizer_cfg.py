import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v21')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Choose the file you want to run on
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
 	'/store/user/wangjian/BptoMMP_bbbar/RunIIAutumn18MiniAOD-102X/200615_142023/0000/BPH-RunIIAutumn18MiniAOD-00158_1.root'
    )
)

process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("histo.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.demo = cms.EDAnalyzer('miniaod_ntuplizer',
    genParticles = cms.InputTag("prunedGenParticles"),
    muons = cms.InputTag("slimmedMuons"),
    #vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    #objects = cms.InputTag("selectedPatTrigger"),
    #trks = cms.InputTag("generalTracks"),
    pfCands = cms.InputTag("packedPFCandidates"),
)

process.p = cms.Path(process.demo)
