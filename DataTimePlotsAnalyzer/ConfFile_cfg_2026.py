import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "file:/afs/cern.ch/work/s/ssantpur/public/ucsb/hlt/menufor2026/CMSSW_16_0_0_pre4/src/SteamRatesEdmWorkflow/Prod/hlt.root"
    )
)

# NOTE: in newer CMSSW, use TryToContinue, not SkipEvent
process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring("ProductNotFound")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out.root")
)

process.demo = cms.EDAnalyzer(
    "DataTimePlotsAnalyzer",
    jets        = cms.InputTag("hltCentralCaloJetptLowPtCollectionProducerSingle", "", "MYHLT"),
    jetTimes    = cms.InputTag("hltCaloJetTimingProducerSingle", "", "MYHLT"),
    jetCells    = cms.InputTag("hltCaloJetTimingProducerSingle", "jetCellsForTiming", "MYHLT"),
    jetEmEnergy = cms.InputTag("hltCaloJetTimingProducerSingle", "jetEcalEtForTiming", "MYHLT"),
    triggerString = cms.string("HLT_L1Tau_DelayedJet40_SingleDelayNeg0p75nsInclusive_BeamHalo"),
)

process.p = cms.Path(process.demo)
