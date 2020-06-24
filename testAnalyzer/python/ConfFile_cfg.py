import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/h/hajheng/private/CMSSW/CMSSW_9_4_13/src/DY1Jet_test.root'
    )
)

process.demo = cms.EDAnalyzer('testAnalyzer',
    outputFile = cms.string('/data4/hjheng/gen/DY1Jet_test_ntuple.root'),
    genParticleSrc = cms.InputTag("genParticles"),
    ak4GenJetSrc = cms.InputTag("ak4GenJets"),
    generatorLabel = cms.InputTag("generator")
)


process.p = cms.Path(process.demo)
