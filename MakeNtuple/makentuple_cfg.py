import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("METSigTuning.MakeNtuple.makentuple_cfi")
#process.load("RecoMET/METProducers.METSignificanceObjects_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/040D9AF7-FB6B-E411-8106-0025907DBA06.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/04442001-036C-E411-9C90-0025901D42C0.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08501832-106C-E411-BCEE-0025904B1420.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0AA8983B-046C-E411-B0D3-0025901D4C3E.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0C924143-0D6C-E411-B48A-0025907BAF70.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CB30B20-186C-E411-A85E-002590AC4BF6.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0CF38497-006C-E411-A713-0025907DCA0C.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E2BE86D-1E6C-E411-9FD4-003048D4DEAE.root",
       #"/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0E50D3F3-056C-E411-8DF6-0025907DCA7E.root"
       #'/store/relval/CMSSW_7_3_0/RelValZMM_13/MINIAODSIM/MCRUN2_73_V7-v1/00000/127CA68E-8981-E411-A524-002590593872.root',
       #'/store/relval/CMSSW_7_3_0/RelValZMM_13/MINIAODSIM/MCRUN2_73_V7-v1/00000/56FE228D-8981-E411-9AD8-0025905A6126.root'
       'file:../../../../../../../../work/s/shtan/private/0432E62A-7A6C-E411-87BB-002590DB92A8.root'
    )
)

process.demo = cms.EDAnalyzer('MakeNtuple',
      src = cms.InputTag("packedPFCandidates"),
      jets = cms.InputTag("slimmedJets"),
      leptons = cms.VInputTag("slimmedElectrons", "slimmedMuons", "slimmedPhotons"),
      met = cms.InputTag("slimmedMETs"),
      muons = cms.InputTag("slimmedMuons"),
      vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)

process.options = cms.untracked.PSet(
      SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.p = cms.Path(
      #process.selectionSequenceForMETSig * 
      process.demo
      )
