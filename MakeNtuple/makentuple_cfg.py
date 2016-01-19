import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("METSigTuning.MakeNtuple.makentuple_cfi")
#process.load("RecoMET/METProducers.METSignificanceObjects_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

##process.source = cms.Source("PoolSource",
##    fileNames = cms.untracked.vstring(
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
##       'file:../../../../../../../../work/s/shtan/private/0432E62A-7A6C-E411-87BB-002590DB92A8.root'
##    )
##)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/B67CF3CA-EAFB-E411-8B79-52540034D38F.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/0A5482D9-D1FB-E411-B007-AC853D9DACE3.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/24ACD574-08FB-E411-A21A-001E67398223.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/0032B047-B7FB-E411-AF23-0025905280BE.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/40A870F8-13FA-E411-BFFD-0025907FD430.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/0E4CEBFE-ECFB-E411-9F0C-842B2B29273C.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_170to250_bcToE_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/020A1AEC-4A02-E511-9877-00259073E456.root'))

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/221B185A-3502-E511-AC72-0025905C4262.root'))

process.demo = cms.EDAnalyzer('MakeNtuple',
      src = cms.InputTag("packedPFCandidates"),
      jets = cms.InputTag("slimmedJets"),
      leptons = cms.VInputTag("slimmedElectrons", "slimmedMuons", "slimmedPhotons"),
      met = cms.InputTag("slimmedMETs"),
      muons = cms.InputTag("slimmedMuons"),
      vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
      genjets = cms.InputTag("slimmedGenJets")
)

process.options = cms.untracked.PSet(
      SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.p = cms.Path(
      #process.selectionSequenceForMETSig * 
      process.demo
      )
