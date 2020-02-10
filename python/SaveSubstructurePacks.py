import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.inputFiles = 'file:N4000_0000_seed1001.root', 'file:N4000_0001_seed1002.root'
options.outputFile = 'substructurepacks.root'
options.maxEvents = -1
options.register(
    'coneSize',
    0.8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Cone size with which to cluster the genJets"
    )
options.parseArguments()

from Configuration.StandardSequences.Eras import eras
process = cms.Process('SoftDropGenJets',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
    )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

# ___________________________________
# Physics

if options.coneSize == 0.8:
    from RecoJets.Configuration.RecoGenJets_cff import ak8GenJetsNoNu
    genJetsNoNu = ak8GenJetsNoNu
else:
    from RecoJets.JetProducers.GenJetParameters_cfi import *
    from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
    genJets = cms.EDProducer(
        "FastjetJetProducer",
        GenJetParameters,
        AnomalousCellParameters,
        jetAlgorithm = cms.string("AntiKt"),
        rParam       = cms.double(options.coneSize)
        )
    genJetsNoNu = genJets.clone(src = cms.InputTag("genParticlesForJetsNoNu"))
    from RecoJets.Configuration.RecoGenJets_cff import *
    recoGenJets        += genJets
    recoAllGenJets     += genJets
    recoAllGenJetsNoNu += genJets

process.genJetsNoNuArea = genJetsNoNu.clone(
    doAreaFastjet = cms.bool(True),
    )

process.genJetsNoNuSoftDrop = genJetsNoNu.clone(
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(0.5),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    doAreaFastjet = cms.bool(True),
    )

process.htProducer = cms.EDProducer("GenHTProducer")

process.substructurePacks = cms.EDProducer(
    "SubstructureProducer",
    jetSrc = cms.InputTag("genJetsNoNuArea"),
    PartTag = cms.InputTag("genParticles"),
    distMax = cms.double(0.8),
    algoTags = cms.VInputTag(
        cms.InputTag("genJetsNoNuSoftDrop"),
        ),
    )

# Output definition
process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
        ),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
    )

process.FEVTDEBUGoutput.outputCommands.extend([
    'drop *',
    'keep *_substructurePacks_*_*',
    'keep *_genParticles_*_*',
    'keep *_genParticlesForJetsNoNu_*_*',
    'keep *_htProducer_*_*'
    ])

# Path and EndPath definitions
process.jet_step = cms.Path(
    process.genJetsNoNuArea
    +
    process.genJetsNoNuSoftDrop
    +
    process.htProducer
    +
    process.substructurePacks
    )

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

# Schedule definition
process.schedule = cms.Schedule(process.jet_step)
process.schedule.extend([process.endjob_step, process.FEVTDEBUGoutput_step])

# #Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)
