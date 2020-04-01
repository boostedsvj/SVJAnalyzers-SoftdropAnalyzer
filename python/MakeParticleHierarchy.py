import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.inputFiles = 'file:susyv1_mmed450_rinv03_0.root'
options.outputFile = 'particlehierarchy.root'
options.maxEvents = 300
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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

process.particleHierarchy = cms.EDProducer(
    'ParticleHierarchyProducer',
    verbose = cms.bool(False)
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
    'keep *_genParticles_*_*',
    'keep *_particleHierarchy_*_*',
    ])

# Path and EndPath definitions
process.hierarchy_step = cms.Path(process.particleHierarchy)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

# Schedule definition
process.schedule = cms.Schedule(process.hierarchy_step)
process.schedule.extend([process.endjob_step, process.FEVTDEBUGoutput_step])

# # #Setup FWK for multithreaded
# process.options.numberOfThreads=cms.untracked.uint32(4)
# process.options.numberOfStreams=cms.untracked.uint32(0)
