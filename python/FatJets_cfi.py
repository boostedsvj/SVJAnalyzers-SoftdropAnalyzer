import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

ak15GenJets = cms.EDProducer(
    "FastjetJetProducer",
    GenJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(1.5)
    )
ak15GenJetsNoNu = ak15GenJets.clone(src = cms.InputTag("genParticlesForJetsNoNu"))

# import RecoJets.Configuration.RecoGenJets_cff
# print RecoJets.Configuration.RecoGenJets_cff
# print dir(RecoJets.Configuration.RecoGenJets_cff)
# print vars(RecoJets.Configuration.RecoGenJets_cff)

from RecoJets.Configuration.RecoGenJets_cff import *

import RecoJets.Configuration.RecoGenJets_cff

# recoGenJetsTask.add(ak15GenJets)
# recoAllGenJetsTask.add(ak15GenJets)
# recoAllGenJetsNoNuTask.add(ak15GenJets)

recoGenJets        += ak15GenJets
recoAllGenJets     += ak15GenJets
recoAllGenJetsNoNu += ak15GenJets
