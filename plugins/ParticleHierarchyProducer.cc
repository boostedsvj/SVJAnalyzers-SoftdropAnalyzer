#include <memory>
#include <vector>
#include <cstdlib>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/SVJFormats/interface/GenParticlePlus.h"

class ParticleHierarchyProducer : public edm::stream::EDProducer<> {
    public:
        explicit ParticleHierarchyProducer(const edm::ParameterSet&);
        ~ParticleHierarchyProducer() {}
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        typedef edm::View<reco::GenParticle> GenParticleView;
        typedef std::vector<GenParticlePlus> GenParticlePlusVector;
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        edm::EDGetTokenT<GenParticleView> particleToken_;
        bool verbose_;
    };

ParticleHierarchyProducer::ParticleHierarchyProducer(const edm::ParameterSet& iConfig) :
    particleToken_(consumes<GenParticleView>(iConfig.getParameter<edm::InputTag>("particleTag"))),
    verbose_(iConfig.getParameter<bool>("verbose"))
    {
    produces<GenParticlePlusVector>();
    }

void ParticleHierarchyProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {  
    std::unique_ptr<GenParticlePlusVector> genParticlesWithHierarchy(new GenParticlePlusVector);
    edm::Handle<GenParticleView> particleHandle;
    iEvent.getByToken(particleToken_, particleHandle);
    for (auto const & particle : particleHandle->ptrs() ) {
        if (particle->pdgId() == 4900023 && particle->isLastCopy() == 1){
            GenParticlePlus zprime(*particle);
            zprime.setQuarksAndFinalProducts();
            genParticlesWithHierarchy->push_back(zprime);
            if (verbose_){zprime.printDecayTree();}
            }
        }    
    iEvent.put(std::move(genParticlesWithHierarchy));
    }

void ParticleHierarchyProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("particleTag", edm::InputTag("genParticles"));
    desc.add<bool>("verbose", false);
    descriptions.add("ParticleHierarchyProducer", desc);
    }

DEFINE_FWK_MODULE(ParticleHierarchyProducer);