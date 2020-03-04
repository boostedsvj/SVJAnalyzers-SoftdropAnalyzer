#include <memory>
#include <vector>
#include <cstdlib>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

/*
This producer takes a complicated GenParticle collection by reading it as a view,
and turns it into a std::vector<reco::GenParticle>.

It is not designed to be cpu, memory, or storage-efficient and should be used only
for small scale studies.
*/
class DumpGenParticlesForJetsNoNu : public edm::stream::EDProducer<> {
    public:
        explicit DumpGenParticlesForJetsNoNu(const edm::ParameterSet&);
        ~DumpGenParticlesForJetsNoNu() {}
        typedef edm::View<reco::Candidate> CandidateView;
        // typedef edm::View<reco::GenParticle> GenParticleView;
        typedef std::vector<reco::GenParticle> GenParticleVector;

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        // edm::EDGetTokenT<GenParticleView> inputParticleToken_;
        edm::EDGetTokenT<CandidateView> inputParticleToken_;
    };

DumpGenParticlesForJetsNoNu::DumpGenParticlesForJetsNoNu(const edm::ParameterSet& iConfig) :
    // inputParticleToken_(consumes<GenParticleView>(iConfig.getParameter<edm::InputTag>("PartTag")))
    inputParticleToken_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("PartTag")))
    {
    produces<GenParticleVector>();
    }

void DumpGenParticlesForJetsNoNu::produce(edm::Event& iEvent, const edm::EventSetup&) {  
    std::unique_ptr<GenParticleVector> outputGenParticles(new GenParticleVector);

    // edm::Handle<GenParticleView> inputParticleHandle;
    edm::Handle<CandidateView> inputParticleHandle;
    iEvent.getByToken(inputParticleToken_, inputParticleHandle);

    for (auto const & candidate : inputParticleHandle->ptrs()){
        std::cout << typeid(candidate).name() << std::endl;
        reco::GenParticle particle(
            candidate->charge(), candidate->p4(), candidate->vertex(),
            candidate->pdgId(), candidate->status(), candidate->charge()
            );
        outputGenParticles->push_back(particle);
        }

    iEvent.put(std::move(outputGenParticles));
    }

DEFINE_FWK_MODULE(DumpGenParticlesForJetsNoNu);