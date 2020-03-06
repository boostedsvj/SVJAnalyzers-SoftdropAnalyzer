#include <memory>
#include <vector>
#include <cstdlib>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"

// #include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        typedef edm::View<reco::Candidate> CandidateView;
        typedef std::vector<reco::GenParticle> GenParticleVector;

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        edm::EDGetTokenT<CandidateView> inputParticleToken_;
        bool onlyFromZPrime_;
    };

DumpGenParticlesForJetsNoNu::DumpGenParticlesForJetsNoNu(const edm::ParameterSet& iConfig) :
    inputParticleToken_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("PartTag"))),
    onlyFromZPrime_(iConfig.getParameter<bool>("onlyFromZPrime"))
    {
    produces<GenParticleVector>();
    }

void DumpGenParticlesForJetsNoNu::produce(edm::Event& iEvent, const edm::EventSetup&) {  
    std::unique_ptr<GenParticleVector> outputGenParticles(new GenParticleVector);
    edm::Handle<CandidateView> inputParticleHandle;
    iEvent.getByToken(inputParticleToken_, inputParticleHandle);

    for (auto const & candidate : inputParticleHandle->ptrs()){
        if (onlyFromZPrime_){
            if (candidate->numberOfMothers() > 0){
                const reco::Candidate * mother = candidate->mother();
                // Recursively find the last mother
                while (mother->numberOfMothers() > 0){
                    if (mother->pdgId() == 4900023){
                        reco::GenParticle particle(
                            candidate->charge(), candidate->p4(), candidate->vertex(),
                            candidate->pdgId(), candidate->status(), candidate->charge()
                            );
                        outputGenParticles->push_back(particle);
                        break;
                        }
                    mother = mother->mother();
                    }
                }
            }
        else {
            // Just push the candidate regardless of what mother it came from
            reco::GenParticle particle(
                candidate->charge(), candidate->p4(), candidate->vertex(),
                candidate->pdgId(), candidate->status(), candidate->charge()
                );
            outputGenParticles->push_back(particle);
            }
        }

    iEvent.put(std::move(outputGenParticles));
    }

void DumpGenParticlesForJetsNoNu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("PartTag", edm::InputTag("genParticlesForJets"));
    desc.add<bool>("onlyFromZPrime", false);    
    descriptions.add("DumpGenParticlesForJetsNoNu", desc);
    }

DEFINE_FWK_MODULE(DumpGenParticlesForJetsNoNu);