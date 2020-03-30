#include <memory>
#include <vector>
#include <string>
#include <cstdlib>
#include <set>
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
class DumpDMParticles : public edm::stream::EDProducer<> {
    public:
        explicit DumpDMParticles(const edm::ParameterSet&);
        ~DumpDMParticles() {}
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        typedef edm::View<reco::GenParticle> GenParticleView;
        typedef std::vector<reco::GenParticle> GenParticleVector;

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        edm::EDGetTokenT<GenParticleView> inputParticleToken_;
        bool onlyFromZPrime_;
        bool verbose_;
    };

DumpDMParticles::DumpDMParticles(const edm::ParameterSet& iConfig) :
    inputParticleToken_(consumes<GenParticleView>(iConfig.getParameter<edm::InputTag>("PartTag"))),
    onlyFromZPrime_(iConfig.getParameter<bool>("onlyFromZPrime")),
    verbose_(iConfig.getParameter<bool>("verbose"))
    {
    produces<GenParticleVector>();
    }

void DumpDMParticles::produce(edm::Event& iEvent, const edm::EventSetup&) {  
    std::unique_ptr<GenParticleVector> outputGenParticles(new GenParticleVector);
    edm::Handle<GenParticleView> inputParticleHandle;
    iEvent.getByToken(inputParticleToken_, inputParticleHandle);
    // Loop over particles, store mothers of pdgid 51-53 in a set
    std::set<const reco::Candidate *> motherSet;
    for (auto const & candidate : inputParticleHandle->ptrs()){
        int abspdgid = abs(candidate->pdgId());
        if (
            candidate->numberOfMothers() < 1
            || !( abspdgid == 51 || abspdgid == 52 || abspdgid == 53 )
            ){
            continue;
            }
        motherSet.insert(candidate->mother());
        }
    // Convert mothers to reco::GenParticle (loses some 
    // information but only care about kinematics anyway)
    // and store the mother in the output vector
    for(auto mother : motherSet) {
        reco::GenParticle motherparticle(
            mother->charge(), mother->p4(), mother->vertex(),
            mother->pdgId(), mother->status(), mother->charge()
            );
        outputGenParticles->push_back(motherparticle);
        }   
    iEvent.put(std::move(outputGenParticles));
    }

void DumpDMParticles::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("PartTag", edm::InputTag("genParticles"));
    desc.add<bool>("onlyFromZPrime", false);
    desc.add<bool>("verbose", false);
    descriptions.add("DumpDMParticles", desc);
    }

DEFINE_FWK_MODULE(DumpDMParticles);