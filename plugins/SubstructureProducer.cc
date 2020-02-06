// system include files
#include <memory>
#include <vector>
#include <cstdlib>

// user include files
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

#include "DataFormats/SVJFormats/interface/SubstructurePack.h"

//
// class declaration
//

class SubstructureProducer : public edm::stream::EDProducer<> {
    public:
        explicit SubstructureProducer(const edm::ParameterSet&);
        ~SubstructureProducer() {}

    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;

        // ----------member data ---------------------------
        // data labels
        float distMax_;
        edm::EDGetTokenT<edm::View<reco::GenJet>> jetToken_;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> particleToken_;
        std::vector<edm::InputTag>                                      algoTags_;
        std::vector< edm::EDGetTokenT< edm::View<reco::BasicJet> > >  algoTokens_;
    };

SubstructureProducer::SubstructureProducer(const edm::ParameterSet& iConfig) :
    distMax_( iConfig.getParameter<double>("distMax") ),
    jetToken_(consumes<edm::View<reco::GenJet>>( iConfig.getParameter<edm::InputTag>("jetSrc") )),
    particleToken_(consumes<edm::View<reco::GenParticle>>( iConfig.getParameter<edm::InputTag>("PartTag") )),
    algoTags_ (iConfig.getParameter<std::vector<edm::InputTag> > ( "algoTags" ))
    {
    algoTokens_ = edm::vector_transform(algoTags_, [this](edm::InputTag const & tag){return consumes< edm::View<reco::BasicJet> >(tag);});
    //register products
    produces<std::vector<SubstructurePack> > ();
    }


// ------------ method called to produce the data  ------------
void SubstructureProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {  
    // Output vector to be put in event
    std::unique_ptr< std::vector< SubstructurePack >> substructurePacks( new std::vector<SubstructurePack> );

    // Get the GenJets and the to-be-ran substructure algorithms using the handeles
    edm::Handle<edm::View<reco::GenJet>> jetHandle;

    std::vector< edm::Handle< edm::View<reco::BasicJet> > > algoHandles;
    iEvent.getByToken( jetToken_, jetHandle );
    algoHandles.resize( algoTags_.size() );
    for ( size_t i = 0; i < algoTags_.size(); ++i ) {
        iEvent.getByToken( algoTokens_[i], algoHandles[i] ); 
        }

    // for ( auto const & genJet : *jetHandle  ) {
    for (auto const & genJet : jetHandle->ptrs() ) {
        // Create a SubstructurePack, a simple collection of associated jets and subjets
        SubstructurePack substructurePack(genJet);
        // Loop over the substructure collections
        for ( auto const & ialgoHandle : algoHandles ) {
            for (auto const & jetFromSubstructureAlgo : ialgoHandle->ptrs() ) {
                // Get the jetFromSubstructureAlgo that matches the jet (by looking at a simple deltaR)
                if ( reco::deltaR( *substructurePack.jet(), *jetFromSubstructureAlgo ) < distMax_ ) {
                    // Add the whole substructure to the substructurePack
                    substructurePack.addSubstructure(jetFromSubstructureAlgo);
                    break;
                    }
                }
            }
        // Probably creates a copy, could be more efficient
        substructurePacks->push_back(substructurePack);
        }

    // Find the Z' genParticle of the event
    edm::Handle<edm::View<reco::GenParticle>> particleHandle;
    iEvent.getByToken(particleToken_, particleHandle);
    for (auto const & particle : particleHandle->ptrs() ) {
        if (particle->pdgId() == 4900023 and particle->isLastCopy() == 1) {
            // Match it with a SubstructurePack
            for(unsigned int i=0; i < substructurePacks->size(); i++){
                SubstructurePack * substructurePack = &(*substructurePacks)[i];
                if ( reco::deltaR( *(substructurePack->jet()), *particle ) < distMax_ ) {
                    substructurePack->addZprime(particle) ;
                    break ;
                    }
                }
            break ;
            }
        }

    // Put the output back in the event
    iEvent.put(std::move(substructurePacks));
    }

//define this as a plug-in
DEFINE_FWK_MODULE(SubstructureProducer);