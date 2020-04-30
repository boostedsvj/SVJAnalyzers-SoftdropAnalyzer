//framework headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

//analysis headers
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "fastjet/config.h"

//ROOT headers
#include <TTree.h>
#include <TLorentzVector.h>
 
//STL headers 
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>
using std::vector;

//user headers
#include "SVJAnalyzers/SoftdropAnalyzer/interface/lester_mt2_bisect.h"
#include "SVJAnalyzers/SoftdropAnalyzer/interface/NjettinessHelper.h"
#include "SVJAnalyzers/SoftdropAnalyzer/interface/ECFHelper.h"
#include "DataFormats/SVJFormats/interface/SubstructurePack.h"
#include "DataFormats/SVJFormats/interface/GenParticlePlus.h"


struct Jet {
    vector<TLorentzVector> p4_;

    void fill(TLorentzVector jet){ p4_.push_back(jet); }
    void fill(edm::Ptr<const reco::Jet> jet){
        fill( &(*jet) );
        }

    void fill(const reco::Jet * jet){
        p4_.push_back( TLorentzVector(jet->px(),jet->py(),jet->pz(),jet->energy()) );
        }

    void linkToTree(TTree* tree, std::string name){
        tree->Branch(name.c_str(), "vector<TLorentzVector>", &p4_, 32000, 0);
        }

    int size(){
        return p4_.size();
        }
    };

struct Subjets{
    /* Much like Jet, but stores an indices vector alongside of it so that subjets per ak15 jet can be retraced */
    vector<TLorentzVector> p4_;
    vector<unsigned int> groupNumber_;
    vector<unsigned int> groupCount_;
    unsigned int currentGroupNumber_ = 0;
    unsigned int currentGroupCount_ = 0;

    void newGroup(){
        groupCount_.push_back(currentGroupCount_);
        currentGroupNumber_++;
        currentGroupCount_ = 0;
        }

    void fill(edm::Ptr<const reco::Candidate> jet){
        fill(&(*jet));
        }

    void fill(const reco::Candidate * jet){
        p4_.push_back( TLorentzVector(jet->px(),jet->py(),jet->pz(),jet->energy()) );
        groupNumber_.push_back(currentGroupNumber_);
        currentGroupCount_++;
        }

    void linkToTree(TTree* tree, std::string name){
        tree->Branch(name.c_str(), "vector<TLorentzVector>", &p4_, 32000, 0);
        tree->Branch((name + "_index").c_str(), "vector<unsigned int>", &groupNumber_, 32000, 0);
        tree->Branch((name + "_count").c_str(), "vector<unsigned int>", &groupCount_, 32000, 0);
        }
    };

struct JetWithSubstructure : Jet {
    static NjettinessHelper * njhelper_;
    static ECFHelper * echelper_;
    static bool doNJ_;
    static bool doECF_;
    static void setNjettinessHelper(NjettinessHelper * fNjhelper){
        njhelper_ = fNjhelper;
        doNJ_ = true;
        }
    static void setECFHelper(ECFHelper * fEchelper){
        echelper_ = fEchelper;
        doECF_ = true;
        }

    vector<double> axisAverage_ ;
    vector<double> axisMajor_ ;
    vector<double> axisMinor_ ;
    vector<double> momentGirth_ ;
    vector<double> momentHalf_ ;
    vector<int>    multiplicity_ ;
    vector<double> ptD_ ;
    vector<double> tau1_ ;
    vector<double> tau2_ ;
    vector<double> tau3_ ;
    vector<double> ECF1_ ;
    vector<double> ECF2_ ;
    vector<double> ECF3_ ;

    void fill(const reco::GenJet * jet){
        Jet::fill(jet);
        double momentGirth = 0.0, momentHalf = 0.0;
        double sumPt = 0.0, sumPt2 = 0.0;
        double sumDeta = 0.0, sumDphi = 0.0, sumDeta2 = 0.0, sumDphi2 = 0.0, sumDetaDphi = 0.0;

        for ( auto const & particle : jet->daughterPtrVector()) {
            if (not(particle.isNonnull() and particle.isAvailable())){ continue; }
            double dphi = reco::deltaPhi(jet->phi(),particle->phi());
            double deta = particle->eta() - jet->eta();
            double dR = reco::deltaR(jet->p4(), particle->p4());
            double pT = particle->pt();
            double pT2 = pT*pT;
            
            momentGirth += pT*dR;
            momentHalf += pT*std::sqrt(dR);
            sumPt += pT;
            sumPt2 += pT2;
            sumDeta += deta*pT2;
            sumDphi += dphi*pT2;
            sumDeta2 += deta*deta*pT2;
            sumDphi2 += dphi*dphi*pT2;
            sumDetaDphi += deta*dphi*pT2;
            }

        //finish axis calculations (eigenvectors)
        sumDeta /= sumPt2;
        sumDphi /= sumPt2;
        sumDeta2 /= sumPt2;
        sumDphi2 /= sumPt2;
        sumDetaDphi /= sumPt2;
        double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
        a = sumDeta2 - sumDeta*sumDeta;
        b = sumDphi2 - sumDphi*sumDphi;
        c = sumDeta*sumDphi - sumDetaDphi;
        d = std::sqrt(std::fabs((a-b)*(a-b)+4*c*c));
        double axisMajor = (a+b+d)>0 ? std::sqrt(0.5*(a+b+d)) : 0.0 ;
        double axisMinor = (a+b-d)>0 ? std::sqrt(0.5*(a+b-d)) : 0.0 ;

        // Finish up class variables
        axisMajor_.push_back(axisMajor);
        axisMinor_.push_back(axisMinor);
        axisAverage_.push_back( std::sqrt(axisMajor*axisMajor + axisMinor*axisMinor) );
        momentGirth_.push_back( momentGirth / jet->pt() );
        momentHalf_.push_back( momentHalf / jet->pt() );
        multiplicity_.push_back( jet->numberOfDaughters() );
        ptD_.push_back( std::sqrt(sumPt2) / sumPt );

        if (doNJ_){
            tau1_.push_back( njhelper_->getTau(1, *jet) );
            tau2_.push_back( njhelper_->getTau(2, *jet) );
            tau3_.push_back( njhelper_->getTau(3, *jet) );
            }
        if (doECF_){
            auto ECFresult = echelper_->getECFs(*jet);
            ECF1_.push_back( (ECFresult[0]) );
            ECF2_.push_back( (ECFresult[1]) );
            ECF3_.push_back( (ECFresult[2]) );
            }
        }

    void linkToTree(TTree* tree, std::string name){
        Jet::linkToTree(tree, name);
        tree->Branch((name + std::string("_AxisAverage")).c_str(),  "vector<double>",         &axisAverage_, 32000, 0);
        tree->Branch((name + std::string("_AxisMajor")).c_str(),    "vector<double>",         &axisMajor_, 32000, 0);
        tree->Branch((name + std::string("_AxisMinor")).c_str(),    "vector<double>",         &axisMinor_, 32000, 0);
        tree->Branch((name + std::string("_ECF1")).c_str(),         "vector<double>",         &ECF1_, 32000, 0);
        tree->Branch((name + std::string("_ECF2")).c_str(),         "vector<double>",         &ECF2_, 32000, 0);
        tree->Branch((name + std::string("_ECF3")).c_str(),         "vector<double>",         &ECF3_, 32000, 0);
        tree->Branch((name + std::string("_MomentGirth")).c_str(),  "vector<double>",         &momentGirth_, 32000, 0);
        tree->Branch((name + std::string("_MomentHalf")).c_str(),   "vector<double>",         &momentHalf_, 32000, 0);
        tree->Branch((name + std::string("_Multiplicity")).c_str(), "vector<int>",            &multiplicity_, 32000, 0);
        tree->Branch((name + std::string("_PtD")).c_str(),          "vector<double>",         &ptD_, 32000, 0);
        tree->Branch((name + std::string("_Tau1")).c_str(),         "vector<double>",         &tau1_, 32000, 0);
        tree->Branch((name + std::string("_Tau2")).c_str(),         "vector<double>",         &tau2_, 32000, 0);
        tree->Branch((name + std::string("_Tau3")).c_str(),         "vector<double>",         &tau3_, 32000, 0);
        }
    };

bool JetWithSubstructure::doNJ_ = false;
bool JetWithSubstructure::doECF_ = false;
NjettinessHelper * JetWithSubstructure::njhelper_ = NULL;
ECFHelper * JetWithSubstructure::echelper_ = NULL;


struct JetWithSubstructure_MTwMET : JetWithSubstructure {
    vector<double> transverseMass_;

    void fill(const reco::GenJet * jet, TLorentzVector & met){
        JetWithSubstructure::fill(jet);
        TLorentzVector jetp4 = p4_.back();
        double Ejet = std::sqrt(
            std::pow(jetp4.Px(),2) + std::pow(jetp4.Py(),2) + std::pow(jetp4.M(),2)
            );
        double Emet = std::sqrt(
            std::pow(met.Px(),2) + std::pow(met.Py(),2) // Assume 0 mass for met
            );
        double MTsq =
              std::pow( Ejet + Emet , 2 )
            - std::pow( jetp4.Px() + met.Px(), 2 )
            - std::pow( jetp4.Py() + met.Py(), 2 )
            ;
        transverseMass_.push_back( std::sqrt( std::max(MTsq, 0.0) ) );        
        }

    void linkToTree(TTree* tree, std::string name){
        JetWithSubstructure::linkToTree(tree, name);
        tree->Branch((name + std::string("_MTwMET")).c_str(), &transverseMass_);
        }
    };


struct GenParticle {
    vector<TLorentzVector> p4_;
    vector<int> pdgId_;

    // template <class GenParticlePtr = const reco::GenParticle *>
    void fill(const reco::GenParticle * particle){
        p4_.push_back(
            TLorentzVector(particle->px(), particle->py(), particle->pz(), particle->energy())
            );
        pdgId_.push_back(particle->pdgId());
        }

    void linkToTree(TTree* tree, std::string prefix){
        tree->Branch((prefix + std::string("_p4")).c_str(),    "vector<TLorentzVector>", &p4_, 32000, 0);
        tree->Branch((prefix + std::string("_pdgId")).c_str(), "vector<int>",            &pdgId_, 32000, 0);
        }
    };


// ____________________________________________________________
// Event format

struct Entry {
    /* Main entry of the flat ntuple */
    int event_number;

    int n_all_ak15jet;
    // shape: n_all_ak15jet
    JetWithSubstructure_MTwMET all_ak15jet;
    Jet all_softdropjet;
    Jet all_summedsubjets;
    Subjets all_subjets;
    vector<bool> is_matched;

    // shape: n_matched_ak15jet
    JetWithSubstructure_MTwMET matched_ak15jet;
    Jet matched_softdropjet;
    Jet matched_summedsubjets;
    Subjets matched_subjets;
    GenParticle zprime;
    GenParticle initial_dark_quarks;
    GenParticle final_dark_quarks;
    GenParticle final_visible_particles;
    GenParticle dark_mesons_decaying_dark;
    GenParticle dark_mesons_decaying_visible;

    // other
    GenParticle genparticles;
    TLorentzVector met;
    double ht;

    void linkToTree(TTree* tree){
        tree->Branch("event_number", &event_number);

        tree->Branch("n_all_ak15jet", &n_all_ak15jet);
        all_ak15jet.linkToTree(tree, "all_ak15jet");
        all_softdropjet.linkToTree(tree, "all_softdropjet");
        all_subjets.linkToTree(tree, "all_subjets");
        all_summedsubjets.linkToTree(tree, "all_summedsubjets");
        tree->Branch("is_matched", "vector<bool>", &is_matched, 32000, 0);

        matched_ak15jet.linkToTree(tree, "matched_ak15jet");
        matched_softdropjet.linkToTree(tree, "matched_softdropjet");
        matched_summedsubjets.linkToTree(tree, "matched_summedsubjets");
        matched_subjets.linkToTree(tree, "matched_subjets");
        zprime.linkToTree(tree, "zprime");
        initial_dark_quarks.linkToTree(tree, "initial_dark_quarks");
        final_dark_quarks.linkToTree(tree, "final_dark_quarks");
        final_visible_particles.linkToTree(tree, "final_visible_particles");
        dark_mesons_decaying_dark.linkToTree(tree, "dark_mesons_decaying_dark");
        dark_mesons_decaying_visible.linkToTree(tree, "dark_mesons_decaying_visible");

        genparticles.linkToTree(tree, "genparticles");
        tree->Branch("met", &met);
        tree->Branch("ht", &ht);
        }
    };


// ____________________________________________________________
// The actual analyzer

class SoftdropAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit SoftdropAnalyzer(const edm::ParameterSet&);
        ~SoftdropAnalyzer() {}
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        void beginJob() override;
        void doBeginRun_(const edm::Run&, const edm::EventSetup&) override {}
        void analyze(const edm::Event&, const edm::EventSetup&) override;
        void doEndRun_(const edm::Run&, const edm::EventSetup&) override {}
        void endJob() override {}

        edm::Service<TFileService> fs;
        TTree* tree_;
        Entry entry;

        edm::EDGetTokenT<vector<reco::GenMET>> tok_met;
        edm::EDGetTokenT<vector<reco::GenJet>> tok_jet;
        edm::EDGetTokenT<vector<reco::GenParticle>> tok_part;
        edm::EDGetTokenT<vector<SubstructurePack>> tok_substructurepacks;
        edm::EDGetTokenT<double> tok_ht;

        NjettinessHelper njhelper_;
        ECFHelper echelper_;
    };

SoftdropAnalyzer::SoftdropAnalyzer(const edm::ParameterSet& iConfig) : 
    tree_(NULL),
    tok_met(consumes<vector<reco::GenMET>>(iConfig.getParameter<edm::InputTag>("METTag"))),
    tok_jet(consumes<vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
    tok_part(consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("PartTag"))),
    tok_substructurepacks(consumes<vector<SubstructurePack>>(iConfig.getParameter<edm::InputTag>("SubstructurePackTag"))),
    tok_ht(consumes<double>(iConfig.getParameter<edm::InputTag>("HTTag"))),
    njhelper_(iConfig.getParameter<edm::ParameterSet>("Nsubjettiness")),
    echelper_(iConfig.getParameter<edm::ParameterSet>("ECF"))
    {
        JetWithSubstructure::setNjettinessHelper(&njhelper_);
        JetWithSubstructure::setECFHelper(&echelper_);
        usesResource("TFileService");
        }

void SoftdropAnalyzer::beginJob() {
    asymm_mt2_lester_bisect::disableCopyrightMessage();
    tree_ = fs->make<TTree>("tree","tree");
    entry.linkToTree(tree_);
    }

void SoftdropAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    entry = Entry(); // Clear the entry

    edm::Handle<vector<reco::GenMET>> h_met;
    iEvent.getByToken(tok_met,h_met);
    edm::Handle<vector<reco::GenJet>> h_jet;
    iEvent.getByToken(tok_jet,h_jet);
    edm::Handle<vector<reco::GenParticle>> h_part;
    iEvent.getByToken(tok_part,h_part);
    edm::Handle<vector<SubstructurePack>> h_subsstructurepacks;
    iEvent.getByToken(tok_substructurepacks, h_subsstructurepacks);

    entry.event_number = iEvent.id().event();

    // Get HT
    edm::Handle<double> h_ht;
    iEvent.getByToken(tok_ht, h_ht);
    entry.ht = *h_ht;

    // Get the MET
    const auto& met = h_met->front();
    entry.met = TLorentzVector(met.px(), met.py(), met.pz(), met.energy());

    // Put all jets in the tree (there is exactly one substructurePack per genJet)
    for(const auto& substructurePack : *(h_subsstructurepacks.product())){
        entry.all_ak15jet.fill(substructurePack.jet(), entry.met);
        entry.all_softdropjet.fill(substructurePack.substructurejet());
        entry.all_summedsubjets.fill(substructurePack.summedsubjets());
        for ( auto const & subjet : substructurePack.subjets()) {
            entry.all_subjets.fill(subjet);
            }
        entry.all_subjets.newGroup();
        entry.is_matched.push_back(substructurePack.hasZprime());
        // Fill matched jets again in a dedicated variable
        if (substructurePack.hasZprime()) {
            entry.matched_ak15jet.fill(substructurePack.jet(), entry.met);
            entry.matched_softdropjet.fill(substructurePack.substructurejet());
            entry.matched_summedsubjets.fill(substructurePack.summedsubjets());
            for ( auto const & subjet : substructurePack.subjets()) {
                entry.matched_subjets.fill(subjet);
                }
            entry.matched_subjets.newGroup();
            // Get the zprime, fill it, look for all its decay products, and fill those too
            const reco::GenParticle * zprime = substructurePack.zprime();
            entry.zprime.fill(zprime);
            GenParticlePlus zprimeWithDecayProducts(*zprime);
            zprimeWithDecayProducts.setQuarksAndFinalProducts();
            for( auto& particle : zprimeWithDecayProducts.initialDarkQuark){
                entry.initial_dark_quarks.fill(&particle);
                }
            for( auto& particle : zprimeWithDecayProducts.finalDarkQuark){
                entry.final_dark_quarks.fill(&particle);
                }
            for( auto& particle : zprimeWithDecayProducts.finalVisibleProduct){
                entry.final_visible_particles.fill(&particle);
                }
            for( auto& particle : zprimeWithDecayProducts.darkMesonDecayingDark){
                entry.dark_mesons_decaying_dark.fill(&particle);
                }
            for( auto& particle : zprimeWithDecayProducts.darkMesonDecayingVisible){
                entry.dark_mesons_decaying_visible.fill(&particle);
                }
            }
        entry.n_all_ak15jet = entry.all_ak15jet.size();
        }
    // Save all genParticles too
    for(const auto& particle : *(h_part.product())){
        entry.genparticles.fill(&particle);
        }
    tree_->Fill();
    }

void SoftdropAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("METTag",edm::InputTag("genMetTrue"));
    desc.add<edm::InputTag>("HTTag", edm::InputTag("htProducer", "genHT"));
    desc.add<edm::InputTag>("JetTag",edm::InputTag("ak8GenJetsNoNu"));
    desc.add<edm::InputTag>("PartTag",edm::InputTag("genParticles"));
    desc.add<edm::InputTag>("SubstructurePackTag", edm::InputTag("packedGenJetsAK8NoNu")); 
    
    edm::ParameterSetDescription desc_nj;
    desc_nj.add<unsigned>("measureDefinition",0);
    desc_nj.add<double>("beta",1.0);
    desc_nj.add<double>("R0",0.8);
    desc_nj.add<double>("Rcutoff",999.0);
    desc_nj.add<unsigned>("axesDefinition",6);
    desc_nj.add<int>("nPass",999);
    desc_nj.add<double>("akAxesR0",999.0);
    desc.add<edm::ParameterSetDescription>("Nsubjettiness", desc_nj);

    edm::ParameterSetDescription desc_ec;
    desc_ec.add<vector<unsigned>>("Njets",{1,2,3});
    desc_ec.add<double>("beta",1.0);
    desc.add<edm::ParameterSetDescription>("ECF", desc_ec);
    
    descriptions.add("SoftdropAnalyzer",desc);
    }

DEFINE_FWK_MODULE(SoftdropAnalyzer);