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


// __________________________________________________
// Class definition

class SingleJetProperties {
    public:
        SingleJetProperties() {}
        ~SingleJetProperties() {}

        static NjettinessHelper * njhelper_;
        static ECFHelper * echelper_;
        static bool doNJ_;
        static bool doECF_;

        // static void setHelpers(const edm::ParameterSet& iConfig){
        //     doNJ_ = true;
        //     doECF_ = true;
        //     njhelper_ = new NjettinessHelper(iConfig.getParameter<edm::ParameterSet>("Nsubjettiness"));
        //     echelper_ = new ECFHelper(iConfig.getParameter<edm::ParameterSet>("ECF"));
        //     }

        static void setNjettinessHelper(NjettinessHelper * fNjhelper){
            njhelper_ = fNjhelper;
            doNJ_ = true;
            }

        static void setECFHelper(ECFHelper * fEchelper){
            echelper_ = fEchelper;
            doECF_ = true;
            }

        vector<TLorentzVector> p4_;
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

        // Take const reco::GenJet * as the default, but needs to be able to take edm::Ptr<reco::GenJet> too
        template <class GenJetPtr = const reco::GenJet *>
        void read(GenJetPtr jet){
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
            p4_.push_back( TLorentzVector(jet->px(),jet->py(),jet->pz(),jet->energy()) );
            momentGirth_.push_back( momentGirth / jet->pt() );
            momentHalf_.push_back( momentHalf / jet->pt() );
            multiplicity_.push_back( jet->numberOfDaughters() );
            ptD_.push_back( sumPt2 / sumPt );

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

        void setTreeAdresses(TTree* tree, std::string prefix){
            tree->Branch((prefix + std::string("p4")).c_str(),           "vector<TLorentzVector>", &p4_, 32000, 0);
            tree->Branch((prefix + std::string("AxisAverage")).c_str(),  "vector<double>",         &axisAverage_, 32000, 0);
            tree->Branch((prefix + std::string("AxisMajor")).c_str(),    "vector<double>",         &axisMajor_, 32000, 0);
            tree->Branch((prefix + std::string("AxisMinor")).c_str(),    "vector<double>",         &axisMinor_, 32000, 0);
            tree->Branch((prefix + std::string("ECF1")).c_str(),         "vector<double>",         &ECF1_, 32000, 0);
            tree->Branch((prefix + std::string("ECF2")).c_str(),         "vector<double>",         &ECF2_, 32000, 0);
            tree->Branch((prefix + std::string("ECF3")).c_str(),         "vector<double>",         &ECF3_, 32000, 0);
            tree->Branch((prefix + std::string("MomentGirth")).c_str(),  "vector<double>",         &momentGirth_, 32000, 0);
            tree->Branch((prefix + std::string("MomentHalf")).c_str(),   "vector<double>",         &momentHalf_, 32000, 0);
            tree->Branch((prefix + std::string("Multiplicity")).c_str(), "vector<int>",            &multiplicity_, 32000, 0);
            tree->Branch((prefix + std::string("PtD")).c_str(),          "vector<double>",         &ptD_, 32000, 0);
            tree->Branch((prefix + std::string("Tau1")).c_str(),         "vector<double>",         &tau1_, 32000, 0);
            tree->Branch((prefix + std::string("Tau2")).c_str(),         "vector<double>",         &tau2_, 32000, 0);
            tree->Branch((prefix + std::string("Tau3")).c_str(),         "vector<double>",         &tau3_, 32000, 0);
            }

        int size(){
            return p4_.size();
            }
    };

bool SingleJetProperties::doNJ_ = false;
bool SingleJetProperties::doECF_ = false;
NjettinessHelper * SingleJetProperties::njhelper_ = NULL;
ECFHelper * SingleJetProperties::echelper_ = NULL;


class SingleParticleProperties {
    public:
        SingleParticleProperties() {}
        ~SingleParticleProperties() {}

        vector<TLorentzVector> p4_;
        vector<int> pdgId_;

        void read(const reco::GenParticle * particle){
            p4_.push_back(
                TLorentzVector(
                    particle->px(),
                    particle->py(),
                    particle->pz(),
                    particle->energy()
                    )
                );
            pdgId_.push_back(particle->pdgId());
            }

        void setTreeAdresses(TTree* tree, std::string prefix){
            tree->Branch((prefix + std::string("p4")).c_str(),    "vector<TLorentzVector>", &p4_, 32000, 0);
            tree->Branch((prefix + std::string("pdgId")).c_str(), "vector<int>",            &pdgId_, 32000, 0);
            }
    };


class SVJSystemProperties {
    public:
        SVJSystemProperties() {}
        ~SVJSystemProperties() {}

        double deltaPhi1_;
        double deltaPhi2_;
        double deltaPhiMin_;

        TLorentzVector dijetVector_;
        vector<TLorentzVector> METSystems_;

        double Mmc_;
        double MT_;
        double MT2_;
        double MAOS_;

        void read(
            TLorentzVector * genJet1,
            TLorentzVector * genJet2,
            TLorentzVector * met,
            TLorentzVector * hvmesonSum
            ){

            //delta phis
            deltaPhi1_ = std::abs(reco::deltaPhi( genJet1->Phi(), met->Phi() ));
            deltaPhi2_ = std::abs(reco::deltaPhi( genJet2->Phi(), met->Phi() ));
            deltaPhiMin_ = std::min(deltaPhi1_, deltaPhi2_);
            
            dijetVector_ = *genJet1 + *genJet2 ;

            // include all jets in MC mass
            Mmc_ = (dijetVector_ + *hvmesonSum).M();

            //assume MET is massless
            MT_ = asymm_mt2_lester_bisect::MT(
                dijetVector_.Px(), met->Px(), dijetVector_.Py(), met->Py(), dijetVector_.M(), 0.0
                );
            MT2_ = asymm_mt2_lester_bisect::get_mT2(
                genJet1->M(), genJet1->Px(), genJet1->Py(),
                genJet2->M(), genJet2->Px(), genJet2->Py(),
                met->Px(), met->Py(), 0.0,
                0.0, 0
                );

            //get invisible systems from MT2 (double check argument ordering)
            auto MET0 = asymm_mt2_lester_bisect::ben_findsols(
                MT2_, 
                genJet1->Px(), genJet1->Py(), genJet1->M(), 0.0,
                genJet2->Px(), genJet2->Py(),
                met->Px(), met->Py(), genJet2->M(), 0.0
                );
            double MET0x = MET0.first;
            double MET0y = MET0.second;
            double MET0t = std::sqrt( std::pow(MET0x,2) + std::pow(MET0y,2) );
            double MET1x = met->Px() - MET0x;
            double MET1y = met->Py() - MET0y;
            double MET1t = std::sqrt( std::pow(MET1x,2) + std::pow(MET1y,2) );
            
            //use MAOS scheme 2 ("modified") to estimate longitudinal momenta of invisible systems
            double MET0z = MET0t * genJet1->Pz() / genJet1->Pt();
            double MET1z = MET1t * genJet2->Pz() / genJet2->Pt();

            TLorentzVector vectorMET0( MET0x, MET0y, MET0z, 0.0 );
            TLorentzVector vectorMET1( MET1x, MET1y, MET1z, 0.0 );
            METSystems_.push_back(vectorMET0);
            METSystems_.push_back(vectorMET1);
            //construct invariant mass of parent
            MAOS_ = (dijetVector_ + vectorMET0 + vectorMET1).M();
            }

        void setTreeAdresses(TTree* tree, std::string prefix){
            tree->Branch((prefix + std::string("deltaPhi1")).c_str(), &deltaPhi1_);
            tree->Branch((prefix + std::string("deltaPhi2")).c_str(), &deltaPhi2_);
            tree->Branch((prefix + std::string("deltaPhiMin")).c_str(), &deltaPhiMin_);

            tree->Branch((prefix + std::string("dijetVector")).c_str(), "TLorentzVector", &dijetVector_, 32000, 0);
            tree->Branch((prefix + std::string("METSystems")).c_str(), "vector<TLorentzVector>", &METSystems_, 32000, 0);

            tree->Branch((prefix + std::string("Mmc")).c_str(), &Mmc_);
            tree->Branch((prefix + std::string("MT")).c_str(), &MT_);
            tree->Branch((prefix + std::string("MT2")).c_str(), &MT2_);
            tree->Branch((prefix + std::string("MAOS")).c_str(), &MAOS_);
            }
    };


class SubstructurePackProperties {
    public:
        SubstructurePackProperties() {}
        ~SubstructurePackProperties() {}

        SingleJetProperties genjet_;
        vector<bool> hasZprime_;
        SingleParticleProperties zprime_;
        SingleJetProperties subjets_;
        vector<int> nSubjets_;
        vector<double> summedSDsubjetmass_;

        void read(const SubstructurePack * substructurePack){
            // Put the GenJetAK8 in the tree
            genjet_.read<const reco::GenJet *>(substructurePack->jet());

            // Note in tree whether the GenJetAK8 was matched to GenParticle Zprime,
            // and if so save the Zprime
            bool hasZprime = substructurePack->hasZprime();
            if (hasZprime){
                zprime_.read(substructurePack->zprime());
                }
            hasZprime_.push_back(hasZprime);

            // Put all the subjets in the tree, and also their summed mass
            TLorentzVector summedSDsubjets;
            for ( auto const subjet : substructurePack->subjets()) {
                subjets_.read<edm::Ptr<reco::GenJet>>(subjet);
                summedSDsubjets += TLorentzVector(subjet->px(),subjet->py(),subjet->pz(),subjet->energy());
                }
            summedSDsubjetmass_.push_back(summedSDsubjets.M());
            nSubjets_.push_back(substructurePack->nSubjets());
            }

        void setTreeAdresses(TTree* tree, std::string prefix){
            genjet_.setTreeAdresses(tree, prefix);
            tree->Branch((prefix + std::string("nSDsubjets")).c_str(), "vector<int>", &nSubjets_, 32000, 0);
            subjets_.setTreeAdresses(tree, prefix + std::string("SDsubjets_"));
            tree->Branch((prefix + std::string("summedSDsubjetmass")).c_str(), "vector<double>", &summedSDsubjetmass_, 32000, 0);
            tree->Branch((prefix + std::string("hasZprime")).c_str(), "vector<bool>", &hasZprime_, 32000, 0);
            zprime_.setTreeAdresses(tree, prefix + std::string("zprime_"));
            }
    };



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

        float distMax_;
        double MET_ = 0.;
        double METphi_ = 0.;
        int nGenJetsAK8_ = 0;
        int eventNum_;

        SubstructurePackProperties subpacks_;
        SingleParticleProperties hvmesons_;
        SingleParticleProperties allDarkParticles_;
        SVJSystemProperties metmasscalculations_;

        edm::Service<TFileService> fs;
        TTree* tree_;

        edm::EDGetTokenT<vector<reco::GenMET>> tok_met;
        edm::EDGetTokenT<vector<reco::GenJet>> tok_jet;
        edm::EDGetTokenT<vector<reco::GenParticle>> tok_part;
        edm::EDGetTokenT<vector<SubstructurePack>> tok_substructurepacks;

        NjettinessHelper njhelper_;
        ECFHelper echelper_;
    };


// __________________________________________________
// constructors and destructor

SoftdropAnalyzer::SoftdropAnalyzer(const edm::ParameterSet& iConfig) : 
    distMax_( iConfig.getParameter<double>("distMax") ),
    tree_(NULL),
    tok_met(consumes<vector<reco::GenMET>>(iConfig.getParameter<edm::InputTag>("METTag"))),
    tok_jet(consumes<vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
    tok_part(consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("PartTag"))),
    tok_substructurepacks(consumes<vector<SubstructurePack>>(iConfig.getParameter<edm::InputTag>("SubstructurePackTag"))),
    njhelper_(iConfig.getParameter<edm::ParameterSet>("Nsubjettiness")),
    echelper_(iConfig.getParameter<edm::ParameterSet>("ECF"))
    {
        // SingleJetProperties::setHelpers(iConfig);
        SingleJetProperties::setNjettinessHelper(&njhelper_);
        SingleJetProperties::setECFHelper(&echelper_);
        // njhelper_ = new NjettinessHelper(iConfig.getParameter<edm::ParameterSet>("Nsubjettiness"));
        // echelper_ = new ECFHelper(iConfig.getParameter<edm::ParameterSet>("ECF"));
        usesResource("TFileService");
        // std::cout << "SoftdropAnalyzer: Warning - Nsubjettiness and ECF variables not available!" << std::endl;
        }


// __________________________________________________
// member functions

void SoftdropAnalyzer::beginJob() {
    asymm_mt2_lester_bisect::disableCopyrightMessage();
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch( "eventNum", &eventNum_ );
    tree_->Branch( "nGenJetsAK8", &nGenJetsAK8_ );
    subpacks_.setTreeAdresses(tree_, "GenJetsAK8_");
    hvmesons_.setTreeAdresses(tree_, "HVMeson_");
    allDarkParticles_.setTreeAdresses(tree_, "allDarkParticles_");
    metmasscalculations_.setTreeAdresses(tree_, "");
    tree_->Branch( "MET", &MET_ );
    tree_->Branch( "METphi", &METphi_ );
    }


void SoftdropAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    eventNum_ = iEvent.id().event();
    edm::Handle<vector<reco::GenMET>> h_met;
    iEvent.getByToken(tok_met,h_met);
    edm::Handle<vector<reco::GenJet>> h_jet;
    iEvent.getByToken(tok_jet,h_jet);
    edm::Handle<vector<reco::GenParticle>> h_part;
    iEvent.getByToken(tok_part,h_part);
    edm::Handle<vector<SubstructurePack>> h_subsstructurepacks;
    iEvent.getByToken(tok_substructurepacks, h_subsstructurepacks);

    subpacks_ = SubstructurePackProperties();
    hvmesons_ = SingleParticleProperties();
    allDarkParticles_ = SingleParticleProperties();
    metmasscalculations_ = SVJSystemProperties();

    // Put all jets in the tree (there is exactly one substructurePack per genJet)
    for(const auto& substructurePack : *(h_subsstructurepacks.product())){
        subpacks_.read(&substructurePack);
        nGenJetsAK8_ = subpacks_.genjet_.size();
        }

    // Get all hvmesons, and calculate their summed 4-vector
    // Also plainly save all dark particles in the tree
    TLorentzVector hvmesonSum;
    for(const auto& particle : *(h_part.product())){
        if( std::abs(particle.pdgId()) == 4900211 ){
            hvmesons_.read(&particle);
            hvmesonSum += TLorentzVector(particle.px(), particle.py(), particle.pz(), particle.energy());
            }
        if(std::abs(particle.pdgId()) - 4900000 > 0){
            allDarkParticles_.read(&particle);
            }
        }

    // Get the MET
    const auto& met = h_met->front();
    MET_ = met.pt();
    METphi_ = met.phi();

    // Calculate variables that are only available if there are at least two genJets
    if( subpacks_.genjet_.size() >= 2 ){
        TLorentzVector met_p4(met.px(), met.py(), met.pz(), met.energy());
        metmasscalculations_.read(
            &(subpacks_.genjet_.p4_[0]),
            &(subpacks_.genjet_.p4_[1]),
            &met_p4,
            &hvmesonSum
            );
        }

    tree_->Fill();
    }


// __________________________________________________
// method fills 'descriptions' with the allowed parameters for the module
void SoftdropAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("METTag",edm::InputTag("genMetTrue"));
    desc.add<edm::InputTag>("JetTag",edm::InputTag("ak8GenJetsNoNu"));
    desc.add<edm::InputTag>("PartTag",edm::InputTag("genParticles"));
    desc.add<edm::InputTag>("SubstructurePackTag", edm::InputTag("packedGenJetsAK8NoNu")); 
    desc.add<double>("distMax",0.8);
    
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

//define this as a plug-in
DEFINE_FWK_MODULE(SoftdropAnalyzer);