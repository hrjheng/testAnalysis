// -*- C++ -*-
//
// Package:    testAnalyzer/testAnalyzer
// Class:      testAnalyzer
//
/**\class testAnalyzer testAnalyzer.cc testAnalyzer/testAnalyzer/plugins/testAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hao-Ren Jheng
//         Created:  Tue, 23 Jun 2020 01:08:25 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <iostream>
#include "testAnalyzer/testAnalyzer/interface/GenParticleParentage.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

void setbit(UShort_t &x, UShort_t bit)
{
    UShort_t a = 1;
    x |= (a << bit);
}

class testAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit testAnalyzer(const edm::ParameterSet &);
    ~testAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void refresh();

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    TTree *tree_out;

    TFile *rootFile_;
    std::string outputFile_;

    // + The variables
    //----------------
    int nGen;
    vector<int> genPID;
    vector<int> genStatus;
    vector<UShort_t> genStatusFlag;
    vector<float> genVtx;
    vector<float> genVty;
    vector<float> genVtz;
    vector<float> genPt;
    vector<float> genMass;
    vector<float> genEta;
    vector<float> genPhi;
    vector<float> genE;
    vector<float> genEt;
    vector<int> genMomPID;
    vector<int> genGMomPID;
    vector<float> genMomPt;
    vector<float> genMomMass;
    vector<float> genMomEta;
    vector<float> genMomPhi;
    vector<int> genParentage;
    // vector<float> genCalIsoDR04;
    int nGenJet;
    vector<float> genJetPt;
    vector<float> genJetEta;
    vector<float> genJetPhi;
    vector<float> genJetE;
    float genWeight;

    edm::EDGetTokenT<vector<reco::GenParticle>> genParticlesCollection_;
    edm::EDGetTokenT<vector<reco::GenJet>> genJetHandle_;
    edm::EDGetTokenT<GenEventInfoProduct> generatorLabel_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
testAnalyzer::testAnalyzer(const edm::ParameterSet &iConfig) : outputFile_(iConfig.getParameter<std::string>("outputFile"))
{
    genParticlesCollection_ = consumes<vector<reco::GenParticle>>(iConfig.getParameter<InputTag>("genParticleSrc"));
    genJetHandle_ = consumes<vector<reco::GenJet>>(iConfig.getParameter<InputTag>("ak4GenJetSrc"));
    generatorLabel_ = consumes<GenEventInfoProduct>(iConfig.getParameter<InputTag>("generatorLabel"));

    //now do what ever initialization is needed
    usesResource("TFileService");
    rootFile_ = TFile::Open(outputFile_.c_str(), "RECREATE"); // open output file to store histograms
}

testAnalyzer::~testAnalyzer()
{
    delete rootFile_;
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void testAnalyzer::beginJob()
{
    rootFile_->cd();
    tree_out = new TTree("EventTree", "GEN");

    // + Create the branches
    //----------------------
    // * Gen level
    tree_out->Branch("nGen", &nGen);
    tree_out->Branch("genPID", &genPID);
    tree_out->Branch("genStatus", &genStatus);
    tree_out->Branch("genStatusFlag", &genStatusFlag);
    tree_out->Branch("genPt", &genPt);
    tree_out->Branch("genEta", &genEta);
    tree_out->Branch("genPhi", &genPhi);
    tree_out->Branch("genE", &genE);
    tree_out->Branch("genEt", &genEt);
    tree_out->Branch("genMomPID", &genMomPID);
    tree_out->Branch("genGMomPID", &genGMomPID);
    tree_out->Branch("genMomPt", &genMomPt);
    tree_out->Branch("genMomMass", &genMomMass);
    tree_out->Branch("genMomEta", &genMomEta);
    tree_out->Branch("genMomPhi", &genMomPhi);
    tree_out->Branch("genParentage", &genParentage);
    // tree_out->Branch("genCalIsoDR04", &genCalIsoDR04);
    // * Generator level jets
    tree_out->Branch("nGenJet", &nGenJet);
    tree_out->Branch("genJetPt", &genJetPt);
    tree_out->Branch("genJetEta", &genJetEta);
    tree_out->Branch("genJetPhi", &genJetPhi);
    tree_out->Branch("genJetE", &genJetE);
    tree_out->Branch("genWeight", &genWeight);
}

// ------------ method called for each event  ------------
void testAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    using namespace edm;
    using namespace std;

    refresh();

    edm::Handle<vector<reco::GenParticle>> genParticlesHandle;
    iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

    for (vector<reco::GenParticle>::const_iterator igen = genParticlesHandle->begin(); igen != genParticlesHandle->end(); igen++)
    {
        // * Get particle status
        // Reference: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
        int status = igen->status();

        bool quarks = abs(igen->pdgId()) < 7;

        // Reference: https://github.com/cmkuo/ggAnalysis/blob/master/ggNtuplizer/plugins/ggNtuplizer_genParticles.cc
        // keep non-FSR photons with pT > 5.0 and all leptons with pT > 3.0; => modify the pt cut to 0 GeV
        // pdgId = 22 : photon
        // pdgId = 11 : electron
        // pdgId = 13 : muon
        // pdgId = 12/14/16 : electron/muon/tau neutrino
        // pdgId = 15 tau
        bool photonOrLepton =
            (igen->pdgId() == 22 && (igen->isPromptFinalState() || igen->isLastCopy())) ||
            (status == 1 && abs(igen->pdgId()) == 11 && (igen->isPromptFinalState() || igen->isLastCopy())) ||
            (status == 1 && abs(igen->pdgId()) == 13 && (igen->isPromptFinalState() || igen->isLastCopy())) ||
            (status == 1 && (abs(igen->pdgId()) == 12 || abs(igen->pdgId()) == 14 || abs(igen->pdgId()) == 16)) ||
            (status == 1 && (abs(igen->pdgId()) >= 11 && abs(igen->pdgId()) <= 16) && igen->pt() > 0.0) ||
            (status < 10 && abs(igen->pdgId()) == 15 && igen->pt() > 0.0);

        // select also Z, W, H, top and b
        bool heavyParticle =
            ((igen->pdgId() == 23 && igen->isHardProcess()) ||
             (abs(igen->pdgId()) == 24 && igen->isHardProcess()) ||
             (igen->pdgId() == 25 && igen->isHardProcess()) ||
             (abs(igen->pdgId()) == 6 && igen->isHardProcess()) ||
             (abs(igen->pdgId()) == 5 && igen->isHardProcess()));

        if (heavyParticle || photonOrLepton || quarks)
        {
            // + Generator particles (non-genJet)
            //--------------------------------
            const reco::Candidate *myParticle = (const reco::Candidate *)&(*igen);
            if (!myParticle->mother())
                continue;

            genPID.push_back(myParticle->pdgId());
            genVtx.push_back(myParticle->vx());
            genVty.push_back(myParticle->vy());
            genVtz.push_back(myParticle->vz());
            genPt.push_back(myParticle->pt());
            genMass.push_back(myParticle->mass());
            genEta.push_back(myParticle->eta());
            genPhi.push_back(myParticle->phi());
            genE.push_back(myParticle->energy());
            genEt.push_back(myParticle->et());
            genStatus.push_back(myParticle->status());

            UShort_t tmpStatusFlag = 0;
            if (igen->fromHardProcessFinalState())
                setbit(tmpStatusFlag, 0);
            if (igen->isPromptFinalState())
                setbit(tmpStatusFlag, 1);
            if (igen->isHardProcess())
                setbit(tmpStatusFlag, 2);

            // if genParticle is W or Z, check its decay type
            if (igen->pdgId() == 23 || abs(igen->pdgId()) == 24)
            {
                for (size_t k = 0; k < myParticle->numberOfDaughters(); ++k)
                {
                    const reco::Candidate *dp = myParticle->daughter(k);
                    if (abs(dp->pdgId()) <= 6)
                        setbit(tmpStatusFlag, 4);
                    else if (abs(dp->pdgId()) == 11 || abs(dp->pdgId()) == 12)
                        setbit(tmpStatusFlag, 5);
                    else if (abs(dp->pdgId()) == 13 || abs(dp->pdgId()) == 14)
                        setbit(tmpStatusFlag, 6);
                    else if (abs(dp->pdgId()) == 15 || abs(dp->pdgId()) == 16)
                        setbit(tmpStatusFlag, 7);
                    else
                        continue;
                }
            }
            if (igen->isLastCopy())
                setbit(tmpStatusFlag, 8);

            if (abs(myParticle->pdgId()) == 13 && igen->isPromptFinalState())
                cout << tmpStatusFlag << " " << ((tmpStatusFlag >> 0) & 1) << " " << ((tmpStatusFlag >> 1) & 1) << endl;
                
            genStatusFlag.push_back(tmpStatusFlag);

            int mcGMomPID_ = -999;
            int mcMomPID_ = -999;
            float mcMomPt_ = -999.;
            float mcMomMass_ = -999.;
            float mcMomEta_ = -999.;
            float mcMomPhi_ = -999.;

            reco::GenParticleRef partRef = reco::GenParticleRef(genParticlesHandle,
                                                                igen - genParticlesHandle->begin());
            genpartparentage::GenParticleParentage particleHistory(partRef);

            genParentage.push_back(particleHistory.hasLeptonParent() * 16 +
                                   particleHistory.hasBosonParent() * 8 +
                                   particleHistory.hasNonPromptParent() * 4 +
                                   particleHistory.hasQCDParent() * 2 +
                                   particleHistory.hasExoticParent());

            if (particleHistory.hasRealParent())
            {
                reco::GenParticleRef momRef = particleHistory.parent();
                if (momRef.isNonnull() && momRef.isAvailable())
                {
                    mcMomPID_ = momRef->pdgId();
                    mcMomPt_ = momRef->pt();
                    mcMomMass_ = momRef->mass();
                    mcMomEta_ = momRef->eta();
                    mcMomPhi_ = momRef->phi();

                    // get Granny
                    genpartparentage::GenParticleParentage motherParticle(momRef);
                    if (motherParticle.hasRealParent())
                    {
                        reco::GenParticleRef granny = motherParticle.parent();
                        mcGMomPID_ = granny->pdgId();
                    }
                }
            }

            // if (abs(myParticle->pdgId()) == 13 && igen->isPromptFinalState())
            //     cout << "I am " << myParticle->pdgId() << " ,my Mom is " << mcMomPID_ << " ,my Gmom is " << mcGMomPID_ << endl;

            genMomPID.push_back(mcMomPID_);
            genGMomPID.push_back(mcGMomPID_);
            genMomPt.push_back(mcMomPt_);
            genMomMass.push_back(mcMomMass_);
            genMomEta.push_back(mcMomEta_);
            genMomPhi.push_back(mcMomPhi_);

            nGen++;
        }
    }

    edm::Handle<vector<reco::GenJet>> genJetHandle;
    iEvent.getByToken(genJetHandle_, genJetHandle);

    for (vector<reco::GenJet>::const_iterator igen = genJetHandle->begin(); igen != genJetHandle->end(); ++igen)
    {
        nGenJet++;
        genJetPt.push_back(igen->pt());
        genJetEta.push_back(igen->eta());
        genJetPhi.push_back(igen->phi());
        genJetE.push_back(igen->energy());
    }

    // + Get event weight
    //-------------------
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    iEvent.getByToken(generatorLabel_, genEventInfoHandle);

    genWeight = 0.;

    if (genEventInfoHandle.isValid())
        genWeight = genEventInfoHandle->weight();

    tree_out->Fill();
}

// ------------ method called once each job just after ending the event loop  ------------
void testAnalyzer::endJob()
{
    // go to *OUR* root file and store histograms
    rootFile_->cd();
    tree_out->Write();
    rootFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void testAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(testAnalyzer);

void testAnalyzer::refresh()
{
    nGen = 0;
    genPID.clear();
    genStatus.clear();
    genStatusFlag.clear();
    genVtx.clear();
    genVty.clear();
    genVtz.clear();
    genPt.clear();
    genMass.clear();
    genEta.clear();
    genPhi.clear();
    genE.clear();
    genEt.clear();
    genMomPID.clear();
    genGMomPID.clear();
    genMomPt.clear();
    genMomMass.clear();
    genMomEta.clear();
    genMomPhi.clear();
    genParentage.clear();
    nGenJet = 0;
    genJetPt.clear();
    genJetEta.clear();
    genJetPhi.clear();
    genJetE.clear();
    genWeight = 0.;
}