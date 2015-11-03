// -*- C++ -*-
//
// Package:    METSigTuning/MakeNtuple
// Class:      MakeNtuple
// 
/**\class MakeNtuple MakeNtuple.cc METSigTuning/MakeNtuple/src/MakeNtuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nathan Mirman
//         Created:  Wed, 26 Nov 2014 14:22:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>

//
// class declaration
//

class MakeNtuple : public edm::EDAnalyzer {
   public:
      explicit MakeNtuple(const edm::ParameterSet&);
      ~MakeNtuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

//      std::vector<reco::Jet> cleanJets(double, double,
      std::vector<pat::Jet> cleanJets(double, double,
//            std::vector<reco::Jet>&, std::vector<reco::Candidate::LorentzVector>&);
            std::vector<pat::Jet>&, std::vector<reco::Candidate::LorentzVector>&);
      edm::EDGetTokenT<edm::View<reco::Candidate> > inputToken_;
//      edm::EDGetTokenT<edm::View<reco::Jet> > jetToken_;
      edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
      edm::EDGetTokenT<edm::View<reco::Jet> > genJetToken_;
      std::vector< edm::EDGetTokenT<edm::View<reco::Candidate> > > lepTokens_;
      edm::InputTag metTag_;
      edm::EDGetTokenT<edm::View<reco::Candidate> > muonToken_;
      edm::InputTag verticesTag_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      TTree *results_tree;
      TFile *OutFile__file;
      std::string OutputFileName_;

      Long64_t run, event, lumi;

      std::vector<double> lep_pt, lep_energy, lep_phi, lep_eta;
      std::vector<double> muon_pt, muon_energy, muon_phi, muon_eta;
      std::vector<double> jet_pt, jet_energy, jet_phi, jet_eta;
      std::vector<double> genjet_pt, genjet_energy, genjet_phi, genjet_eta;
      std::vector<double> jet_sigmapt, jet_sigmaphi;
      std::vector<int> jet_num_constituents;
      std::vector<double> jet_chargedEmEnergy, jet_chargedHadronEnergy, jet_neutralEmEnergy, jet_neutralHadronEnergy;
      std::vector<double> jet_emEnergyInHF, jet_hadEnergyInHF;
      std::vector<bool> jet_isCaloJet, jet_isPFJet, jet_isJPTJet, jet_isBasicJet;
      double met_pt, met_energy, met_phi, met_eta, met_sumpt;
      int nvertices;

      double jetThreshold;

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
MakeNtuple::MakeNtuple(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   inputToken_ = consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("src"));
//   jetToken_ = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
   jetToken_ = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
   genJetToken_ = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("genjets"));
   std::vector<edm::InputTag> srcLeptonsTags = iConfig.getParameter< std::vector<edm::InputTag> >("leptons");
   for(std::vector<edm::InputTag>::const_iterator it=srcLeptonsTags.begin();it!=srcLeptonsTags.end();it++) {
      lepTokens_.push_back( consumes<edm::View<reco::Candidate> >( *it ) );
   }
   metTag_ = iConfig.getParameter<edm::InputTag>("met");
   muonToken_ = consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("muons"));
   verticesTag_ = iConfig.getParameter<edm::InputTag>("vertices");

   jetThreshold = 20;

   OutputFileName_ = "ntuple.root";
}


MakeNtuple::~MakeNtuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MakeNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // clear all vectors
   lep_pt.clear();
   lep_energy.clear();
   lep_phi.clear();
   lep_eta.clear();
   muon_pt.clear();
   muon_energy.clear();
   muon_phi.clear();
   muon_eta.clear();
   jet_pt.clear();
   jet_energy.clear();
   jet_phi.clear();
   jet_eta.clear();
   genjet_pt.clear();
   genjet_energy.clear();
   genjet_phi.clear();
   genjet_eta.clear();
   jet_sigmapt.clear();
   jet_sigmaphi.clear();
   jet_num_constituents.clear();
   jet_chargedEmEnergy.clear();
   jet_chargedHadronEnergy.clear();
   jet_neutralEmEnergy.clear();
   jet_neutralHadronEnergy.clear();
   jet_emEnergyInHF.clear();
   jet_hadEnergyInHF.clear();
   jet_isCaloJet.clear();
   jet_isPFJet.clear();
   jet_isBasicJet.clear();
   jet_isJPTJet.clear();

   run = iEvent.id().run();
   event = iEvent.id().event();
   lumi = iEvent.id().luminosityBlock();

   Handle<View<reco::Candidate> > input;
   iEvent.getByToken(inputToken_, input);

   //std::cout << "Leptons: ";
   // leptons
   std::vector<reco::CandidatePtr> footprint;
   std::vector<reco::Candidate::LorentzVector> leptons;
   for ( std::vector<EDGetTokenT<View<reco::Candidate> > >::const_iterator srcLeptons_i = lepTokens_.begin();
         srcLeptons_i != lepTokens_.end(); ++srcLeptons_i ) {
      Handle<reco::CandidateView> leptons_i;
      iEvent.getByToken(*srcLeptons_i, leptons_i);
      for ( reco::CandidateView::const_iterator lepton = leptons_i->begin();
            lepton != leptons_i->end(); ++lepton ) {
         // cut on lepton pt
         //std::cout << lepton->pt() << " ";
         if( lepton->pt() > 10 ){
            leptons.push_back(lepton->p4());
            for( unsigned int n=0; n < lepton->numberOfSourceCandidatePtrs(); n++){
               if( lepton->sourceCandidatePtr(n).isNonnull() and lepton->sourceCandidatePtr(n).isAvailable() ){
                  footprint.push_back(lepton->sourceCandidatePtr(n));
               }
            }
         }
      }
   }
   //std::cout << std::endl;

   // muons (for event selection)
   Handle<reco::CandidateView> muons;
   iEvent.getByToken(muonToken_, muons);
   double nmuons = 0;
   for ( reco::CandidateView::const_iterator muon = muons->begin();
         muon != muons->end(); ++muon ) {
      if( muon->pt() > 20 and fabs(muon->eta()) < 2.4 ){
         muon_pt.push_back( muon->pt() );
         muon_energy.push_back( muon->energy() );
         muon_phi.push_back( muon->phi() );
         muon_eta.push_back( muon->eta() );
         nmuons++;
      }
   }

/*   double dimuon_mass = 0;
   if( nmuons == 2 ){
      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiE( muon_pt[0], muon_eta[0], muon_phi[0], muon_energy[0] );
      mu2.SetPtEtaPhiE( muon_pt[1], muon_eta[1], muon_phi[1], muon_energy[1] );
      dimuon_mass = (mu1+mu2).M();
   }
*/
   // jets
//   Handle<View<reco::Jet>> inputJets;
   Handle<View<pat::Jet>> inputJets;
   iEvent.getByToken( jetToken_, inputJets );
//   std::vector<reco::Jet> jets;
   std::vector<pat::Jet> jets;
//   for(View<reco::Jet>::const_iterator jet = inputJets->begin(); jet != inputJets->end(); ++jet) {
   for(View<pat::Jet>::const_iterator jet = inputJets->begin(); jet != inputJets->end(); ++jet) {
      jets.push_back( *jet );
   }

   // genjets
   Handle<View<reco::Jet>> inputGenJets;

   iEvent.getByToken( genJetToken_, inputGenJets );
   std::vector<reco::Jet> genjets;
   for(View<reco::Jet>::const_iterator genjet = inputGenJets->begin(); genjet != inputGenJets->end(); ++genjet) {
      genjets.push_back( *genjet );
   }

   // disambiguate jets and leptons
//   std::vector<reco::Jet> cleanjets = cleanJets(jetThreshold, 0.4, jets, leptons);
   std::vector<pat::Jet> cleanjets = cleanJets(jetThreshold, 0.4, jets, leptons);

   // loop over jets to disambiguate candidates
//   for(std::vector<reco::Jet>::const_iterator jet = cleanjets.begin(); jet != cleanjets.end(); ++jet) {
   for(std::vector<pat::Jet>::const_iterator jet = cleanjets.begin(); jet != cleanjets.end(); ++jet) {
      for( unsigned int n=0; n < jet->numberOfSourceCandidatePtrs(); n++){
         if( jet->sourceCandidatePtr(n).isNonnull() and jet->sourceCandidatePtr(n).isAvailable() ){
            footprint.push_back(jet->sourceCandidatePtr(n));
         }
      }
   }

   // met
   Handle<View<reco::MET> > metHandle;
   iEvent.getByLabel(metTag_, metHandle);
   reco::MET met = (*metHandle)[0];

   met_pt = met.pt();
   met_energy = met.energy();
   met_phi = met.phi();
   met_eta = met.eta();

   // candidates
   std::vector<reco::Candidate::LorentzVector> candidates;
   for(View<reco::Candidate>::const_iterator cand = input->begin();
         cand != input->end(); ++cand) {
      unsigned int iter = cand - input->begin();
      if (std::find(footprint.begin(), footprint.end(),
               reco::CandidatePtr(input,iter)) != footprint.end()) {
         continue;
      }
      candidates.push_back( cand->p4() );
   }

   // resolutions
   std::string path = "CondFormats/JetMETObjects/data";
   std::string resEra_ = "Spring10";
   std::string resAlg_ = "AK5PF";
   std::string ptResFileName  = path + "/" + resEra_ + "_PtResolution_" +resAlg_+".txt";
   std::string phiResFileName = path + "/" + resEra_ + "_PhiResolution_"+resAlg_+".txt";

   FileInPath fpt(ptResFileName);
   FileInPath fphi(phiResFileName);

   JetResolution *ptRes_  = new JetResolution(fpt.fullPath().c_str(),false);
   JetResolution *phiRes_ = new JetResolution(fphi.fullPath().c_str(),false);

   //
   // begin ttree variables
   //

   // calculate met_sumpt
   met_sumpt = 0;
   for( std::vector<reco::Candidate::LorentzVector>::const_iterator cand = candidates.begin();
         cand != candidates.end(); ++cand){
      met_sumpt += cand->Pt();
   }

   // loop over leptons
   for ( std::vector<reco::Candidate::LorentzVector>::const_iterator lepton = leptons.begin();
         lepton != leptons.end(); ++lepton ) {
      lep_pt.push_back( lepton->Pt() );
      lep_energy.push_back( lepton->E() );
      lep_phi.push_back( lepton->Phi() );
      lep_eta.push_back( lepton->Eta() );
   }

/*   //std::cout << "Jets: ";
   // loop over jets
   for(std::vector<reco::Jet>::const_iterator jet = cleanjets.begin(); jet != cleanjets.end(); ++jet) {
      double jpt  = jet->pt();
      double jeta = jet->eta();
      //std::cout << jpt << " ";

      // jet energy resolutions
      double jeta_res = (fabs(jeta) < 9.9) ? jeta : 9.89; // JetResolutions defined for |eta|<9.9
      TF1* fPtEta    = ptRes_ -> parameterEta("sigma",jeta_res);
      TF1* fPhiEta   = phiRes_-> parameterEta("sigma",jeta_res);
      double sigmapt = fPtEta->Eval(jpt);
      double sigmaphi = fPhiEta->Eval(jpt);
      delete fPtEta;
      delete fPhiEta;

      // split into high-pt and low-pt sector
      if( jpt > jetThreshold ){
         // high-pt jets enter into the covariance matrix via JER

         // subtract the pf constituents in each jet out of the met_sumpt
         for(unsigned int i=0; i < jet->numberOfDaughters(); i++){
            //met_sumpt -= jet->daughter(i)->pt();
         }

         jet_pt.push_back( jet->pt() );
         jet_energy.push_back( jet->energy() );
         jet_phi.push_back( jet->phi() );
         jet_eta.push_back( jet->eta() );
         jet_sigmapt.push_back( sigmapt );
         jet_sigmaphi.push_back( sigmaphi );

      }else{

         // subtract the pf constituents in each jet out of the met_sumpt
         for(unsigned int i=0; i < jet->numberOfDaughters(); i++){
            //met_sumpt -= jet->daughter(i)->pt();
         }
         // add the (corrected) jet to the met_sumpt
         met_sumpt += jpt;

      }
   }
   //std::cout << std::endl;
*/
//   for(std::vector<reco::Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet) {
   for(std::vector<pat::Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet) {
       jet_pt.push_back( jet->pt() );
       jet_energy.push_back( jet->energy() );
       jet_phi.push_back( jet->phi() );
       jet_eta.push_back( jet->eta() );
       jet_num_constituents.push_back( jet->getJetConstituents().size() );
       jet_chargedEmEnergy.push_back( jet->chargedEmEnergy() );
       jet_neutralEmEnergy.push_back( jet->neutralEmEnergy() );
       jet_chargedHadronEnergy.push_back( jet->chargedHadronEnergy() );
       jet_neutralHadronEnergy.push_back( jet->neutralHadronEnergy() );

       if ( jet->isCaloJet()){
            jet_emEnergyInHF.push_back( jet->emEnergyInHF() );
            jet_hadEnergyInHF.push_back( jet->hadEnergyInHF() );
       } else {
            jet_emEnergyInHF.push_back( -99999 );
            jet_hadEnergyInHF.push_back( -99999 );
        }
       //std::cout<< "bool calojet = " << jet->isCaloJet()<<std::endl;
       jet_isCaloJet.push_back( jet->isCaloJet() );
       jet_isPFJet.push_back( jet->isPFJet() );
       jet_isJPTJet.push_back( jet->isJPTJet() );
       jet_isBasicJet.push_back( jet->isBasicJet() );

   }

   for(std::vector<reco::Jet>::const_iterator genjet = genjets.begin(); genjet != genjets.end(); ++genjet) {
       genjet_pt.push_back( genjet->pt() );
       genjet_energy.push_back( genjet->energy() );
       genjet_phi.push_back( genjet->phi() );
       genjet_eta.push_back( genjet->eta() );
   }


   //std::cout << jets.size() << std::endl;
//   for(std::vector<reco::Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet) {
   for(std::vector<pat::Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet) {
       //std::cout << jet->eta() << std::endl;
   }
   //std::cout << cleanjets.size() << std::endl;
   //std::cout << genjets.size()<< std::endl;
   for(std::vector<reco::Jet>::const_iterator jet = genjets.begin(); jet != genjets.end(); ++jet) {
       //std::cout << jet->eta() << std::endl;
   }
   //std::cout << std::endl;

   //std::cout << "sumPt: " << met_sumpt << std::endl;

   // offline primary vertices
   edm::Handle<edm::View<reco::Vertex> > vertices;
   iEvent.getByLabel(verticesTag_, vertices);
   nvertices = int(vertices->size());
   
   //bool pass_selection = (nmuons == 2) and (dimuon_mass > 60) and (dimuon_mass < 120);
   //if( pass_selection ){
      results_tree -> Fill();
   //}

   delete ptRes_;
   delete phiRes_;

}

//   std::vector<reco::Jet>
   std::vector<pat::Jet>
MakeNtuple::cleanJets(double ptThreshold, double dRmatch,
//      std::vector<reco::Jet>& jets, std::vector<reco::Candidate::LorentzVector>& leptons)
      std::vector<pat::Jet>& jets, std::vector<reco::Candidate::LorentzVector>& leptons)
{
   double dR2match = dRmatch*dRmatch;
//   std::vector<reco::Jet> retVal;
   std::vector<pat::Jet> retVal;
//   for ( std::vector<reco::Jet>::const_iterator jet = jets.begin();
   for ( std::vector<pat::Jet>::const_iterator jet = jets.begin();
         jet != jets.end(); ++jet ) {
      bool isOverlap = false;
      for ( std::vector<reco::Candidate::LorentzVector>::const_iterator lepton = leptons.begin();
            lepton != leptons.end(); ++lepton ) {
         TLorentzVector ljet, llep;
         ljet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
         llep.SetPtEtaPhiE( lepton->pt(), lepton->eta(), lepton->phi(), lepton->energy() );
         if ( pow(ljet.DeltaR( llep ),2) < dR2match ) isOverlap = true;
      }
      //if ( jet->pt() > ptThreshold && !isOverlap ){
      if ( !isOverlap ){
         retVal.push_back(*jet);
      }
   }

   return retVal;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeNtuple::beginJob()
{
   OutFile__file  = new TFile( OutputFileName_.c_str(), "RECREATE" );

   results_tree = new TTree("events", "events");
   results_tree -> Branch("run", &run, "run/I");
   results_tree -> Branch("lumi", &lumi, "lumi/I");
   results_tree -> Branch("event", &event, "event/I");

   results_tree -> Branch("muon_pt", &muon_pt);
   results_tree -> Branch("muon_energy", &muon_energy);
   results_tree -> Branch("muon_phi", &muon_phi);
   results_tree -> Branch("muon_eta", &muon_eta);

   results_tree -> Branch("lep_pt", &lep_pt);
   results_tree -> Branch("lep_energy", &lep_energy);
   results_tree -> Branch("lep_phi", &lep_phi);
   results_tree -> Branch("lep_eta", &lep_eta);

   results_tree -> Branch("jet_pt", &jet_pt);
   results_tree -> Branch("jet_energy", &jet_energy);
   results_tree -> Branch("jet_phi", &jet_phi);
   results_tree -> Branch("jet_eta", &jet_eta);
   results_tree -> Branch("jet_sigmapt", &jet_sigmapt);
   results_tree -> Branch("jet_sigmaphi", &jet_sigmaphi);
   results_tree -> Branch("jet_num_constituents", &jet_num_constituents);
   results_tree -> Branch("jet_chargedEmEnergy", &jet_chargedEmEnergy);
   results_tree -> Branch("jet_neutralEmEnergy", &jet_neutralEmEnergy);
   results_tree -> Branch("jet_chargedHadronEnergy", &jet_chargedHadronEnergy);
   results_tree -> Branch("jet_neutralHadronEnergy", &jet_neutralHadronEnergy);
   results_tree -> Branch("jet_emEnergyInHF", &jet_emEnergyInHF);
   results_tree -> Branch("jet_hadEnergyInHF", &jet_hadEnergyInHF);
   results_tree -> Branch("jet_isCaloJet", &jet_isCaloJet);
   results_tree -> Branch("jet_isPFJet", &jet_isPFJet);
   results_tree -> Branch("jet_isJPTJet", &jet_isJPTJet);
   results_tree -> Branch("jet_isBasicJet", &jet_isBasicJet);


   results_tree -> Branch("genjet_pt", &genjet_pt);
   results_tree -> Branch("genjet_energy", &genjet_energy);
   results_tree -> Branch("genjet_phi", &genjet_phi);
   results_tree -> Branch("genjet_eta", &genjet_eta);

   results_tree -> Branch("met_pt", &met_pt);
   results_tree -> Branch("met_energy", &met_energy);
   results_tree -> Branch("met_phi", &met_phi);
   results_tree -> Branch("met_eta", &met_eta);
   results_tree -> Branch("met_sumpt", &met_sumpt);

   results_tree -> Branch("nvertices", &nvertices);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MakeNtuple::endJob() 
{
    std::cout<<"write"<<std::endl;
   OutFile__file -> Write();
   OutFile__file -> Close();
   std::cout<<"written"<<std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MakeNtuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MakeNtuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MakeNtuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MakeNtuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeNtuple);
