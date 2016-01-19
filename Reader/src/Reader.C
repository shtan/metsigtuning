#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TString.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TVirtualFFT.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <list>

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"

#include "Reader.h"

#define array_size(array) (sizeof((array))/sizeof((array[0])))
using namespace std;

//
// constructor and destructor
//

Plotter::Plotter(){

}

Plotter::~Plotter(){
}

void Plotter::ReadNtuple(string filename, vector<event>& eventref_temp){

   std::vector<double> *lep_pt=0, *lep_energy=0, *lep_phi=0, *lep_eta=0;
   std::vector<double> *jet_pt=0, *jet_energy=0, *jet_phi=0, *jet_eta=0;
   std::vector<double> *genjet_pt=0, *genjet_energy=0, *genjet_phi=0, *genjet_eta=0;
   std::vector<double> *jet_sigmapt=0, *jet_sigmaphi=0;
   std::vector<double> *jet_chargedEmEnergy=0, *jet_neutralEmEnergy=0, *jet_chargedHadronEnergy=0, *jet_neutralHadronEnergy=0;
   std::vector<double> *jet_emEnergyInHF=0, *jet_hadEnergyInHF=0;
   std::vector<bool> *jet_isCaloJet=0, *jet_isPFJet=0, *jet_isBasicJet=0, *jet_isJPTJet=0;
   std::vector<int> *jet_num_constituents=0;
   double met_pt=0, met_energy=0, met_phi=0, met_eta=0, met_sumpt=0;
   int nvertices=0;

   TFile *file = new TFile(filename.c_str());
   TTree *tree = (TTree*)file->Get("events");

   tree->SetBranchAddress("lep_pt", &lep_pt);
   tree->SetBranchAddress("lep_energy", &lep_energy);
   tree->SetBranchAddress("lep_phi", &lep_phi);
   tree->SetBranchAddress("lep_eta", &lep_eta);

   tree->SetBranchAddress("jet_pt", &jet_pt);
   tree->SetBranchAddress("jet_energy", &jet_energy);
   tree->SetBranchAddress("jet_phi", &jet_phi);
   tree->SetBranchAddress("jet_eta", &jet_eta);
   tree->SetBranchAddress("jet_sigmapt", &jet_sigmapt);
   tree->SetBranchAddress("jet_sigmaphi", &jet_sigmaphi);
   tree -> SetBranchAddress("jet_num_constituents", &jet_num_constituents);
   tree -> SetBranchAddress("jet_chargedEmEnergy", &jet_chargedEmEnergy);
   tree -> SetBranchAddress("jet_neutralEmEnergy", &jet_neutralEmEnergy);
   tree -> SetBranchAddress("jet_chargedHadronEnergy", &jet_chargedHadronEnergy);
   tree -> SetBranchAddress("jet_neutralHadronEnergy", &jet_neutralHadronEnergy);
   tree -> SetBranchAddress("jet_emEnergyInHF", &jet_emEnergyInHF);
   tree -> SetBranchAddress("jet_hadEnergyInHF", &jet_hadEnergyInHF);
   tree -> SetBranchAddress("jet_isCaloJet", &jet_isCaloJet);
   tree -> SetBranchAddress("jet_isPFJet", &jet_isPFJet);
   tree -> SetBranchAddress("jet_isJPTJet", &jet_isJPTJet);
   tree -> SetBranchAddress("jet_isBasicJet", &jet_isBasicJet);

   tree->SetBranchAddress("genjet_pt", &genjet_pt);
   tree->SetBranchAddress("genjet_energy", &genjet_energy);
   tree->SetBranchAddress("genjet_phi", &genjet_phi);
   tree->SetBranchAddress("genjet_eta", &genjet_eta);

   tree->SetBranchAddress("met_pt", &met_pt);
   tree->SetBranchAddress("met_energy", &met_energy);
   tree->SetBranchAddress("met_phi", &met_phi);
   tree->SetBranchAddress("met_eta", &met_eta);
   tree->SetBranchAddress("met_sumpt", &met_sumpt);

   tree->SetBranchAddress("nvertices", &nvertices);

   for( int ev = 0; ev < tree->GetEntries(); ev++){
   //for( int ev = 0; ev < 10; ev++){
       tree->GetEntry(ev);
       
       event evtemp;
       evtemp.weight = 1.0;
       
       // nvertices
       evtemp.nvertices = nvertices;
        
       // leptons
       evtemp.lepton_pt = *lep_pt;
       evtemp.lepton_phi = *lep_phi;
        
       // jets
       evtemp.genjet_pt = *genjet_pt;
       evtemp.genjet_phi = *genjet_phi;
       evtemp.genjet_eta = *genjet_eta;

       // jets
       evtemp.jet_pt = *jet_pt;
       evtemp.jet_phi = *jet_phi;
       evtemp.jet_eta = *jet_eta;
       evtemp.jet_energy = *jet_energy;
       evtemp.jet_sigmapt = *jet_sigmapt;
       evtemp.jet_sigmaphi = *jet_sigmaphi;
       evtemp.jet_num_constituents = *jet_num_constituents;
       evtemp.jet_chargedEmEnergy = *jet_chargedEmEnergy;
       evtemp.jet_neutralEmEnergy = *jet_neutralEmEnergy;
       evtemp.jet_chargedHadronEnergy = *jet_chargedHadronEnergy;
       evtemp.jet_neutralHadronEnergy = *jet_neutralHadronEnergy;
       
       // met
       evtemp.met_pt = met_pt;
       evtemp.met_phi = met_phi;
       
       // pseudo-jet
       evtemp.pjet_scalpt = met_sumpt;
       
       eventref_temp.push_back( evtemp );
   }

   return;

}

string Plotter::num2string(float i){
    stringstream ss;
    ss << i;
    string strnum = ss.str();
    return strnum;
}

void Plotter::Analyse(vector<event>& eventref_temp){

    hists_["all"]["all"] = new TH1D("hall", "hall", 50, 0, 2);
    for( int j=0; j<(int)( array_size(etaBins) - 1 ); j++){
        //stringstream ss;
        //ss << etaBins.at(j);
        //string strnum = ss.str();
        string strnum = num2string(etaBins[j]);
        //std::cout<<" ini j = " << j << std::endl;
        //std::cout<<" etaBinsj = "<< etaBins[j]<<std::endl;
        //std::cout<<"strnum = "<< strnum<<std::endl;
        hists_["eta"][strnum] = new TH1D( ("heta"+strnum).c_str(), ("heta"+strnum).c_str(), 50, 0, 2);

    }


    for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){
        Matcher( *ev );

        for( int i=0; i<(int)( ev->matchedijets.size()); i++){

            double ptratio = ev->jet_pt.at(ev->matchedijets.at(i)) / ev->genjet_pt.at(ev->matchedigenjets.at(i));
            hists_["all"]["all"] -> Fill( ptratio );
            //std::cout<<"fill all" << std::endl;

            for( int j=0; j<(int)( array_size(etaBins) - 1); j++){
                double eta = ev->jet_eta.at(ev->matchedijets.at(i));
                if( eta >= etaBins[j] and eta < etaBins[j+1] ){

                    hists_["eta"][num2string(etaBins[j])] -> Fill( ptratio );
                    
                    for( int jpt=0; jpt<(int)( array_size(ptBins) - 1); jpt++){
                        double pt = ev->jet_pt.at(ev->matchedijets.at(i));
                        if( pt >= ptBins[jpt] and pt < ptBins[jpt+1] ){
                            histss_["res_pt" + num2string(ptBins[jpt]) + "_eta" + num2string(etaBins[j])] -> Fill( ptratio );
                        }
                    }//end pt loop
                }
                //std::cout<< "j=" << j <<std::endl;
            }//end eta loop

            for( int imass=0; imass<(int)( array_size(massBins) - 1); imass++){
                TLorentzVector jetlorentz;
                int matchedi = ev->matchedijets.at(i);
                jetlorentz.SetPtEtaPhiE(ev->jet_pt.at(matchedi), ev->jet_eta.at(matchedi), ev->jet_phi.at(matchedi), ev->jet_energy.at(matchedi) );
                double jetmass = jetlorentz.M();

                if( jetmass >= massBins[imass] and jetmass < massBins[imass+1] ){
                    histss_["res_jetMass_index" + num2string( (double)(imass)) + "_" + num2string( (double)(massBins[imass]) ) + "to" + num2string( (double)(massBins[imass+1]) )] -> Fill( ptratio);
                    std::cout<<"mass"<<std::endl;
                }
            }


            for( int inumconst=0; inumconst<(int)( array_size(numConstBins) - 1); inumconst++){
                int numconst = ev->jet_num_constituents.at(ev->matchedijets.at(i));
                if( numconst >= numConstBins[inumconst] and numconst < numConstBins[inumconst+1] ){
                    
                    histss_["res_numConstituents_index" + num2string( (double)(inumconst)) + "_" + num2string( (double)(numConstBins[inumconst]) ) + "to" + num2string( (double)(numConstBins[inumconst+1]) )] -> Fill( ptratio);
                    std::cout<<"numconst"<<std::endl;
                }
            }//end numconst loop

            for( int ineutralChargedRatio=0; ineutralChargedRatio<(int)( array_size(neutralChargedRatioBins) - 1); ineutralChargedRatio++){

                double chargedEmEnergy = ev->jet_chargedEmEnergy.at(ev->matchedijets.at(i));
                double chargedHadronEnergy = ev->jet_chargedHadronEnergy.at(ev->matchedijets.at(i));
                double neutralEmEnergy = ev->jet_neutralEmEnergy.at(ev->matchedijets.at(i));
                double neutralHadronEnergy = ev->jet_neutralHadronEnergy.at(ev->matchedijets.at(i));
                double neutralChargedRatio = (neutralEmEnergy + neutralHadronEnergy)/(chargedEmEnergy + chargedHadronEnergy);
                if( neutralChargedRatio >= neutralChargedRatioBins[ineutralChargedRatio] and neutralChargedRatio < neutralChargedRatioBins[ineutralChargedRatio+1] ){

                    if (ev->jet_eta.at(ev->matchedijets.at(i)) > 2.5 ) continue;

                    histss_["res_neutralChargedRatio_index" + num2string( (double)(ineutralChargedRatio) ) +"_" + num2string( (double)(neutralChargedRatioBins[ineutralChargedRatio]) ) + "to" + num2string( (double)(neutralChargedRatioBins[ineutralChargedRatio + 1])) ]  -> Fill(ptratio);
                    std::cout<<"neutralcharged"<<std::endl;
                }
            }//end neutralChargedRatio loop

            for( int ihadronicEmRatio=0; ihadronicEmRatio<(int)( array_size(hadronicEmRatioBins) - 1); ihadronicEmRatio++){
                double chargedEmEnergy = ev->jet_chargedEmEnergy.at(ev->matchedijets.at(i));
                double chargedHadronEnergy = ev->jet_chargedHadronEnergy.at(ev->matchedijets.at(i));
                double neutralEmEnergy = ev->jet_neutralEmEnergy.at(ev->matchedijets.at(i));
                double neutralHadronEnergy = ev->jet_neutralHadronEnergy.at(ev->matchedijets.at(i));
                double hadronicEmRatio = (chargedHadronEnergy + neutralHadronEnergy)/(chargedEmEnergy + neutralEmEnergy);

                if( hadronicEmRatio >= hadronicEmRatioBins[ihadronicEmRatio] and hadronicEmRatio < hadronicEmRatioBins[ihadronicEmRatio+1] ){
                    
                    histss_["res_hadronicEmRatio_index" + num2string((double)(ihadronicEmRatio)) +"_" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio]) ) + "to" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio+1]) ) ] -> Fill(ptratio);

                    if (ev->jet_eta.at(ev->matchedijets.at(i)) <= 2.5){
                         histss_["res_hadronicEmRatio_etalt2p5_index" + num2string((double)(ihadronicEmRatio)) +"_" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio]) ) + "to" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio+1]) ) ] -> Fill(ptratio);
                    }

                    if (ev->jet_eta.at(ev->matchedijets.at(i)) >= 2.5){
                         histss_["res_hadronicEmRatio_etagt3p5_index" + num2string((double)(ihadronicEmRatio)) +"_" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio]) ) + "to" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio+1]) ) ] -> Fill(ptratio);
                    }




                    std::cout<<"hadronicem"<<std::endl;
                }
                std::cout<<"outloop"<<std::endl;
            }//end hadronicEmRatio loop

                
        }
    }

    TFile *outfile = new TFile("plots.root","RECREATE");
    outfile->cd();

    typedef map<string, TH1D*> tmap;
    typedef map<string, tmap> hmap;
    for(hmap::iterator h = hists_.begin(); h != hists_.end(); h++){
        int w = 0;
        for(tmap::iterator hin = h->second.begin(); hin != h->second.end(); hin++){
            hin->second->Write();
            std::cout << "w = " << w << std::endl;
            w++;
        }
        //h->second.second -> Write();
    }

    std::vector<int> colourArray = {1,4,9,6,2,5,3,7,8};
    std::vector<double> maxRangeHadronicEmRatioArray, maxRangeNeutralChargedRatioArray, maxRangeHadronicEmRatio2p5Array, maxRangeHadronicEmRatio3p5Array, maxRangeJetMassArray;
    for( tmap::iterator h = histss_.begin(); h != histss_.end(); h++){
    //for( int iRange=0; iRange<(int)(hadronicEmRatioBins.size()); iRange++){
        std::size_t foundHadronicEmRatio = (h->first).find("hadronicEmRatio_index");
        if( foundHadronicEmRatio != string::npos ){
            h->second->Scale(10000/h->second->Integral());
        }
        std::size_t foundHadronicEmRatio2p5 = (h->first).find("hadronicEmRatio_etalt2p5");
        if( foundHadronicEmRatio2p5 != string::npos ){
            h->second->Scale(10000/h->second->Integral());
        }
        std::size_t foundHadronicEmRatio3p5 = (h->first).find("hadronicEmRatio_etagt3p5");
        if( foundHadronicEmRatio3p5 != string::npos ){
            h->second->Scale(10000/h->second->Integral());
        }
        std::size_t foundNeutralChargedRatio = (h->first).find("neutralChargedRatio");
        if( foundNeutralChargedRatio != string::npos){
            h->second->Scale(10000/h->second->Integral());
        }
        std::size_t foundJetMass = (h->first).find("jetMass");
        if( foundJetMass != string::npos){
            h->second->Scale(10000/h->second->Integral());
        }
    }
    for( tmap::iterator h = histss_.begin(); h != histss_.end(); h++){
    //for( int iRange=0; iRange<(int)(hadronicEmRatioBins.size()); iRange++){
        std::size_t foundHadronicEmRatio = (h->first).find("hadronicEmRatio_index");
        if( foundHadronicEmRatio != string::npos ){
            maxRangeHadronicEmRatioArray.push_back( h->second->GetMaximum() );
        }
        std::size_t foundHadronicEmRatio2p5 = (h->first).find("hadronicEmRatio_etalt2p5");
        if( foundHadronicEmRatio2p5 != string::npos ){
            maxRangeHadronicEmRatio2p5Array.push_back( h->second->GetMaximum() );
        }
        std::size_t foundHadronicEmRatio3p5 = (h->first).find("hadronicEmRatio_etagt3p5");
        if( foundHadronicEmRatio3p5 != string::npos ){
            maxRangeHadronicEmRatio3p5Array.push_back( h->second->GetMaximum() );
        }
        std::size_t foundNeutralChargedRatio = (h->first).find("neutralChargedRatio");
        if( foundNeutralChargedRatio != string::npos){
            maxRangeNeutralChargedRatioArray.push_back( h->second->GetMaximum() );
        }
        std::size_t foundJetMass = (h->first).find("jetMass");
        if( foundJetMass != string::npos){
            maxRangeJetMassArray.push_back( h->second->GetMaximum() );
        }
    }
    auto maxRangeHadronicEmRatioIt = std::max_element(std::begin(maxRangeHadronicEmRatioArray), std::end(maxRangeHadronicEmRatioArray));
    double maxRangeHadronicEmRatio = *maxRangeHadronicEmRatioIt;
    auto maxRangeHadronicEmRatio2p5It = std::max_element(std::begin(maxRangeHadronicEmRatio2p5Array), std::end(maxRangeHadronicEmRatio2p5Array));
    double maxRangeHadronicEmRatio2p5 = *maxRangeHadronicEmRatio2p5It;
    auto maxRangeHadronicEmRatio3p5It = std::max_element(std::begin(maxRangeHadronicEmRatio3p5Array), std::end(maxRangeHadronicEmRatio3p5Array));
    double maxRangeHadronicEmRatio3p5 = *maxRangeHadronicEmRatio3p5It;
    std::cout<<"LALALA"<<std::endl;
    auto maxRangeNeutralChargedRatioIt = std::max_element(std::begin(maxRangeNeutralChargedRatioArray), std::end(maxRangeNeutralChargedRatioArray));
    std::cout<<"0"<<std::endl;
    double maxRangeNeutralChargedRatio = *maxRangeNeutralChargedRatioIt;
    auto maxRangeJetMassIt = std::max_element(std::begin(maxRangeJetMassArray), std::end(maxRangeJetMassArray));
    std::cout<<"0"<<std::endl;
    double maxRangeJetMass = *maxRangeJetMassIt;
    std::cout<<"1"<<std::endl;
    TCanvas *chadronicEmRatio = new TCanvas("chadronicEmRatio","chadronicEmRatio",700,700);
    TCanvas *chadronicEmRatio2p5 = new TCanvas("chadronicEmRatio2p5","chadronicEmRatio2p5",700,700);
    TCanvas *chadronicEmRatio3p5 = new TCanvas("chadronicEmRatio3p5","chadronicEmRatio3p5",700,700);
    TCanvas *cneutralChargedRatio = new TCanvas("cneutralChargedRatio","cneutralChargedRatio",700,700);
    TCanvas *cjetMass = new TCanvas("cjetMass","cjetMass",700,700);
    chadronicEmRatio->cd();

    std::cout<<"2"<<std::endl;
    int colourRotate = 0;
    for(tmap::iterator h = histss_.begin(); h != histss_.end(); h++){
        std::cout<<"3"<<std::endl;
        h->second->Write();
        std::size_t foundHadronicEmRatio = (h->first).find("hadronicEmRatio_index");
        if( foundHadronicEmRatio != string::npos ){
            //chadronicEmRatio->cd();
            std::cout<<"4"<<std::endl;
            h->second->SetLineColor(colourArray.at(colourRotate % (int)(colourArray.size())));
            if( colourRotate==0){
                //h->second->SetMaximum(maxRangeHadronicEmRatio + 100);
                //h->second->Scale(1000/h->second->Integral());
                h->second->SetMaximum(maxRangeHadronicEmRatio + 100);
                h->second->DrawCopy();
            } else {
                //h->second->Scale(1000/h->second->Integral());
                h->second->DrawCopy("same");
            }
            colourRotate++;
            //chadronicEmRatio->ls();
        }
    }
    std::cout<<"HAHAHA"<<std::endl;
    chadronicEmRatio->Write();
    chadronicEmRatio->ls();

    chadronicEmRatio2p5->cd();

    std::cout<<"2"<<std::endl;
    colourRotate = 0;
    for(tmap::iterator h = histss_.begin(); h != histss_.end(); h++){
        std::cout<<"3"<<std::endl;
        //h->second->Write();
        std::size_t foundHadronicEmRatio2p5 = (h->first).find("hadronicEmRatio_etalt2p5");
        if( foundHadronicEmRatio2p5 != string::npos ){
            //chadronicEmRatio->cd();
            std::cout<<"4"<<std::endl;
            h->second->SetLineColor(colourArray.at(colourRotate % (int)(colourArray.size())));
            if( colourRotate==0){
                //h->second->SetMaximum(maxRangeHadronicEmRatio + 100);
                //h->second->Scale(1000/h->second->Integral());
                h->second->SetMaximum(maxRangeHadronicEmRatio2p5 + 100);
                h->second->DrawCopy();
            } else {
                //h->second->Scale(1000/h->second->Integral());
                h->second->DrawCopy("same");
            }
            colourRotate++;
            //chadronicEmRatio->ls();
        }
    }
    std::cout<<"HAHAHA"<<std::endl;
    chadronicEmRatio2p5->Write();
    chadronicEmRatio2p5->ls();

    chadronicEmRatio3p5->cd();

    std::cout<<"2"<<std::endl;
    colourRotate = 0;
    for(tmap::iterator h = histss_.begin(); h != histss_.end(); h++){
        std::cout<<"3"<<std::endl;
        //h->second->Write();
        std::size_t foundHadronicEmRatio3p5 = (h->first).find("hadronicEmRatio_etagt3p5");
        if( foundHadronicEmRatio3p5 != string::npos ){
            //chadronicEmRatio->cd();
            std::cout<<"4"<<std::endl;
            h->second->SetLineColor(colourArray.at(colourRotate % (int)(colourArray.size())));
            if( colourRotate==0){
                //h->second->SetMaximum(maxRangeHadronicEmRatio + 100);
                //h->second->Scale(1000/h->second->Integral());
                h->second->SetMaximum(maxRangeHadronicEmRatio3p5 + 100);
                h->second->DrawCopy();
            } else {
                //h->second->Scale(1000/h->second->Integral());
                h->second->DrawCopy("same");
            }
            colourRotate++;
            //chadronicEmRatio->ls();
        }
    }
    std::cout<<"HAHAHA"<<std::endl;
    chadronicEmRatio3p5->Write();
    chadronicEmRatio3p5->ls();

    cneutralChargedRatio->cd();

    colourRotate = 0;
    for(tmap::iterator h = histss_.begin(); h != histss_.end(); h++){

        std::size_t foundNeutralChargedRatio = (h->first).find("neutralChargedRatio");
        if( foundNeutralChargedRatio != string::npos ){
            //cneutralChargedRatio->cd();
            h->second->SetLineColor(colourArray.at(colourRotate % (int)(colourArray.size())));
            if( colourRotate==0){
                //h->second->SetMaximum(maxRangeNeutralChargedRatio + 100);
                //h->second->Scale(1000/h->second->Integral());
                h->second->SetMaximum(maxRangeNeutralChargedRatio + 100);
                h->second->DrawCopy();
            } else {
                //h->second->Scale(1000/h->second->Integral());
                h->second->DrawCopy("same");
            }
            colourRotate++;
            //cneutralChargedRatio->ls();
        }
            

    }
    std::cout<<"GAGAGA"<<std::endl;
    //chadronicEmRatio->Write();
    cneutralChargedRatio->Write();
    //chadronicEmRatio->ls();
    cneutralChargedRatio->ls();
    //hists_["all"]["all"] -> Write();

    cjetMass->cd();

    colourRotate = 0;
    for(tmap::iterator h = histss_.begin(); h != histss_.end(); h++){

        std::size_t foundJetMass = (h->first).find("jetMass");
        if( foundJetMass != string::npos ){
            //cneutralChargedRatio->cd();
            h->second->SetLineColor(colourArray.at(colourRotate % (int)(colourArray.size())));
            if( colourRotate==0){
                //h->second->SetMaximum(maxRangeNeutralChargedRatio + 100);
                //h->second->Scale(1000/h->second->Integral());
                h->second->SetMaximum(maxRangeJetMass + 100);
                h->second->DrawCopy();
            } else {
                //h->second->Scale(1000/h->second->Integral());
                h->second->DrawCopy("same");
            }
            colourRotate++;
            //cneutralChargedRatio->ls();
        }
            

    }
    std::cout<<"GAGAGA"<<std::endl;
    //chadronicEmRatio->Write();
    cjetMass->Write();
    //chadronicEmRatio->ls();
    cjetMass->ls();
 
    
    outfile->Write();
    outfile-> Close();

    for(tmap::iterator h = histss_.begin(); h != histss_.end(); h++){
        //std::size_t foundHadronicEmRatio = (h->first).find("hadronicEmRatio");
        //if( foundHadronicEmRatio != string::npos ){
            double mean = h->second->GetMean(1);
            double rms = h->second->GetRMS();
            std::cout << h->first << ": Mean = " << mean << "; RMS = " << rms << std::endl;
        //}
 /*           
        std::size_t foundNeutralChargedRatio = (h->first).find("neutralChargedRatio");
        if( foundNeutralChargedRatio != string::npos ){
            double mean = h->second->GetMean(1);
            double rms = h->second->GetRMS();
            std::cout << h->first << ": Mean = " << mean << "; RMS = " << rms << std::endl;
        }
 */           

    }


}

void Plotter::DeclareHists(){

    for( int iptbins=0; iptbins < (int)(array_size(ptBins)); iptbins++ ){
        for( int ietabins=0; ietabins < (int)(array_size(etaBins) - 1 ); ietabins++){
            
            string histName = "res_pt" + num2string(ptBins[iptbins]) + "_eta" + num2string(etaBins[ietabins]);
            histss_[histName] = new TH1D(histName.c_str(), histName.c_str(), 50, 0, 2);

        }
    }

    for( int imass=0; imass<(int)( array_size(massBins) - 1); imass++){
            
        string histName = "res_jetMass_index" + num2string((double)(imass)) + "_" + num2string( (double)(massBins[imass])) + "to" + num2string( (double)(massBins[imass+1] )) ;
        histss_[histName] = new TH1D(histName.c_str(), histName.c_str(), 50, 0, 2);
    }

    for( int inumconst=0; inumconst<(int)( array_size(numConstBins) - 1); inumconst++){
            
        string histName = "res_numConstituents_index" + num2string((double)(inumconst)) + "_" + num2string( (double)(numConstBins[inumconst])) + "to" + num2string( (double)(numConstBins[inumconst+1] )) ;
        histss_[histName] = new TH1D(histName.c_str(), histName.c_str(), 50, 0, 2);
    }
    std::cout<<"hoho"<<std::endl;
    for( int ineutralChargedRatio=0; ineutralChargedRatio<(int)( array_size(neutralChargedRatioBins) - 1); ineutralChargedRatio++){
                     
        string histName = "res_neutralChargedRatio_index" + num2string((double)(ineutralChargedRatio)) +"_" + num2string( (double)(neutralChargedRatioBins[ineutralChargedRatio])) + "to" + num2string( (double)(neutralChargedRatioBins[ineutralChargedRatio+1])) ;

        histss_[histName] = new TH1D(histName.c_str(), histName.c_str(), 50, 0, 2);
     }//end neutralChargedRatio loop
    std::cout<<"hehe"<<std::endl;
     for( int ihadronicEmRatio=0; ihadronicEmRatio<(int)( array_size(hadronicEmRatioBins) - 1); ihadronicEmRatio++){
        string histName = "res_hadronicEmRatio_index" + num2string((double)(ihadronicEmRatio)) +"_" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio])) + "to" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio+1]) ) ;

        string histName2 = "res_hadronicEmRatio_etalt2p5_index" + num2string((double)(ihadronicEmRatio)) +"_" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio])) + "to" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio+1]) ) ;

        string histName3 = "res_hadronicEmRatio_etagt3p5_index" + num2string((double)(ihadronicEmRatio)) +"_" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio])) + "to" + num2string( (double)(hadronicEmRatioBins[ihadronicEmRatio+1]) ) ;

        histss_[histName] = new TH1D(histName.c_str(), histName.c_str(), 50, 0, 2);
        histss_[histName2] = new TH1D(histName2.c_str(), histName2.c_str(), 50, 0, 2);
        histss_[histName3] = new TH1D(histName3.c_str(), histName3.c_str(), 50, 0, 2);
     }//end hadronicEmRatio loop

     std::cout<<"hihi"<<std::endl;

}

void Plotter::Matcher( event &evt ){

//   for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){
       
    event* ev = &evt;

        vector<int> ijets, igenjets;
        ijets.clear();
        igenjets.clear();

       // Match jets with genjets
       for( int ijet=0; ijet < (int)(ev->jet_pt.size()); ijet++){
           //std::cout<< "lulu"<< ev->jet_eta.at(ijet) << std::endl;
          for( int igenjet=0; igenjet < (int)(ev->genjet_pt.size()); igenjet++ ){
              //std::cout<< ev->genjet_eta.at(igenjet) << std::endl;
             
             double deltaEta = ev->jet_eta.at(ijet) - ev->genjet_eta.at(igenjet);
             double deltaPhi = ev->jet_phi.at(ijet) - ev->genjet_phi.at(igenjet);
             double deltaR2 = pow(deltaEta,2) + pow(deltaPhi,2);
             double deltaR = sqrt(deltaR2);
             //std::cout<< "deltaR = " << deltaR << std::endl;

             if( deltaR < deltaRmax ){
                 std::cout <<"ha" << std::endl;
                 std::cout << ijet <<" "<< igenjet << std::endl;
                ijets.push_back(ijet);
                igenjets.push_back(igenjet);
             }

           }
        }

       std::cout<<std::endl;
       for(int i=0; i<(int)(igenjets.size() ); i++){
           std::cout<< igenjets.at(i) <<" ";
    }
       std::cout<<std::endl;

        vector<int> iijets_todelete;

        for( int i=0; i<(int)(ijets.size()); i++){
            std::cout << ijets.at(i) << " ";
            for (int j=0; j<i; j++){
                if( i != j and ijets.at(i) == ijets.at(j) ) {
                    iijets_todelete.push_back(j);
                    iijets_todelete.push_back(i);
                }
            }
        }
        std::cout<<std::endl;
                    
        std::sort( iijets_todelete.begin(), iijets_todelete.end() );
        std::unique( iijets_todelete.begin(), iijets_todelete.end() );
        std::sort( iijets_todelete.begin(), iijets_todelete.end(), greater<int>() );

        for ( int i=0; i<(int)(iijets_todelete.size()); i++){
            ijets.erase( ijets.begin() + iijets_todelete.at(i) );
            igenjets.erase( igenjets.begin() + iijets_todelete.at(i) );
        }

        for( int i=0; i<(int)(ijets.size()); i++){
            std::cout << ijets.at(i) << " ";
        }
        std::cout<<std::endl;




        vector<int> iigenjets_todelete;

        for( int i=0; i<(int)(igenjets.size()); i++){
            std::cout << igenjets.at(i) << " ";
            for (int j=0; j<i; j++){
                if( i != j and igenjets.at(i) == igenjets.at(j) ) {
                    iigenjets_todelete.push_back(j);
                    iigenjets_todelete.push_back(i);
                }
            }
        }
        std::cout<<std::endl;
                    
        std::sort( iigenjets_todelete.begin(), iigenjets_todelete.end() );
        std::unique( iigenjets_todelete.begin(), iigenjets_todelete.end() );
        std::sort( iigenjets_todelete.begin(), iigenjets_todelete.end(), greater<int>() );

        for ( int i=0; i<(int)(iigenjets_todelete.size()); i++){
            ijets.erase( ijets.begin() + iigenjets_todelete.at(i) );
            igenjets.erase( igenjets.begin() + iigenjets_todelete.at(i) );
        }

        for( int i=0; i<(int)(ijets.size()); i++){
            std::cout << ijets.at(i) << " ";
        }
        std::cout<<std::endl;


        for( int i=0; i<(int)(igenjets.size()); i++){
            std::cout << igenjets.at(i) << " ";
        }
        std::cout<<std::endl;
 
        ev->matchedijets = ijets;
        ev->matchedigenjets = igenjets;


    //}

}





