#ifndef READER_H
#define READER_H

#include <vector>
#include <cmath>
#include <map>
#include <TH1.h>
#include <TH2.h>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
using namespace std;

struct event {

   string process;

   int nvertices;
   double weight;

   double met;


   // pseudojet
   double pjet_pt;
   double pjet_phi;
   double pjet_scalpt;

   // variables for ROC
   double metsig2011;

   // met
   double met_pt;
   double met_phi;

   // leptons
   vector<double> lepton_pt;
   vector<double> lepton_phi;

   // genjets
   vector<double> genjet_phi;
   vector<double> genjet_eta;
   vector<double> genjet_pt;

   // jets
   vector<double> jet_phi;
   vector<double> jet_eta;
   vector<double> jet_pt;
   vector<double> jet_energy;
   vector<double> jet_sigmapt;
   vector<double> jet_sigmaphi;
   vector<int> jet_num_constituents;
   vector<double> jet_chargedEmEnergy, jet_neutralEmEnergy, jet_chargedHadronEnergy, jet_neutralHadronEnergy;

   // matched jets
   vector<double> matchedjet_phi;
   vector<double> matchedjet_eta;
   vector<double> matchedjet_pt;
 
   // matched genjets
   vector<double> matchedgenjet_phi;
   vector<double> matchedgenjet_eta;
   vector<double> matchedgenjet_pt;
 
   //matched id's
   vector<int> matchedijets;
   vector<int> matchedigenjets;

   // ratio ptjet to ptgenjet; should be same size as matchedijets and matchedigenjets
   vector<double> ptratios;

   event(){
      process = "";

      nvertices = 0;
      weight = 0;

      met = 0;

      // pseudojet
      pjet_pt = 0;
      pjet_phi = 0;
      pjet_scalpt = 0;

      // met
      met_pt = 0;
      met_phi = 0;
   }

}; 

class Plotter{
   public:
      Plotter();
      ~Plotter();

      void ReadNtuple(string, vector<event>&);
      void Matcher( event& );
      void Analyse(vector<event>&);
      void DeclareHists();
      string num2string( float );

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;


   private:
      double Min2LL(const double*);

      ROOT::Math::IMultiGenFunction* fFunc;

      vector<event>* eventvecPnt;

      double deltaRmax = 0.1;
      map< string, map<string, TH1D*> > hists_;
      map< string, TH1D* > histss_;
      
      double ptBins[17] = {15.,20.,30.,45.,60.,80.,120.,170.,230.,300.,380.,470.,600.,800.,1000.,1400.,1800.};
      double etaBins[11] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.};
      int numConstBins[7] = {0,8,12,16,20,30,50};
      double neutralChargedRatioBins[9] = {0,0.5,1,2,4,6,8,10,20};
      double hadronicEmRatioBins[9] = {0,0.5,1,2,4,6,8,10,20};
      double massBins[10] = {15.,30.,60.,120.,230.,380.,600.,1000.,1800.,3200.};
 

};
#endif

