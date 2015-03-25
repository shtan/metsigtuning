#ifndef METSIG_FIT_H
#define METSIG_FIT_H

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
   double sig;
   double sig_init;
   double det;
   double cov_xx;
   double cov_xy;
   double cov_yy;

   double cov_xx_highpt;
   double cov_xx_pjet;
   double cov_dtt;
   double cov_dff;

   double qt;
   double qx;
   double ut;
   double ut_par;
   double ut_perp;

   double resp_correction;

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

   // jets
   vector<double> jet_phi;
   vector<double> jet_eta;
   vector<double> jet_pt;
   vector<double> jet_sigmapt;
   vector<double> jet_sigmaphi;

   event(){
      process = "";

      nvertices = 0;
      weight = 0;

      met = 0;
      sig = 0;
      sig_init = 0;
      det = 0;
      cov_xx = 0;
      cov_xy = 0;
      cov_yy = 0;

      cov_xx_highpt = 0;
      cov_xx_pjet = 0;
      cov_dtt = 0;
      cov_dff = 0;

      qt = 0;
      qx = 0;
      ut = 0;
      ut_par = 0;
      ut_perp = 0;

      resp_correction = 0;

      // pseudojet
      pjet_pt = 0;
      pjet_phi = 0;
      pjet_scalpt = 0;

      // met
      met_pt = 0;
      met_phi = 0;
   }

}; 

class Fitter{
   public:
      Fitter();
      ~Fitter();

      void ReadNtuple(string, vector<event>&);
      void RunMinimizer(vector<event>&);
      void FindSignificance(const double*, vector<event>&);
      void MakePlots(vector<event>&);

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      bool significance_cut;

   private:
      double Min2LL(const double*);

      ROOT::Math::IMultiGenFunction* fFunc;

      vector<event>* eventvecPnt;
};
#endif

