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

#include "METSigFit.h"
using namespace std;

//
// constructor and destructor
//

Fitter::Fitter(){
   // MINUIT variables
   gMinuit = 0;
   fFunc = 0;

   // significance cut for minimization
   significance_cut = false;
}

Fitter::~Fitter(){
   if (gMinuit) delete gMinuit;
   if (fFunc) delete fFunc;
}

void Fitter::ReadNtuple(string filename, vector<event>& eventref_temp){

   std::vector<double> *lep_pt=0, *lep_energy=0, *lep_phi=0, *lep_eta=0;
   std::vector<double> *jet_pt=0, *jet_energy=0, *jet_phi=0, *jet_eta=0;
   std::vector<double> *jet_sigmapt=0, *jet_sigmaphi=0;
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

   tree->SetBranchAddress("met_pt", &met_pt);
   tree->SetBranchAddress("met_energy", &met_energy);
   tree->SetBranchAddress("met_phi", &met_phi);
   tree->SetBranchAddress("met_eta", &met_eta);
   tree->SetBranchAddress("met_sumpt", &met_sumpt);

   tree->SetBranchAddress("nvertices", &nvertices);

   for( int ev = 0; ev < tree->GetEntries(); ev++){
      tree->GetEntry(ev);

      event evtemp;
      evtemp.weight = 1.0;

      // nvertices
      evtemp.nvertices = nvertices;

      // leptons
      evtemp.lepton_pt = *lep_pt;
      evtemp.lepton_phi = *lep_phi;

      // jets
      evtemp.jet_pt = *jet_pt;
      evtemp.jet_phi = *jet_phi;
      evtemp.jet_eta = *jet_eta;
      evtemp.jet_sigmapt = *jet_sigmapt;
      evtemp.jet_sigmaphi = *jet_sigmaphi;

      // met
      evtemp.met_pt = met_pt;
      evtemp.met_phi = met_phi;

      // pseudo-jet
      evtemp.pjet_scalpt = met_sumpt;

      eventref_temp.push_back( evtemp );
   }

   return;
}

void Fitter::FindSignificance(const double *x, vector<event>& eventref_temp){

   for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){

      // metsig covariance
      double cov_xx = 0;
      double cov_xy = 0;
      double cov_yy = 0;

      // jets
      for(unsigned int i=0; i < ev->jet_pt.size(); i++){
         double jpt = ev->jet_pt[i];
         double jeta = ev->jet_eta[i];
         double feta = fabs(jeta);
         double c = cos(ev->jet_phi[i]);
         double s = sin(ev->jet_phi[i]);

         int index=-1;
         if(feta<0.5) index=0;
         else if(feta<1.1) index=1;
         else if(feta<1.7) index=2;
         else if(feta<2.3) index=3;
         else{
            index=4;
         }

         double dpt = x[index]*jpt*ev->jet_sigmapt[i];
         double dph =          jpt*ev->jet_sigmaphi[i];

         double dtt = dpt*dpt;
         double dff = dph*dph;
         cov_xx += dtt*c*c + dff*s*s;
         cov_xy += (dtt-dff)*c*s;
         cov_yy += dff*c*c + dtt*s*s;
      }

      // pseudo-jet
      double ctt = x[5]*x[5] + x[6]*x[6]*(ev->pjet_scalpt);
      cov_xx += ctt;
      cov_yy += ctt;

      // compute significance
      double met_x = ev->met_pt * cos(ev->met_phi);
      double met_y = ev->met_pt * sin(ev->met_phi);

      double det = cov_xx*cov_yy - cov_xy*cov_xy;
      double ncov_xx = cov_yy / det;
      double ncov_xy = -cov_xy / det;
      double ncov_yy = cov_xx / det;

      double sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy;

      // load into eventvec
      ev->sig = sig;
      ev->det = det;

      ev->cov_xx = cov_xx;
      ev->cov_xy = cov_xy;
      ev->cov_yy = cov_yy;
      std::cout << "SIG: " << sig << std::endl;

      if( ev->pjet_scalpt < 0 ) cout << ev->pjet_scalpt << endl;
   }

   return;
}

void Fitter::RunMinimizer(vector<event>& eventref_temp){

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );

   gMinuit->SetStrategy(0);
   gMinuit->SetPrintLevel(2);
   gMinuit->SetTolerance(0.1);

   fFunc = new ROOT::Math::Functor ( this, &Fitter::Min2LL, 7 );
   gMinuit->SetFunction( *fFunc );
   gMinuit->SetVariable(0, "a1", 1.0, 0.05);
   gMinuit->SetVariable(1, "a2", 1.0, 0.05);
   gMinuit->SetVariable(2, "a3", 1.0, 0.05);
   gMinuit->SetVariable(3, "a4", 1.0, 0.05);
   gMinuit->SetVariable(4, "a5", 1.0, 0.05);
   gMinuit->SetVariable(5, "N1", 0.0, 0.05);
   gMinuit->SetVariable(6, "S1", 0.5, 0.05);

   // set event vector and minimize
   cout << " -----> minimize, first pass" << endl;
   eventvecPnt = &eventref_temp;
   gMinuit->Minimize();

   // significance for cut
   for( vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){
      ev->sig_init = ev->sig;
   }

   // minimize core sig
   cout << " -----> minimize, core sig" << endl;
   significance_cut = true;
   gMinuit->SetStrategy(1);
   gMinuit->Minimize();

   significance_cut = false;

   // load best-fit significance values
   cout << " -----> fill event vec with best-fit significance" << endl;
   const double *xmin = gMinuit->X();
   FindSignificance(xmin, eventref_temp);

}

double Fitter::Min2LL(const double *x){

   // load significance values into eventvec
   FindSignificance(x, *eventvecPnt);

   // event loop
   double m2ll = 0;
   for( vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){
      if( !(significance_cut and ev->sig_init > 9) ){
         m2ll += ev->weight*(ev->sig + log(ev->det));
      }
   }

   return m2ll;
}

void Fitter::MakePlots(vector<event>& eventref){
   
   // histograms
   map<string, TH1D*> hists_;
   map<string, TH2D*> hists2d_;
   map<string, TProfile*> profs_;

   hists_["sig"] = new TH1D("hsig", "Significance", 50, 0, 50);
   hists_["pchi2"] = new TH1D("hpchi2", "#chi^{2} Probability", 50, 0, 1);
   hists_["nvert"] = new TH1D("hnvert", "Number of Vertices", 50, 0, 50);
   hists2d_["sig_nvert"] = new TH2D("hsig_nvert", "Significance vs. N Vertices", 50, 0, 50, 50, 0, 50);
   profs_["sig_nvert"] = new TProfile("psig_nvert", "<Significance> vs. N Vertices", 50, 0, 50, 0, 10);
   hists_["njets"] = new TH1D("njets", "Number of jets", 10, 0, 10);
   hists_["sumEt"] = new TH1D("sumEt", "sumEt object", 50, 0, 2000);

   hists_["pchi2_0jets"] = new TH1D("hpchi2_0jets", "#chi^{2} Probability, 0 jet events", 50, 0, 1);
   hists_["pchi2_mt0jets"] = new TH1D("hpchi2_mt0jets", "#chi^{2} Probability, >0 jet events", 50, 0, 1);
   profs_["sig_nvert_0jets"] = new TProfile("psig_nvert_0jets", "<Significance> vs. N Vertices, 0 jet events", 50, 0, 50, 0, 10);
   profs_["sig_nvert_mt0jets"] = new TProfile("psig_nvert_mt0jets", "<Significance> vs. N Vertices, >0 jet events", 50, 0, 50, 0, 10);
   hists_["pchi2_sumEtl500"] = new TH1D("hpchi2_sumEtl500", "#chi^{2} Probability, sumEt < 500", 50, 0, 1);
   profs_["sig_nvert_sumEtl500"] = new TProfile("psig_nvert_sumEtl500", "<Significance> vs. N Vertices, sumEt < 500", 50, 0, 50, 0, 10);
   hists_["pchi2_sumEtl700"] = new TH1D("hpchi2_sumEtl700", "#chi^{2} Probability, sumEt < 700", 50, 0, 1);
   profs_["sig_nvert_sumEtl700"] = new TProfile("psig_nvert_sumEtl700", "<Significance> vs. N Vertices, sumEt < 700", 50, 0, 50, 0, 10);
   hists_["pchi2_sumEtm700"] = new TH1D("hpchi2_sumEtm700", "#chi^{2} Probability, sumEt > 700", 50, 0, 1);
   profs_["sig_nvert_sumEtm700"] = new TProfile("psig_nvert_sumEtm700", "<Significance> vs. N Vertices, sumEt > 700", 50, 0, 50, 0, 10);
   hists_["pchi2_sumEtm900"] = new TH1D("hpchi2_sumEtm900", "#chi^{2} Probability, sumEt > 900", 50, 0, 1);
   profs_["sig_nvert_sumEtm900"] = new TProfile("psig_nvert_sumEtm900", "<Significance> vs. N Vertices, sumEt > 900", 50, 0, 50, 0, 10);

   // event loop, fill histograms
   for(vector<event>::iterator ev = eventref.begin(); ev < eventref.end(); ev++){

      // significance
      hists_["sig"]->Fill( ev->sig );
      hists_["pchi2"]->Fill( TMath::Prob(ev->sig,2) );

      // pile up
      hists_["nvert"]->Fill( ev->nvertices );

      // significance vs. pile up
      hists2d_["sig_nvert"]->Fill( ev->nvertices, ev->sig );
      if( ev->sig < 9 ){
         profs_["sig_nvert"]->Fill( ev->nvertices, ev->sig );
      }

      // hadronic activity
      hists_["njets"]->Fill( ev->jet_pt.size() );
      hists_["sumEt"]->Fill( ev->pjet_scalpt );

      // bins of hadronic activity
      if( ev->jet_pt.size() == 0 ){
         hists_["pchi2_0jets"]->Fill( TMath::Prob(ev->sig,2) );
         if( ev->sig < 9 ){
            profs_["sig_nvert_0jets"]->Fill( ev->nvertices, ev->sig );
         }
      }else{
         hists_["pchi2_mt0jets"]->Fill( TMath::Prob(ev->sig,2) );
         if( ev->sig < 9 ){
            profs_["sig_nvert_mt0jets"]->Fill( ev->nvertices, ev->sig );
         }
      }
      if( ev->pjet_scalpt < 500 ){
         hists_["pchi2_sumEtl500"]->Fill( TMath::Prob(ev->sig,2) );
         if( ev->sig < 9 ){
            profs_["sig_nvert_sumEtl500"]->Fill( ev->nvertices, ev->sig );
         }
      }
      if( ev->pjet_scalpt < 700 ){
         hists_["pchi2_sumEtl700"]->Fill( TMath::Prob(ev->sig,2) );
         if( ev->sig < 9 ){
            profs_["sig_nvert_sumEtl700"]->Fill( ev->nvertices, ev->sig );
         }
      }
      if( ev->pjet_scalpt >= 700 ){
         hists_["pchi2_sumEtm700"]->Fill( TMath::Prob(ev->sig,2) );
         if( ev->sig < 9 ){
            profs_["sig_nvert_sumEtm700"]->Fill( ev->nvertices, ev->sig );
         }
      }
      if( ev->pjet_scalpt >= 900 ){
         hists_["pchi2_sumEtm900"]->Fill( TMath::Prob(ev->sig,2) );
         if( ev->sig < 9 ){
            profs_["sig_nvert_sumEtm900"]->Fill( ev->nvertices, ev->sig );
         }
      }


   }

   TFile *file = new TFile("metsigplots.root","RECREATE");
   file->cd();

   for(map<string,TH1D*>::const_iterator it = hists_.begin(); it != hists_.end(); it++){
      it->second->Write();
   }
   for(map<string,TH2D*>::const_iterator it = hists2d_.begin(); it != hists2d_.end(); it++){
      it->second->Write();
   }
   for(map<string,TProfile*>::const_iterator it = profs_.begin(); it != profs_.end(); it++){
      it->second->Write();
   }

   file->Write();

}
