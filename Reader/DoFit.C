#include "Reader.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
using namespace std;

int main(){

   // declarations
   Plotter plotter;
   vector<event> eventvec;

   //fitter.ReadNtuple( "ntuple_DYminiaod_PU20bx25_20150206.root", eventvec );
   //fitter.ReadNtuple( "ntuple_DYminiaod_PU4bx50_20150223.root", eventvec );
   //plotter.ReadNtuple( "../MakeNtuple/ntuples/ntuple_oct20.root", eventvec );
   //plotter.ReadNtuple( "../MakeNtuple/ntuple.root", eventvec );
   plotter.ReadNtuple( "../MakeNtuple/ntuples/ntuple_dec16_combined.root", eventvec );


   //double x [] = {1.15,1.08,1.04,1.13,1.56,0.0,0.55}; // 2012 data
   //double x [] = {1.42,1.29,1.41,1.40,2.52,0.0,0.673}; // 20bx25
   //double x [] = {1.39,1.31,1.36,1.34,2.27,-2.23,0.630}; // 4bx50
//   double x [] = {1.41,1.29,1.41,1.40,2.53,0.0,0.674}; // 20bx25
//   fitter.FindSignificance( x, eventvec );
   //fitter.RunMinimizer( eventvec );

//   fitter.MakePlots( eventvec );

   plotter.DeclareHists();
   plotter.Analyse( eventvec );
   return 0;
}
