#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace RooFit;

int MINOS() {
    // Setting the observables 
    RooRealVar reco_en ("reco_en", "reco_en", 0.5, 14.);
    RooRealVar mc_noosc ("mc_noosc", "mc_noosc", 0.5, 14.);
    RooRealVar mixing ("mixing", "sin^2(2theta)", 0., 1.);
    RooRealVar dm2 ("dm2", "mass_splitting", 0., 1.); 
    RooRealVar distance ("distance", "distance", 730.); //constant
    
    // Computing the oscillation probability with observed data
    RooFormulaVar osc_prob ("osc_prob", "1 - @1 * (TMath::Sin(1.267 * @2 * @3 / @0))^2", RooArgList(reco_en, mixing, dm2, distance));

    // Retrieving measured and simulated MC data from file
    RooDataSet* ds_reco_en = RooDataSet::read("minos_2013_data.dat", reco_en, "reco_en");  //tutti puntatori perchè o tutti o nessuno
    RooDataSet* ds_mc_noosc = RooDataSet::read("minos_2013_mc.dat", mc_noosc, "mc_noosc"); 

    // Getting a copy of data_mc_noosc with only "mc_noosc" column and binning it
    RooDataSet* ds_mc_reduced = (RooDataSet *)ds_mc_noosc->reduce(RooArgSet(mc_noosc)); //questo giro è assolutamente inutile, ha solo una colonna di partenza
    //RooDataHist* dh_mc_noosc = ds_mc_reduced->binnedClone();
    RooDataHist* dh_mc_noosc = ds_mc_reduced->binnedClone();

    // Creating from MC data a histo-based function for non-osc neutrinos
    RooHistFunc func_noosc ("func_noosc", "No oscillation", mc_noosc, *dh_mc_noosc, 2); //2 è interpolation order, ovvero ordine polinomio che va fittato

    // Defining fitting function
    auto model = RooGenericPdf("model", "Oscillation probability * MC expected oscillation", "@0 * @1", RooArgList(osc_prob, func_noosc));

    // Fitting and plotting
    model.fitTo(*ds_reco_en);
    RooPlot *osc = reco_en.frame();
    ds_reco_en->plotOn(osc, Name("Data"));
    model.plotOn(osc, Name("Model Fit"));

    TCanvas *c1 = new TCanvas("c1", "B0_decay", 800, 400);
    osc->Draw();

    return 0;
}