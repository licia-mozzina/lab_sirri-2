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
    //RooRealVar mc_noosc ("mc_noosc", "mc_noosc", 0.5, 14.); // non serve per l'altro data set, anzi, devono essere la stessa variabile perchè altrimenti non fitta
    RooRealVar mixing ("mixing", "sin^{2}(2#theta)", 0.8, 1.);
    RooRealVar dm2 ("dm2", "(#Deltam^{2})", 2.2e-3, 2.6e-3);  
    RooRealVar distance ("distance", "distance", 730); //constant
    
    // Computing the oscillation probability with observed data
    RooFormulaVar osc_prob ("osc_prob", "1 - @1 * (TMath::Sin(1.267 * @2 * @3 / @0))^2", RooArgList(reco_en, mixing, dm2, distance));

    // Retrieving measured and simulated MC data from file
    RooDataSet* ds_reco_en = RooDataSet::read("minos_2013_data.dat", reco_en, "v");  //tutti puntatori perchè o tutti o nessuno
    RooDataSet* ds_mc_noosc = RooDataSet::read("minos_2013_mc.dat", reco_en, "v"); 

    // Getting a copy of data_mc_noosc with only "mc_noosc" column and binning it
    RooDataSet* ds_mc_reduced = (RooDataSet *)ds_mc_noosc->reduce(RooArgSet(reco_en)); //questo giro è assolutamente inutile, ha solo una colonna di partenza
    RooDataHist* dh_mc_noosc = ds_mc_reduced->binnedClone();

    // Creating from MC data a histo-based function for non-osc neutrinos
    RooHistFunc func_noosc ("func_noosc", "No oscillation", reco_en, *dh_mc_noosc, 2); //2 è interpolation order, ovvero ordine polinomio che va fittato

    // Defining fitting function
    auto model = RooGenericPdf("model", "Oscillation probability * MC expected oscillation", "@0 * @1", RooArgList(osc_prob, func_noosc));

    // Fitting and plotting
    model.fitTo(*ds_reco_en);
    RooPlot *osc = reco_en.frame();
    ds_reco_en->plotOn(osc, Name("Data"));
    model.plotOn(osc, Name("Model Fit"));

    TCanvas *c1 = new TCanvas("c1", "MINOS", 1600, 800);
    osc->Draw();
    c1->Print("minos_data.png");

    // Defining the negative log likelihood and the MINUIT interface object
    RooAbsReal* nll = model.RooAbsPdf::createNLL(*ds_reco_en);
    RooMinuit m(*nll);

    // Minimizing via MINUIT with MIGRAD and printing the results
    m.setVerbose(kTRUE);
    m.migrad();
    dm2.Print();
    mixing.Print();

    RooFitResult * result = m.save("m", "MIGRAD");
    result->Print("v"); // con v stampa anche distance e le global corrections 

    ofstream myfile;
    myfile.open("minos_fit.txt");
    myfile << "MIGRAD results: " << *result << '\n'; 

    RooPlot* contour = m.contour(mixing, dm2, 1, 2, 0); // up to 3 sigmas
    TCanvas *c2 = new TCanvas("c2", "MINOS_contour", 1600, 800);
    contour->Draw();    
    c2->Print("minos_likelihood.png");

    // Running HESSE for the parameters and printing the results
    m.setVerbose(kFALSE);
    m.hesse();               
    dm2.Print();
    mixing.Print();

    result = m.save("m1", "HESSE");
    result->Print("v"); 
    myfile << "HESSE results: " << *result << '\n'; 

    // Running MINOS only on mass splitting parameter, printing the result
    m.minos(dm2);
    dm2.Print();

    result = m.save("m2", "MINOS");
    result->Print("v"); 
    myfile << "MINOS results: " << *result << '\n'; 

    myfile.close();

  
    return 0;
}