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
#include "TMath.h"

using namespace RooFit;

int MINOS() {
    // Setting the observables 
    RooRealVar reco_en ("reco_en", "reco_en", 0.5, 14.);
    //RooRealVar mc_noosc ("mc_noosc", "mc_noosc", 0.5, 14.); // non serve per l'altro data set, anzi, devono essere la stessa variabile perchè altrimenti non fitta
    RooRealVar mixing ("mixing", "sin^{2}(2#theta)", 0.8, 1.);
    RooRealVar dm2 ("dm2", "(#Deltam^{2})", 2.0e-3, 2.6e-3);  
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
    osc->SetTitle("#nu_{#mu} oscillations in MINOS");
    osc->GetXaxis()->SetTitle("#nu_{#mu} energy (GeV)");
    ds_reco_en->plotOn(osc, Name("Data"));
    model.plotOn(osc, Name("Model Fit"));

    TCanvas *c1 = new TCanvas("c1", "MINOS", 1600, 800);
    osc->Draw();

    TLegend *leg1 = new TLegend(0.70, 0.70, 0.85, 0.87);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kBlack);
    leg1->AddEntry("Data", "Data", "P");
    leg1->AddEntry("Model Fit", "Fit function", "LP");
   
    leg1->Draw("SAME");

    c1->Print("minos_data.png");

    // Defining the negative log likelihood and the MINUIT interface object
    RooAbsReal* nll = model.RooAbsPdf::createNLL(*ds_reco_en);
    RooMinuit m(*nll);

    // Creating txt file to host fit results
    ofstream myfile;
    myfile.open("minos_fit.txt");

    // Minimizing via MINUIT with MIGRAD and printing the results to txt file
    m.setVerbose(kTRUE);
    m.migrad();
    dm2.Print();
    mixing.Print();

    //RooFitResult * result = m.save("m", "MIGRAD"); // forse lui lo vuole solo alla fine, quindi su file di testo scrivo solo gli altri risultati
    //result->Print("v"); // con v stampa anche distance e le global corrections 

    myfile << "MIGRAD results: \n(Delta m)^2 = " << dm2 << " +/- " << dm2.getError() << " eV^2" << '\n'; 
//    myfile << "\t \t \t (sin(2theta))^2 = " << mixing << " +/- " << mixing.getError() << '\n'; 
    myfile << "(sin(2theta))^2 = " << mixing << " +/- " << mixing.getError() << '\n'; 
    //result->printArgs(myfile); 
    myfile << "\n*****************************\n"; 

    RooPlot* contour = m.contour(mixing, dm2, 1, 2, 3); // up to 3 sigmas
    contour->SetTitle("Contour plot");
    contour->GetXaxis()->SetTitle("sin^{2}(2#theta)");
    contour->GetYaxis()->SetTitle("#Deltam^{2}");
    TCanvas *c2 = new TCanvas("c2", "MINOS_contour", 1600, 800);
    contour->Draw();    
    c2->Print("minos_likelihood.png");

    // Running HESSE for the parameters and printing the results to txt file
    m.setVerbose(kFALSE);
    m.hesse();               
    dm2.Print();
    mixing.Print();

    myfile << "HESSE results: \n(Delta m)^2 = " << dm2 << " +/- " << dm2.getError() << " eV^2" << '\n'; 
//    myfile << "\t \t \t (sin(2theta))^2 = " << mixing << " +/- " << mixing.getError() << '\n';
    myfile << "(sin(2theta))^2 = " << mixing << " +/- " << mixing.getError() << '\n';
    myfile << "\n*****************************\n"; 


    // result = m.save("m1", "HESSE");
    // result->Print("v"); 
    // myfile << "HESSE results: " << *result << '\n'; 

    // Running MINOS only on mass splitting parameter and printing the result to txt file
    m.minos(dm2);
    dm2.Print();
    myfile << "MINOS results: \n(Delta m)^2 = " << dm2 << " +/- " << dm2.getError() << " eV^2" << '\n';     
    //myfile << "\t \t \t";

    // Saving all results from MINOS fit on txt file
    RooFitResult * result = m.save("m", "MINOS"); 
    result = m.save("m2", "MINOS");
    result->Print("v"); 
    //myfile << "MINOS results: " << *result << '\n'; 

//  result->floatParsFinal().printMultiline(myfile, 1111, kTRUE);
    result->printMultiline(myfile, 1111, kTRUE);
    result->printValue(myfile);
    
    myfile.close();

  
    return 0;
}