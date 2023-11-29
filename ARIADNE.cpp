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

int ARIADNE(){

    RooRealVar drift_time ("drift_time", "drift_time", 460, 530);
    drift_time.setBins(70);

    RooRealVar mean ("mean", "mean", 485., 500.);
    RooRealVar sigma ("sigma", "sigma", 0., 3.);
    RooLandau landau ("landau", "landau", drift_time, mean, sigma);

    RooUniform unif ("unif", "unif", drift_time);

    RooRealVar fs("fs", "fraction of signal" , 0.6, -0.5, 1.5);

    RooAddPdf model("non_extended model", "non_extended model", RooArgList(landau, unif), RooArgList(fs));
    
    RooDataHist dh_drift_time("histo_drift", "histo_drift", drift_time); //non va con il puntatore

    ifstream myfile("ariadne_g006_plus_400.dat");
    Double_t val, weight;
    while(!myfile.eof()) {
        myfile >> val >> weight;
        drift_time.setVal(val);
        dh_drift_time.set(drift_time, weight);
    }

    model.fitTo(dh_drift_time);
    RooPlot *ari = drift_time.frame();
    dh_drift_time.plotOn(ari, Name("Data"));
    model.plotOn(ari, Name("Model Fit")); 

    TCanvas *c1 = new TCanvas("c1", "ARIADNE", 1600, 800);
    ari->Draw();
    c1->Print("ariadne.png");       

    return 0;
}