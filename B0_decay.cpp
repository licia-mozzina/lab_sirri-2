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

int B0_decay() {
  // first declare the observable
  RooRealVar x("x", "m(p-pbar)", 5090., 5590.);

  // background model
  RooRealVar c("c", "exponential par", -0.001, -1.,
               0.); // mancavano i punti decimali, era confuso
  RooExponential bkg("bkg", "Comb. bkg.", x, c);

  // first gaussian
  RooRealVar m0("m0", "B0 mass", 5279., 5220., 5320.);
  RooRealVar s0("s0", "B0 width", 10., 2., 50.);
  RooGaussian gaus0("gaus0", "B0 peak", x, m0, s0);

  // second gaussian
  RooRealVar m1("m1", "B0s mass", 5380., 5320., 5420.);
  RooRealVar s1("s1", "B0s width", 10., 2., 50.);
  RooGaussian gaus1("gaus1", "B0s peak", x, m1, s1);

  /*
  // composite model
  // not extended model  3 p.d.s. => 2 fractions
  //   f0 * gaus0(x) + f1 * gaus1(x) + (1 -f0 -f1)* bkg(x)

  RooRealVar f0("f0", "fraction of B0" , 0.2, 0.00001, 1);
  RooRealVar f1("f1", "fraction of B0s", 0.2, 0.00001, 1);

  RooAddPdf model("model", "model",
     RooArgList(gaus0, gaus1, bkg), RooArgList(f0, f1) );
  */

  // composite model
  // extended model  3 p.d.s. => 3 expected n. of events
  //   N0 * gaus0(x) + N1 * gaus1(x) + Nb* bkg(x)

  RooRealVar n0("n0", "number of B0 events", 100, 0, 10000);
  RooRealVar n1("n1", "number of B0s events", 10, 0, 10000);
  RooRealVar nb("nb", "number of background events", 1000, 0, 100000);

  RooAddPdf model("extended model", "extended model", RooArgList(gaus0, gaus1, bkg), RooArgList(n0, n1, nb));

  // read data
  // RooDataHist for binned dataset, RooDataSet for unbinned

  RooDataSet *data = RooDataSet::read("rarest_b0_decay.dat", x, "x");

  // RooDataHist data = *RooDataHist::read("rarest_b0_decay.dat", x, "x");

  // fit the model to the data

  model.fitTo(*data);
  
  // draw data and model
  RooPlot *xframe = x.frame();

  data->plotOn(xframe, Name("Data"));
  model.plotOn(xframe, Components(gaus0), LineColor(kViolet+1), Name("gaus0"));
  model.plotOn(xframe, Components(gaus1), LineColor(kViolet-1), Name("gaus1"));
  model.plotOn(xframe, Components(bkg), LineColor(kMagenta), Name("exp_bkg"));
  model.plotOn(xframe, Components(gaus0, gaus1, bkg), LineColor(kBlue), Name("fit"));

  // Creating histograms from residual and pulls distributions wrt the model
  RooHist *hresid = xframe->residHist();
  RooHist *hpull = xframe->pullHist();

  // Frames and plots of residual and pulls distributions
  RooPlot *res = x.frame(Title("Residual Distribution"));
  res->addPlotable(hresid, "P");
  RooPlot *pull = x.frame(Title("Pulls Distribution"));
  pull->addPlotable(hpull, "P");

  TCanvas *c1 = new TCanvas("c1", "B0_decay", 800, 400);
  c1->Divide(3);
  
  c1->cd(1);
  xframe->Draw();

  TLegend *leg1 = new TLegend(0.65, 0.73, 0.86, 0.87);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kBlack);
  leg1->AddEntry("Data", "Data", "P");
  leg1->AddEntry("gaus0", "First Gaussian Peak", "LP");
  leg1->AddEntry("gaus1", "Second Gaussian Peak", "LP");
  leg1->AddEntry("exp_bkg", "Exponential Background", "LP");
  leg1->AddEntry("fit", "Composite Fit", "LP");

  leg1->Draw("SAME");

  c1->cd(2);
  res->Draw();

  c1->cd(3);
  pull->Draw();


  return 0;
}