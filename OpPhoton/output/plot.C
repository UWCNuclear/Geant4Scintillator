#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

#include <math.h>
#include <numeric>
void plot()
{
  gStyle->SetOptStat(00000000);

  TCanvas *c1 = new TCanvas ("c1", "Histograms", 800, 500);
  c1->Divide(1, 1);
  
  c1->cd(1);
  c1->SetGrid();
  
  double x[3] = {5., 10., 20.};
  double y900[3] = {20.5, 6.1, 0.};
  double y925[3] = {24.4, 18.1, 0.};
  double y950[3] = {27.9, 33.4, 8.};
  double y975[3] = {31.1, 45.1, 45.2};

  TGraph* gr900 = new TGraph(3, x, y900);
  TGraph* gr925 = new TGraph(3, x, y925);
  TGraph* gr950 = new TGraph(3, x, y950);
  TGraph* gr975 = new TGraph(3, x, y975);

  TMultiGraph *mg  = new TMultiGraph();

  mg->Add(gr900);
  gr900->SetLineColor(kRed);

  mg->GetXaxis()->SetRangeUser(5., 20.);
  mg->GetYaxis()->SetRangeUser(0., 55.);
  mg->GetXaxis()->SetTitle("Stick length (cm)");
  mg->GetYaxis()->SetTitle("Efficiency (%)");

  mg->GetXaxis()->SetNdivisions(4);

  mg->Add(gr925);
  gr925->SetLineColor(kGreen);
  mg->Add(gr950);
  gr950->SetLineColor(kBlue);
  mg->Add(gr975);
  gr975->SetLineColor(kBlack);
  mg->Draw("APL");

  TLegend *legCount = new TLegend(0.1, 0.75, .3, .9);
  legCount->AddEntry(gr900, "Reflectivity = 0.900", "l");
  legCount->AddEntry(gr925, "Reflectivity = 0.925", "l");
  legCount->AddEntry(gr950, "Reflectivity = 0.950", "l");
  legCount->AddEntry(gr975, "Reflectivity = 0.975", "l");
  legCount->Draw();

  c1->Print("efficiency.png");
}