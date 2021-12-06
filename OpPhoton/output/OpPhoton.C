void OpPhoton()
{
  gStyle->SetOptStat(0);
  // gStyle->SetOptStat(1001111);

  // For the binning
  int nbBins = 50;
  // in centimeters (for consistency better to keep 1 bin = 2 mm)
  double stickLength = 10.;

  // Reading the input file
  int nbFiles = 4;

  TFile *input[nbFiles];

  // input[0] = TFile::Open("OpPhoton_EJ200_900_5cm.root");
  // input[1] = TFile::Open("OpPhoton_EJ200_925_5cm.root");
  // input[2] = TFile::Open("OpPhoton_EJ200_950_5cm.root");
  // input[3] = TFile::Open("OpPhoton_EJ200_975_5cm.root");

  input[0] = TFile::Open("OpPhoton_EJ200_900_10cm.root");
  input[1] = TFile::Open("OpPhoton_EJ200_925_10cm.root");
  input[2] = TFile::Open("OpPhoton_EJ200_950_10cm.root");
  input[3] = TFile::Open("OpPhoton_EJ200_975_10cm.root");

  // input[0] = TFile::Open("OpPhoton_EJ200_900_20cm.root");
  // input[1] = TFile::Open("OpPhoton_EJ200_925_20cm.root");
  // input[2] = TFile::Open("OpPhoton_EJ200_950_20cm.root");
  // input[3] = TFile::Open("OpPhoton_EJ200_975_20cm.root");

  // input[0] = TFile::Open("OpPhoton_EJ230_900_10cm.root");
  // input[1] = TFile::Open("OpPhoton_EJ230_925_10cm.root");
  // input[2] = TFile::Open("OpPhoton_EJ230_950_10cm.root");
  // input[3] = TFile::Open("OpPhoton_EJ230_975_10cm.root");

  char imageCount[64] = "count_EJ230_5cm.png";
  char imageTime[64] = "time_EJ230_5cm.png";
    
  // Reading the data (in TTree format)
  TTree *OpPhoton[nbFiles];
  bool rowWise;
  TBranch* eventBranch[nbFiles];

  // Setting up the variables to be read from the data
  double ComptonE[nbFiles];
  double ComptonDepth[nbFiles];
  double nbFrontPhot[nbFiles];
  double timeFront20[nbFiles];
  double timeFront100[nbFiles];
  double nbBackPhot[nbFiles];
  double timeBack20[nbFiles];
  double timeBack100[nbFiles];

  int entries[nbFiles];

  TH1F *ComptonProfile[nbFiles];
  TH1F *ComptonProfile_noCut[nbFiles];
  TH1F *ElostProfile[nbFiles];
  TH1F *ComptonFront[nbFiles];
  TH1F *ComptonBack[nbFiles];
  TH1F *timeFront[nbFiles];
  TH1F *timeBack[nbFiles];

  // Loop over the files to extract their data
  for (int i = 0; i < nbFiles; i++)
  {
    OpPhoton[i] = (TTree *) input[i]->Get("OpPhoton;1");
    
    rowWise = true;
    eventBranch[i] = OpPhoton[i]->FindBranch("row_wise");
    if (!eventBranch[i]) rowWise = false;

    // Connect these variables to that of the TTree data
    if (!rowWise)
    {
      OpPhoton[i]->SetBranchAddress("ComptonE", &ComptonE[i]);
      OpPhoton[i]->SetBranchAddress("ComptonDepth", &ComptonDepth[i]);
      OpPhoton[i]->SetBranchAddress("nbFrontPhot", &nbFrontPhot[i]);
      OpPhoton[i]->SetBranchAddress("timeFront20", &timeFront20[i]);
      OpPhoton[i]->SetBranchAddress("timeFront100", &timeFront100[i]);
      OpPhoton[i]->SetBranchAddress("nbBackPhot", &nbBackPhot[i]);
      OpPhoton[i]->SetBranchAddress("timeBack20", &timeBack20[i]);
      OpPhoton[i]->SetBranchAddress("timeBack100", &timeBack100[i]);
    }

    else cout << "WARNING: row_wiseWise = true" << endl;

    entries[i] = OpPhoton[i]->GetEntries();
    // cout << "Numbre of entries = " << entries[i] << endl;

    // Histograms
    // Depth-dose profile
    char param0[64], param1[64], param2[64], param3[64], param4[64], param5[64],
         param6[64];
  
    sprintf(param0, "ComptonProfile%i", i);
    sprintf(param1, "ComptonProfile_noCut%i", i);
    sprintf(param2, "ElostProfile%i", i);
    sprintf(param3, "ComptonFront%i", i);
    sprintf(param5, "timeFront%i", i);
    sprintf(param4, "ComptonBack%i", i);
    sprintf(param6, "timeBack%i", i);

    ComptonProfile[i] = new TH1F(param0, "", nbBins, 0., stickLength);

    ComptonProfile_noCut[i] = new TH1F(param1, "", nbBins, 0., stickLength);


    ElostProfile[i] = new TH1F(param2,
                               "Elost vs depth - 1st Compton scattering",
                               nbBins, 0., stickLength);

    ComptonFront[i] = new TH1F(param3, "", nbBins, 0., stickLength);
    timeFront[i] = new TH1F(param4, "", nbBins, 0., stickLength);
    
    ComptonBack[i] = new TH1F(param5, "", nbBins, 0., stickLength);
    timeBack[i] = new TH1F(param6, "", nbBins, 0., stickLength);

    // int nb20 = 0;
    int nb100 = 0;

    // Loop over the entries (one entry = annihilation gamma)
    for (int j = 0; j < entries[i]; j++)
    {
      // Read the single entry inside the data file (Now all the variables are
      // linked to the values at this entry)
      OpPhoton[i]->GetEntry(j);

      // if (timeFront20[i] > 0. && timeBack20[i] > 0.)
      // {
      //   ComptonProfile[i]->Fill(ComptonDepth[i]);
      //   ElostProfile[i]->Fill(ComptonDepth[i], ComptonE[i]);
      //   ComptonFront[i]->Fill(ComptonDepth[i], nbFrontPhot[i]);
      //   ComptonBack[i]->Fill(ComptonDepth[i], nbBackPhot[i]);
      //   timeFront[i]->Fill(ComptonDepth[i], timeFront20[i]);
      //   timeBack[i]->Fill(ComptonDepth[i], timeBack20[i]);

      //   // nb20++;
      // }

      if (timeFront100[i] > 0. && timeBack100[i] > 0.)
      {
        ComptonProfile[i]->Fill(ComptonDepth[i]);
        ElostProfile[i]->Fill(ComptonDepth[i], ComptonE[i]);
        ComptonFront[i]->Fill(ComptonDepth[i], nbFrontPhot[i]);
        ComptonBack[i]->Fill(ComptonDepth[i], nbBackPhot[i]);
        timeFront[i]->Fill(ComptonDepth[i], timeFront100[i]);
        timeBack[i]->Fill(ComptonDepth[i], timeBack100[i]);

        nb100++;
      }

      // ComptonProfile_noCut[i]->Fill(ComptonDepth[i]);
    }

    cout << "Efficiency (%) = " << 100.*nb100/200000 << endl;
  }

  // Optical photon count analysis
  TCanvas *C1 = new TCanvas ("C1", "Histograms", 800, 500);
  C1->Divide(2, 2);

  TH1F *frontNorm[nbFiles];
  TH1F *backNorm[nbFiles];
  TH1F *sum[nbFiles];
  TH1F *ratio[nbFiles];

  C1->cd(1);
  frontNorm[0] = (TH1F*) ComptonFront[0]->Clone();
  frontNorm[0]->Divide(ElostProfile[0]);
  frontNorm[0]->GetYaxis()->SetRangeUser(0., 6.);
  frontNorm[0]->SetTitle(
    "Nb of opt photons vs Compton scattering depth");
  frontNorm[0]->GetXaxis()->SetTitle("Compton scattering depth (cm)");
  frontNorm[0]->GetYaxis()->SetTitle("Front end count / Compton keV");
  frontNorm[0]->SetLineColor(kRed);
  frontNorm[0]->Draw("HIST");

  for (int i = 1; i < nbFiles; i++)
  {
    frontNorm[i] = (TH1F*) ComptonFront[i]->Clone();
    frontNorm[i]->Divide(ElostProfile[i]);
    if (i == 1) frontNorm[i]->SetLineColor(kGreen);
    if (i == 2) frontNorm[i]->SetLineColor(kBlue);
    if (i == 3) frontNorm[i]->SetLineColor(kBlack);
    frontNorm[i]->Draw("HIST, SAME");
  }

  TLegend *legCount = new TLegend(0.5, 0.7, .9, 0.9);
  legCount->AddEntry(frontNorm[0], "R = 0.900", "l");
  if (nbFiles > 1) legCount->AddEntry(frontNorm[1], "R = 0.925", "l");
  if (nbFiles > 2) legCount->AddEntry(frontNorm[2], "R = 0.950", "l");
  if (nbFiles > 3) legCount->AddEntry(frontNorm[3], "R = 0.975", "l");
  legCount->Draw();

  C1->cd(2);
  backNorm[0] = (TH1F*) ComptonBack[0]->Clone();
  backNorm[0]->Divide(ElostProfile[0]);
  backNorm[0]->GetYaxis()->SetRangeUser(0., 6.);
  backNorm[0]->SetTitle(
    "Nb of opt photons reaching the back end vs Compton scattering depth");
  backNorm[0]->GetXaxis()->SetTitle("Compton scattering depth (cm)");
  backNorm[0]->GetYaxis()->SetTitle("Back end count / Compton keV");
  backNorm[0]->SetLineColor(kRed);
  backNorm[0]->Draw("HIST");

  for (int i = 1; i < nbFiles; i++)
  {
    backNorm[i] = (TH1F*) ComptonBack[i]->Clone();
    backNorm[i]->Divide(ElostProfile[i]);
    if (i == 1) backNorm[i]->SetLineColor(kGreen);
    if (i == 2) backNorm[i]->SetLineColor(kBlue);
    if (i == 3) backNorm[i]->SetLineColor(kBlack);
    backNorm[i]->Draw("HIST, SAME");
  }

  C1->cd(3);
  sum[0] = (TH1F*) frontNorm[0]->Clone();
  sum[0]->Add(backNorm[0]);
  sum[0]->GetYaxis()->SetRangeUser(0., 6.);
  sum[0]->SetTitle("Sum");
  sum[0]->GetXaxis()->SetTitle("Compton scattering depth (cm)");
  sum[0]->GetYaxis()->SetTitle("Sum front end + back end counts / Compton keV");
  sum[0]->SetLineColor(kRed);
  sum[0]->Draw("HIST");

  for (int i = 1; i < nbFiles; i++)
  {
    sum[i] = (TH1F*) frontNorm[i]->Clone();
    sum[i]->Add(backNorm[i]);
    if (i == 1) sum[i]->SetLineColor(kGreen);
    if (i == 2) sum[i]->SetLineColor(kBlue);
    if (i == 3) sum[i]->SetLineColor(kBlack);
    sum[i]->Draw("HIST, SAME");
  }

  C1->cd(4);
  ratio[0] = (TH1F*) frontNorm[0]->Clone();
  ratio[0]->Divide(sum[0]);
  ratio[0]->GetYaxis()->SetRangeUser(0., 1.);
  ratio[0]->SetTitle("Ratio");
  ratio[0]->GetXaxis()->SetTitle("Compton scattering depth (cm)");
  ratio[0]->GetYaxis()->SetTitle(
    "Front end count / Sum front end + back end counts");
  ratio[0]->SetLineColor(kRed);
  ratio[0]->Draw("HIST");

  for (int i = 1; i < nbFiles; i++)
  {
    ratio[i] = (TH1F*) frontNorm[i]->Clone();
    ratio[i]->Divide(sum[i]);
    if (i == 1) ratio[i]->SetLineColor(kGreen);
    if (i == 2) ratio[i]->SetLineColor(kBlue);
    if (i == 3) ratio[i]->SetLineColor(kBlack);
    ratio[i]->Draw("HIST, SAME");
  }

  C1->Print(imageCount);

  TCanvas *C2 = new TCanvas ("C2", "Histograms", 800, 500);
  C2->Divide(2, 2);

  TH1F *timeFrontNorm[nbFiles];
  TH1F *timeBackNorm[nbFiles];
  TH1F *deltaTime[nbFiles];

  C2->cd(1);
  timeFrontNorm[0] = (TH1F*) timeFront[0]->Clone();
  timeFrontNorm[0]->Divide(ComptonProfile[0]);
  timeFrontNorm[0]->GetYaxis()->SetRangeUser(0., 8.);
  timeFrontNorm[0]->GetXaxis()->SetTitle("Depth (cm)");
  timeFrontNorm[0]->GetYaxis()->SetTitle("Time front (ns)");
  timeFrontNorm[0]->SetLineColor(kRed);
  timeFrontNorm[0]->Draw("HIST");

  for (int i = 1; i < nbFiles; i++)
  {
    timeFrontNorm[i] = (TH1F*) timeFront[i]->Clone();
    timeFrontNorm[i]->Divide(ComptonProfile[i]);
    if (i == 1) timeFrontNorm[i]->SetLineColor(kGreen);
    if (i == 2) timeFrontNorm[i]->SetLineColor(kBlue);
    if (i == 3) timeFrontNorm[i]->SetLineColor(kBlack);
    timeFrontNorm[i]->Draw("HIST, SAME");
  }

  TLegend *legTime = new TLegend(0.1, 0.7, .3, 0.9);
  legTime->AddEntry(timeFrontNorm[0], "Reflectivity = 0.900", "l");
  if (nbFiles > 1) legTime->AddEntry(timeFrontNorm[1], "R = 0.925", "l");
  if (nbFiles > 2) legTime->AddEntry(timeFrontNorm[2], "R = 0.950", "l");
  if (nbFiles > 3) legTime->AddEntry(timeFrontNorm[3], "R = 0.975", "l");
  legTime->Draw();

  C2->cd(2);
  timeBackNorm[0] = (TH1F*) timeBack[0]->Clone();
  timeBackNorm[0]->Divide(ComptonProfile[0]);
  timeBackNorm[0]->GetYaxis()->SetRangeUser(0., 8.);
  timeBackNorm[0]->GetXaxis()->SetTitle("Depth (cm)");
  timeBackNorm[0]->GetYaxis()->SetTitle("Time back (ns)");
  timeBackNorm[0]->SetLineColor(kRed);
  timeBackNorm[0]->Draw("HIST");

  for (int i = 1; i < nbFiles; i++)
  {
    timeBackNorm[i] = (TH1F*) timeBack[i]->Clone();
    timeBackNorm[i]->Divide(ComptonProfile[i]);
    if (i == 1) timeBackNorm[i]->SetLineColor(kGreen);
    if (i == 2) timeBackNorm[i]->SetLineColor(kBlue);
    if (i == 3) timeBackNorm[i]->SetLineColor(kBlack);
    timeBackNorm[i]->Draw("HIST, SAME");
  }

  C2->cd(3);
  deltaTime[0] = (TH1F*) timeBackNorm[0]->Clone();
  deltaTime[0]->Add(timeFrontNorm[0], -1);
  deltaTime[0]->GetYaxis()->SetRangeUser(-8., 8.);
  deltaTime[0]->GetXaxis()->SetTitle("Depth (cm)");
  deltaTime[0]->GetYaxis()->SetTitle("Delta time front - time back (ns)");
  deltaTime[0]->SetLineColor(kRed);
  deltaTime[0]->Draw("HIST");

  for (int i = 1; i < nbFiles; i++)
  {
    deltaTime[i] = (TH1F*) timeBackNorm[i]->Clone();
    deltaTime[i]->Add(timeFrontNorm[i], -1);
    if (i == 1) deltaTime[i]->SetLineColor(kGreen);
    if (i == 2) deltaTime[i]->SetLineColor(kBlue);
    if (i == 3) deltaTime[i]->SetLineColor(kBlack);
    deltaTime[i]->Draw("HIST, SAME");
  }

  C2->Print(imageTime);
}