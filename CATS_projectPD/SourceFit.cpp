#include"SketchRealisticPot.h"
void SetStyle(bool graypalette) {
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  //gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if (graypalette)
    gStyle->SetPalette(8, 0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.065);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelFont(43, "xyz");
  gStyle->SetTitleFont(43, "xyz");
  gStyle->SetLabelSize(28, "xyz");
  gStyle->SetTitleSize(28, "xyz");
  gStyle->SetLabelOffset(0.0, "xy");
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetErrorX(0.5);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.5);
  gStyle->SetPalette(kCividis);
}
void EffectiveGaussianPd() {
   SetStyle(kTRUE);
  const double CoreSize = 0.906;  // we have from pp source paper ~0.94 fm

  //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
  DLM_CleverMcLevyResoTM MagicSource;

  //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
  MagicSource.InitStability(1, 2 - 1e-6, 2 + 1e-6);
  MagicSource.InitScale(38, 0.15, 2.0);
  MagicSource.InitRad(257 * 2, 0, 64);
  MagicSource.InitType(2);
  ///////////////////

  //for p-Xi, set up the amount of secondaries
  //first for the protons (64.22%)
  MagicSource.SetUpReso(0, 0.6422);
  //than for the Xis, here its 0% (we have ONLY primordials)
  MagicSource.SetUpReso(1, 0.0);

  //the cut off scale in k*, for which the angular distributions from EPOS
  //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
  const double k_CutOff = 200;

  //to be used for the NTuple later on
  Float_t k_D;
  Float_t fP1;
  Float_t fP2;
  Float_t fM1;
  Float_t fM2;
  Float_t Tau1;
  Float_t Tau2;
  Float_t AngleRcP1;
  Float_t AngleRcP2;
  Float_t AngleP1P2;
  //random generator dimi style. The input is incompatible with the ROOT random generator,
  //do not mix and match, do not ask me how I know this. Ask Bernie.
  //11 is the seed, you can change that to you favorite number
  DLM_Random RanGen(11);
  //dummies to save random shit
  double RanVal1;
  double RanVal2;
  double RanVal3;

  //open the magic file from dimi with the angular distributions.
  TFile *F_EposDisto_pReso_Xim = new TFile(
      "/home/sbhawani/Desktop/CATS_CF_Pd/RootFiles/eposdisto-preso-omega.root");
  //set up the ntuple, do not change anything unless told so by dimi
  TNtuple *T_EposDisto_pReso_Xim = (TNtuple*) F_EposDisto_pReso_Xim->Get(
      "InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xim->GetEntries();
  T_EposDisto_pReso_Xim->SetBranchAddress("k_D", &k_D);
  T_EposDisto_pReso_Xim->SetBranchAddress("P1", &fP1);
  T_EposDisto_pReso_Xim->SetBranchAddress("P2", &fP2);
  T_EposDisto_pReso_Xim->SetBranchAddress("M1", &fM1);
  T_EposDisto_pReso_Xim->SetBranchAddress("M2", &fM2);
  T_EposDisto_pReso_Xim->SetBranchAddress("Tau1", &Tau1);
  T_EposDisto_pReso_Xim->SetBranchAddress("Tau2", &Tau2);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP1", &AngleRcP1);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP2", &AngleRcP2);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleP1P2", &AngleP1P2);
  //iterate over the ntuple
  for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Xim; uEntry++) {
    //get each entry
    T_EposDisto_pReso_Xim->GetEntry(uEntry);
    //disregard the entry of you are outside the desired k*
    if (k_D > k_CutOff)
      continue;
    //overwrite the value for the lifetime. This is computed from the
    //stat. hadronization model (Vale) or thermal fist (Max)
    //this is the value for the secondary protons
    Tau1 = 1.65;
    //for primoridials (the Xis) we put 0
    Tau2 = 0;
    //put in the average mass of the resonances (again from SHM or TF)
    //this is the value for protons
    fM1 = 1362;
    //generate a random path length for the propagation of the resonances
    //nothing to change!
    RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
    //adds a single entry into the PDF for the angular distribution to be used
    MagicSource.AddBGT_RP(RanVal1, cos(AngleRcP1));

  }
  delete F_EposDisto_pReso_Xim;
  //if you have resonances contributing to both particles, we need to repeat the above procedure
  //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

  const unsigned NumSourceBins = 128;
  const double rMin = 0;
  const double rMax = 16;
  TFile *fOutput = new TFile("fOutput.root", "recreate");
  TH1F *hSource = new TH1F("hSource", "hSource", NumSourceBins, rMin, rMax);

  //fill the histo fro the source
  for (unsigned uBin = 0; uBin < NumSourceBins; uBin++) {
    //get the x-axis (r value) of the current bin
    double xaxis = hSource->GetBinCenter(uBin + 1);
    //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
    double parameters[2];
    parameters[0] = CoreSize;
    parameters[1] = 2.0;
    double SourceValue = MagicSource.RootEval(&xaxis, parameters);
    hSource->SetBinContent(uBin + 1, SourceValue);
    //infinite errors for now
    hSource->SetBinError(uBin + 1, 1000.);
  }

  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName =
      Form(
          "%s/PD_FullPTResonancePdFor2p0pT.pdf",
          "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD");  //plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");

  //idea: fit the source distribution only in a range around its peak
  //to do this: silly idea: put very large uncertainties in the bins outside of this range
  //we can get this range automatically, by evaluating the central (median) integral of the source distribution
  //with this set up, we fit the 68% most central yield of the source distribution
  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
  unsigned lowerbin = hSource->FindBin(lowerlimit);
  unsigned upperbin = hSource->FindBin(upperlimit);
  for (unsigned uBin = lowerbin; uBin <= upperbin; uBin++) {
    hSource->SetBinError(uBin + 1, 0.01);
  }

  printf("Core size of %.3f fm\n", CoreSize);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n", lowerlimit,
         upperlimit);
  //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
  TF1 *fSource = new TF1("fSource", GaussSourceTF1, rMin, rMax, 1);
  fSource->SetParameter(0, CoreSize);
  fSource->SetParLimits(0, CoreSize * 0.5, CoreSize * 2.0);
  hSource->Fit(fSource, "S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",
         fSource->GetParameter(0), fSource->GetParError(0));

  //get rid of weird plotting
  for (unsigned uBin = 0; uBin < NumSourceBins; uBin++) {
    hSource->SetBinError(uBin + 1, 0.01);
  }
  //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  hSource->SetTitle(0);
  //hSource->GetXaxis()->SetTitle("r (fm)");
 // hSource->GetXaxis()->SetTitleSize(0.045);
  //hSource->GetXaxis()->SetTitleOffset(1.35);
  //hSource->GetYaxis()->SetTitleSize(0.045);
  //hSource->GetYaxis()->SetTitleOffset(1.35);
  //hSource->GetYaxis()->SetTitle("S_{4#pi}(r)(fm^{-1})");
  hSource->SetTitle(";r (fm) ;S_{4#pi}(r)(fm^{-1})");
  //fSource->SetTitle(";r (fm) ;S_{4#pi}(r)(fm^{-1})");
  hSource->GetXaxis()->SetRangeUser(0, 10);
  hSource->GetYaxis()->SetRangeUser(-0.1, 0.5);
  //hSource->GetXaxis()->SetLabelSize(40);
 // hSource->GetYaxis()->SetLabelSize(40);
  //hSource->GetXaxis()->SetLabelOffset(.015);
  hSource->GetYaxis()->SetLabelOffset(.008);
  hSource->SetLineColor(kBlue);
  hSource->SetLineWidth(2);
  //gSource->MarkerStyle(20);
  //gSource->MarkerSize(1.5);
  TLegend *leg1 = new TLegend(0.43, 0.5, 0.7, 0.65);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.04);
  leg1->SetLineColor(0);
  leg1->SetNColumns(1);
  leg1->AddEntry(hSource, Form(" Source dist.(Reso. with r_{core}= %0.3f)",CoreSize));
  leg1->AddEntry(fSource, Form(" Gaussian Fit (%0.3f +/- %0.3f)",fSource->GetParameter(0),fSource->GetParError(0)));

  hSource->Draw(" ");
  fSource->Draw("same");
  leg1->Draw("same");
  printf("reached -3!");
  Plot->Print(pdfName);

  hSource->Write();
  fSource->Write();
  delete Plot;
  delete hSource;
  delete fSource;
  delete fOutput;
}
