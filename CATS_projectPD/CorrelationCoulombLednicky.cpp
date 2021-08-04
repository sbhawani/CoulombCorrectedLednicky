//#ifdef __CINT__
#include "TGraph.h"
#include "TH1F.h"
#include "TGraph2D.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
//#endif
#include "CorrelationCoulombLednicky.h"

void Plot_Source_Pd(double rad) {
  printf("this is working!");
  const unsigned NumRadBins = 100;
  const double MinRad = 0;
  const double MaxRad = 10;

  TGraph *gSource = new TGraph();
  gSource->SetName("gSource");
  unsigned COUNTER = 0;
  for (double RAD = 0.1; RAD < MaxRad; RAD += 0.1) {
    gSource->SetPoint(COUNTER, RAD, SOURCE_4PI(RAD, rad));
    // printf("%f  %f", RAD, SOURCE_4PI(RAD, rad));
    COUNTER++;
  }
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName =
      Form(
          "%sSourcePD.png",
          "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD");  //plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");
  //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  gSource->GetXaxis()->SetTitle("r (fm)");
  gSource->GetXaxis()->SetTitleSize(0.045);
  gSource->GetXaxis()->SetTitleOffset(0.75);
  gSource->GetYaxis()->SetTitleSize(0.045);
  gSource->GetYaxis()->SetTitleOffset(0.95);
  gSource->GetYaxis()->SetTitle("S_{4#pi}(r)(fm^{-1})");

  gSource->GetXaxis()->SetRangeUser(0, 10);
  gSource->GetYaxis()->SetRangeUser(0, 0.5);

  gSource->SetLineColor(kRed);
  gSource->SetLineWidth(2);
  //gSource->MarkerStyle(20);
  //gSource->MarkerSize(1.5);
  TLegend *leg1 = new TLegend(0.5, 0.5, 0.7, 0.65);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.04);
  leg1->SetLineColor(0);
  leg1->SetNColumns(1);
  leg1->AddEntry(gSource, "  Gaussian Source");

  gSource->Draw("ACP");
  leg1->Draw("same");
  printf("reached -3!");
  Plot->Print(pdfName);
  delete Plot;
  delete gSource;
}

void CHECK_FOR_FUNCTIONS() {
  printf(
      "Plotting for Ac(), chi(), B(),P(),G() for pi pi tp check calculation is consistent!");
  const unsigned NumRadBins = 100;
  const double MinRad = 0;
  const double Maxk = 30;

  TGraph *gAc_repul = new TGraph();
  gAc_repul->SetName("gAc_repul");

  TGraph *gAc_att = new TGraph();
  gAc_att->SetName("gAc_att");

  TGraph *ReChi = new TGraph();
  ReChi->SetName("ReChi");

  TGraph *ImChi_repul = new TGraph();
  ImChi_repul->SetName("ImChi_repul");

  TGraph *ImChi_att = new TGraph();
  ImChi_att->SetName("ImChi");

  //TGraph *gB = new TGraph();
  // gB->SetName("gB");

  TGraph *gB1 = new TGraph();
  gB1->SetName("gB1");
  TGraph *gB2 = new TGraph();
  gB2->SetName("gB2");
  TGraph *gB3 = new TGraph();
  gB3->SetName("gB3");

  TGraph *gB12 = new TGraph();
  gB12->SetName("gB12");
  TGraph *gB22 = new TGraph();
  gB22->SetName("gB22");
  TGraph *gB32 = new TGraph();
  gB32->SetName("gB32");

  TGraph *gP1 = new TGraph();
  gP1->SetName("gP1");
  TGraph *gP2 = new TGraph();
  gP2->SetName("gP2");
  TGraph *gP3 = new TGraph();
  gP3->SetName("gP3");

  TGraph *gP12 = new TGraph();
  gP12->SetName("gP12");
  TGraph *gP22 = new TGraph();
  gP22->SetName("gP22");
  TGraph *gP32 = new TGraph();
  gP32->SetName("gP32");

  TGraph *gReG1 = new TGraph();
  gReG1->SetName("gReG1");
  TGraph *gReG2 = new TGraph();
  gReG2->SetName("gReG2");
  TGraph *gReG3 = new TGraph();
  gReG3->SetName("gReG3");

  TGraph *gImG1 = new TGraph();
  gImG1->SetName("gImG1");
  TGraph *gImG2 = new TGraph();
  gImG2->SetName("gImG2");
  TGraph *gImG3 = new TGraph();
  gImG3->SetName("gImG3");

  unsigned COUNTER = 0;
  for (double Inveta = 2; Inveta < Maxk; Inveta += 1) {
    gAc_repul->SetPoint(COUNTER, Inveta, Ac(1 / Inveta));
    gAc_att->SetPoint(COUNTER, Inveta, Ac(-1 / Inveta));
    ReChi->SetPoint(COUNTER, Inveta, h(1 / Inveta));
    ImChi_repul->SetPoint(COUNTER, Inveta, Ac(1 / Inveta) / 2 * Inveta);
    ImChi_att->SetPoint(COUNTER, Inveta, -Ac(-1 / Inveta) / 2 * Inveta);
    //printf("%f  %f\n", k, GET_CORRELATION(k, ScatLend, EffecRange, ChargeRad, SourceRad));
    COUNTER++;
  }

  unsigned counter = 0;
  for (double Inveta = 0.001; Inveta < Maxk; Inveta += 1) {
    ReChi->SetPoint(counter, Inveta, h(1 / Inveta));
    ImChi_repul->SetPoint(counter, Inveta, Ac(1 / Inveta) / 2 * Inveta);
    ImChi_att->SetPoint(counter, Inveta, -Ac(-1 / Inveta) / 2 * Inveta);
    //printf("%f  %f\n", k, GET_CORRELATION(k, ScatLend, EffecRange, ChargeRad, SourceRad));
    counter++;
  }

  unsigned counter2 = 0;
  const double rval1 = 5 * FmToNu;  // fm to nu
  const double rval2 = 15 * FmToNu;  // fm to nu
  const double rval3 = 50 * FmToNu;  // fm to nu
  const double ac = -387.5 * FmToNu;  // fm to nu charged radius

  //double B(double rho, double eta);  //recursive sum on bn
  // double P(double rho, double eta);//recursive sum on pn
  //TILDE_G(double rho, double eta);//

  for (double Q = 0.01; Q < 50; Q += 1) {
    gB1->SetPoint(counter2, Q, B(Q / 2 * rval1, 2 / Q / ac));
    gB2->SetPoint(counter2, Q, B(Q / 2 * rval2, 2 / Q / ac));
    gB3->SetPoint(counter2, Q, B(Q / 2 * rval3, 2 / Q / ac));

    gB12->SetPoint(counter2, Q, sin(Q / 2 * rval1) / (Q / 2 * rval1));
    gB22->SetPoint(counter2, Q, sin(Q / 2 * rval2) / (Q / 2 * rval2));
    gB32->SetPoint(counter2, Q, sin(Q / 2 * rval3) / (Q / 2 * rval3));

    gP1->SetPoint(counter2, Q, P(Q / 2 * rval1, 2 / Q / ac));
    gP2->SetPoint(counter2, Q, P(Q / 2 * rval2, 2 / Q / ac));
    gP3->SetPoint(counter2, Q, P(Q / 2 * rval3, 2 / Q / ac));

    gP12->SetPoint(counter2, Q, cos(Q / 2 * rval1));
    gP22->SetPoint(counter2, Q, cos(Q / 2 * rval2));
    gP32->SetPoint(counter2, Q, cos(Q / 2 * rval3));

    gReG1->SetPoint(counter2, Q, real(TILDE_G(Q / 2 * rval1, 2 / Q / ac)));
    gReG2->SetPoint(counter2, Q, real(TILDE_G(Q / 2 * rval2, 2 / Q / ac)));
    gReG3->SetPoint(counter2, Q, real(TILDE_G(Q / 2 * rval3, 2 / Q / ac)));

    gImG1->SetPoint(counter2, Q, imag(TILDE_G(Q / 2 * rval1, 2 / Q / ac)));
    gImG2->SetPoint(counter2, Q, imag(TILDE_G(Q / 2 * rval2, 2 / Q / ac)));
    gImG3->SetPoint(counter2, Q, imag(TILDE_G(Q / 2 * rval3, 2 / Q / ac)));

    counter2++;
  }

  TFile *OutputFile =
      new TFile(
          TString::Format(
              "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD/CHECK_FOR_FUNCTION_2Approx.root"),
          "recreate");
  printf("File Created\n");
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName =
      Form(
          "%sCkPD.png",
          "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD");  //plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");
  //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  /*gCk->GetXaxis()->SetTitle("k (MeV/c)");
   gCk->GetXaxis()->SetTitleSize(0.045);
   gCk->GetXaxis()->SetTitleOffset(0.75);
   gCk->GetYaxis()->SetTitleSize(0.045);
   gCk->GetYaxis()->SetTitleOffset(0.95);
   gCk->GetYaxis()->SetTitle("#it{C}(k)");

   gCk->GetXaxis()->SetRangeUser(0, 200);
   gCk->GetYaxis()->SetRangeUser(0, 10);

   gCk->SetLineColor(kRed);
   gCk->SetLineWidth(2);
   //gSource->MarkerStyle(20);
   //gSource->MarkerSize(1.5);
   TLegend *leg1 = new TLegend(0.5, 0.5, 0.7, 0.65);
   leg1->SetFillStyle(0);
   leg1->SetTextSize(0.04);
   leg1->SetLineColor(0);
   leg1->SetNColumns(1);
   leg1->AddEntry(gCk, " Ck pd Lednicky Coulomb");

   gCk->Draw("ACP");
   leg1->Draw("same");
   printf("reached -3!");*/
  gAc_repul->Write();
  gAc_att->Write();
  ReChi->Write();
  ImChi_repul->Write();
  ImChi_att->Write();

  //gB->Write();
  gB1->Write();
  gB2->Write();
  gB3->Write();
  gB12->Write();
  gB22->Write();
  gB32->Write();

  //gP->Write();

  gP1->Write();
  gP2->Write();
  gP3->Write();
  gP12->Write();
  gP22->Write();
  gP32->Write();

  gReG1->Write();
  gReG2->Write();
  gReG3->Write();
  gImG1->Write();
  gImG2->Write();
  gImG3->Write();

  Plot->Print(pdfName);
  delete Plot;
  delete gAc_repul;
  delete gAc_att;
  delete ReChi;
  delete ImChi_repul;
  delete ImChi_att;

  delete gB1;
  delete gB2;
  delete gB3;
  delete gB12;
  delete gB22;
  delete gB32;

  delete gP1;
  delete gP2;
  delete gP3;
  delete gP12;
  delete gP22;
  delete gP32;

  delete gReG1;
  delete gReG2;
  delete gReG3;
  delete gImG1;
  delete gImG2;
  delete gImG3;

  delete OutputFile;
}

void Plot_CoulombLed_WF_ForPd(double k, double t, double ScatLend,
                              double EffecRange, double ChargeRad) {

  printf("Plotting WF from Coulomb Lednicky!");
  const double MaxRad = 100;

  TGraph *gReWF = new TGraph();
  TGraph *gImWF = new TGraph();
  gReWF->SetName("gReWF");
  gImWF->SetName("gImWF");
  unsigned COUNTER = 0;
  for (double RAD = 0.001; RAD < MaxRad; RAD += 1) {
    gReWF->SetPoint(
        COUNTER, RAD,
        real(LEDNICKY_COULOMB_WF(k, RAD, t, ScatLend, EffecRange, ChargeRad)));
    gImWF->SetPoint(
        COUNTER, RAD,
        imag(LEDNICKY_COULOMB_WF(k, RAD, t, ScatLend, EffecRange, ChargeRad)));
    printf(
        "%f  %f  %f\n", RAD,
        real(LEDNICKY_COULOMB_WF(k, RAD, t, ScatLend, EffecRange, ChargeRad)),
        imag(LEDNICKY_COULOMB_WF(k, RAD, t, ScatLend, EffecRange, ChargeRad)));
    COUNTER++;
  }
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName =
      Form(
          "%sWF_k=.png",
          "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD");  //plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");
  //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  gReWF->GetXaxis()->SetTitle("r (fm)");
  gReWF->GetXaxis()->SetTitleSize(0.045);
  gReWF->GetXaxis()->SetTitleOffset(0.75);
  gReWF->GetYaxis()->SetTitleSize(0.045);
  gReWF->GetYaxis()->SetTitleOffset(0.95);
  gReWF->GetYaxis()->SetTitle("#psi(r)");

  gReWF->GetXaxis()->SetRangeUser(0, 100);
  gReWF->GetYaxis()->SetRangeUser(-10, 10);

  gImWF->SetLineColor(kBlue);
  gImWF->SetLineWidth(2);

  gReWF->SetLineColor(kRed);
  gReWF->SetLineWidth(2);
  //gSource->MarkerStyle(20);
  //gSource->MarkerSize(1.5);
  TLegend *leg1 = new TLegend(0.65, 0.6, 0.75, 0.75);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.06);
  leg1->SetLineColor(0);
  leg1->SetNColumns(1);

  leg1->AddEntry(gReWF, "#Rgothic [#psi]", "pef");
  leg1->AddEntry(gImWF, "#Jgothic [#psi]", "pef");

  gReWF->Draw("ACP");
  gImWF->Draw("ACPsame");
  leg1->Draw("same");
  printf("reached -2!");
  Plot->Print(pdfName);
  delete Plot;
  delete gReWF;
  delete gImWF;
}

void Plot_CoulombLed_Density_ForPd(double k, double t, double ScatLend,
                                   double EffecRange, double ChargeRad,
                                   bool draw2D) {

  printf("Plotting WF from Coulomb Lednicky!");
  const double MaxRad = 30;
  const double kMax = 300;
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName =
      Form(
          "%sDensityWF_k=.png",
          "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD");  //plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");

  if (draw2D) {

    TGraph2D *g2DDensity = new TGraph2D();
    g2DDensity->SetName("gDensity");
    unsigned COUNTER = 0;
    for (double kval = 0.0001; kval < kMax; kval += 2) {
      for (double RAD = 0.0001; RAD < MaxRad; RAD += 2) {
        g2DDensity->SetPoint(
            COUNTER,
            RAD,
            kval,
            PROB_LEDNICKY_COULOMB_WF(kval, RAD, t, ScatLend, EffecRange,
                                     ChargeRad));
        printf(
            "%f  %f\n",
            RAD,
            PROB_LEDNICKY_COULOMB_WF(kval, RAD, t, ScatLend, EffecRange,
                                     ChargeRad));
        COUNTER++;
      }
    }

    //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
    g2DDensity->GetXaxis()->SetTitle("r(fm)");
    g2DDensity->GetXaxis()->SetTitleSize(0.045);
    g2DDensity->GetXaxis()->SetTitleOffset(1.0);
    g2DDensity->GetYaxis()->SetTitleSize(0.045);
    g2DDensity->GetYaxis()->SetTitleOffset(1.2);
    //g2DDensity->GetYaxis()->SetTitle("k*(MeV/c)");
    g2DDensity->GetZaxis()->SetTitle("|#psi_{k}^{2}(k)|)");
    g2DDensity->SetTitle("|#psi_{k}^{2}(k)|); r(fm); k*(MeV/c)");
    g2DDensity->GetXaxis()->SetRangeUser(0, 80);
    g2DDensity->GetYaxis()->SetRangeUser(0, 200);
    //g2DDensity->SetLineColor(kRed);
    g2DDensity->SetLineWidth(1);
    //g2DDensity->MarkerStyle(20);
    //gSource->MarkerSize(1.5);
    TLegend *leg1 = new TLegend(0.65, 0.6, 0.75, 0.75);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);
    leg1->SetLineColor(0);
    leg1->SetNColumns(1);

    leg1->AddEntry(g2DDensity, "|#psi_{k}^{2}(r)|", "pef");

    g2DDensity->Draw("surf1");
    //leg1->Draw("same");
    printf("reached -2!");
    Plot->Print(pdfName);
    delete Plot;
    delete g2DDensity;

  } else {

    TGraph *gDensity = new TGraph();
    gDensity->SetName("gDensity");
    unsigned COUNTER = 0;
    for (double RAD = 0.001; RAD < MaxRad; RAD += 1) {
      gDensity->SetPoint(
          COUNTER, RAD,
          PROB_LEDNICKY_COULOMB_WF(k, RAD, t, ScatLend, EffecRange, ChargeRad));
      // printf(
      //    "%f  %f\n", RAD,
      //    PROB_LEDNICKY_COULOMB_WF(k, RAD, t, ScatLend, EffecRange, ChargeRad));
      COUNTER++;
    }

    //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
    gDensity->GetXaxis()->SetTitle("r(fm)");
    gDensity->GetXaxis()->SetTitleSize(0.045);
    gDensity->GetXaxis()->SetTitleOffset(0.75);
    gDensity->GetYaxis()->SetTitleSize(0.045);
    gDensity->GetYaxis()->SetTitleOffset(0.95);
    gDensity->GetYaxis()->SetTitle("|#psi_{k}^{2}(k)|)");

    gDensity->GetXaxis()->SetRangeUser(0, 100);
    gDensity->GetYaxis()->SetRangeUser(-10, 10);
    gDensity->SetLineColor(kRed);
    gDensity->SetLineWidth(2);
    //gSource->MarkerStyle(20);
    //gSource->MarkerSize(1.5);
    TLegend *leg1 = new TLegend(0.65, 0.6, 0.75, 0.75);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);
    leg1->SetLineColor(0);
    leg1->SetNColumns(1);

    leg1->AddEntry(gDensity, "|#psi_{k}^{2}(r)|", "pef");

    gDensity->Draw("ACP");
    leg1->Draw("same");
    printf("reached -2!");
    Plot->Print(pdfName);
    delete Plot;
    delete gDensity;
  }
}

void Plot_CoulombLed_Correlation_ForPd(double ScatLend1, double EffecRange1,
                                       double ScatLend2, double EffecRange2,
                                       double ChargeRad, bool SourceRad) {
  printf("this is working!");
  const unsigned NumRadBins = 100;
  const double MinRad = 0;
  const double Maxk = 140;

  TGraph *gCk = new TGraph();
  gCk->SetName("gCk");
  unsigned COUNTER = 0;
  for (double k = 0.1; k < Maxk; k += 10) {
    gCk->SetPoint(
        COUNTER,
        k,
        0.666 * GET_CORRELATION(k, ScatLend1, EffecRange1, ChargeRad, SourceRad)
            + 0.333
                * GET_CORRELATION(k, ScatLend2, EffecRange2, ChargeRad,
                                  SourceRad));
    /* printf(
     "%f  %f\n",
     k,
     2 / 3 * GET_CORRELATION(k, ScatLend1, EffecRange1, ChargeRad, SourceRad)
     + 1 / 3
     * GET_CORRELATION(k, ScatLend2, EffecRange2, ChargeRad,
     SourceRad));*/
    COUNTER++;
  }

  TFile *OutputFile =
      new TFile(
          TString::Format(
              "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD/ProtonDeuteronCkLedCoulombUpdatedScatterParAprrox.root"),
          "recreate");
  printf("File Created\n");
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName =
      Form(
          "%sCkPD.png",
          "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD");  //plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");
  //h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  gCk->GetXaxis()->SetTitle("k (MeV/c)");
  gCk->GetXaxis()->SetTitleSize(0.045);
  gCk->GetXaxis()->SetTitleOffset(0.75);
  gCk->GetYaxis()->SetTitleSize(0.045);
  gCk->GetYaxis()->SetTitleOffset(0.95);
  gCk->GetYaxis()->SetTitle("#it{C}(k)");

  gCk->GetXaxis()->SetRangeUser(0, 200);
  gCk->GetYaxis()->SetRangeUser(0, 10);

  gCk->SetLineColor(kRed);
  gCk->SetLineWidth(2);
  //gSource->MarkerStyle(20);
  //gSource->MarkerSize(1.5);
  TLegend *leg1 = new TLegend(0.5, 0.5, 0.7, 0.65);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.04);
  leg1->SetLineColor(0);
  leg1->SetNColumns(1);
  leg1->AddEntry(gCk, " Ck pd Lednicky Coulomb");

  gCk->Draw("ACP");
  leg1->Draw("same");
  printf("reached -3!");
  gCk->Write();
  Plot->Print(pdfName);
  delete Plot;
  delete gCk;
  delete OutputFile;
}
void PLOT_CORRELATION_SIMPSON_METHOD(double ScatLend1, double EffecRange1,
                                     double ScatLend2, double EffecRange2,
                                     double ChargeRad, double SourceRad,
                                     bool plotall, bool QS) {
  bool doublet = false;
  bool quartet = false;
  printf("this is working!");
  const unsigned NumRadBins = 100;
  double MinRad = 0;
  double Maxk = 400;
  double Scat1_Brokman = 11.4;  //quartet
  double EffRange1_Brokman = 2.05;
  double Scat2_Brokman = 1.3;  //doublet
  double EffRange2_Brokman = 0.0;

  double Scat1_Arvieux = 11.88;  //quartet
  double EffRange1_Arvieux = 2.63;
  double Scat2_Arvieux = 2.73;  //doublet
  double EffRange2_Arvieux = 0.0;

  double Scat1_Huttel = 11.1;  //quartet
  double EffRange1_Huttel = 0.0;
  double Scat2_Huttel = 4.0;  //doublet
  double EffRange2_Huttel = 0.0;

  double Scat1_Keivsky = 13.8;  //quartet
  double EffRange1_Keivsky = 0.0;
  double Scat2_Keivsky = 0.024;  //doublet
  double EffRange2_Keivsky = 0.0;

  double Scat1_Black = 14.7;  //quartet
  double EffRange1_Black = 0.0;
  double Scat2_Black = -0.13;   //doublet
  double EffRange2_Black = 0.0;
  float w1 = 2.0 / 3.0;
  float w2 = 1.0 / 3.0;
  double EREquartet = 2.63;
  double EREdoublet = 2.27;

 // double EREquartet = 0.0;
 //  double EREdoublet = 0.0;
  if (doublet) {
    Scat1_Brokman = 0.0;  //quartet
    EffRange1_Brokman = 0.0;
    Scat1_Arvieux = 0.0;  //quartet
    EffRange1_Arvieux = 0.0;
    Scat1_Huttel = 0.0;  //quartet
    EffRange1_Huttel = 0.0;
    Scat1_Keivsky = 0.0;  //quartet
    EffRange1_Keivsky = 0.0;
    Scat1_Black = 0.0;  //quartet
    EffRange1_Black = 0.0;
  } else if (quartet) {
    Scat2_Brokman = 0.0;  //doublet
    EffRange2_Brokman = 0.0;
    Scat2_Arvieux = 0.0;  //doublet
    EffRange2_Arvieux = 0.0;
    Scat2_Huttel = 0.0;  //doublet
    EffRange2_Huttel = 0.0;
    Scat2_Keivsky = 0.0;  //doublet
    EffRange2_Keivsky = 0.0;
    Scat2_Black = 0.0;   //doublet
    EffRange2_Black = 0.0;
  } else {
    Scat1_Brokman = 11.4;  //quartet
    EffRange1_Brokman = EREquartet;
    Scat2_Brokman = 1.3;  //doublet
    EffRange2_Brokman = EREdoublet;

    Scat1_Arvieux = 11.88;  //quartet
    EffRange1_Arvieux = EREquartet;
    Scat2_Arvieux = 2.73;  //doublet
    EffRange2_Arvieux = EREdoublet;

    Scat1_Huttel = 11.1;  //quartet
    EffRange1_Huttel = EREquartet;
    Scat2_Huttel = 4.0;  //doublet
    EffRange2_Huttel = EREdoublet;

    Scat1_Keivsky = 13.8;  //quartet
    EffRange1_Keivsky = EREquartet;
    Scat2_Keivsky = 0.024;  //doublet
    EffRange2_Keivsky = EREdoublet;

    Scat1_Black = 14.7;  //quartet
    EffRange1_Black = EREquartet;
    Scat2_Black = -0.13;   //doublet
    EffRange2_Black = EREdoublet;
  }

  TGraph *gCk = new TGraph();
  TGraph *gCk_Brockman = new TGraph();
  TGraph *gCk_Arveux = new TGraph();
  TGraph *gCk_Huttel = new TGraph();
  TGraph *gCk_Keivsky = new TGraph();
  TGraph *gCk_Black = new TGraph();

  TH1F *hCk = new TH1F("hCk", "hCk", NumRadBins, 0.0, Maxk);
  TH1F *hCk_Brockman = new TH1F("hCk_Brockman", "hCk_Brockman", NumRadBins, 0.0,
                                Maxk);
  TH1F *hCk_Arveux = new TH1F("hCk_Arveux", "hCk_Arveux", NumRadBins, 0.0,
                              Maxk);
  TH1F *hCk_Huttel = new TH1F("hCk_Huttel", "hCk_Huttel", NumRadBins, 0.0,
                              Maxk);
  TH1F *hCk_Keivsky = new TH1F("hCk_Keivsky", "hCk_Keivsky", NumRadBins, 0.0,
                               Maxk);
  TH1F *hCk_Black = new TH1F("hCk_Black", "hCk_Black", NumRadBins, 0.0, Maxk);

  gCk->SetName("TestBlack ");
  gCk_Brockman->SetName("gCk_Brockman");
  gCk_Arveux->SetName("gCk_Arveux");
  gCk_Huttel->SetName("gCk_Huttel");
  gCk_Keivsky->SetName("gCk_Keivsky");
  gCk_Black->SetName("gCk_Black");

  float CkValues = 0.0;
  float CkBrockmanValues = 0.0;
  float CkArveuxValues = 0.0;
  float CkHuttelValues = 0.0;
  float CkKeivskyValues = 0.0;
  float CkBlackValues = 0.0;
  unsigned COUNTER = 0;
  for (double k = 4; k < Maxk; k += 4) {
    if (plotall) {

      CkValues = 0.0;
      if (doublet) {
        CkBrockmanValues =  GET_CORRELATION_SIMPSON2D(k, Scat2_Brokman, EffRange2_Brokman,ChargeRad, SourceRad);
        CkArveuxValues = GET_CORRELATION_SIMPSON2D(k, Scat2_Arvieux, EffRange2_Arvieux,ChargeRad, SourceRad);
        CkHuttelValues =  GET_CORRELATION_SIMPSON2D(k, Scat2_Huttel, EffRange2_Huttel,ChargeRad, SourceRad);
        CkKeivskyValues = GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,ChargeRad, SourceRad);
        CkBlackValues =  GET_CORRELATION_SIMPSON2D(k, Scat2_Black, EffRange2_Black,ChargeRad, SourceRad);

      } else if (quartet) {
        CkBrockmanValues = GET_CORRELATION_SIMPSON2D(k, Scat1_Brokman, EffRange1_Brokman,ChargeRad, SourceRad);
        CkArveuxValues =  GET_CORRELATION_SIMPSON2D(k, Scat1_Arvieux, EffRange1_Arvieux,ChargeRad, SourceRad);
        CkHuttelValues = GET_CORRELATION_SIMPSON2D(k, Scat1_Huttel, EffRange1_Huttel,ChargeRad, SourceRad);
        CkKeivskyValues =  GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,ChargeRad, SourceRad);
        CkBlackValues =  GET_CORRELATION_SIMPSON2D(k, Scat1_Black, EffRange1_Black,ChargeRad, SourceRad);
        gCk_Brockman->SetPoint(COUNTER, k, CkBrockmanValues);
      } else {
        CkBrockmanValues = w1
            * GET_CORRELATION_SIMPSON2D(k, Scat1_Brokman, EffRange1_Brokman,
                                        ChargeRad, SourceRad)
            + w2
                * GET_CORRELATION_SIMPSON2D(k, Scat2_Brokman, EffRange2_Brokman,
                                            ChargeRad, SourceRad);
        CkArveuxValues = w1
            * GET_CORRELATION_SIMPSON2D(k, Scat1_Arvieux, EffRange1_Arvieux,
                                        ChargeRad, SourceRad)
            + w2
                * GET_CORRELATION_SIMPSON2D(k, Scat2_Arvieux, EffRange2_Arvieux,
                                            ChargeRad, SourceRad);
        CkHuttelValues = w1
            * GET_CORRELATION_SIMPSON2D(k, Scat1_Huttel, EffRange1_Huttel,
                                        ChargeRad, SourceRad)
            + w2
                * GET_CORRELATION_SIMPSON2D(k, Scat2_Huttel, EffRange2_Huttel,
                                            ChargeRad, SourceRad);
        CkKeivskyValues = w1
            * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                        ChargeRad, SourceRad)
            + w2
                * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                            ChargeRad, SourceRad);
        CkBlackValues = w1
            * GET_CORRELATION_SIMPSON2D(k, Scat1_Black, EffRange1_Black,
                                        ChargeRad, SourceRad)
            + w2
                * GET_CORRELATION_SIMPSON2D(k, Scat2_Black, EffRange2_Black,
                                            ChargeRad, SourceRad);
      }
      gCk_Brockman->SetPoint(COUNTER, k, CkBrockmanValues);
      hCk_Brockman->SetBinContent(COUNTER + 1, CkBrockmanValues);
      gCk_Arveux->SetPoint(COUNTER, k, CkArveuxValues);
      hCk_Arveux->SetBinContent(COUNTER + 1, CkArveuxValues);
      gCk_Huttel->SetPoint(COUNTER, k, CkHuttelValues);
      hCk_Huttel->SetBinContent(COUNTER + 1, CkHuttelValues);
      gCk_Keivsky->SetPoint(COUNTER, k, CkKeivskyValues);
      hCk_Keivsky->SetBinContent(COUNTER + 1, CkKeivskyValues);
      gCk_Black->SetPoint(COUNTER, k, CkBlackValues);
      hCk_Black->SetBinContent(COUNTER + 1, CkBlackValues);
    } else {
      if (QS) {
        gCk->SetPoint(
            COUNTER,
            k,
            (1.0 / 4.0
                * GET_CORRELATION_SIMPSON2D(k, ScatLend1, EffecRange1,
                                            ChargeRad, SourceRad)
                + 3.0 / 4.0
                    * GET_CORRELATION_SIMPSON2D(k, ScatLend2, EffecRange2,
                                                ChargeRad, SourceRad))
                - 0.5 * exp(-SourceRad * SourceRad * 4. * k * k));
      } else {

        CkValues = 2.0 / 3.0
            * GET_CORRELATION_SIMPSON2D(k, ScatLend1, EffecRange1, ChargeRad,
                                        SourceRad)
            + 1.0 / 3.0
                * GET_CORRELATION_SIMPSON2D(k, ScatLend2, EffecRange2,
                                            ChargeRad, SourceRad);
        gCk->SetPoint(COUNTER, k, CkValues);
        hCk->SetBinContent(COUNTER + 1, CkValues);

      }
    }
    COUNTER++;
  }
  TString OutputFileName = "";
  if (plotall) {

    if (doublet) {
      OutputFileName =
          OutputFileName
              + "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/2Doublet_AllScatteringCoreReso.root";
    } else if (quartet) {
      OutputFileName =
          OutputFileName
              + "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/2Quartet_AllScatteringCoreReso.root";
    } else {
      OutputFileName =
          OutputFileName
              + "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/UpEffectiveRangeAllScatteringCoreReso.root";
    }

  } else {
    OutputFileName =
        OutputFileName
            + "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/AllScatteringCoreReso.root";
  }
  TFile *OutputFile = new TFile(OutputFileName, "recreate");
  printf("File Created\n");
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName = "";
  if (plotall) {

    if (doublet) {
      pdfName =
          pdfName
              + Form(
                  "%s2Doublet_AllScatteringCoreReso.pdf",
                  "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/");
    } else if (quartet) {
      pdfName =
          pdfName
              + Form(
                  "%s2Quartet_AllScatteringCoreReso.pdf",
                  "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/");
    } else {
      pdfName =
          pdfName
              + Form(
                  "%sAllScatteringCoreReso.pdf",
                  "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/");
    }

  } else {
    pdfName =
        pdfName
            + Form(
                "%sKieveskyDoublet4.pdf",
                "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/");
  }  //plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");
//h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  if (plotall) {
    gCk_Brockman->GetXaxis()->SetTitle("k (MeV/c)");
    gCk_Brockman->GetXaxis()->SetTitleSize(0.045);
    gCk_Brockman->GetXaxis()->SetTitleOffset(0.5);
    gCk_Brockman->GetYaxis()->SetTitleSize(0.045);
    gCk_Brockman->GetYaxis()->SetTitleOffset(0.65);
    gCk_Brockman->GetYaxis()->SetTitle("#it{C}(k)");

    gCk_Brockman->GetXaxis()->SetRangeUser(0, 400);
    gCk_Brockman->GetYaxis()->SetRangeUser(0, 50);

    gCk_Arveux->SetLineColor(kRed);
    gCk_Arveux->SetLineWidth(1.5);

    gCk_Huttel->SetLineColor(kBlack);
    gCk_Huttel->SetLineWidth(1.5);

    gCk_Keivsky->SetLineColor(kGreen);
    gCk_Keivsky->SetLineWidth(1.5);

    gCk_Black->SetLineColor(kBlue);
    gCk_Black->SetLineWidth(1.5);

    gCk_Arveux->SetLineStyle(10);
    gCk_Keivsky->SetLineStyle(8);
    gCk_Black->SetLineStyle(5);
    gCk_Huttel->SetLineStyle(4);
    gCk_Brockman->SetLineStyle(3);
    gCk_Arveux->SetMarkerStyle(25);
    gCk_Keivsky->SetMarkerStyle(26);
    gCk_Black->SetMarkerStyle(27);
    gCk_Huttel->SetMarkerStyle(28);
    gCk_Brockman->SetMarkerStyle(29);

    hCk_Brockman->GetXaxis()->SetTitle("k (MeV/c)");
    hCk_Brockman->GetXaxis()->SetTitleSize(0.045);
    hCk_Brockman->GetXaxis()->SetTitleOffset(0.5);
    hCk_Brockman->GetYaxis()->SetTitleSize(0.045);
    hCk_Brockman->GetYaxis()->SetTitleOffset(0.65);
    hCk_Brockman->GetYaxis()->SetTitle("#it{C}(k)");

    hCk_Brockman->GetXaxis()->SetRangeUser(0, 400);
    hCk_Brockman->GetYaxis()->SetRangeUser(0, 50);

    hCk_Arveux->SetLineColor(kRed);
    hCk_Arveux->SetLineWidth(1.5);

    hCk_Huttel->SetLineColor(kBlack);
    hCk_Huttel->SetLineWidth(1.5);

    hCk_Keivsky->SetLineColor(kGreen);
    hCk_Keivsky->SetLineWidth(1.5);

    hCk_Black->SetLineColor(kBlue);
    hCk_Black->SetLineWidth(1.5);

    hCk_Arveux->SetLineStyle(10);
    hCk_Keivsky->SetLineStyle(8);
    hCk_Black->SetLineStyle(5);
    hCk_Huttel->SetLineStyle(4);
    hCk_Brockman->SetLineStyle(3);
    hCk_Arveux->SetMarkerStyle(25);
    hCk_Keivsky->SetMarkerStyle(26);
    hCk_Black->SetMarkerStyle(27);
    hCk_Huttel->SetMarkerStyle(28);
    hCk_Brockman->SetMarkerStyle(29);

    //gSource->MarkerStyle(20);
    //gSource->MarkerSize(1.5);
    TLegend *leg1 = new TLegend(0.5, 0.5, 0.7, 0.65);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(0);
    leg1->SetNColumns(1);
    leg1->AddEntry(gCk_Brockman, "Ck van Oers,Brockman (1967)");
    leg1->AddEntry(gCk_Arveux, "Ck Arvieux (1973)");
    leg1->AddEntry(gCk_Huttel, "Ck Huttel et al.(1983)");
    leg1->AddEntry(gCk_Keivsky, "Ck Kievsky et al.(1997)");
    leg1->AddEntry(gCk_Black, " Ck Black et al.(1999)");
    gCk_Brockman->Draw("ACP");
    gCk_Arveux->Draw("same");
    gCk_Huttel->Draw("same ");
    gCk_Keivsky->Draw("same");
    gCk_Black->Draw("same");
    leg1->Draw("same");
    printf("reached -3!");
    gCk_Brockman->Write();
    gCk_Arveux->Write();
    gCk_Huttel->Write();
    gCk_Keivsky->Write();
    gCk_Black->Write();

    hCk_Brockman->Write();
    hCk_Arveux->Write();
    hCk_Huttel->Write();
    hCk_Keivsky->Write();
    hCk_Black->Write();
    Plot->Print(pdfName);
  } else {
    gCk->GetXaxis()->SetTitle("k (MeV/c)");
    gCk->GetXaxis()->SetTitleSize(0.045);
    gCk->GetXaxis()->SetTitleOffset(0.75);
    gCk->GetYaxis()->SetTitleSize(0.045);
    gCk->GetYaxis()->SetTitleOffset(0.95);
    gCk->GetYaxis()->SetTitle("#it{C}(k)");

    gCk->GetXaxis()->SetRangeUser(0, 400);
    gCk->GetYaxis()->SetRangeUser(0, 50);

    gCk->SetLineColor(kRed);
    gCk->SetLineWidth(2);

    hCk->GetXaxis()->SetTitle("k (MeV/c)");
    hCk->GetXaxis()->SetTitleSize(0.045);
    hCk->GetXaxis()->SetTitleOffset(0.75);
    hCk->GetYaxis()->SetTitleSize(0.045);
    hCk->GetYaxis()->SetTitleOffset(0.95);
    hCk->GetYaxis()->SetTitle("#it{C}(k)");

    hCk->GetXaxis()->SetRangeUser(0, 400);
    hCk->GetYaxis()->SetRangeUser(0, 50);

    hCk->SetLineColor(kRed);
    hCk->SetLineWidth(2);
    //gSource->MarkerStyle(20);
    //gSource->MarkerSize(1.5);
    TLegend *leg1 = new TLegend(0.2, 0.5, 0.7, 0.65);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(0);
    leg1->SetNColumns(1);
    leg1->AddEntry(gCk, " Ck D^{-}-p [a_{0}= 0.07, a_{1} = -0.45]");

    gCk->Draw("ACP");
    //hCk->Draw(" ");
    leg1->Draw("same");
    printf("reached -3!");
    gCk->Write();
    hCk->Write();
    Plot->Print(pdfName);
  }
  delete Plot;
  delete gCk;
  delete gCk_Brockman;
  delete gCk_Arveux;
  delete gCk_Huttel;
  delete gCk_Keivsky;
  delete gCk_Black;
  delete hCk;
  delete hCk_Brockman;
  delete hCk_Arveux;
  delete hCk_Huttel;
  delete hCk_Keivsky;
  delete hCk_Black;
  delete OutputFile;
}

void PLOT_RadiusVsKeivesky(double ScatLend1, double EffecRange1,
                           double ScatLend2, double EffecRange2,
                           double ChargeRad) {

  printf("this is working!");
  const unsigned NumRadBins = 80;
  double MinRad = 0;
  double Maxk = 320;
  double EREquartet = 2.63;
  double EREdoublet = 2.27;

 // double EREquartet = 0.0;
 //  double EREdoublet = 0.0;
  double Scat1_Keivsky = 13.8;  //quartet
  double EffRange1_Keivsky = EREquartet;
  double Scat2_Keivsky = 0.024;  //doublet
  double EffRange2_Keivsky = EREdoublet;

  TGraph *gR1 = new TGraph();
  TGraph *gR2 = new TGraph();
  TGraph *gR3 = new TGraph();
  TGraph *gR4 = new TGraph();
  TGraph *gR5 = new TGraph();
  TGraph *gR6 = new TGraph();
  TGraph *gR7 = new TGraph();
  TGraph *gR8 = new TGraph();
  TGraph *gR9 = new TGraph();
  TGraph *gR10 = new TGraph();
  TGraph *gR11 = new TGraph();
  TGraph *gR12 = new TGraph();
  TGraph *gR13 = new TGraph();
  TGraph *gR14 = new TGraph();
  TGraph *gR15 = new TGraph();
  TGraph *gR16 = new TGraph();
  TGraph *gR17 = new TGraph();

  TH1F *hR1 = new TH1F("hR1", "hR1", NumRadBins, 0.0, Maxk);
  TH1F *hR2 = new TH1F("hR2", "hR2", NumRadBins, 0.0, Maxk);
  TH1F *hR3 = new TH1F("hR3", "hR3", NumRadBins, 0.0, Maxk);
  TH1F *hR4 = new TH1F("hR4", "hR4", NumRadBins, 0.0, Maxk);
  TH1F *hR5 = new TH1F("hR5", "hR5", NumRadBins, 0.0, Maxk);
  TH1F *hR6 = new TH1F("hR6", "hR6", NumRadBins, 0.0, Maxk);
  TH1F *hR7 = new TH1F("hR7", "hR7", NumRadBins, 0.0, Maxk);
  TH1F *hR8 = new TH1F("hR8", "hR8", NumRadBins, 0.0, Maxk);
  TH1F *hR9 = new TH1F("hR9", "hR9", NumRadBins, 0.0, Maxk);
  TH1F *hR10 = new TH1F("hR10", "hR10", NumRadBins, 0.0, Maxk);
  TH1F *hR11 = new TH1F("hR11", "hR11", NumRadBins, 0.0, Maxk);
  TH1F *hR12 = new TH1F("hR12", "hR12", NumRadBins, 0.0, Maxk);

  TH1F *hR13 = new TH1F("hR13", "hR13", NumRadBins, 0.0, Maxk);
  TH1F *hR14 = new TH1F("hR14", "hR14", NumRadBins, 0.0, Maxk);
  TH1F *hR15 = new TH1F("hR15", "hR15", NumRadBins, 0.0, Maxk);
  TH1F *hR16 = new TH1F("hR16", "hR16", NumRadBins, 0.0, Maxk);
  TH1F *hR17 = new TH1F("hR17", "hR17", NumRadBins, 0.0, Maxk);
  gR1->SetName("gR1");
  gR2->SetName("gR2");
  gR3->SetName("gR3");
  gR4->SetName("gR4");
  gR5->SetName("gR5");
  gR6->SetName("gR6");
  gR7->SetName("gR7");
  gR8->SetName("gR8");
  gR9->SetName("gR9");
  gR10->SetName("gR10");
  gR11->SetName("gR11");
  gR12->SetName("gR12");
  gR13->SetName("gR13");
  gR14->SetName("gR14");
  gR15->SetName("gR15");
  gR16->SetName("gR16");
  gR17->SetName("gR17");

  float CkR1 = 0.0;
  float CkR2 = 0.0;
  float CkR3 = 0.0;
  float CkR4 = 0.0;
  float CkR5 = 0.0;
  float CkR6 = 0.0;
  float CkR7 = 0.0;
  float CkR8 = 0.0;
  float CkR9 = 0.0;
  float CkR10 = 0.0;
  float CkR11 = 0.0;
  float CkR12 = 0.0;
  float CkR13 = 0.0;
  float CkR14 = 0.0;
  float CkR15 = 0.0;
  float CkR16 = 0.0;
  float CkR17 = 0.0;

  float R1 = 1.0;
  float R2 = 1.2;
  float R3 = 1.5;
  float R4 = 1.8;
  float R5 = 2.0;
  float R6 = 2.5;
  float R7 = 3.0;
  float R8 = 4.0;
  float R9 = 5.0;
  float R10 = 6.0;
  float R11 = 8.0;
  float R12 = 10.0;
  float R13 = 12.0;
  float R14 = 15.0;
  float R15 = 20.0;
  float R16 = 40.0;
  float R17 = 70.0;
  unsigned COUNTER = 0;
  for (double k = 4; k < Maxk; k += 4) {

    CkR1 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R1)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R1);
    CkR2 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R2)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R2);
    CkR3 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R3)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R3);
    CkR4 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R4)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R4);
    CkR5 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R5)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R5);
    CkR6 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R6)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R6);
    CkR7 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R7)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R7);
    CkR8 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R8)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R8);
    CkR9 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R9)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R9);
    CkR10 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R10)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R10);
    CkR11 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R11)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R11);
    CkR12 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R12)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R12);
    CkR13 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R13)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R13);
    CkR14 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R14)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R14);
    CkR15 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R15)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R15);
    CkR16 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R16)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R16);
    CkR17 = 0.666
        * GET_CORRELATION_SIMPSON2D(k, Scat1_Keivsky, EffRange1_Keivsky,
                                    ChargeRad, R17)
        + 0.333
            * GET_CORRELATION_SIMPSON2D(k, Scat2_Keivsky, EffRange2_Keivsky,
                                        ChargeRad, R17);

    gR1->SetPoint(COUNTER, k, CkR1);
    gR2->SetPoint(COUNTER, k, CkR2);
    gR3->SetPoint(COUNTER, k, CkR3);
    gR4->SetPoint(COUNTER, k, CkR4);
    gR5->SetPoint(COUNTER, k, CkR5);
    gR6->SetPoint(COUNTER, k, CkR6);
    gR7->SetPoint(COUNTER, k, CkR7);
    gR8->SetPoint(COUNTER, k, CkR8);
    gR9->SetPoint(COUNTER, k, CkR9);
    gR10->SetPoint(COUNTER, k, CkR10);
    gR11->SetPoint(COUNTER, k, CkR11);
    gR12->SetPoint(COUNTER, k, CkR12);
    gR13->SetPoint(COUNTER, k, CkR13);
    gR14->SetPoint(COUNTER, k, CkR14);
    gR15->SetPoint(COUNTER, k, CkR15);
    gR16->SetPoint(COUNTER, k, CkR16);
    gR17->SetPoint(COUNTER, k, CkR17);

    hR1->SetBinContent(COUNTER + 1, CkR1);
    hR2->SetBinContent(COUNTER + 1, CkR2);
    hR3->SetBinContent(COUNTER + 1, CkR3);
    hR4->SetBinContent(COUNTER + 1, CkR4);
    hR5->SetBinContent(COUNTER + 1, CkR5);
    hR6->SetBinContent(COUNTER + 1, CkR6);
    hR7->SetBinContent(COUNTER + 1, CkR7);
    hR8->SetBinContent(COUNTER + 1, CkR8);
    hR9->SetBinContent(COUNTER + 1, CkR9);
    hR10->SetBinContent(COUNTER + 1, CkR10);
    hR11->SetBinContent(COUNTER + 1, CkR11);
    hR12->SetBinContent(COUNTER + 1, CkR12);
    hR13->SetBinContent(COUNTER + 1, CkR13);
    hR14->SetBinContent(COUNTER + 1, CkR14);
    hR15->SetBinContent(COUNTER + 1, CkR15);
    hR16->SetBinContent(COUNTER + 1, CkR16);
    hR17->SetBinContent(COUNTER + 1, CkR17);
    COUNTER++;
  }
  TString OutputFileName = "";
  OutputFileName =
      OutputFileName
          + "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/LowKeiveskyVsRadPlots.root";

  TFile *OutputFile = new TFile(OutputFileName, "recreate");
  printf("File Created\n");
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName = "";
  pdfName =
      pdfName
          + Form(
              "%sLowKeiveskyVsRad.pdf",
              "/home/sbhawani/cernbox/ProtonDeuteron/Outputs/AODs/Correlations/3HeComparision/SourceFiles/");
//plots are output as .pdf If you prefer other formats simply change the ending
  printf("reached -2!");
//h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  gR1->GetXaxis()->SetTitle("k (MeV/c)");
  gR1->GetXaxis()->SetTitleSize(0.045);
  gR1->GetXaxis()->SetTitleOffset(0.5);
  gR1->GetYaxis()->SetTitleSize(0.045);
  gR1->GetYaxis()->SetTitleOffset(0.65);
  gR1->GetYaxis()->SetTitle("#it{C}(k)");

  gR1->GetXaxis()->SetRangeUser(0, 400);
  gR1->GetYaxis()->SetRangeUser(0, 50);

  hR1->GetXaxis()->SetTitle("k (MeV/c)");
  hR1->GetXaxis()->SetTitleSize(0.045);
  hR1->GetXaxis()->SetTitleOffset(0.5);
  hR1->GetYaxis()->SetTitleSize(0.045);
  hR1->GetYaxis()->SetTitleOffset(0.65);
  hR1->GetYaxis()->SetTitle("#it{C}(k)");

  hR1->GetXaxis()->SetRangeUser(0, 400);
  hR1->GetYaxis()->SetRangeUser(0, 50);

  gR1->Draw("ACP");
  gR1->Write();
  gR2->Write();
  gR3->Write();
  gR4->Write();
  gR5->Write();
  gR6->Write();
  gR7->Write();
  gR8->Write();
  gR9->Write();
  gR10->Write();
  gR11->Write();
  gR12->Write();
  gR13->Write();
  gR14->Write();
  gR15->Write();
  gR16->Write();
  gR17->Write();

  hR1->Write();
  hR2->Write();
  hR3->Write();
  hR4->Write();
  hR5->Write();
  ;
  hR6->Write();
  hR7->Write();
  hR8->Write();
  hR9->Write();
  hR10->Write();
  hR11->Write();
  hR12->Write();
  hR13->Write();
  hR14->Write();
  hR15->Write();
  hR16->Write();
  hR17->Write();
  printf("reached -3!");
  Plot->Print(pdfName);

  delete gR1;
  delete gR2;
  delete gR3;
  delete gR4;
  delete gR5;
  delete gR6;
  delete gR7;
  delete gR8;
  delete gR9;
  delete gR10;
  delete gR11;
  delete gR12;
  delete gR13;
  delete gR14;
  delete gR15;
  delete gR16;
  delete gR17;

  delete hR1;
  delete hR2;
  delete hR3;
  delete hR4;
  delete hR5;
  delete hR6;
  delete hR7;
  delete hR8;
  delete hR9;
  delete hR10;
  delete hR11;
  delete hR12;
  delete hR13;
  delete hR14;
  delete hR15;
  delete hR16;
  delete hR17;
}

void PLOT_CORRELATION_PROTON_DMESON(double ChargeRad, double SourceRad,
                                    bool plotall) {
  printf("this is working!");
  const unsigned NumRadBins = 80;
  const double MinRad = 0;
  const double Maxk = 400;

//quark gluon exchange Model
  const double aI0 = -0.13;  //ISOSpin 0
  const double rI0 = 0.0;
  const double aI1 = -0.29;  //ISOSpin 1
  const double rI1 = 0.0;
//Full Model
  const double aI0full = 0.07;  //ISOSpin 0
  const double rI0full = 0.0;
  const double aI1full = -0.45;  //ISOSpin 1
  const double rI1full = 0.0;

  const double wChannel1 = 0.5;  // ISO Spin 0
  const double wChannel2 = 0.5;  // ISO Spin 1

  TGraph *gCk = new TGraph();
  TGraph *gCk_PDqgModel = new TGraph();
  TGraph *gCk_FullModel = new TGraph();

  gCk->SetName("TestBlack ");
  gCk_PDqgModel->SetName("gCk_PDqgModel");
  gCk_FullModel->SetName("gCk_FullModel");

  unsigned COUNTER = 0;
  for (double k = 5; k < Maxk; k += 5) {
    if (plotall) {
      gCk_PDqgModel->SetPoint(
          COUNTER,
          k,
          wChannel1
              * GET_CORRELATION_SIMPSON2D(k, aI0, rI0, ChargeRad, SourceRad)
              + wChannel1
                  * GET_CORRELATION_SIMPSON2D(k, aI1, rI1, ChargeRad,
                                              SourceRad));
      gCk_FullModel->SetPoint(
          COUNTER,
          k,
          wChannel1
              * GET_CORRELATION_SIMPSON2D(k, aI0full, rI0full, ChargeRad,
                                          SourceRad)
              + wChannel1
                  * GET_CORRELATION_SIMPSON2D(k, aI1full, rI1full, ChargeRad,
                                              SourceRad));
    } else {
      gCk->SetPoint(
          COUNTER,
          k,
          (0.5 * GET_CORRELATION_SIMPSON2D(k, 0.00, 0.0, ChargeRad, SourceRad)
              + 0.5
                  * GET_CORRELATION_SIMPSON2D(k, 0.0, 0.0, ChargeRad, SourceRad)));

    }
    COUNTER++;
  }
  TString OutputFileName = "";
  if (plotall) {
    OutputFileName =
        OutputFileName

            + "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD/ProtonDMeson/CkProtonDMesonStrongCoulomb.root";
  } else {
    OutputFileName =
        OutputFileName
            + "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD/ProtonDMeson/CkProtonDMesonCoulomb.root";
  }
  TFile *OutputFile = new TFile(OutputFileName, "recreate");
  printf("File Created\n");
  TCanvas *Plot = new TCanvas("Plot", "Plot", 0, 0, 800, 600);
  TString pdfName = "";
  if (plotall) {
    pdfName =
        pdfName
            + Form(
                "%sCkProtonDMeson.pdf",
                "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD/ProtonDMeson/");
  } else {
    pdfName =
        pdfName
            + Form(
                "%sCkProtonDMesonCoulombOnly.pdf",
                "/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD/ProtonDMeson/");
  }
  printf("reached -2!");
//h1->SetTitle(Form("%s Channel [#Lambda = %i MeV, #it{m}_{#gamma} = %0.1f MeV]", Name.Data(), Cutoff, photonmass));
  if (plotall) {
    gCk_PDqgModel->GetXaxis()->SetTitle("k (MeV/c)");
    gCk_PDqgModel->GetXaxis()->SetTitleSize(0.045);
    gCk_PDqgModel->GetXaxis()->SetTitleOffset(0.5);
    gCk_PDqgModel->GetYaxis()->SetTitleSize(0.045);
    gCk_PDqgModel->GetYaxis()->SetTitleOffset(0.65);
    gCk_PDqgModel->GetYaxis()->SetTitle("#it{C}(k)");

    gCk_PDqgModel->GetXaxis()->SetRangeUser(0, 410);
    gCk_PDqgModel->GetYaxis()->SetRangeUser(0, 10);

    gCk_FullModel->SetLineColor(kRed);
    gCk_FullModel->SetLineWidth(1.5);

    gCk_FullModel->SetLineStyle(10);
    gCk_PDqgModel->SetLineStyle(3);
    gCk_FullModel->SetMarkerStyle(25);
    gCk_PDqgModel->SetMarkerStyle(29);

    TLegend *leg1 = new TLegend(0.5, 0.5, 0.7, 0.65);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(0);
    leg1->SetNColumns(1);
    leg1->AddEntry(gCk_PDqgModel, "p-D^{-} q-g Model");
    leg1->AddEntry(gCk_FullModel, "p-D^{-} Full Model");
    gCk_PDqgModel->Draw("ACP");
    gCk_FullModel->Draw("same");

    leg1->Draw("same");
    printf("reached -3!");
    gCk_PDqgModel->Write();
    gCk_FullModel->Write();
    Plot->Print(pdfName);
  } else {
    gCk->GetXaxis()->SetTitle("k (MeV/c)");
    gCk->GetXaxis()->SetTitleSize(0.045);
    gCk->GetXaxis()->SetTitleOffset(0.75);
    gCk->GetYaxis()->SetTitleSize(0.045);
    gCk->GetYaxis()->SetTitleOffset(0.95);
    gCk->GetYaxis()->SetTitle("#it{C}(k)");

    gCk->GetXaxis()->SetRangeUser(0, 400);
    gCk->GetYaxis()->SetRangeUser(0, 10);

    gCk->SetLineColor(kRed);
    gCk->SetLineWidth(2);
    //gSource->MarkerStyle(20);
    //gSource->MarkerSize(1.5);
    TLegend *leg1 = new TLegend(0.2, 0.5, 0.7, 0.65);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(0);
    leg1->SetNColumns(1);
    leg1->AddEntry(gCk, " Ck D^{-}-p [a_{0}= 0.07, a_{1} = -0.45]");

    gCk->Draw("ACP");
    leg1->Draw("same");
    printf("reached -3!");
    gCk->Write();
    Plot->Print(pdfName);
  }
  delete Plot;
  delete gCk;
  delete gCk_PDqgModel;
  delete gCk_FullModel;
  delete OutputFile;
}
