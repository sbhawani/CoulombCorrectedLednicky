#include"SketchRealisticPot.h"

void Ledni_SmallRad(TString PotentialName) {

	double Mass_p = 940;
	double Mass_Kch  = 1116;
	double Mass_L  = 1116;
	int NumMomBins = 60;
	//const TString PotentialName = "NSC97b";
	//const TString PotentialName = "NF48";
	//const TString PotentialName = "emma";
	//const TString PotentialName = "custom";

	const double kMin = 0;
	const double kMax = 300;
	const double kStep = 3;
	const unsigned nMom = TMath::Nint(kMax / kStep);
	const double Radius = 1.0;

	CATSparameters sPars(CATSparameters::tSource, 1, true);
	sPars.SetParameter(0, Radius);
	CATS Kitty_SE;
	Kitty_SE.SetMomBins(nMom, kMin, kMax);
	Kitty_SE.SetAnaSource(GaussSource, sPars);
	Kitty_SE.SetUseAnalyticSource(true);
	Kitty_SE.SetQ1Q2(0);
	Kitty_SE.SetQuantumStatistics(false);
	printf(" Hello reached -1 \n");
	Kitty_SE.SetRedMass(Mass_L * 0.5);
	//Kitty_SE.SetRedMass((Mass_p * Mass_Kch) / (Mass_p + Mass_Kch));
	Kitty_SE.SetNumChannels(1);
	Kitty_SE.SetNumPW(0, 1);
	Kitty_SE.SetSpin(0, 0);
	Kitty_SE.SetChannelWeight(0, 1.);
	CATSparameters pPars(CATSparameters::tPotential, 7, true);
	double c_f0, c_d0;
	if (PotentialName == "NSC97b") {
		pPars.SetParameter(0, -78.42);
		pPars.SetParameter(1, 1.0);
		pPars.SetParameter(2, 741.76);
		pPars.SetParameter(3, 0.45);
		c_f0 = 0.397;
		c_d0 = 10.360;
	} else if (PotentialName == "NF48") {
		pPars.SetParameter(0, -1647.40);
		pPars.SetParameter(1, 0.6);
		pPars.SetParameter(2, 3888.96);
		pPars.SetParameter(3, 0.45);
		c_f0 = 1.511;
		c_d0 = 2.549;
	} else if (PotentialName == "emma") {
		//for emma -> NSC97f basis
		pPars.SetParameter(0, -106.53 * 0.85);
		pPars.SetParameter(1, 1.0 * 1.18);
		pPars.SetParameter(2, 1469.33);
		pPars.SetParameter(3, 0.45 * 1.1);
		c_f0 = 0.350;
		c_d0 = 16.330;
	} else if (PotentialName == "Toy1") {
		pPars.SetParameter(0, -144.5);
		pPars.SetParameter(1, 2.11);
		pPars.SetParameter(2, 520.0);
		pPars.SetParameter(3, 0.54);
		c_f0 = -0.73;
		c_d0 = 7.72;
	} else if (PotentialName == "ND46") {
		pPars.SetParameter(0, -144.89);
		pPars.SetParameter(1, 1.0);
		pPars.SetParameter(2, 127.87);
		pPars.SetParameter(3, 0.45);
		c_f0 = -4.621;
		c_d0 = 1.3;
	} else if (PotentialName == "Yukawa1") {
		pPars.SetParameter(0, 1.0);
		pPars.SetParameter(1, 1.0);
		pPars.SetParameter(2, 100.0);
		pPars.SetParameter(3, 1.0);
		pPars.SetParameter(4, 0.4);
		c_f0 = 0;
		c_d0 = 1;
	} else if (PotentialName == "pKplusI0") {
		pPars.SetParameter(0, 0.0);
		c_f0 = 0.03;
		c_d0 = 0.0;
	} else if (PotentialName == "pKplusI1") {
		pPars.SetParameter(0, 1.0);
		c_f0 = -0.3;
		c_d0 = 0.0;
	} else if (PotentialName == "pKplusYuki") {
		pPars.SetParameter(0, 0.376); //0.376;0.335
		pPars.SetParameter(1, sqrt(200.*(Mass_p * Mass_Kch) / (Mass_p + Mass_Kch)));
		pPars.SetParameter(2, 3);
		pPars.SetParameter(3, 2084);
		pPars.SetParameter(4, 50.81);
		pPars.SetParameter(5, 18.34);
		pPars.SetParameter(6, -1.752);
		c_f0 = -0.3;
		c_d0 = 0.0;
	} else {
		pPars.SetParameter(0, -5.50337);
		pPars.SetParameter(1, 1. / sqrt(2.148805));
		pPars.SetParameter(2, 0);
		pPars.SetParameter(3, 1);
		//pPars.SetParameter(0,-78.42*0.39*4.5);//0.39,0.4
		//pPars.SetParameter(1,1.0*1.35);
		//pPars.SetParameter(2,741.76*4.5);
		//pPars.SetParameter(3,0.45*1.4);
		//NF46 as a stariting point
		//pPars.SetParameter(0,-1327.26*1.0);
		//pPars.SetParameter(1,0.6);
		//pPars.SetParameter(2,2561.56);
		//pPars.SetParameter(3,0.45);
		c_f0 = 0.02;
		c_d0 = 30.0;
		printf(" Hello\n");

	}
	Kitty_SE.SetEpsilonConv(5e-9);
	Kitty_SE.SetEpsilonProp(5e-9);
	//if (PotentialName.Contains("Yukawa"))
	//	Kitty_SE.SetShortRangePotential(0, 0, YukawaDimiCore, pPars);
	if (PotentialName.Contains("pKplusI")) {
		Kitty_SE.SetShortRangePotential(0, 0, KpProtonEquivalentPotential, pPars);
		Kitty_SE.SetEpsilonConv(1e-9);
		Kitty_SE.SetEpsilonProp(1e-9);
		//} else if (PotentialName.Contains("pKplusYuki")) {
		///	Kitty_SE.SetShortRangePotential(0, 0, SingleGaussDynamic, pPars);
		//Kitty_SE.SetEpsilonConv(1e-9);
		//Kitty_SE.SetEpsilonProp(1e-9);
	} else Kitty_SE.SetShortRangePotential(0, 0, DoubleGaussSum, pPars);
	printf(" Hello reached -2 \n");
	Kitty_SE.KillTheCat();
	printf(" Hello reached -3 \n");
	TFile* OutputFile = new TFile(
	    TString::Format("/home/sbhawani/ProtonDeuteron/Outputs/OutputCATSProjects/OutputCATSProjectPD/RealisticPotDrawNF48.root"), "recreate");
	printf("File Created\n");
	TGraph gKitty;
	gKitty.SetName(TString::Format("gKitty"));
	gKitty.Set(NumMomBins);
	for (unsigned uBin = 0; uBin < nMom; uBin++) {
		printf("C(%.2f) = %.2f\n", Kitty_SE.GetMomentum(uBin), Kitty_SE.GetCorrFun(uBin));
		gKitty.SetPoint(uBin, Kitty_SE.GetMomentum(uBin), Kitty_SE.GetCorrFun(uBin));
	} 
	printf(" Hello reached -4 \n");
	gKitty.Write();
	delete OutputFile;
	//delete cPars;
}