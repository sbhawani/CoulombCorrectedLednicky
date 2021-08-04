#include "GenerateDecayMatrix.h"

void GetDecayMatrix(){
    DLM_DecayMatrix decmat;
    decmat.SetFileName("/home/sbhawani/cernbox/ProtonDeuteron/Outputs/CATSOutput/DecayMatrices/Test2020_pL_pp.root");
    decmat.SetHistoName("Test2020_pL_pp");
    decmat.SetBins(1000,0,1000);
    decmat.SetNumDaughters1(1);
    decmat.SetNumDaughters2(2);
    decmat.SetDaughterMass1(0,Mass_p);
    decmat.SetDaughterMass2(0,Mass_p);
    decmat.SetDaughterMass2(1,Mass_Pic);
    decmat.SetMotherMass1(Mass_p);
    decmat.SetMotherMass2(Mass_L);
    decmat.SetMeanMomentum(0);
    decmat.SetMomentumSpread(350);
    decmat.Run(1,1250);
}
