#ifndef CORRELATIONCOULOMBLEDNICKY_H
#define CORRELATIONCOULOMBLEDNICKY_H
#include"CoulombLednicky.h"
#endif

void Plot_Source_Pd(double rad);
void Plot_CoulombLed_WF_ForPd(double k, double t, double ScatLend,
                              double EffecRange, double ChargeRad);
void Plot_CoulombLed_Density_ForPd(double k, double t, double ScatLend,
                                   double EffecRange, double ChargeRad,
                                   bool draw2D);

void Plot_CoulombLed_Correlation_ForPd(double ScatLend1, double EffecRange1,
                                       double ScatLend2, double EffecRange2,
                                       double ChargeRad, bool SourceRad);
void PLOT_CORRELATION_SIMPSON_METHOD(double ScatLend1, double EffecRange1,
                                     double ScatLend2, double EffecRange2,
                                     double ChargeRad, double SourceRad,bool plotall,bool QS);
void CHECK_FOR_FUNCTIONS();
void PLOT_CORRELATION_PROTON_DMESON(double ChargeRad, double SourceRad,bool plotall);
void PLOT_RadiusVsKeivesky(double ScatLend1, double EffecRange1,
                                     double ScatLend2, double EffecRange2,
                                     double ChargeRad);
