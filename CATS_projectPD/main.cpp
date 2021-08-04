#include <iostream>
#include <stdio.h>
#include <string.h>
//#include <omp.h>
#include <complex>

#include"SourceFit.h"
#include"SketchRealisticPot.h"
#include"CoulombLednicky.h"
#include"GenerateDecayMatrix.h"
#include"CorrelationCoulombLednicky.h"

using namespace std;

int main(int argc, char *argv[]) {
  printf("Hello there! This project is based on CATS and ROOT \n");
  DLM_Timer TIMER;

  //EffectiveGaussianPd();
  // printf("B = %f \n", B(1.0,2.0));
  //  printf("P = %f \n",P(0.1,1.0));
  //double val = abs(conj(CONFLUENT_HYPG_1F1(-1,1))*CONFLUENT_HYPG_1F1(-1,1));

  //LEDNICKY_COULOMB_WF(5,20, 1,  10.5, 0, 55);
  //PROB_LEDNICKY_COULOMB_WF(5,20, 1,  10.5, 0, 55);
  //Plot_Source_Pd(0.94);
  //Plot_CoulombLed_WF_ForPd(200,0,10.1,2.3,55);
  // Plot_CoulombLed_Density_ForPd(60,-0.5,10.1,2.3,55,true);
  //Plot_CoulombLed_Correlation_ForPd(14.7, 0.0, -1.13, 0.0, 43.207150, 0.94);
  // CHECK_FOR_FUNCTIONS();

  //printf("Ck 2/3 = %f\n", 2.0 / 3.0 * GET_CORRELATION(60, 14.7, 0.0,43.207150, 0.94));

  //printf("Ck 1/2 = %f\n", 1 / 3 * GET_CORRELATION(60, -1.13, 0.0, 43.207150, 0.94));

  //printf("charge radius for pi+pi+ = %f \n",BOHR_RADIUS());
  // printf("%f\n", doubleIntegral(0.01, 0.015, 2.3, 12.5, 3.7, 14.3));
  // printf("Ck 2/3 from ROOT = %f\n", 2.0 / 3.0 * GET_CORRELATION(60, 14.7, 0.0,43.207150, 0.94));
  // printf("Ck 2/3 from SimpsonMethod %f\n", 2.0 / 3.0 *GET_CORRELATION_SIMPSON2D(60, 14.7, 0.0,43.207150, 0.94));
 /* printf(
      "Ck(%f) from SimpsonMethod %f\n",
      22.0,
      1.0 / 4.0 * GET_CORRELATION_SIMPSON2D(22.0, 7.8064, 2.788, 56.207150, 1.2)
          + 3.0 / 4.0
              * GET_CORRELATION_SIMPSON2D(22.0, 0.0, 0.0, 56.207150, 1.2));

  printf(
      "Ck(%f) from SimpsonMethod %f\n",
      28.0,
      1.0 / 4.0 * GET_CORRELATION_SIMPSON2D(28.0, 7.8064, 2.788, 56.207150, 1.2)
          + 3.0 / 4.0
              * GET_CORRELATION_SIMPSON2D(28.0, 0.0, 0.0, 56.207150, 1.2));
  printf(
      "Ck(%f) from SimpsonMethod %f\n",
      20.0,
      1.0 / 4.0 * GET_CORRELATION_SIMPSON2D(20.0, 7.8064, 2.788, 56.207150, 1.2)
          + 3.0 / 4.0
              * GET_CORRELATION_SIMPSON2D(20.0, 0.0, 0.0, 43.207150, 0.99));*/

  //GetDecayMatrix();
  //PLOT_CORRELATION_SIMPSON_METHOD(-14.7, 0.0,-0.024,0.0, 43.207150, 1.073,true,false);
  //PLOT_RadiusVsKeivesky(-14.7, 0.0,-0.024,0.0, 43.207150);

  /*printf(
       "Ck(%f) from SimpsonMethod %f\n",
       50.0,0.5 * GET_CORRELATION_SIMPSON2D(50.0, -0.13,0.0, -56.207150, 1.0)
           + 0.5
               * GET_CORRELATION_SIMPSON2D(50.0, -0.29, 0.0, -56.207150, 1.0));*/
  //PLOT_CORRELATION_PROTON_DMESON(-43.21,1.0,false);
 // printf("Bohr radius of p-D^{-} =  %f",BOHR_RADIUS());



  Ledni_SmallRad("NF48");


  long long ExeTime = TIMER.Stop() / 1000.;
  char *strtime = new char[128];
  ShowTime(ExeTime, strtime, 0, true, 6);
  printf("The script terminated after: %s\n", strtime);
  delete[] strtime;
  return 0;
}

