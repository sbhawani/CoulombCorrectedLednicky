//********************************************************//
//  References                                            //
//  1) Phy Letters B 802(2020)135223                      //
//  2) ISSN 1063-7796 Physics of Particle and Nuclie      //
//     2009, Vol.40,No.3,pp.3-7-352                       //
//********************************************************//

#ifndef COULOMBLEDNICKY_H
#define COULOMBLEDNICKY_H

#include "CATSconstants.h"
#include "TF2.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <complex>
#include <iomanip>
#include <complex.h>
#include <tgmath.h>
#include <arb.h>
#include <acb.h>
#include <acb_hypgeom.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#endif

using namespace std;

//double givenFunction(float x, float y);
//double doubleIntegral(float h, float k,
         //            float lx, float ux,
          //           float ly, float uy);
double BOHR_RADIUS();
double Ac(double eta);  //Gamow factor
double h(double eta);  // h function
// strong amplitude modified by the Coulomb interaction
complex<double> SCAT_AMP(double f0, double d0, double ac, double eta);
double B(double rho, double eta);  //recursive sum on bn
double P(double rho, double eta);  //recursive sum on pn
complex<double> TILDE_G(double rho, double eta);  //s-wave coulomb function
void ARB_HYPG_1F1(double *re, double *im, double a, double b, double z);
complex<double> CONFLUENT_HYPG_1F1(double eta, double zei);
complex<double> LEDNICKY_COULOMB_WF(double k, double r, double t,
                                    double ScatLend, double EffecRange,
                                    double ChargeRad);
double PROB_LEDNICKY_COULOMB_WF(double k, double r, double t, double ScatLend,
                                double EffecRange, double ChargeRad);
double PROB_LEDNICKY_COULOMB_WF_DEF_2(double k, double r, double t,
                                      double ScatLend, double EffecRange,
                                      double ChargeRad);
double SOURCE(double r, double R);
double SOURCE_4PI(double r, double R);
double DIFF_CORRELATION(double *arg, double *par);
double DIFF_CORRELATION_SIMPSON(double kval, double rval, double tval, double ScatLend,
                                double EffecRange, double ChargeRad,double Rval);
double GET_CORRELATION(double k, double ScatLend, double EffecRange,
                       double ChargeRad,double SourceRad);
double GET_CORRELATION_SIMPSON2D(double kval, double ScatLend,
                                 double EffecRange, double ChargeRad,
                                 double SourceRad);
