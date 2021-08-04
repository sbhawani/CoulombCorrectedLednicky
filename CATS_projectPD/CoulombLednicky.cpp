#include"CoulombLednicky.h"
#include <bits/stdc++.h>
using namespace std;

double BOHR_RADIUS() {
  const double me = 0.510;//d mass
  const double mp = 938.27208816;  //p mass
  const double md = 1869.62;//d mass
  double mu_pp = mp / 2.0;
  double mu_pd = (mp * md) / (mp + md);
  double ac = 0.0;
  ac = me / mu_pd * 5.2917*pow(10,4.0);
  //ac = mu_pp / mu_pd * 57.6;
  return ac;
}

double Ac(double eta) {
  return 2 * M_PI * eta * 1 / (exp(2 * M_PI * eta) - 1);
}

double h(double eta) {
  const int N = 1000;  // check for numerical precision
  double t[N];
  double sum_accel, err;
  double sum = 0;
  int n;
  double h = 0;
  //if (abs(eta) < 0.3) {
  //  h = 1.2 * eta * eta - log(abs(eta)) - M_EULER;
 // } else {
    for (n = 1; n < N; n++) {
      t[n - 1] = (1.0 / (n * (n * n + eta * eta)));
      sum += t[n - 1];
   }
    h = eta * eta * sum - M_EULER - log(abs(eta));
  //}

  return h;
}

complex<double> SCAT_AMP(double k, double f0, double d0, double ac,
                         double eta) {
  complex<double> ScattAmpl = pow(1.0 / f0 + 0.5 * d0 * k * k - 2.0 / ac * h(eta)- i * k * Ac(eta),-1.0);
  return ScattAmpl;
}

double B(double rho, double eta) {
  const int N = 1000;  // check for numerical precision
  double b[N + 2];
  double sum = 0;
  int n;
  b[0] = 1.0;
  b[1] = rho * eta;
 // if (abs(rho * eta) < 0.0001) {
  //  sum = sin(rho) / rho;  // kind of bessel approximation
 // } else {
    for (n = 1; n <= N; n++) {
      b[n + 1] = (2.0 * eta * rho * b[n] - rho * rho * b[n - 1])
          / ((n + 1) * (n + 2));
      sum += b[n - 1];
    }
  //}
  return sum;
}

double P(double rho, double eta) {
  const int N = 1000;  // check for numerical precision
  double b[N + 2];
  double p[N + 2];
  double sum = 0;
  int n;
  b[0] = 1;
  b[1] = rho * eta;
  p[0] = 1;
  p[1] = 0;
 // if (abs(rho * eta) < 0.0001) {
  //  sum = cos(rho);  // kind of bessel approximation
 // } else {
    for (n = 1; n <= N; n++) {
      b[n + 1] = (2.0 * eta * rho * b[n] - rho * rho * b[n - 1])
          / ((n + 1) * (n + 2));
      p[n + 1] = (2.0 * eta * rho * p[n] - rho * rho * p[n - 1]
          - (2.0 * n + 1) * 2.0 * eta * rho * b[n]) / (n * (n + 1));
      sum += p[n - 1];
      //printf("p[%i] = %f \n",n-1, p[n-1]);
    }
 // }
  return sum;
}

complex<double> TILDE_G(double rho, double eta) {
  return P(rho, eta)
      + 2.0 * rho * eta
          * (log(abs(2.0 * rho * eta)) + 2.0 * M_EULER - 1.0 + h(eta)
              + i * Ac(eta) / (2.0 * eta)) * B(rho, eta);
}

void ARB_HYPG_1F1(double *re, double *im, double a, double b, double z) {

  acb_t aa, bb, zz, rr;
  acb_init(aa);
  acb_init(bb);
  acb_init(zz);
  acb_init(rr);

  acb_set_d(aa, a);
  acb_set_d(bb, b);
  acb_set_d(zz, z);
  acb_mul_onei(aa, aa);
  acb_mul_onei(zz, zz);
  slong prec = 64;

  acb_hypgeom_m(rr, aa, bb, zz, 0, prec);

  *re = arf_get_d(arb_midref(acb_realref(rr)), ARF_RND_DOWN);
  *im = arf_get_d(arb_midref(acb_imagref(rr)), ARF_RND_DOWN);
  acb_clear(aa);
  acb_clear(bb);
  acb_clear(zz);
  acb_clear(zz);
}

complex<double> CONFLUENT_HYPG_1F1(double eta, double zei) {
  double re, im;
  ARB_HYPG_1F1(&re, &im, eta, 1, zei);
  return re + i * im;
}

complex<double> LEDNICKY_COULOMB_WF(double k, double r, double t,
                                    double ScatLend, double EffecRange,
                                    double ChargeRad) {

  const double eta = 1.0 / (k * ChargeRad) / FmToNu;
  const double rho = k * r * FmToNu;
  const double zei = rho * (1 + t);
  const double f0 = ScatLend * FmToNu;
  const double d0 = EffecRange * FmToNu;
  const double ac = ChargeRad * FmToNu;
  complex<double> WF;
  WF = pow(Ac(eta), 0.5)
      * (exp(-i * k * r * t) * CONFLUENT_HYPG_1F1(-eta, zei)
          + SCAT_AMP(k, f0, d0, ac, eta) * TILDE_G(rho, eta) / r);
  return WF;
}

double PROB_LEDNICKY_COULOMB_WF(double k, double r, double t, double ScatLend,
                                double EffecRange, double ChargeRad) {

  const double eta = 1.0 / (k * ChargeRad) / FmToNu;
  const double rho = k * r * FmToNu;
  const double zei = rho * (1.0 + t);
  const double f0 = ScatLend * FmToNu;
  const double d0 = EffecRange * FmToNu;
  const double ac = ChargeRad * FmToNu;

  double Prob = 0.0;
  complex<double> WF;
  WF = pow(Ac(eta), 0.5)
      * (exp(-i * k * r * t) * CONFLUENT_HYPG_1F1(-eta, zei)
          + SCAT_AMP(k, f0, d0, ac, eta) * TILDE_G(rho, eta) / r);
  Prob = abs(conj(WF) * WF);
  return Prob;
}
double PROB_LEDNICKY_COULOMB_WF_DEF_2(double k, double r, double t,
                                      double ScatLend, double EffecRange,
                                      double ChargeRad) {
  const double rval = r * FmToNu;
  const double ac = ChargeRad * FmToNu;
  const double eta = 1.0 / (k * ac);
  const double rho = k * rval;
  const double zei = rho * (1.0 + t);
  const double f0 = ScatLend * FmToNu;
  const double d0 = EffecRange * FmToNu;

  double PDF = 0.0;
  PDF = Ac(eta)*(abs(conj(CONFLUENT_HYPG_1F1(-eta, zei))*CONFLUENT_HYPG_1F1(-eta, zei))
          + abs(
              conj(SCAT_AMP(k, f0, d0, ac, eta) * TILDE_G(rho, eta) / rval)
                  * SCAT_AMP(k, f0, d0, ac, eta) * TILDE_G(rho, eta) / rval)
          + real(
              conj(exp(-i * k * rval * t) * CONFLUENT_HYPG_1F1(-eta, zei))* SCAT_AMP(k, f0, d0, ac, eta) * TILDE_G(rho, eta) / rval)  + real(
                          conj(SCAT_AMP(k, f0, d0, ac, eta) * TILDE_G(rho, eta) / rval) *exp(-i * k * rval * t) * CONFLUENT_HYPG_1F1(-eta, zei)));

  return PDF;
}
double SOURCE(double r, double R) {
  const double Rad = R * FmToNu;  //source size
  const double rad = r * FmToNu;
  return 1.0 / pow((4.0 * M_PI * Rad * Rad), 3.0 / 2.0)
      * exp(-rad * rad / (4.0 * Rad * Rad));
}

double SOURCE_4PI(double r, double R) {
  const double Rad = R;  //source size
  const double rad = r;
  return 4.0 * M_PI * rad * rad / pow((4 * M_PI * Rad * Rad), 3.0 / 2.0)
      * exp(-rad * rad / (4.0 * Rad * Rad));
}

//Used root to perform 2D integration
double DIFF_CORRELATION(double *arg, double *par) {

  double tval = arg[0];
  double rval = arg[1];

  double kval = par[0];
  double ScatLend = par[1];
  double EffecRange = par[2];
  double ChargeRad = par[3];
  double Rval = par[4];
  return PROB_LEDNICKY_COULOMB_WF(kval, rval, tval, ScatLend, EffecRange,
                                  ChargeRad) * SOURCE(rval, Rval) * rval * rval
      * FmToNu * FmToNu * FmToNu * 2.0 * M_PI;

}

double DIFF_CORRELATION_SIMPSON(double kval, double rval, double tval,
                                double ScatLend, double EffecRange,
                                double ChargeRad, double Rval) {
  return PROB_LEDNICKY_COULOMB_WF_DEF_2(kval, rval, tval, ScatLend, EffecRange,
                                  ChargeRad) * SOURCE(rval, Rval) * rval * rval
      * FmToNu * FmToNu * FmToNu * 2.0 * M_PI;

}

double GET_CORRELATION(double k, double ScatLend, double EffecRange,
                       double ChargeRad, double SourceRad) {
  double CkVal = 0.0;
  TF2 *C_k = new TF2("C_k", DIFF_CORRELATION, -1, 1, 0.00, 100, 5);
  //C_k->SetParameter(0, k);
  // C_k->SetParameter(1, ScatLend);
  // C_k->SetParameter(2, EffecRange);
  // C_k->SetParameter(3, ChargeRad);
  // C_k->SetParameter(4, SourceRad);
  C_k->FixParameter(0, k);
  C_k->FixParameter(1, ScatLend);
  C_k->FixParameter(2, EffecRange);
  C_k->FixParameter(3, ChargeRad);
  C_k->FixParameter(4, SourceRad);

  CkVal = C_k->Integral(-1, 1, 0.0, 100.0);
  return CkVal;
}

double GET_CORRELATION_SIMPSON2D(double kval, double ScatLend,
                                 double EffecRange, double ChargeRad,
                                 double SourceRad) {
  float h = 0.1;
  float k = 0.10;
  float lx = -1.0;
  float ux = 1.0;
  float ly = 0.1;
  float uy = 30.0;
  //double DIFF_CORRELATION_SIMPSON(double kval, double rval, double tval, double ScatLend,
  //        double EffecRange, double ChargeRad,double Rval)
  int nx, ny;

  // z stores the table
  // ax[] stores the integral wrt y
  // for all x points considered
  float answer;

  // Calculating the numner of points
  // in x and y integral
  nx = (ux - lx) / h + 1;
  ny = (uy - ly) / k + 1;
  float z[nx + 4][ny + 5], ax[nx + 5];
  // Calculating the values of the table
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      z[i][j] = DIFF_CORRELATION_SIMPSON(kval, ly + j * k, lx + i * h, ScatLend,
                                         EffecRange, ChargeRad, SourceRad);
      //givenFunction(lx + i * h, ly + j * k);
    }
  }

  // Calculating the integral value
  // wrt y at each point for x
  for (int i = 0; i < nx; ++i) {
    ax[i] = 0;
    for (int j = 0; j < ny; ++j) {
      if (j == 0 || j == ny - 1)
        ax[i] += z[i][j];
      else if (j % 2 == 0)
        ax[i] += 2 * z[i][j];
      else
        ax[i] += 4 * z[i][j];
    }
    ax[i] *= (k / 3);
  }

  answer = 0;

  // Calculating the final integral value
  // using the integral obtained in the above step
  for (int i = 0; i < nx; ++i) {
    if (i == 0 || i == nx - 1)
      answer += ax[i];
    else if (i % 2 == 0)
      answer += 2 * ax[i];
    else
      answer += 4 * ax[i];
  }
  answer *= (h / 3);

  return answer;

}
