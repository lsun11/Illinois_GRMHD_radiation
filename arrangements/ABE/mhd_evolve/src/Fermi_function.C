#include "math.h"
#define SQR(x) ((x) * (x))
#define eta0 1.0e-3

double Fermi3p(double &eta){
  double K31=4.9348;
  double K32=11.3644;
  double K33=1.9039;

  double X = 1.0/eta;
  return (0.25 + K31*SQR(X) + K32*SQR(X)*SQR(X))/(1.0 - exp(-K33*eta));
}


double Fermi3m(double &eta){
  double A31=0.0559;
  double A32=0.9069;

  return 6.0/(1.0 + A31*exp(A32*eta));
}

double Fermi3(double &eta){
  if (eta > eta0){
    return pow(eta,4.0)*Fermi3p(eta);
  }
  else{
    return exp(eta)*Fermi3m(eta);
  }
}




double Fermi4p(double &eta){
  double K41=6.5797;
  double K42=45.4576;
  double K43=1.9484;

  double X = 1.0/eta;
  return (0.2 + K41*SQR(X) + K42*SQR(X)*SQR(X))/(1.0 - exp(-K43*eta));
}


double Fermi4m(double &eta){
  double A41=0.0287;
  double A42=0.9257;

  return 24.0/(1.0 + A41*exp(A42*eta));
}

double Fermi4(double &eta){
  if (eta > eta0){
    return pow(eta,5.0)*Fermi4p(eta);
  }
  else{
    return exp(eta)*Fermi4m(eta);
  }
}


double Fermi5p(double &eta){
  double K51=8.2247;
  double K52=113.6439;
  double K53=236.5323;
  double K54=1.9727;

  double X = 1.0/eta;
  return (1.0/6.0 + K51*SQR(X) + K52*SQR(X)*SQR(X) + K53*SQR(X)*SQR(X)*SQR(X))/(1.0 + exp(-K54*eta));
}


double Fermi5m(double &eta){
  double A51=0.0147;
  double A52=0.9431;

  return 120.0/(1.0 + A51*exp(A52*eta));
}

double Fermi5(double &eta){
  if (eta > eta0){
    return pow(eta,6.0)*Fermi5p(eta);
  }
  else{
    return exp(eta)*Fermi5m(eta);
  }
}
