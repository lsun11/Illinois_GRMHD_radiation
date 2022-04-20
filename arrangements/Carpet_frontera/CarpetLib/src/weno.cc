// 
// WENO Interpolation 
//
// Interpolate a function f at point x using stencils x_{j-2}, x_{j-1}, x_j, 
// x_{j+1} and x_{j+2}. The spacings in x, h_k = x_k-x_{k-1} (k=j-1,j,j+1,j+2), 
// are assumed to be a constant (independent of k). 
//
// In function weno_5pt, the value of f at point x is interpolated using 
// the point values of f at x_{j+k} (k=-2,-1,...,2). The interpolation 
// is 5th order accurate if f is at least C^4 in (x_{j-2},x_{j+2}), but reduces 
// to 3rd order if f is C^n with n<4.
//
// In function weno_5avg1, fbarh2(x), the average value of f at x
// fbarh2(x) = \int_{x-h/4}^{x+h/4} f(x') dx' , 
// is computed using values of fbar(x),
// fbarh(x) = \int_{x-h/2}^{x+h/2} f(x') dx' ,
// at x=x_{j+k} (k=-2,-1,0,1,2). The interpolation
// is 5th order accurate if f is at least C^4 in (x_{j-2},x_{j+2}), but reduces
// to 3rd order if f is C^n with n<4.
//
// Below, f is an array of size 5, storing 
// f[0] = f(x_{j-2}), f[1] = f(x_{j-1}), f[2] = f(x_j), 
// f[3] = f(x_{j+1}), f[4] = f(x_{j+2}); 
// xj = x_j.
//

#include "math.h"

#define SQR(x) ((x) * (x))

double weno_5pt(double *f, double x, double xj, double h) {

  double xbar = (x-xj)/h;

  // NOTE: f[0] = f_{j-2}, f[1]=f_{j-1}, f[2]=f_j, f[3]=f_{j+1}, f[4]=f_{j+2}
  double pj = f[2] + 0.5*(3.0*f[2]-4.0*f[1]+f[0])*xbar
                   + 0.5*(f[2]-2.0*f[1]+f[0])*SQR(xbar);
  double pj1 = f[2] + 0.5*(f[3]-f[1])*xbar + 0.5*(f[3]-2.0*f[2]+f[1])*SQR(xbar);
  double pj2 = f[2] - 0.5*(f[4]-4.0*f[3]+3.0*f[2])*xbar
                    + 0.5*(f[4]-2.0*f[3]+f[2])*SQR(xbar);

  double C0 = (xbar-1.0)*(xbar-2.0)/12.0;
  double C1 = (4.0-xbar*xbar)/6.0;
  double C2 = (xbar+1.0)*(xbar+2.0)/12.0;

  double IS0 = 13.0/12.0*SQR(f[0]-2.0*f[1]+f[2]) + 0.25*SQR(f[0]-4.0*f[1]+3.0*f[2]);
  double IS1 = 13.0/12.0*SQR(f[1]-2.0*f[2]+f[3]) + 0.25*SQR(f[1]-f[3]);
  double IS2 = 13.0/12.0*SQR(f[2]-2.0*f[3]+f[4]) + 0.25*SQR(3.0*f[2]-4.0*f[3]+f[4]);
  // Normalize ISk
  double ISk_norm = (fabs(IS0)+fabs(IS1)+fabs(IS2))/3.0;
  if (ISk_norm > 0) {
    double inv_ISk_norm = 1.0/ISk_norm;
    IS0 *= inv_ISk_norm; IS1 *= inv_ISk_norm; IS2 *= inv_ISk_norm;
  } else {
    IS0 = 1.0; IS1=1.0; IS2 = 1.0;
  }

  double eps = 1.e-6;
  double alpha0 = C0/SQR(eps + IS0);
  double alpha1 = C1/SQR(eps + IS1);
  double alpha2 = C2/SQR(eps + IS2);
  double inv_alpha_sum = 1.0/(alpha0 + alpha1 + alpha2);
  double omega0 = alpha0 * inv_alpha_sum;
  double omega1 = alpha1 * inv_alpha_sum;
  double omega2 = alpha2 * inv_alpha_sum;

  double f_int = omega0*pj + omega1*pj1 + omega2*pj2;

  return f_int;
}


double weno_5avg1(double *f, double x, double xj, double h) {

  double xbar = (x-xj)/h;
  double xbar2 = SQR(xbar);
  double xbar3 = xbar*SQR(xbar);
  double xbar4 = SQR(xbar2);

  // NOTE: f[0] = f_{j-2}, f[1]=f_{j-1}, f[2]=f_j, f[3]=f_{j+1}, f[4]=f_{j+2}
  double pj = (31.0*f[2]+2.0*f[1]-f[0])/32.0
             + 0.5*(3.0*f[2]-4.0*f[1]+f[0])*xbar
             + 0.5*(f[2]-2.0*f[1]+f[0])*xbar2;
  double pj1 = (34.0*f[2]-f[3]-f[1])/32.0
              + 0.5*(f[3]-f[1])*xbar
              + 0.5*(f[3]+f[1]-2.0*f[2])*xbar2;
  double pj2 = (31.0*f[2]+2.0*f[3]-f[4])/32.0
              + 0.5*(4.0*f[3]-3.0*f[2]-f[4])*xbar
              + 0.5*(f[2]-2.0*f[3]+f[4])*xbar2;

  double C0 = (21.0+608.0*xbar-352.0*xbar2-512.0*xbar3+256.0*xbar4)/
                (192.0* (16.0*xbar2 + 16.0*xbar - 1.0) );
  double C2 = (21.0-608.0*xbar-352.0*xbar2+512.0*xbar3+256.0*xbar4)/
                (192.0* (16.0*xbar2 - 16.0*xbar - 1.0) );
  double C1 = 1.0-C0-C2;

  double IS0 = 13.0/12.0*SQR(f[0]-2.0*f[1]+f[2]) + 0.25*SQR(f[0]-4.0*f[1]+3.0*f[2]);
  double IS1 = 13.0/12.0*SQR(f[1]-2.0*f[2]+f[3]) + 0.25*SQR(f[1]-f[3]);
  double IS2 = 13.0/12.0*SQR(f[2]-2.0*f[3]+f[4]) + 0.25*SQR(3.0*f[2]-4.0*f[3]+f[4]);
  // Normalize ISk
  double ISk_norm = (fabs(IS0)+fabs(IS1)+fabs(IS2))/3.0;
  if (ISk_norm > 0) {
    double inv_ISk_norm = 1.0/ISk_norm;
    IS0 *= inv_ISk_norm; IS1 *= inv_ISk_norm; IS2 *= inv_ISk_norm;
  } else {
    IS0 = 1.0; IS1=1.0; IS2 = 1.0;
  }

  double eps = 1.e-6;
  double alpha0 = C0/SQR(eps + IS0);
  double alpha1 = C1/SQR(eps + IS1);
  double alpha2 = C2/SQR(eps + IS2);
  double inv_alpha_sum = 1.0/(alpha0 + alpha1 + alpha2);
  double omega0 = alpha0 * inv_alpha_sum;
  double omega1 = alpha1 * inv_alpha_sum;
  double omega2 = alpha2 * inv_alpha_sum;

  double f_int = omega0*pj + omega1*pj1 + omega2*pj2;

  return f_int;
}
