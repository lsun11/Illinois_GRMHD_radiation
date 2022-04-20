//----------------------------------------------------------------------------
//
// $Id: $
//
//----------------------------------------------------------------------------
//
// Contains methods for integrating the OV equations
//
//----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <string.h>
#include <fstream.h>
#include "Symmetry.h"

//#define PI acos(-1.0) // <-- BAD.  Use M_PI.
#define N_INT 1000

void Setup_Hybrid_EOS(double K_poly, double n, char * eos_file,int neos,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab);

static char *rcsid = "$Meow...$";
CCTK_FILEVERSION(disk_initialdata_readfiles)

extern "C" void disk_initialdata_readfiles(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Setup hybrid EOS table
  double n = 3.0;
  Setup_Hybrid_EOS(K_poly, n, "eos.tab",neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
  
}


// void Setup_Hybrid_EOS(double K_poly, double n, char * eos_file,
// 		      int neos,double *rho_tab,double *P_tab,
// 		      double *eps_tab,double *k_tab,double *gamma_tab) {
//   ifstream eosdat;
//   eosdat.open(eos_file);
//   if (!eosdat) {
//     neos = 1;
//     cout << "Can't open " << eos_file << endl;
//     cout << "Assuming polytropic EOS..." << endl;
//     cout << "Polytropic index = " << n << endl;
//     cout << "Polytropic constant K = " << K_poly << endl;
//     //rho_tab = new double[1]; P_tab = new double[1]; eps_tab = new double[1];
//     //gamma_tab = new double[2]; k_tab = new double[2];
//     rho_tab[0] = 0.03;
//     P_tab[0] = K_poly * pow(rho_tab[0],1.0+1.0/n);
//     eps_tab[0] = n*P_tab[0]/rho_tab[0];
//     gamma_tab[0] = 1.0 + 1.0/n; k_tab[0] = K_poly;
//     gamma_tab[1] = 1.0 + 1.0/n; k_tab[1] = K_poly;
//   } else {
//     char buf[100],c;
//     eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> neos;
//     //rho_tab = new double[neos]; P_tab = new double[neos];
//     //eps_tab = new double[neos]; gamma_tab = new double[neos+1]; k_tab = new double[neos+1];
//     eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> gamma_tab[0];
//     eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> gamma_tab[neos];
//     for (int i=0; i<neos; i++) {
//         eosdat >> rho_tab[i] >> P_tab[i];
//     }
//     k_tab[0] = P_tab[0]/pow(rho_tab[0],gamma_tab[0]);
//     k_tab[neos] = P_tab[neos-1]/pow(rho_tab[neos-1],gamma_tab[neos]);
//     for (int i=0; i<neos-2; i++) {
//         gamma_tab[i+1] = log(P_tab[i+1]/P_tab[i])/log(rho_tab[i+1]/rho_tab[i]);
//         k_tab[i+1] = P_tab[i]/pow(rho_tab[i],gamma_tab[i+1]);
//     } 
//     eps_tab[0] = P_tab[0]/rho_tab[0]/(gamma_tab[0]-1.0);
//     for (int i=1; i<neos; i++) {
//         eps_tab[i] = eps_tab[i-1] +
//           (P_tab[i]/rho_tab[i] - P_tab[i-1]/rho_tab[i-1])/(gamma_tab[i]-1.0);
//     }
//   } 
//   eosdat.close();

// }

