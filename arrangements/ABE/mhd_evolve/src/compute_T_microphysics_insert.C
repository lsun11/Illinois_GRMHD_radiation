#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sys/time.h>
#include "cctk.h"

#include "primitives_solver_header.h"
#define SQR(x) ((x) * (x))
double max_insert(double val1,double val2);
double max_insert(double val1, double val2) {
  if(val1>val2) return val1;
  return val2;
}
double min_insert(double val1,double val2);
double min_insert(double val1, double val2) {
  if(val1<val2) return val1;
  return val2;
}

void compute_P_T_microphys_insert(double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b);

double find_T_insert(double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
void solve_T_insert(double &Ta, double &Tb, double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);



extern "C" void CCTK_FCALL CCTK_FNAME(compute_P_T_microphys_insert)
  (double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b);


void compute_P_T_microphys_insert(double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b){

  //|||| Now include thermal and nuclear parts   ||||//                                                          
  // Note that the code unit quantities are powers of km, we need to convert them to powers of cm first!!!///                                                                                           
  bool recom = false;
  double eps_thermal = eps - eps_cold;

  if (rho_b <1.0e-9 && rho_b > 1.45e-13 &&eps_thermal < 1.0e-9)
    {                                 
      eps_thermal = eps_cold/rho_b*149.0;                                                                  
      //Need to recompute eps_tot if eps_thermal is revalued.                                                          
      eps = eps_thermal + eps_cold;                                                        
    }                                              

  double eng_thermal = rho_b*eps_thermal;
  double eng_thermal_cgs = fasterpow_prim(c_cgs,4.0) / G_cgs * eng_thermal * fasterpow_prim(k2c, -2.0);
  double rho_b_cgs = fasterpow_prim(c_cgs,2.0) / G_cgs * rho_b * fasterpow_prim(k2c, -2.0);
  double n_nucleon = rho_b_cgs/m_n_cgs;
  double N_nu = 3.0;
  double T_fluid_cgs;
  double T_scal= 1.0e0;
  //double T_scal = 1.0; 
  //  printf("Begin compute_P_T_microphys_insert, T_fluid!!! T_fluid=%e \n", T_fluid);                                                      

  double T_fluid_cgs_floor = 1.0e1;

  //check if the thermal part if positive:                                                                                                                                                                                                                                             
  if (eps_thermal <= 0.0){
    //printf("eps_thermal = %e \n", eps_thermal);                                                                                                                                      
    T_fluid_cgs = T_fluid_cgs_floor;
  }
  else{
    // If we use this part rad_const_cgs is set to be 7.5646e-15 erg cm^-3 K^-4            
    // Solve for T_fluid                                                                                                                                                                                                                                                               
    // As 1.5nkBT >> 11/4*aT^4 in most cases, make the initial guess of Tcgs = eng_thermal_cgs/(1.5*n*kB) and then solve for Tcgs.                                                                                                                                                     
    // If T_fluid is postive (calculated in the last iteration), use it as the guess, otherwise  use Tcgs = eng_thermal_cgs/(1.5*n*kB).                                                                                                                                                
    // T guess needs to be adjusted for high temperature, low density region (rare).                                                                                                                                                                                                   
    //double threshold = fasterpow_prim(eng_thermal_cgs/(kb_cgs*0.21512), 0.75);                                                                                                                                                                                                       
    double threshold = fasterpow_prim(eng_thermal_cgs,0.75)*0.8585*fasterpow_prim(rad_const_cgs,0.25)/kb_cgs;
    double T_fluid_cgs_guess;
    if (n_nucleon < threshold){
      //T_fluid_cgs_guess= eng_thermal_cgs/(1.5 * kb_cgs * n_nucleon);                                                                                                                                                                                                                 
      T_fluid_cgs = fasterpow_prim(4.0*eng_thermal_cgs/rad_const_cgs/11.0, 0.25);
      T_fluid_cgs_guess = fasterpow_prim(4.0*eng_thermal_cgs/rad_const_cgs/11.0, 0.25);
    }
    else{
      //T_fluid_cgs = fasterpow_prim(4.0*eng_thermal_cgs/rad_const_cgs/11.0, 0.25);                                                                                                                                                                                                    
      T_fluid_cgs = eng_thermal_cgs/(1.5 * kb_cgs * n_nucleon);
      T_fluid_cgs_guess= eng_thermal_cgs/(1.5 * kb_cgs * n_nucleon);
    }

    //The final T should be smaller than T_guess since T_guess only consists of one of the two parts of eng_thermal.                                                                                                                                                                   
    // Therefore, no matter we set F1(T_guess) or F2(T_guess) to eng_thermal, F1(T_guess) + F2(T_guess) > eng_thermal and thus real T is smaller.                                                                                                                                      
    double Ta = T_fluid_cgs*0.8;
    double Tb = T_fluid_cgs*1.3;

    double rad_constl = rad_const_cgs;
    //  if (T_fluid_cgs == 0.0) printf("inside compute_P_T_microphys_insert, T_fluid_cgs is 0!!! T_fluid=%e \n", T_fluid);                                                                                                                                                                    
    solve_T_insert(Ta, Tb, T_fluid_cgs, eng_thermal_cgs, rad_constl, n_nucleon, N_nu, T_scal);

    // WE DO A CHECK AFTER GETTING THE TEMPERATURE!!!!                                                                                                                                                                                                                                 
    double E1 = (11.0/4.0 * rad_const_cgs)* fasterpow_prim(T_fluid_cgs, 4.0);
    double E2 = 1.5*n_nucleon*T_fluid_cgs*kb_cgs;
    double E_ratio = (E1 + E2)/eng_thermal_cgs;
    if (T_fluid_cgs > T_fluid_cgs_floor && E_ratio > 1.05){
      printf("after solve_T_insert, E_ratio too large, is %.16e, E1=%.16e, E2=%.16e, eng_thermal_cgs=%.16e, T_fluid_cgs=%.16e, T_fluid_cgs_guess=%.16e, n_nucleon=%.16e \n", E_ratio, E1, E2, eng_thermal_cgs, T_fluid_cgs, T_fluid_cgs_guess, n_nucleon);
    }

    //set T_fluid_cgs floor (10^6 K)                                                                                                                                                                                                                                
    if (T_fluid_cgs < T_fluid_cgs_floor){
      T_fluid_cgs = T_fluid_cgs_floor;
    }

  }

  //  Here fasterpow_prim(c2k, 1.0) means converting from cm^1 to km^1 (Temperature scales as [L]);                                                                                                                                                                                    
  T_fluid = T_fluid_cgs * kb_cgs * G_cgs / fasterpow_prim(c_cgs,4.0) * fasterpow_prim(c2k, 1.0);

  //  double Fac2 = min_insert_val(1.0,  2.0*fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*m_n_cgs/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(n_nucleon, -2.0/3.0)*kb_cgs*T_fluid_cgs  );                                                                                                             
  double Fac2 = 1.0;
  double P_nucl_cgs = (5.0/3.0-1.0)*(1.5*n_nucleon*kb_cgs*T_fluid_cgs*Fac2);
  double P_nucl = fasterpow_prim(c_cgs, -4.0) * fasterpow_prim(G_cgs, 1.0) * P_nucl_cgs * fasterpow_prim(c2k, -2.0);

  //double P_rad_m_cgs = (3.0 + 7.0 * N_nu/8.0) * rad_const_cgs * fasterpow_prim(T_fluid_cgs, 4.0)/3.0;                                                                                                                                                                                
  double P_rad_ph_cgs = (11.0/12.0)*rad_const_cgs * fasterpow_prim(T_fluid_cgs, 4.0);
  double P_rad_ph = fasterpow_prim(c_cgs, -4.0) * fasterpow_prim(G_cgs, 1.0) * P_rad_ph_cgs * fasterpow_prim(c2k, -2.0);

  //  double P_thermal =  E_rad/3.0 + P_nucl;                                                                                                                                                                                                                                          
  double P_thermal = P_rad_ph + P_nucl;
  P = P_cold + P_thermal;


}


extern "C" void CCTK_FCALL CCTK_FNAME(compute_P_T_microphys_insert)
  (double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b) {
  compute_P_T_microphys_insert(P,T_fluid,P_cold,eps,eps_cold,rho_b);
}



double find_T_insert(double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal){
  //  double Fac = min_insert_val(1.0,  2.0*fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*m_n_cgs/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(n_nucleon, -2.0/3.0)*kb_cgs*T );               
  double Fac = 1.0;
  double diff = (11.0/4.0 * rad_const)* fasterpow_prim(T*T_scal, 4.0) + 1.5*n_nucleon*T*T_scal*kb_cgs*Fac - eng_thermal_cgs;
  //  if (diff == 0) {printf("Inside find_T_insert, diff is 0!!!, T=%.16e, (11.0/4.0 * rad_const)* fasterpow_prim(T*T_scal, 4.0)=%.16e, 1.5*n_nucleon*T*T_scal*kb_cgs*Fac=%.16e,  eng_thermal_cgs=%.16e \n", T, (11.0/4.0 * rad_const)* fasterpow_prim(T*T_scal, 4.0), 1.5*n_nucleon*T*T_scal*kb_cgs, eng_thermal_cgs);}                                                                                                                                                                                                                                                                                
return diff/eng_thermal_cgs;
}


void solve_T_insert(double &Ta, double &Tb,
             double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal)
{
  double fa, fb, fc;
  int bisec = 1;
  double T_init =T;
  double tiny = 1.0;
  int iter=0;
  int iter1=0;

  //If T_cgs_guess is smaller than 1K, set T_cgs to 1K, don't do the bisection!                          
  if (T_init < tiny){
    T = tiny;
    bisec = 0;
  }
  else{
    Ta /=T_scal;
    Tb /=T_scal;

    // Obtain the result fa = F1(T_guess_a) + F2(T_guess_a) - eng_thermal_cgs & fb = F1(T_guess_b) + F2(T_guess_b) - eng_thermal_cgs             
    // In principle, fb >  F1(T_guess) + F2(T_guess) - eng_thermal_cgs > 0, and we want fa < 0                                                                                          
    fa = find_T_insert(Ta, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
    fb = find_T_insert(Tb, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);

    while (fa *fb >= 0)
      {
        if (iter > 20) {//We fail to converge to find a proper Ta -- Tb range!!!!                                                   
          T = T_init;
          bisec = 0;
          printf("break loop: Ta = %e, Tb =%e, Ta_init= %e, Tb_init=%e, eng_thermal_cgs = %e, fa=%e, fb=%e \n", Ta, Tb, T_init*0.8, T_init*1.05, eng_thermal_cgs, fa, fb);
          break;
        }
        Ta *= 0.5;
        Tb *= 1.5;
        fa = find_T_insert(Ta, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
        fb = find_T_insert(Tb, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
        iter += 1;
      }
    if (bisec){
      T = Ta;
      while ( fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
        {
          // Find middle point           
          T = (Ta+Tb)/2;
          fc = find_T_insert(T, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
          // Check if middle point is root                                                                                                   
          if (fc == 0.0)
            break;
          // Decide the side to repeat the steps                                                                                             
          else if (fc*fa < 0)
            Tb = T;
          else
            Ta = T;

          iter1 +=1;

          if (iter1 > 50){
            // Failed to converge after 50 iteration, do 50 more!!!!!!!               
            //   printf("inside solve_T, over iteration, fabs(Tb-Ta)/(Tb+Ta)=%e \n", fabs(Tb-Ta)/(Tb+Ta));          
            int iter2 = 0;
            while (fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
              {
                // Find middle point                                                                                                                                
                T = (Ta+Tb)/2;
                fc = find_T_insert(T, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
                // Check if middle point is root          
                if (fc == 0.0)
                  break;
                // Decide the side to repeat the steps                                                                                                                            
                else if (fc*fa < 0)
                  Tb = T;
                else
                  Ta = T;

                iter2 +=1;
              }
            if (iter2 > 50){
              printf("inside solve_T, over iteration AGAIN, fabs(Tb-Ta)/(Tb+Ta)=%e \n", fabs(Tb-Ta)/(Tb+Ta));
              break;
            }
          }

        }

      T *=T_scal;
      Ta *=T_scal;
      Tb *=T_scal;

    } // End of bisection 
  } // End of T > T_tiny                                                                                                                                                                                                                                                     
}
