#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "primitives_solver_header.h"

/*
extern "C" void CCTK_FCALL compute_t_fluid_cpp_
  (double &X, double &Y, double &Z, int &num_CO,
   double &xNS1, double &yNS1, double &xNS2, double &yNS2, double &rhob, double &rho_star, double &P,
   double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor,
   double &rho_b_atm, double &T_fluid);
*/

void compute_T_fluid_cpp(double &X, double &Y, double &Z, int&num_CO,
                         double &xNS1, double &yNS1, double &xNS2, double &yNS2, double &rhob, double &rho_star, double &P,
                         double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor,
                         double &rho_b_atm, double &T_fluid);


void compute_T_fluid_cpp(double &X, double &Y, double &Z, int&num_CO,
                         double &xNS1, double &yNS1, double &xNS2, double &yNS2, double &rhob, double &rho_star, double &P,
                         double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor,
			 double &rho_b_atm, double &T_fluid){
  
        double xc1 = X-xNS1;
        double yc1 = Y-yNS1;
        double zc1 = Z;
        double rad1 = sqrt(xc1*xc1+yc1*yc1+zc1*zc1);

        double xc2 = X-xNS2;
        double yc2 = Y-yNS2;
        double zc2 = Z;
        double rad2 = sqrt(xc2*xc2+yc2*yc2+zc2*zc2);


        if (num_CO == 1){ 
         if (rad1 <= rad_T_cutoff){
            T_fluid = P*M_B/rho_star ;
          }
          else{
            T_fluid = P*M_B/rho_star * pow( rad_T_fac/rad1, rad_T_pow);
          }
        }
	else
	  {
	    double rad = fmin(rad1, rad2);
	    if (rad <= rad_T_cutoff){
	      T_fluid = P*M_B/rho_star ;
	    }
	    else{
	      T_fluid = P*M_B/rho_star * (pow( rad_T_fac/rad1, rad_T_pow) +  pow( rad_T_fac/rad2, rad_T_pow));
	    }
	  }
        //set floor temperature                                                                                                      
        if (T_fluid < rad_T_floor){
	  T_fluid = rad_T_floor;
        }	

}


/*
extern "C" void CCTK_FCALL compute_t_fluid_cpp_
  (double &X, double &Y, double &Z, int &num_CO,
   double &xNS1, double &yNS1, double &xNS2, double &yNS2, double &rhob, double &rho_star, double &P,
   double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor, 
   double &rho_b_atm, double &T_fluid)
{
  compute_T_fluid_cpp(X,Y,Z, num_CO, xNS1, yNS1, xNS2, yNS2, rhob, rho_star, P,
 		      M_B, rad_T_fac, rad_T_cutoff, rad_T_pow, rad_T_floor, rho_b_atm, T_fluid);

}
*/


//A much simpler function for shocktest.
//extern "C" void CCTK_FCALL CCTK_FNAME(compute_T_fluid_shock_cpp)                                                                                                                                            (double &rho_star, double &P,double &M_B,double &rad_T_floor, double &T_fluid);   



void compute_T_fluid_shock_cpp(double &rho_star, double &P,double &M_B, 
			       double &rad_T_floor,double &T_fluid);


void compute_T_fluid_shock_cpp(double &rho_star, double &P,double &M_B,
			       double &rad_T_floor,double &T_fluid){
  T_fluid = P*M_B/rho_star;
  //set floor temperature                                                                                                                     
   if (T_fluid < rad_T_floor)
   {
    T_fluid = rad_T_floor;
   }
}


//Use E_rad to compute temperature for LTE.

void compute_T_fluid_OS(double &E_rad, double &rad_const, double &T_fluid);

void compute_T_fluid_OS(double &E_rad, double &rad_const, double &T_fluid){

  T_fluid = pow((fabs(E_rad)/rad_const), 0.25);

}


void compute_T_fluid_shock_cpp2(double &P, double &rho, double &T_fluid);

void compute_T_fluid_shock_cpp2(double &P, double &rho, double &T_fluid){
  T_fluid = P/rho;
}
