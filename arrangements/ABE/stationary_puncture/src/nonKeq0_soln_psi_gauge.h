#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include "nr.h"

//==========================================================
// class  prclSoln()
//       (abbreviation of precollapsedSolution) 
//
// Last modified Nov 14 2007
//
// *Class uses Numerical Recipes Vec_DP class
//
// *Class constructor requires the name of a compatible data
//  file only
// 
// *Class methods calculate and return the value of the
//  conformal factor, lapse and shift for an input radial
//  coordinate greater than that of the event horizon 
//
// *File 'prclSoln26.txt' (attached) may be used to find 
//  these values for the equilibrium solution of a single 
//  black hole in vacuum with 1+log slicing; this file
//  was produced using spectral methods with 26 collocation
//  points, and is believed to be accurate to order   
//  less than E-10.  r_{event horizon} = 0.8304040892.
//
// *To store the value of the conformal factor for radial
//  coordinate r in variable psi:
//
//    bool legal;
//    double psi = confFact(r,legal);
//
//  'legal' now holds true if r>r_{event horizon}, false 
//  otherwise.  'psi' holds the value of the conformal
//  factor at r if legal is true, else it holds the error
//  value -117.  
//
// *The lapse and shift values can be attained similarly,
//  using the functions
//    
//    lapse(double r, bool & legal)
//    shift(double r, bool & legal)
//
//==========================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// trans()
//
// Returns the s coordinate of input r, and 
// returns by reference tag 't' recording which
// domain the r value falls into
// (t=0 : within the event horizon)
// (t=1 : between the event horizon and shell)
// (t=2 : outside of the shell)
//
// @ r = radial coordinate to transform
// @ t = tag variable
// @ return = related s coordinate
//          = tag 't' assigned by reference
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double trans(double r, double r_EH, double r_outr_bdry,int & t){
  
  if(r < r_EH){                           // Inside event horizon
    t = 0;
    return -117;
  }
  
  if(r<=r_outr_bdry){                            // Between EH and shell
    t=1;
    return (2.0*r_EH*r_outr_bdry/r-r_EH-r_outr_bdry)/(r_EH-r_outr_bdry);
  }
  
  t=2;                                  // Outisde of shell
  return -2.0*r_outr_bdry/r+1.0;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// makeT()
//
// Construct vector of Chebychev polynomial
// values T for input s
//
// @ s = s coordinate to construct T with
// @ return = T constructed
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void makeT(double s, int N, Vec_DP * TA){
  Vec_DP &T         = *TA;              
  for(int n=0;n<=N;n++){
    T[n]     = cos(n*acos(s));
  }  
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// vecMul()
// 
// Multiplies two vectors, returns their 
// scalar product
//
// @ a = first vector to be multiplied
// @ b = second vector to be multiplied
// @ return = vector product of a and b
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double vecMul(Vec_DP a, Vec_DP b){
  int aS = a.size();
  int bS = b.size();
  if(aS!=bS){
    cout << "Cannot multiply vectors; not of the same size" << endl;
    return -117;
  } 
  double c=0.0;
  for(int i=0;i<aS;i++){
    c+=a[i]*b[i];
  }
  return c;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// confFact(), lapse(), shift()
//
// Return value of conformal factor, lapse,
// and/or shift for input r coordinate.
// Returns by reference whether or not
// the r is 'legal' (ie falls outside of the
// event horizon), and error value -1337 for
// an illegal r.
//
// @ r = radial coordinate to find value for
// @ legal = variable to store legality of r
// @ return = value of field variable at r
//          = legality tag 'legal' set by 
//            reference
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double output_1var_1radius(double r, int N, double r_EH, double r_outr_bdry, Vec_DP * TA, Vec_DP *varIA, Vec_DP *varOA, bool & legal){
  int tag;
  legal = true;
  double s = trans(r, r_EH, r_outr_bdry, tag);
  if(tag==0){                 // If r falls within horizon,
    legal = false;            // return error values
    return -117;
  }
  makeT(s,N,TA);
  if(tag==1){ return vecMul(*TA,*varIA); }
  return vecMul(*TA,*varOA);
}







