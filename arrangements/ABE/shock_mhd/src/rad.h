// Auxiliary variables
struct auxvars {
  // n: polytropic index, am4 = a_R m^4, fR: asymptoptic right-state F^x
  double n, kap_a, kap_s, am4, U1, U2, U3, U4, U5, fR;
  // current value of ux (for initial guess in the recovery of primitive variables)
  double ux;
  // These parameters are needed for shooting
  double xL, xR, rho0L, PL, uxL, EL;
  int nx;
};

//void GetInput (char *filename, auxvars &u, double &rho0, double &P, double &ux,
//               double &E, double &fL, double &fR, double &xL, double &xR, int &n,
//               int &flag, char *outfile);


void getUfromP(auxvars &u, double rho0, double P, double ux, double E, double f);
void computeRightState(auxvars &u, double &rho0r, double &Pr, double &uxr, double &Er);
void odeIntegrate(double x1, double x2, int n, auxvars &u);
void odeIntegrateOut(double x1, double x2, int n, auxvars &u, double *xsave, double *rhosave,
		  double *Psave, double *uxsave, double *Esave, double *Fxsave, double *tausave);
void outputShootingData(double fL, double fR, auxvars u, double *xsave, double *rhosave,
			double *Psave, double *uxsave, double *Esave, double *Fxsave, double *tausave);
void uxFromU(double ux, double &f, double &dfdux, double &fnorm, auxvars &u);
void uxFromUandLTE(double ux, double &f, double &dfdux, double &fnorm, auxvars &u);
double newt(double xg, auxvars &u, 
            void (*f)(double, double &, double &, double &, auxvars &), 
            double abserr, double relerr);
void getPfromU(auxvars &u, double &rho0, double &P, double &ux, double &u0, 
               double &E, double &Fx, double &F0);
void shootingFun(double fL, double fR, double *f, double *fnorm, auxvars u);
void shootingFJ(double *x, double *f, double *fnorm, double *fjac, auxvars u);
void odeRHS(double *U, double w, auxvars u, double *derivs);
void newt2d(double *x, double abserr, double relerr, auxvars &u,
            void (*FJac)(double *, double *, double *, double *, auxvars));


template <typename BidirectionalIterator, typename T>
BidirectionalIterator getClosest(BidirectionalIterator first,
				   BidirectionalIterator last,
				   const T & value)
{
  BidirectionalIterator before = std::lower_bound(first, last, value);

  if (before == first) return first;
  if (before == last)  return --last; // iterator must be bidirectional                                                                                                                                   

  BidirectionalIterator after = before;
  --before;

  return (*after - value) < (value - *before) ? after : before;
};


template <typename BidirectionalIterator, typename T>
  std::size_t getClosestIndex(BidirectionalIterator first,
			      BidirectionalIterator last,
			      const T & value)
{
  return std::distance(first, getClosest(first, last, value));
};
