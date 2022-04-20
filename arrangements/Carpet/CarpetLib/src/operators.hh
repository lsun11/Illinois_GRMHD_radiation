#ifndef OPERATORS_HH
#define OPERATORS_HH

// Transport (i.e., prolongation and restriction) operator types

enum operator_type
{
  op_error,                     // illegal operator type
  op_none,                      // do not transport
  op_sync,                      // transport only on the same level
                                // (error if called between levels)
  op_restrict,                  // restrict only, do not prolongate
  op_copy,                      // use simple copying for prolongation
                                // (needs only one time level)
  op_Lagrange,                  // Lagrange interpolation (standard)
  op_ENO,                       // use ENO stencils (for hydro)
  op_WENO,                      // use WENO stencils (for hydro)
  op_Lagrange_monotone,         // monotone Lagrange interpolation (for hydro)
  op_Lag3,
  op_STAGGER011,                // use stagger grid of type 011 (store gf at i,j+1/2,k+1/2)
  op_STAGGER101,                // use stagger grid of type 011 (store gf at i+1/2,j,k+1/2)
  op_STAGGER110,                // use stagger grid of type 011 (store gf at i+1/2,j+1/2,k)
  op_STAGGER111,                // use stagger grid of type 011 (store gf at i+1/2,j+1/2,k+1/2)
  op_FAKE                       // Keep check of refinement boundaries
};

#endif // OPERATORS_HH
