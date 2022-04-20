 /*@@
   @file      Complex.c
   @date      Tue Dec 14 12:09:43 1999
   @author    Tom Goodale
   @desc
              Complex variable stuff
   @enddesc
   @version   $Id: Complex.c 4418 2007-04-27 19:41:24Z schnetter $
 @@*/

#include <math.h>

#include "cctk_Flesh.h"

#ifndef DEFINE_CCTK_COMPLEX_INLINE_FUNCTIONS
#  define DEFINE_CCTK_COMPLEX_EXTERN_FUNCTIONS
#  include "cctk_Complex.h"
#  undef DEFINE_CCTK_COMPLEX_EXTERN_FUNCTIONS
#endif


#ifndef DEFINE_CCTK_COMPLEX_INLINE_FUNCTIONS

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_Complex_c);

#endif


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTK_Cmplx
   @date       Tue Dec 14 12:16:01 1999
   @author     Tom Goodale
   @desc
               Turns two reals into a complex number
   @enddesc
 
   @var        Re
   @vdesc      Real part
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        Im
   @vdesc      Imaginary part
   @vtype      CCTK_REAL
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX(CCTK_Cmplx, cctk_real, cctk_complex)                \
cctk_complex CCTK_Cmplx (cctk_real Re, cctk_real Im)                          \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = Re;                                                             \
  result.Im = Im;                                                             \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxReal
   @date       Tue Dec 14 12:22:41 1999
   @author     Tom Goodale
   @desc
               Returns the real part of a complex number.
   @enddesc

   @var        complex_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The real part
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_REAL(CCTK_Cmplx, cctk_real, cctk_complex)           \
cctk_real CCTK_Cmplx##Real (cctk_complex complex_number)                      \
{                                                                             \
  return (complex_number.Re);                                                 \
}


 /*@@
   @routine    CCTK_CmplxImag
   @date       Tue Dec 14 12:22:41 1999
   @author     Tom Goodale
   @desc
               Returns the imaginary part of a complex number.
   @enddesc

   @var        complex_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The imaginary part
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_IMAG(CCTK_Cmplx, cctk_real, cctk_complex)           \
cctk_real CCTK_Cmplx##Imag (cctk_complex complex_number)                      \
{                                                                             \
  return (complex_number.Im);                                                 \
}


 /*@@
   @routine    CCTK_CmplxNeg
   @date       2006-06-07
   @author     Erik Schnetter
   @desc
               Returns the negative of a complex number.
   @enddesc
 
   @var       in
   @vdesc     The complex number
   @vtype     CCTK_COMPLEX
   @vio       in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The negative
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_NEG(CCTK_Cmplx, cctk_real, cctk_complex)            \
cctk_complex CCTK_Cmplx##Neg (cctk_complex complex_number)                    \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = -complex_number.Re;                                             \
  result.Im = -complex_number.Im;                                             \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxConjg
   @date       Tue Dec 14 12:22:41 1999
   @author     Tom Goodale
   @desc
               Returns the complex conjugate of a complex number.
   @enddesc
 
   @var       in
   @vdesc     The complex number
   @vtype     CCTK_COMPLEX
   @vio       in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The complex conjugate
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_CONJG(CCTK_Cmplx, cctk_real, cctk_complex)          \
cctk_complex CCTK_Cmplx##Conjg (cctk_complex complex_number)                  \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re =  complex_number.Re;                                             \
  result.Im = -complex_number.Im;                                             \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxAbs
   @date       Tue Dec 14 12:26:33 1999
   @author     Tom Goodale
   @desc
               Return the absolute value of a complex number.
   @enddesc
 
   @var        in
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The absolute value of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_ABS(CCTK_Cmplx, cctk_real, cctk_complex, type)      \
cctk_real CCTK_Cmplx##Abs (cctk_complex complex_number)                       \
{                                                                             \
  return ((cctk_real)hypot##type (complex_number.Re, complex_number.Im));     \
}


 /*@@
   @routine    CCTK_CmplxNorm
   @date       2006-07-06
   @author     Erik Schnetter
   @desc
               Return the absolute value squared of a complex number.
   @enddesc
 
   @var        in
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The absolute value squared of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_NORM(CCTK_Cmplx, cctk_real, cctk_complex, type)     \
cctk_real CCTK_Cmplx##Norm (cctk_complex complex_number)                      \
{                                                                             \
  return (pow##type (complex_number.Re, 2) + pow##type (complex_number.Im, 2));\
}


 /*@@
   @routine    CCTK_CmplxArg
   @date       2005-11-17
   @author     Erik Schnetter
   @desc
               Return the argument of a complex number.
   @enddesc
 
   @var        in
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The argument of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_ARG(CCTK_Cmplx, cctk_real, cctk_complex, type)      \
cctk_real CCTK_Cmplx##Arg (cctk_complex complex_number)                       \
{                                                                             \
  return (atan2##type (complex_number.Im, complex_number.Re));                \
}


 /*@@
   @routine    CCTK_CmplxAdd
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Adds two complex numbers
   @enddesc
 
   @var        a
   @vdesc      First summand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second summand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The sum of a and b
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_ADD(CCTK_Cmplx, cctk_real, cctk_complex)            \
cctk_complex CCTK_Cmplx##Add (cctk_complex a, cctk_complex b)                 \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = a.Re + b.Re;                                                    \
  result.Im = a.Im + b.Im;                                                    \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxSub
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Subtracts two complex numbers
   @enddesc
 
   @var        a
   @vdesc      First operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The difference
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_SUB(CCTK_Cmplx, cctk_real, cctk_complex)            \
cctk_complex CCTK_Cmplx##Sub (cctk_complex a, cctk_complex b)                 \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = a.Re - b.Re;                                                    \
  result.Im = a.Im - b.Im;                                                    \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxMul
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Multiplies two complex numbers
   @enddesc
 
   @var        a
   @vdesc      First operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The product
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_MUL(CCTK_Cmplx, cctk_real, cctk_complex)            \
cctk_complex CCTK_Cmplx##Mul (cctk_complex a, cctk_complex b)                 \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = a.Re*b.Re - a.Im*b.Im;                                          \
  result.Im = a.Im*b.Re + a.Re*b.Im;                                          \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxDiv
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Divide two complex numbers.
               Make sure that the intermediate values do not overflow.
               See e.g. the C99 compliant implementation in libgcc
               (which ships with gcc).
               This implementation here does not handle nan and inf.
   @enddesc
 
   @var        a
   @vdesc      First operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The quotient
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_DIV(CCTK_Cmplx, cctk_real, cctk_complex, type)      \
cctk_complex CCTK_Cmplx##Div (cctk_complex a, cctk_complex b)                 \
{                                                                             \
  cctk_real afact, bfact, factor;                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  afact = (cctk_real)(fabs##type(a.Re) + fabs##type(a.Im));                   \
  a.Re /= afact;                                                              \
  a.Im /= afact;                                                              \
  bfact = (cctk_real)(fabs##type(b.Re) + fabs##type(b.Im));                   \
  b.Re /= bfact;                                                              \
  b.Im /= bfact;                                                              \
  factor = b.Re*b.Re + b.Im*b.Im;                                             \
  result.Re = afact / bfact * (a.Re*b.Re + a.Im*b.Im) / factor;               \
  result.Im = afact / bfact * (a.Im*b.Re - a.Re*b.Im) / factor;               \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxCPow
   @date       2005-11-17
   @author     Erik Schnetter
   @desc
               Raises a complex number to a given power
               This algorithm was taken from glibc 2.3.5,
               file sysdeps/generic/s_cpow.c.
   @enddesc
 
   @var        a
   @vdesc      The base
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      The exponent
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               a to the power of b
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_CPOW(CCTK_Cmplx, cctk_real, cctk_complex)           \
cctk_complex CCTK_Cmplx##CPow (cctk_complex a, cctk_complex b)                \
{                                                                             \
  return CCTK_Cmplx##Exp (CCTK_Cmplx##Mul (b, CCTK_Cmplx##Log (a)));          \
}


 /*@@
   @routine    CCTK_CmplxSin
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the sine of a complex number.
   @enddesc
 
   @var        complex_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The sine
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_SIN(CCTK_Cmplx, cctk_real, cctk_complex, type)      \
cctk_complex CCTK_Cmplx##Sin (cctk_complex complex_number)                    \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = sin##type (complex_number.Re) * cosh##type (complex_number.Im); \
  result.Im = cos##type (complex_number.Re) * sinh##type (complex_number.Im); \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxCos
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the cosine of a complex number.
   @enddesc
 
   @var        complex_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The cosine
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_COS(CCTK_Cmplx, cctk_real, cctk_complex, type)      \
cctk_complex CCTK_Cmplx##Cos (cctk_complex complex_number)                    \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = cos##type (complex_number.Re) * cosh##type (complex_number.Im); \
  result.Im = sin##type (complex_number.Re) * sinh##type (complex_number.Im); \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxExp
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the exponential of a complex number.
   @enddesc
 
   @var        complex_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The exponential
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_EXP(CCTK_Cmplx, cctk_real, cctk_complex, type)      \
cctk_complex CCTK_Cmplx##Exp (cctk_complex complex_number)                    \
{                                                                             \
  cctk_real rho, theta;                                                       \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  rho = (cctk_real)exp##type (complex_number.Re);                             \
  theta = complex_number.Im;                                                  \
  result.Re = rho * (cctk_real)cos##type (theta);                             \
  result.Im = rho * (cctk_real)sin##type (theta);                             \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxLog
   @date       2005-11-17
   @author     Erik Schnetter
   @desc
               Returns the logarithm of a complex number.
               This algorithm was taken from glibc 2.3.5,
               file sysdeps/generic/s_clog.c.
   @enddesc
 
   @var        complex_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The logarithm
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_LOG(CCTK_Cmplx, cctk_real, cctk_complex, type)      \
cctk_complex CCTK_Cmplx##Log (cctk_complex complex_number)                    \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result.Re = log##type (hypot##type (complex_number.Re, complex_number.Im)); \
  result.Im = atan2##type (complex_number.Im, complex_number.Re);             \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxSqrt
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the square root of a complex number.
               This algorithm was taken from glibc 2.3.5,
               file sysdeps/generic/s_csqrt.c.
   @enddesc
   @history
               This code bears a resemblance to that in the GSL.  That
	       code contains references.
   @endhistory

   @var        complex_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The square root
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_SQRT(CCTK_Cmplx, cctk_real, cctk_complex, type)     \
cctk_complex CCTK_Cmplx##Sqrt (cctk_complex complex_number)                   \
{                                                                             \
  cctk_real d, r, s;                                                          \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  d = hypot##type (complex_number.Re, complex_number.Im);                     \
  /* Use the identity   2  Re res  Im res = Im x                              \
     to avoid cancellation error in  d +/- Re x.  */                          \
  if (complex_number.Re > 0)                                                  \
  {                                                                           \
    r = sqrt##type (0.5##type * d + 0.5##type * complex_number.Re);           \
    s = (0.5##type * complex_number.Im) / r;                                  \
  }                                                                           \
  else                                                                        \
  {                                                                           \
    s = sqrt##type (0.5##type * d - 0.5##type * complex_number.Re);           \
    r = fabs##type ((0.5##type * complex_number.Im) / s);                     \
  }                                                                           \
                                                                              \
  result.Re = r;                                                              \
  result.Im = copysign##type (s, complex_number.Im);                          \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxPow
   @date       
   @author     Yaakoub Y El Khamra
   @desc
               Raises a complex number to a given power
               This algorithm was taken from glibc 2.3.5,
               file sysdeps/generic/s_cpow.c.
   @enddesc
 
   @var        cx_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        w
   @vdesc      The exponent
   @vtype      CCTK_REAL
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The power of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_POW(CCTK_Cmplx, cctk_real, cctk_complex)            \
cctk_complex CCTK_Cmplx##Pow (cctk_complex complex_number, cctk_real w)       \
{                                                                             \
  cctk_complex result;                                                        \
                                                                              \
                                                                              \
  result = CCTK_Cmplx##Log (complex_number);                                  \
  result.Re *= w;                                                             \
  result.Im *= w;                                                             \
  result = CCTK_Cmplx##Exp (result);                                          \
                                                                              \
  return (result);                                                            \
}


 /*@@
   @routine    CCTK_CmplxIPow
   @date       
   @author     Erik Schnetter
   @desc
               Raises a complex number to a given power
               This algorithm was taken from glibc 2.3.5,
               file sysdeps/generic/s_cpow.c.
   @enddesc
 
   @var        cx_number
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        w
   @vdesc      The exponent
   @vtype      int
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The power of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_IPOW(CCTK_Cmplx, cctk_real, cctk_complex)     \
cctk_complex CCTK_Cmplx##IPow (cctk_complex complex_number, int w)      \
{                                                                       \
  cctk_complex result;                                                  \
                                                                        \
                                                                        \
  result = CCTK_Cmplx (1, 0);                                           \
  if (w < 0)                                                            \
  {                                                                     \
    result = CCTK_Cmplx##Div (result,                                   \
                              CCTK_Cmplx##IPow (complex_number, - w));  \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    while (w)                                                           \
    {                                                                   \
      if (w % 2)                                                        \
      {                                                                 \
        result = CCTK_Cmplx##Mul (result, complex_number);              \
      }                                                                 \
      w /= 2;                                                           \
      complex_number = CCTK_Cmplx##Mul (complex_number, complex_number); \
    }                                                                   \
  }                                                                     \
                                                                        \
  return (result);                                                      \
}



#ifdef DEFINE_CCTK_COMPLEX_INLINE_FUNCTIONS
#  define PREFIX static inline
#else
#  define PREFIX
#endif


#ifdef __cplusplus
  /* define C++ overloaded operators */
#  define DEFINE_CMPLX_CXX_OPERATORS(CCTK_Cmplx, cctk_real, cctk_complex) \
PREFIX cctk_complex operator+  (cctk_complex const & a)                         { return a; } \
PREFIX cctk_complex operator-  (cctk_complex const & a)                         { return CCTK_Cmplx##Neg (a); } \
PREFIX cctk_complex conj       (cctk_complex const & a)                         { return CCTK_Cmplx##Conjg (a); } \
PREFIX cctk_real    abs        (cctk_complex const & a)                         { return CCTK_Cmplx##Abs (a); } \
PREFIX cctk_real    arg        (cctk_complex const & a)                         { return CCTK_Cmplx##Arg (a); } \
PREFIX cctk_real    norm       (cctk_complex const & a)                         { return CCTK_Cmplx##Norm (a); } \
PREFIX cctk_complex operator+  (cctk_complex const & a, cctk_complex const & b) { return CCTK_Cmplx##Add (a, b); } \
PREFIX cctk_complex operator+  (cctk_real    const & a, cctk_complex const & b) { return CCTK_Cmplx##Add (CCTK_Cmplx (a, 0), b); } \
PREFIX cctk_complex operator+  (cctk_complex const & a, cctk_real    const & b) { return CCTK_Cmplx##Add (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex operator-  (cctk_complex const & a, cctk_complex const & b) { return CCTK_Cmplx##Sub (a, b); } \
PREFIX cctk_complex operator-  (cctk_real    const & a, cctk_complex const & b) { return CCTK_Cmplx##Sub (CCTK_Cmplx (a, 0), b); } \
PREFIX cctk_complex operator-  (cctk_complex const & a, cctk_real    const & b) { return CCTK_Cmplx##Sub (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex operator*  (cctk_complex const & a, cctk_complex const & b) { return CCTK_Cmplx##Mul (a, b); } \
PREFIX cctk_complex operator*  (cctk_real    const & a, cctk_complex const & b) { return CCTK_Cmplx##Mul (CCTK_Cmplx (a, 0), b); } \
PREFIX cctk_complex operator*  (cctk_complex const & a, cctk_real    const & b) { return CCTK_Cmplx##Mul (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex operator/  (cctk_complex const & a, cctk_complex const & b) { return CCTK_Cmplx##Div (a, b); } \
PREFIX cctk_complex operator/  (cctk_real    const & a, cctk_complex const & b) { return CCTK_Cmplx##Div (CCTK_Cmplx (a, 0), b); } \
PREFIX cctk_complex operator/  (cctk_complex const & a, cctk_real    const & b) { return CCTK_Cmplx##Div (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex operator+= (cctk_complex       & a, cctk_complex const & b) { return a = CCTK_Cmplx##Add (a, b); } \
PREFIX cctk_complex operator+= (cctk_complex       & a, cctk_real    const & b) { return a = CCTK_Cmplx##Add (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex operator-= (cctk_complex       & a, cctk_complex const & b) { return a = CCTK_Cmplx##Sub (a, b); } \
PREFIX cctk_complex operator-= (cctk_complex       & a, cctk_real    const & b) { return a = CCTK_Cmplx##Sub (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex operator*= (cctk_complex       & a, cctk_complex const & b) { return a = CCTK_Cmplx##Mul (a, b); } \
PREFIX cctk_complex operator*= (cctk_complex       & a, cctk_real    const & b) { return a = CCTK_Cmplx##Mul (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex operator/= (cctk_complex       & a, cctk_complex const & b) { return a = CCTK_Cmplx##Div (a, b); } \
PREFIX cctk_complex operator/= (cctk_complex       & a, cctk_real    const & b) { return a = CCTK_Cmplx##Div (a, CCTK_Cmplx (b, 0)); } \
PREFIX cctk_complex pow        (cctk_complex const & a, int                  i) { return CCTK_Cmplx##IPow (a, i); } \
PREFIX cctk_complex pow        (cctk_complex const & a, cctk_real    const & r) { return CCTK_Cmplx##Pow (a, r); } \
PREFIX cctk_complex pow        (cctk_complex const & a, cctk_complex const & b) { return CCTK_Cmplx##CPow (a, b); } \
PREFIX cctk_complex sin        (cctk_complex const & a)                         { return CCTK_Cmplx##Sin (a); } \
PREFIX cctk_complex cos        (cctk_complex const & a)                         { return CCTK_Cmplx##Cos (a); } \
PREFIX cctk_complex exp        (cctk_complex const & a)                         { return CCTK_Cmplx##Exp (a); } \
PREFIX cctk_complex log        (cctk_complex const & a)                         { return CCTK_Cmplx##Log (a); } \
PREFIX cctk_complex sqrt       (cctk_complex const & a)                         { return CCTK_Cmplx##Sqrt (a); } \
PREFIX std::ostream & operator << (std::ostream & os, cctk_complex const & a) { return os << CCTK_Cmplx##Real (a) << " " << CCTK_Cmplx##Imag (a); }
#else
  /* define no C++ overloaded operators */
#  define DEFINE_CMPLX_CXX_OPERATORS(CCTK_Cmplx, cctk_real, cctk_complex)
#endif


/* macro to define a set of complex functions for a given precision */
#define DEFINE_CMPLX_FUNCTIONS(CCTK_Cmplx, cctk_real, cctk_complex, type)     \
PREFIX DEFINE_CCTK_CMPLX       (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_REAL  (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_IMAG  (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_NEG   (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_CONJG (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_ABS   (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_NORM  (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_ARG   (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_ADD   (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_SUB   (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_MUL   (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_DIV   (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_CPOW  (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_SIN   (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_COS   (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_EXP   (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_LOG   (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_SQRT  (CCTK_Cmplx, cctk_real, cctk_complex, type)    \
PREFIX DEFINE_CCTK_CMPLX_POW   (CCTK_Cmplx, cctk_real, cctk_complex)    \
PREFIX DEFINE_CCTK_CMPLX_IPOW  (CCTK_Cmplx, cctk_real, cctk_complex)    \
DEFINE_CMPLX_CXX_OPERATORS(CCTK_Cmplx, cctk_real, cctk_complex)


/* define complex functions for all available precisions */
/* NOTE: The code in this file will function correctly with C float
 * and double types, but not with long double.
*/
#ifdef HAVE_CCTK_REAL4
  #define _REAL_TYPE   f
  DEFINE_CMPLX_FUNCTIONS (CCTK_Cmplx8, CCTK_REAL4, CCTK_COMPLEX8, _REAL_TYPE)
  #undef _REAL_TYPE
#endif

#ifdef HAVE_CCTK_REAL8
  #define _REAL_TYPE
  DEFINE_CMPLX_FUNCTIONS (CCTK_Cmplx16, CCTK_REAL8, CCTK_COMPLEX16, _REAL_TYPE)
  #undef _REAL_TYPE
#endif

#ifdef HAVE_CCTK_REAL16
  #define _REAL_TYPE   l
  DEFINE_CMPLX_FUNCTIONS (CCTK_Cmplx32, CCTK_REAL16, CCTK_COMPLEX32, _REAL_TYPE)
  #undef _REAL_TYPE
#endif

#undef PREFIX
