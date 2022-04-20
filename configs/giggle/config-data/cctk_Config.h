/* cctk_Config.h.  Generated automatically by configure.  */
#ifndef _CCTK_CONFIG_H_
#define _CCTK_CONFIG_H_

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* Define if standard C headers are available. */
#define STDC_HEADERS 1

/* Modifier for Fortran function definitions. */
#define CCTK_FCALL 

/* What debugging options to use */
/* #undef CCTK_DEBUG */
/* #undef CCTK_TRACEMEMORY */

/* Various library functions */
#define HAVE_GETHOSTBYNAME 1
#define HAVE_GETOPT_LONG_ONLY 1
#define HAVE_CRYPT 1
#define HAVE_FINITE 1
#define HAVE_ISNAN 1
#define HAVE_ISINF 1
#define HAVE_MKSTEMP 1
#define HAVE_VA_COPY 1

/* Do we have mode_t ? */
#define HAVE_MODE_T 1

/* Do we have SOCKET ? */
/* #undef HAVE_SOCKET_TYPE */

/* Do we have socklen_t ? Default to 'int' if not. */
#define HAVE_SOCKLEN_T 1
#ifdef HAVE_SOCKLEN_T
#  define CCTK_SOCKLEN_T socklen_t
#else
#  define CCTK_SOCKLEN_T int
#endif

/* Do we have hrtime_t ? */
/* #undef HAVE_HRTIME_T */

/* Some timing functions */
/* #undef HAVE_GETHRTIME */
/* #undef HAVE_READ_REAL_TIME */
/* #undef HAVE_TIME_BASE_TO_TIME */
#define HAVE_CLOCK_GETTIME 1
/* #undef HAVE_MACH_ABSOLUTE_TIME */

/* Cray UNICOS _rtc() (real-time clock) intrinsic */
/* #undef HAVE__RTC */

/* Some include things */
#define HAVE_TIME_H 1
/* #undef HAVE_SYS_FILIO_H */
#define HAVE_SYS_IOCTL_H 1
#define HAVE_SYS_SOCKET_H 1
#define HAVE_SYS_TIME_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_UNISTD_H 1
#define HAVE_STRING_H 1
#define HAVE_ASSERT_H 1
#define HAVE_TGMATH_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_GETOPT_H 1
#define HAVE_REGEX_H 1
#define HAVE_NETINET_IN_H 1
#define HAVE_NETDB_H 1
#define HAVE_ARPA_INET_H 1
/* #undef HAVE_WINSOCK2_H */
#define HAVE_CRYPT_H 1
#define HAVE_DIRENT_H 1
/* #undef HAVE_C_ASM_H */
/* #undef HAVE_INTRINSICS_H */
/* #undef HAVE_MACH_MACH_TIME_H */
#define HAVE_SIGNAL_H 1
#define HAVE_MALLOC_H 1
#define HAVE_MALLINFO 1
#define HAVE_MALLOPT 1
#define HAVE_M_MMAP_THRESHOLD_VALUE 1

#define TIME_WITH_SYS_TIME 1

/* Which format does the C++ STL on this machine provide. */
#define HAVE_VECTOR 1
/* #undef HAVE_VECTOR_H */

/* Timing stuff */
#define HAVE_TIME_GETTIMEOFDAY 1
#define GETTIMEOFDAY_NEEDS_TIMEZONE 1
#define HAVE_TIME_GETRUSAGE 1
/* #undef HAVE_TIME__FTIME */

/* Cache stuff */
#define CCTK_CACHELINE_BYTES 64
#define CCTK_CACHE_SIZE 2048*1024

/* The name of the NULL device for redirecting things to */
#define NULL_DEVICE "/dev/null"

/* Define the machine architecture for the build */
#define CCTK_BUILD_OS "linux-gnu"
#define CCTK_BUILD_CPU "x86_64"
#define CCTK_BUILD_VENDOR "unknown"

/******************************************************************************/

/* Size info for various datatypes */
#define SIZEOF_SHORT_INT 2
#define SIZEOF_INT 4
#define SIZEOF_LONG_INT 8
#define SIZEOF_LONG_LONG 8
#define SIZEOF_LONG_DOUBLE 16
#define SIZEOF_DOUBLE 8
#define SIZEOF_FLOAT 4
#define SIZEOF_CHAR_P 8

/* The chosen CCTK precision */

/* Floating point precision */
/* #undef CCTK_REAL_PRECISION_16 */
#define CCTK_REAL_PRECISION_8 1
/* #undef CCTK_REAL_PRECISION_4 */

/* Integer precision */
/* #undef CCTK_INTEGER_PRECISION_8 */
#define CCTK_INTEGER_PRECISION_4 1
/* #undef CCTK_INTEGER_PRECISION_2 */
/* #undef CCTK_INTEGER_PRECISION_1 */

/* Integer sizes */
#define HAVE_CCTK_INT8 1
#define HAVE_CCTK_INT4 1
#define HAVE_CCTK_INT2 1
#define HAVE_CCTK_INT1 1

/* Float sizes */
#define HAVE_CCTK_REAL16 1
#define HAVE_CCTK_REAL8 1
#define HAVE_CCTK_REAL4 1

/******************************************************************************/

#ifdef CCODE

/* CCTK C/C++ Integer datatypes */
#define CCTK_INT8 long int
#define CCTK_INT4 int
#define CCTK_INT2 short int
#define CCTK_INT1 signed char

/* CCTK C/C++ Float datatypes */
#define CCTK_REAL16 long double
#define CCTK_REAL8 double
#define CCTK_REAL4 float

/****************************************************************************/
/* C specific stuff */
/****************************************************************************/
#ifndef __cplusplus

/* Define to empty if the 'inline' keyword does not work. */
/* #undef inline */

/* Define to empty if the 'const' keyword does not work. */
/* #undef const */

/* Define to empty if the 'restrict' keyword does not work. */
#define CCTK_C_RESTRICT __restrict__
#ifdef CCTK_C_RESTRICT
#define restrict CCTK_C_RESTRICT
#endif

/* Allow the use of CCTK_RESTRICT as a qualifier always. */
#ifdef CCTK_C_RESTRICT
#define CCTK_RESTRICT CCTK_C_RESTRICT
#else
#define CCTK_RESTRICT restrict
#endif

/* Some C compilers (e.g. c99 ones) define bool */
/* #undef HAVE_CCTK_C_BOOL */

/* deprecated in beta15 */
#ifdef HAVE_CCTK_C_BOOL
#define CCTK_HAVE_C_BOOL
#endif

/* Since this is non-standard leave commented out for the moment */
#if 0
#ifndef HAVE_CCTK_C_BOOL
typedef enum {false, true} bool;
#endif /* HAVE_CCTK_CXX_BOOL */
#endif

#endif /* ! defined __cplusplus */
/****************************************************************************/

/****************************************************************************/
/* C++ specific stuff */
/****************************************************************************/
#ifdef __cplusplus

/* Some C++ compilers don't have bool ! */
#define HAVE_CCTK_CXX_BOOL 1

#ifndef HAVE_CCTK_CXX_BOOL
typedef enum {false, true} bool;
#else
/* deprecated in beta15 */
#define CCTK_HAVE_CXX_BOOL
#endif

/* Some C++ compilers recognise the restrict keyword */
#define CCTK_CXX_RESTRICT __restrict__

/* Since this is non-standard leave commented out for the moment */
#if 0
/* Define to empty if the keyword does not work. */
#ifdef CCTK_CXX_RESTRICT
#define restrict CCTK_CXX_RESTRICT
#endif
#endif

/* Allow the use of CCTK_RESTRICT as a qualifier always. */
#ifdef CCTK_CXX_RESTRICT
#define CCTK_RESTRICT CCTK_CXX_RESTRICT
#else
#define CCTK_RESTRICT restrict
#endif

#endif /* __cplusplus */
/****************************************************************************/

#endif /* CCODE */

#ifdef FCODE

#define HAVE_CCTK_FORTRAN_REAL4 1
#define HAVE_CCTK_FORTRAN_REAL8 1
#define HAVE_CCTK_FORTRAN_REAL16 1

#define HAVE_CCTK_FORTRAN_COMPLEX8 1
#define HAVE_CCTK_FORTRAN_COMPLEX16 1
#define HAVE_CCTK_FORTRAN_COMPLEX32 1

#endif /* FCODE */

/* Now include the code to pick an appropriate precison for reals and ints */
#include "cctk_Types.h"

/* Include any other stuff which is specific to this architecture */
#include "cctk_Archdefs.h"

/* Include any extra stuff from optional extra packages. */
#include "cctk_Extradefs.h"

#endif /* _CCTK_CONFIG_H_ */
