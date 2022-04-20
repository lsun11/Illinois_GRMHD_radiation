#ifndef _CCTK_LOOP_H_
#define _CCTK_LOOP_H_

#include "cctk_Config.h"



#ifdef CCODE
/* C and C++ variants of iterators */

#include "cctk_WarnLevel.h"
#include "cGH.h"



/* 1D */

#define CCTK_LOOP1(name, i,                     \
                   imin, imax, ilsh)            \
  do {                                          \
    typedef int lc_loop1_##name;                \
    int const lc_imin = (imin);                 \
    int const lc_imax = (imax);                 \
    for (int i=lc_imin; i<lc_imax; ++i) {
#define CCTK_ENDLOOP1(name)                             \
    }                                                   \
    typedef lc_loop1_##name lc_ensure_proper_nesting;   \
  } while (0)

#define CCTK_LOOP1_ALL(name, cctkGH, i)                                 \
  do {                                                                  \
    typedef int lc_loop1_all_##name;                                    \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 1) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP1_ALL can only be used in 1 dimension"); \
    }                                                                   \
    CCTK_LOOP1(name##_all,                                              \
               i,                                                       \
               0,                                                       \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)],                \
               lc_cctkGH->cctk_lsh[0]) {
#define CCTK_ENDLOOP1_ALL(name)                                 \
    } CCTK_ENDLOOP1(name##_all);                                \
    typedef lc_loop1_all_##name lc_ensure_proper_nesting;       \
  } while (0)

#define CCTK_LOOP1_INTERIOR(name, cctkGH, i,                            \
                            iblo, ibhi)                                 \
  do {                                                                  \
    typedef int lc_loop1_interior_##name;                               \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 1) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP1_INTERIOR can only be used in 1 dimension"); \
    }                                                                   \
    CCTK_LOOP1(name##_interior,                                         \
               i,                                                       \
               (iblo),                                                  \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)]-(ibhi),         \
               lc_cctkGH->cctk_lsh[0]) {
#define CCTK_ENDLOOP1_INTERIOR(name)                            \
    } CCTK_ENDLOOP1(name##_interior);                           \
    typedef lc_loop1_interior_##name lc_ensure_proper_nesting;  \
  } while(0)

#define CCTK_LOOP1_BOUNDARIES(name, cctkGH, i,                          \
                              iblo, ibhi)                               \
  do {                                                                  \
    typedef int lc_loop1_boundaries_##name;                             \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 1) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP1_BOUNDARIES can only be used in 1 dimension"); \
    }                                                                   \
    int const lc_blo[] = { (iblo) };                                    \
    int const lc_bhi[] = { (ibhi) };                                    \
    for (int lc_dir=0; lc_dir<1; ++lc_dir) {                            \
      for (int lc_face=0; lc_face<2; ++lc_face) {                       \
        int lc_bmin[1], lc_bmax[1];                                     \
        for (int lc_d=0; lc_d<1; ++lc_d) {                              \
          lc_bmin[lc_d] = 0;                                            \
          lc_bmax[lc_d] = lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,lc_d)];  \
          if (lc_d<lc_dir) {                                            \
            lc_bmin[lc_d] += lc_blo[lc_d];                              \
            lc_bmax[lc_d] -= lc_bhi[lc_d];                              \
          }                                                             \
        }                                                               \
        if (lc_face==0) {                                               \
          lc_bmax[lc_dir] = lc_bmin[lc_dir]+lc_blo[lc_dir];             \
        } else {                                                        \
          lc_bmin[lc_dir] = lc_bmax[lc_dir]-lc_bhi[lc_dir];             \
        }                                                               \
        CCTK_LOOP1(name##_boundaries,                                   \
                   i,                                                   \
                   lc_bmin[0],                                          \
                   lc_bmax[0],                                          \
                   lc_cctkGH->cctk_lsh[0]) {
#define CCTK_ENDLOOP1_BOUNDARIES(name)                                  \
        } CCTK_ENDLOOP1(name##_boundaries);                             \
      } /* for face */                                                  \
    }  /* for dir */                                                    \
    typedef lc_loop1_boundaries_##name lc_ensure_proper_nesting;        \
  } while (0)



/* 2D */

#define CCTK_LOOP2(name, i,j,                           \
                   imin,jmin, imax,jmax, ilsh,jlsh)     \
  do {                                                  \
    typedef int lc_loop2_##name;                        \
    int const lc_imin = (imin);                         \
    int const lc_jmin = (jmin);                         \
    int const lc_imax = (imax);                         \
    int const lc_jmax = (jmax);                         \
    for (int j=lc_jmin; j<lc_jmax; ++j) {               \
      for (int i=lc_imin; i<lc_imax; ++i) {
#define CCTK_ENDLOOP2(name)                             \
      }                                                 \
    }                                                   \
    typedef lc_loop2_##name lc_ensure_proper_nesting;   \
  } while (0)

#define CCTK_LOOP2_ALL(name, cctkGH, i,j)                               \
  do {                                                                  \
    typedef int lc_loop2_all_##name;                                    \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 2) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP2_ALL can only be used in 2 dimensions"); \
    }                                                                   \
    CCTK_LOOP2(name##_all,                                              \
               i,j,                                                     \
               0,0,                                                     \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)],                \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,1)],                \
               lc_cctkGH->cctk_lsh[0],                                  \
               lc_cctkGH->cctk_lsh[1]) {
#define CCTK_ENDLOOP2_ALL(name)                                 \
    } CCTK_ENDLOOP2(name##_all);                                \
    typedef lc_loop2_all_##name lc_ensure_proper_nesting;       \
  } while (0)

#define CCTK_LOOP2_INTERIOR(name, cctkGH, i,j,                          \
                            iblo,jblo, ibhi,jbhi)                       \
  do {                                                                  \
    typedef int lc_loop2_interior_##name;                               \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 2) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP2_INTERIOR can only be used in 2 dimensions"); \
    }                                                                   \
    CCTK_LOOP2(name##_interior,                                         \
               i,j,                                                     \
               (iblo),(jblo),                                           \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)]-(ibhi),         \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,1)]-(jbhi),         \
               lc_cctkGH->cctk_lsh[0],                                  \
               lc_cctkGH->cctk_lsh[1]) {
#define CCTK_ENDLOOP2_INTERIOR(name)                            \
    } CCTK_ENDLOOP2(name##_interior);                           \
    typedef lc_loop2_interior_##name lc_ensure_proper_nesting;  \
  } while(0)

#define CCTK_LOOP2_BOUNDARIES(name, cctkGH, i,j,                        \
                              iblo,jblo, ibhi,jbhi)                     \
  do {                                                                  \
    typedef int lc_loop2_boundaries_##name;                             \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 2) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP2_BOUNDARIES can only be used in 2 dimensions"); \
    }                                                                   \
    int const lc_blo[] = { (iblo),(jblo) };                             \
    int const lc_bhi[] = { (ibhi),(jbhi) };                             \
    for (int lc_dir=0; lc_dir<2; ++lc_dir) {                            \
      for (int lc_face=0; lc_face<2; ++lc_face) {                       \
        int lc_bmin[2], lc_bmax[2];                                     \
        for (int lc_d=0; lc_d<2; ++lc_d) {                              \
          lc_bmin[lc_d] = 0;                                            \
          lc_bmax[lc_d] = lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,lc_d)];  \
          if (lc_d<lc_dir) {                                            \
            lc_bmin[lc_d] += lc_blo[lc_d];                              \
            lc_bmax[lc_d] -= lc_bhi[lc_d];                              \
          }                                                             \
        }                                                               \
        if (lc_face==0) {                                               \
          lc_bmax[lc_dir] = lc_bmin[lc_dir]+lc_blo[lc_dir];             \
        } else {                                                        \
          lc_bmin[lc_dir] = lc_bmax[lc_dir]-lc_bhi[lc_dir];             \
        }                                                               \
        CCTK_LOOP2(name##_boundaries,                                   \
                   i,j,                                                 \
                   lc_bmin[0],lc_bmin[1],                               \
                   lc_bmax[0],lc_bmax[1],                               \
                   lc_cctkGH->cctk_lsh[0],                              \
                   lc_cctkGH->cctk_lsh[1]) {
#define CCTK_ENDLOOP2_BOUNDARIES(name)                                  \
        } CCTK_ENDLOOP2(name##_boundaries);                             \
      } /* for face */                                                  \
    }  /* for dir */                                                    \
    typedef lc_loop2_boundaries_##name lc_ensure_proper_nesting;        \
  } while (0)



/* 3D */

#define CCTK_LOOP3(name, i,j,k,                                         \
                   imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh)      \
  do {                                                                  \
    typedef int lc_loop3_##name;                                        \
    int const lc_imin = (imin);                                         \
    int const lc_jmin = (jmin);                                         \
    int const lc_kmin = (kmin);                                         \
    int const lc_imax = (imax);                                         \
    int const lc_jmax = (jmax);                                         \
    int const lc_kmax = (kmax);                                         \
    for (int k=lc_kmin; k<lc_kmax; ++k) {                               \
      for (int j=lc_jmin; j<lc_jmax; ++j) {                             \
        for (int i=lc_imin; i<lc_imax; ++i) {
#define CCTK_ENDLOOP3(name)                             \
        }                                               \
      }                                                 \
    }                                                   \
    typedef lc_loop3_##name lc_ensure_proper_nesting;   \
  } while (0)

#define CCTK_LOOP3_ALL(name, cctkGH, i,j,k)                             \
  do {                                                                  \
    typedef int lc_loop3_all_##name;                                    \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 3) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP3_ALL can only be used in 3 dimensions"); \
    }                                                                   \
    CCTK_LOOP3(name##_all,                                              \
               i,j,k,                                                   \
               0,0,0,                                                   \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)],                \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,1)],                \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,2)],                \
               lc_cctkGH->cctk_lsh[0],                                  \
               lc_cctkGH->cctk_lsh[1],                                  \
               lc_cctkGH->cctk_lsh[2]) {
#define CCTK_ENDLOOP3_ALL(name)                                 \
    } CCTK_ENDLOOP3(name##_all);                                \
    typedef lc_loop3_all_##name lc_ensure_proper_nesting;       \
  } while (0)

#define CCTK_LOOP3_INTERIOR(name, cctkGH, i,j,k,                        \
                            iblo,jblo,kblo, ibhi,jbhi,kbhi)             \
  do {                                                                  \
    typedef int lc_loop3_interior_##name;                               \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 3) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP3_INTERIOR can only be used in 3 dimensions"); \
    }                                                                   \
    CCTK_LOOP3(name##_interior,                                         \
               i,j,k,                                                   \
               (iblo),(jblo),(kblo),                                    \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)]-(ibhi),         \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,1)]-(jbhi),         \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,2)]-(kbhi),         \
               lc_cctkGH->cctk_lsh[0],                                  \
               lc_cctkGH->cctk_lsh[1],                                  \
               lc_cctkGH->cctk_lsh[2]) {
#define CCTK_ENDLOOP3_INTERIOR(name)                            \
    } CCTK_ENDLOOP3(name##_interior);                           \
    typedef lc_loop3_interior_##name lc_ensure_proper_nesting;  \
  } while(0)

#define CCTK_LOOP3_BOUNDARIES(name, cctkGH, i,j,k,                      \
                              iblo,jblo,kblo, ibhi,jbhi,kbhi)           \
  do {                                                                  \
    typedef int lc_loop3_boundaries_##name;                             \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 3) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP3_BOUNDARIES can only be used in 3 dimensions"); \
    }                                                                   \
    int const lc_blo[] = { (iblo),(jblo),(kblo) };                      \
    int const lc_bhi[] = { (ibhi),(jbhi),(kbhi) };                      \
    for (int lc_dir=0; lc_dir<3; ++lc_dir) {                            \
      for (int lc_face=0; lc_face<2; ++lc_face) {                       \
        int lc_bmin[3], lc_bmax[3];                                     \
        for (int lc_d=0; lc_d<3; ++lc_d) {                              \
          lc_bmin[lc_d] = 0;                                            \
          lc_bmax[lc_d] = lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,lc_d)];  \
          if (lc_d<lc_dir) {                                            \
            lc_bmin[lc_d] += lc_blo[lc_d];                              \
            lc_bmax[lc_d] -= lc_bhi[lc_d];                              \
          }                                                             \
        }                                                               \
        if (lc_face==0) {                                               \
          lc_bmax[lc_dir] = lc_bmin[lc_dir]+lc_blo[lc_dir];             \
        } else {                                                        \
          lc_bmin[lc_dir] = lc_bmax[lc_dir]-lc_bhi[lc_dir];             \
        }                                                               \
        CCTK_LOOP3(name##_boundaries,                                   \
                   i,j,k,                                               \
                   lc_bmin[0],lc_bmin[1],lc_bmin[2],                    \
                   lc_bmax[0],lc_bmax[1],lc_bmax[2],                    \
                   lc_cctkGH->cctk_lsh[0],                              \
                   lc_cctkGH->cctk_lsh[1],                              \
                   lc_cctkGH->cctk_lsh[2]) {
#define CCTK_ENDLOOP3_BOUNDARIES(name)                                  \
        } CCTK_ENDLOOP3(name##_boundaries);                             \
      } /* for face */                                                  \
    }  /* for dir */                                                    \
    typedef lc_loop3_boundaries_##name lc_ensure_proper_nesting;        \
  } while (0)



/* 4D */

#define CCTK_LOOP4(name, i,j,k,l,                                       \
                   imin,jmin,kmin,lmin, imax,jmax,kmax,lmax, ilsh,jlsh,klsh,llsh) \
  do {                                                                  \
    typedef int lc_loop4_##name;                                        \
    int const lc_imin = (imin);                                         \
    int const lc_jmin = (jmin);                                         \
    int const lc_kmin = (kmin);                                         \
    int const lc_lmin = (lmin);                                         \
    int const lc_imax = (imax);                                         \
    int const lc_jmax = (jmax);                                         \
    int const lc_kmax = (kmax);                                         \
    int const lc_lmax = (lmax);                                         \
    for (int l=lc_lmin; l<lc_lmax; ++l) {                               \
     for (int k=lc_kmin; k<lc_kmax; ++k) {                              \
       for (int j=lc_jmin; j<lc_jmax; ++j) {                            \
         for (int i=lc_imin; i<lc_imax; ++i) {
#define CCTK_ENDLOOP4(name)                             \
          }                                             \
        }                                               \
      }                                                 \
    }                                                   \
    typedef lc_loop4_##name lc_ensure_proper_nesting;   \
  } while (0)

#define CCTK_LOOP4_ALL(name, cctkGH, i,j,k,l)                           \
  do {                                                                  \
    typedef int lc_loop4_all_##name;                                    \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 4) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP4_ALL can only be used in 4 dimensions"); \
    }                                                                   \
    CCTK_LOOP4(name##_all,                                              \
               i,j,k,l,                                                 \
               0,0,0,0,                                                 \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)],                \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,1)],                \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,2)],                \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,3)],                \
               lc_cctkGH->cctk_lsh[0],                                  \
               lc_cctkGH->cctk_lsh[1],                                  \
               lc_cctkGH->cctk_lsh[2],                                  \
               lc_cctkGH->cctk_lsh[3]) {
#define CCTK_ENDLOOP4_ALL(name)                                 \
    } CCTK_ENDLOOP4(name##_all);                                \
    typedef lc_loop4_all_##name lc_ensure_proper_nesting;       \
  } while (0)

#define CCTK_LOOP4_INTERIOR(name, cctkGH, i,j,k,l,                      \
                            iblo,jblo,kblo,lblo, ibhi,jbhi,kbhi,lbhi)   \
  do {                                                                  \
    typedef int lc_loop4_interior_##name;                               \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 4) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP4_INTERIOR can only be used in 4 dimensions"); \
    }                                                                   \
    CCTK_LOOP4(name##_interior,                                         \
               i,j,k,l,                                                 \
               (iblo),(jblo),(kblo),(lblo),                             \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,0)]-(ibhi),         \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,1)]-(jbhi),         \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,2)]-(kbhi),         \
               lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,3)]-(lbhi),         \
               lc_cctkGH->cctk_lsh[0],                                  \
               lc_cctkGH->cctk_lsh[1],                                  \
               lc_cctkGH->cctk_lsh[2],                                  \
               lc_cctkGH->cctk_lsh[3]) {
#define CCTK_ENDLOOP4_INTERIOR(name)                            \
    } CCTK_ENDLOOP4(name##_interior);                           \
    typedef lc_loop4_interior_##name lc_ensure_proper_nesting;  \
  } while(0)

#define CCTK_LOOP4_BOUNDARIES(name, cctkGH, i,j,k,l,                    \
                              iblo,jblo,kblo,lblo, ibhi,jbhi,kbhi,lbhi) \
  do {                                                                  \
    typedef int lc_loop4_boundaries_##name;                             \
    cGH const *CCTK_RESTRICT const lc_cctkGH = (cctkGH);                \
    if (lc_cctkGH->cctk_dim != 4) {                                     \
      CCTK_WARN (CCTK_WARN_ABORT,                                       \
                 "Macro CCTK_LOOP4_BOUNDARIES can only be used in 4 dimensions"); \
    }                                                                   \
    int const lc_blo[] = { (iblo),(jblo),(kblo),(lblo) };               \
    int const lc_bhi[] = { (ibhi),(jbhi),(kbhi),(lbhi) };               \
    for (int lc_dir=0; lc_dir<4; ++lc_dir) {                            \
      for (int lc_face=0; lc_face<2; ++lc_face) {                       \
        int lc_bmin[4], lc_bmax[4];                                     \
        for (int lc_d=0; lc_d<4; ++lc_d) {                              \
          lc_bmin[lc_d] = 0;                                            \
          lc_bmax[lc_d] = lc_cctkGH->cctk_lssh[CCTK_LSSH_IDX(0,lc_d)];  \
          if (lc_d<lc_dir) {                                            \
            lc_bmin[lc_d] += lc_blo[lc_d];                              \
            lc_bmax[lc_d] -= lc_bhi[lc_d];                              \
          }                                                             \
        }                                                               \
        if (lc_face==0) {                                               \
          lc_bmax[lc_dir] = lc_bmin[lc_dir]+lc_blo[lc_dir];             \
        } else {                                                        \
          lc_bmin[lc_dir] = lc_bmax[lc_dir]-lc_bhi[lc_dir];             \
        }                                                               \
        CCTK_LOOP4(name##_boundaries,                                   \
                   i,j,k,l,                                             \
                   lc_bmin[0],lc_bmin[1],lc_bmin[2],lc_bmin[3],         \
                   lc_bmax[0],lc_bmax[1],lc_bmax[2],lc_bmax[3],         \
                   lc_cctkGH->cctk_lsh[0],                              \
                   lc_cctkGH->cctk_lsh[1],                              \
                   lc_cctkGH->cctk_lsh[2],                              \
                   lc_cctkGH->cctk_lsh[3]) {
#define CCTK_ENDLOOP4_BOUNDARIES(name)                                  \
        } CCTK_ENDLOOP4(name##_boundaries);                             \
      } /* for face */                                                  \
    }  /* for dir */                                                    \
    typedef lc_loop4_boundaries_##name lc_ensure_proper_nesting;        \
  } while (0)



#endif



#ifdef FCODE



/* 3D */

#define CCTK_LOOP3_DECLARE(name)                                \
   integer :: name/**/_imin, name/**/_jmin, name/**/_kmin    && \
   integer :: name/**/_imax, name/**/_jmax, name/**/_kmax

#define CCTK_LOOP3_OMP_PRIVATE(name)

#define CCTK_LOOP3(name, i,j,k,                                         \
                   imin,jmin,kmin, imax,jmax,kmax, ilsh,jlsh,klsh)      \
   name/**/_imin = imin                                              && \
   name/**/_jmin = jmin                                              && \
   name/**/_kmin = kmin                                              && \
   name/**/_imax = imax                                              && \
   name/**/_jmax = jmax                                              && \
   name/**/_kmax = kmax                                              && \
   do k = name/**/_kmin, name/**/_kmax                               && \
      do j = name/**/_jmin, name/**/_jmax                            && \
         do i = name/**/_imin, name/**/_imax
#define CCTK_ENDLOOP3(name)                     \
         end do                              && \
      end do                                 && \
   end do

#define CCTK_LOOP3_ALL_DECLARE(name)            \
   CCTK_LOOP3_DECLARE(name)
#define CCTK_LOOP3_ALL_OMP_PRIVATE(name)        \
   CCTK_LOOP3_OMP_PRIVATE(name)
#define CCTK_LOOP3_ALL(name, i,j,k)                     \
   CCTK_LOOP3(name,                                     \
              i,j,k,                                    \
              1,1,1,                                    \
              CCTK_LSSH(0,1),                           \
              CCTK_LSSH(0,2),                           \
              CCTK_LSSH(0,3),                           \
              cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
#define CCTK_ENDLOOP3_ALL(name)                 \
   CCTK_ENDLOOP3(name)

#define CCTK_LOOP3_INTERIOR_DECLARE(name)       \
   CCTK_LOOP3_DECLARE(name)
#define CCTK_LOOP3_INTERIOR_OMP_PRIVATE(name)   \
   CCTK_LOOP3_OMP_PRIVATE(name)
#define CCTK_LOOP3_INTERIOR(name, i,j,k,                        \
                            iblo,jblo,kblo, ibhi,jbhi,kbhi)     \
   CCTK_LOOP3(name,                                             \
              i,j,k,                                            \
              (iblo),(jblo),(kblo),                             \
              CCTK_LSSH(0,1)-(ibhi),                            \
              CCTK_LSSH(0,2)-(jbhi),                            \
              CCTK_LSSH(0,3)-(kbhi),                            \
              cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
#define CCTK_ENDLOOP3_INTERIOR(name)            \
   CCTK_ENDLOOP3(name)

#define CCTK_LOOP3_BOUNDARIES_DECLARE(name)     \
   CCTK_LOOP3_DECLARE(name)                  && \
   integer :: lc_bmin(3), lc_bmax(3)         && \
   integer :: lc_blo(3), lc_bhi(3)           && \
   integer :: lc_dir, lc_face                && \
   integer :: lc_d
#define CCTK_LOOP3_BOUNDARIES_OMP_PRIVATE(name) \
   CCTK_LOOP3_OMP_PRIVATE(name)
#define CCTK_LOOP3_BOUNDARIES(name, i,j,k,                      \
                              iblo,jblo,kblo, ibhi,jbhi,kbhi)   \
   lc_blo = (/ (iblo),(jblo),(kblo) /)                       && \
   lc_bhi = (/ (ibhi),(jbhi),(kbhi) /)                       && \
   do lc_dir=1,3                                             && \
      do lc_face=1,2                                         && \
         do lc_d=1,3                                         && \
            lc_bmin(lc_d) = 0                                && \
            lc_bmax(lc_d) = CCTK_LSSH(0,lc_d)                && \
            if (lc_d<lc_dir) then                            && \
               lc_bmin(lc_d) = lc_bmin(lc_d)+lc_blo(lc_d)    && \
               lc_bmax(lc_d) = lc_bmax(lc_d)-lc_bhi(lc_d)    && \
            end if                                           && \
         end do                                              && \
         if (lc_face==0) then                                && \
            lc_bmax(lc_dir) = lc_bmin(lc_dir)+lc_blo(lc_dir) && \
         else                                                && \
            lc_bmin(lc_dir) = lc_bmax(lc_dir)-lc_bhi(lc_dir) && \
         end if                                              && \
         CCTK_LOOP3(name,                                       \
                    i,j,k,                                      \
                    lc_bmin(1),lc_bmin(2),lc_bmin(3),           \
                    lc_bmax(1),lc_bmax(2),lc_bmax(3),           \
                    cctk_lsh(1),cctk_lsh(2),cctk_lsh(3))
#define CCTK_ENDLOOP3_BOUNDARIES(name)          \
         CCTK_ENDLOOP3(name)                 && \
      end do /* face */                      && \
   end do /* dir */
  


#endif

#endif /* #ifndef _CCTK_LOOP_H_ */
