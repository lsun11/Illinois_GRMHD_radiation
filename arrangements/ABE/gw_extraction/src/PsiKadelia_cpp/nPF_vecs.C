// From brandt@aei-potsdam.mpg.de Thu May 28 09:31:22 1998
// To: Hisaaki Shinkai <shinkai@wurel.wustl.edu>
//
// There are three vectors, v1, v2, v3.  Each has an x, y, and z
// component.  However, the epsilon gives you v3 once you have v1
// and v2, so you don not need to supply it.  To get basis vectors
// aligned with the x, y, and z directions you could set v1x = 1,
// v1yL = v1zL = 0, v2xL = v2zL = 0, v2yL = 1.  To get spherical basis
// vectors one needs v1x = x, v1yL = y, v1zL = z, v2xL = -y, v2yL = x, 
// v2zL = 0. Or something like that.  This gives an r and a phi basis 
// vector (I believe it is the default -- the problem you had occurs 
// when x = y = z = 0.0).

// Intiialize vector 
//      "psif_vec" contains "radial", or as a default

if(whichvec==NEWAGEPSIFRIENDS_RADIAL_VEC) {
  
  
  //        Note to Steve:  These vectors are orthogonal, but are
  //        not normalized to 1.  Should they be?  Also, the vectors
  //        vanish at the origin and this will produce NaNs. Maybe
  //        a warning is in order if someone choses this option.  Miguel.

  //         v1x = xx
  //         v1yL = yy
  //         v1zL = zz
  
  //         v2xL = -yy
  //         v2yL = xx
  //         v2zL = 0.0D0
  
  v2xL = xx;
  v2yL = yy;
  v2zL = zz;
  
  v1xL = -yy;
  v1yL = xx;
  v1zL = 0.0;
  
 } else if(whichvec==NEWAGEPSIFRIENDS_CARTESIAN_VEC) {

  //      "psif_vec" contains "cartesian"
  
  //        Note to Steve:  Here I changed the original expressions to make
  //        the vectors orthonormal and avoid NaNs along the axis.  Miguel.
  //        write(*,*) "which=2 in pnd_nPF"
  
  v1xL = 1.0;
  v1yL = 0.0;
  v1zL = 0.0;
  
  v2xL = 0.0;
  v2yL = 1.0;
  v2zL = 0.0;
  
 } else if(whichvec==NEWAGEPSIFRIENDS_METDIAG_VEC) {
  
                 // 980611 9806018
  v1xL = 0.0;
  v1yL = 0.0;
  v1zL = 1.0/sqrt(ggzz);
  v2xL = 1.0/sqrt(ggxx);
  v2yL = 0.0;
  v2zL = 0.0;
  
 } else if(whichvec==NEWAGEPSIFRIENDS_SHOCK_VEC) {
  v1xL = 1.0/sqrt(ggxx);
  v1yL = 0.0;
  v1zL = 0.0;
  v2xL = 0.0;
  v2yL = 1.0/sqrt(ggyy);
  v2zL = 0.0;
  
 } else {
  printf("Psi Vector Choice improper\n");
  exit(1);
 }

//      write(*,*) "zxe",sqrt(xx*xx+yy*yy+zz*zz),sqrt(v1xL*v1xL*ggxx + v1yL*v1yL*ggyy + v1zL*v1zL*ggzz + 2.0*v1xL*v1yL*ggxy + 2.0*v1xL*v1zL*ggxz + 2.0*v1yL*v1zL*ggyz)

//#include "set_basis_vectors.x"

o1 = SQR(v1xL);
o2 = o1*ggxx;
o3 = 2.0*ggxy*v1xL*v1yL;
o4 = SQR(v1yL);
o5 = o4*ggyy;
o6 = 2.0*ggxz*v1xL*v1zL;
o7 = 2.0*ggyz*v1yL*v1zL;
o8 = SQR(v1zL);
o9 = o8*ggzz;
o10 = o2 + o3 + o5 + o6 + o7 + o9;
o11 = sqrt(o10);
o12 = 1.0/o11;
v1xL = o12*v1xL;
v1yL = o12*v1yL;
v1zL = o12*v1zL;
o1 = ggxx*v2xL;
o2 = ggxy*v2yL;
o3 = ggxz*v2zL;
o4 = o1 + o2 + o3;
o5 = o4*v1xL;
o6 = ggxy*v2xL;
o7 = ggyy*v2yL;
o8 = ggyz*v2zL;
o9 = o6 + o7 + o8;
o10 = o9*v1yL;
o11 = ggxz*v2xL;
o12 = ggyz*v2yL;
o13 = ggzz*v2zL;
o14 = o11 + o12 + o13;
o15 = o14*v1zL;
o16 = o5 + o10 + o15;
v2xL = -(o16*v1xL) + v2xL;
v2yL = -(o16*v1yL) + v2yL;
v2zL = -(o16*v1zL) + v2zL;
o1 = SQR(v2xL);
o2 = o1*ggxx;
o3 = 2.0*ggxy*v2xL*v2yL;
o4 = SQR(v2yL);
o5 = o4*ggyy;
o6 = 2.0*ggxz*v2xL*v2zL;
o7 = 2.0*ggyz*v2yL*v2zL;
o8 = SQR(v2zL);
o9 = o8*ggzz;
o10 = o2 + o3 + o5 + o6 + o7 + o9;
o11 = sqrt(o10);
o12 = 1.0/o11;
v2xL = o12*v2xL;
v2yL = o12*v2yL;
v2zL = o12*v2zL;
o1 = -(ggxz*ggyy);
o2 = ggxy*ggyz;
o3 = o1 + o2;
o4 = SQR(ggxz);
o5 = -(o4*ggyy);
o6 = 2.0*ggxy*ggxz*ggyz;
o7 = SQR(ggyz);
o8 = -(o7*ggxx);
o9 = SQR(ggxy);
o10 = -(o9*ggzz);
o11 = ggxx*ggyy*ggzz;
o12 = o5 + o6 + o8 + o10 + o11;
o13 = sqrt(o12);
o14 = 1.0/o13;
o15 = ggxz*ggyz;
o16 = -(ggxy*ggzz);
o17 = o15 + o16;
o18 = -o7;
o19 = ggyy*ggzz;
o20 = o18 + o19;
o21 = ggxy*ggxz;
o22 = -(ggxx*ggyz);
o23 = o21 + o22;
o24 = -o4;
o25 = ggxx*ggzz;
o26 = o24 + o25;
o27 = -o9;
o28 = ggxx*ggyy;
o29 = o27 + o28;
v3xL = -(o3*o14*v2xL*v1yL) + o3*o14*v1xL*v2yL + o14*o17*v2xL*v1zL - o14*o20*v2yL*v1zL - o14*o17*v1xL*v2zL + o14*o20*v1yL*v2zL;
v3yL = -(o14*o23*v2xL*v1yL) + o14*o23*v1xL*v2yL + o14*o26*v2xL*v1zL - o14*o17*v2yL*v1zL - o14*o26*v1xL*v2zL + o14*o17*v1yL*v2zL;
v3zL = -(o14*o29*v2xL*v1yL) + o14*o29*v1xL*v2yL + o14*o23*v2xL*v1zL - o3*o14*v2yL*v1zL - o14*o23*v1xL*v2zL + o3*o14*v1yL*v2zL;
//---------------------------------------------------

v1xtmp = v1xL;
v1ytmp = v1yL;
v1ztmp = v1zL;

v1xL = v2xL;
v1yL = v2yL;
v1zL = v2zL;

v2xL = v1xtmp;
v2yL = v1ytmp;
v2zL = v1ztmp;

v3xL = -v3xL;
v3yL = -v3yL;
v3zL = -v3zL;
