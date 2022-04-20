o1 = SQR(v2xL);
o2 = -0.5*o1*magxx;
o3 = SQR(v3xL);
o4 = 0.5*o3*magxx;
o5 = -(magxy*v2xL*v2yL);
o6 = SQR(v2yL);
o7 = -0.5*o6*magyy;
o8 = magxy*v3xL*v3yL;
o9 = SQR(v3yL);
o10 = 0.5*o9*magyy;
o11 = -(magxz*v2xL*v2zL);
o12 = -(magyz*v2yL*v2zL);
o13 = SQR(v2zL);
o14 = -0.5*o13*magzz;
o15 = magxz*v3xL*v3zL;
o16 = magyz*v3yL*v3zL;
o17 = SQR(v3zL);
o18 = 0.5*o17*magzz;
o19 = -0.5*o1*elecxx;
o20 = 0.5*o3*elecxx;
o21 = -(elecxy*v2xL*v2yL);
o22 = -0.5*o6*elecyy;
o23 = elecxy*v3xL*v3yL;
o24 = 0.5*o9*elecyy;
o25 = -(elecxz*v2xL*v2zL);
o26 = -(elecyz*v2yL*v2zL);
o27 = -0.5*o13*eleczz;
o28 = elecxz*v3xL*v3zL;
o29 = elecyz*v3yL*v3zL;
o30 = 0.5*o17*eleczz;
o31 = 0.5*elecxx*v1xL*v3xL;
o32 = 0.5*elecxy*v3xL*v1yL;
o33 = 0.5*elecxy*v1xL*v3yL;
o34 = 0.5*elecyy*v1yL*v3yL;
o35 = 0.5*elecxz*v3xL*v1zL;
o36 = 0.5*elecyz*v3yL*v1zL;
o37 = 0.5*elecxz*v1xL*v3zL;
o38 = 0.5*elecyz*v1yL*v3zL;
o39 = 0.5*eleczz*v1zL*v3zL;
o40 = -0.5*magxx*v1xL*v3xL;
o41 = -0.5*magxy*v3xL*v1yL;
o42 = -0.5*magxy*v1xL*v3yL;
o43 = -0.5*magyy*v1yL*v3yL;
o44 = -0.5*magxz*v3xL*v1zL;
o45 = -0.5*magyz*v3yL*v1zL;
o46 = -0.5*magxz*v1xL*v3zL;
o47 = -0.5*magyz*v1yL*v3zL;
o48 = -0.5*magzz*v1zL*v3zL;
o49 = SQR(v1xL);
o50 = SQR(v1yL);
o51 = SQR(v1zL);
psi0im_loc = o2 + o4 + o5 + o7 + o8 + o10 + o11 + o12 + o14 + o15 + o16 + o18 + elecxx*v2xL*v3xL + elecxy*v3xL*v2yL + elecxy*v2xL*v3yL + elecyy*v2yL*v3yL + elecxz*v3xL*v2zL + elecyz*v3yL*v2zL + elecxz*v2xL*v3zL + elecyz*v2yL*v3zL + eleczz*v2zL*v3zL;
psi0re_loc = o19 + o20 + o21 + o22 + o23 + o24 + o25 + o26 + o27 + o28 + o29 + o30 - magxx*v2xL*v3xL - magxy*v3xL*v2yL - magxy*v2xL*v3yL - magyy*v2yL*v3yL - magxz*v3xL*v2zL - magyz*v3yL*v2zL - magxz*v2xL*v3zL - magyz*v2yL*v3zL - magzz*v2zL*v3zL;
//                 psi1im_loc = o31 + o32 + o33 + o34 + o35 + o36 + o37 + o38 + o39 - 5&
  //                      &.d-1*magxx*v1xL*v2xL - 0.5*magxy*v2xL*v1yL - 0.5*magxy*v1xL*v2yL -&
//                      & 0.5*magyy*v1yL*v2yL - 0.5*magxz*v2xL*v1zL - 0.5*magyz*v2yL*v1zL&
//                      & - 0.5*magxz*v1xL*v2zL - 0.5*magyz*v1yL*v2zL - 0.5*magzz*v1zL*v&
//                      &2zL
//                 psi1re_loc = o40 + o41 + o42 + o43 + o44 + o45 + o46 + o47 + o48 - 5&
//                      &.d-1*elecxx*v1xL*v2xL - 0.5*elecxy*v2xL*v1yL - 0.5*elecxy*v1xL*v2&
//                      &yL - 0.5*elecyy*v1yL*v2yL - 0.5*elecxz*v2xL*v1zL - 0.5*elecyz*v&
//                      &2yL*v1zL - 0.5*elecxz*v1xL*v2zL - 0.5*elecyz*v1yL*v2zL - 0.5*ele&
//                      &czz*v1zL*v2zL
//                 psi2im_loc = -0.5*o49*magxx - 0.5*o50*magyy - 0.5*o51*magzz - &
//                      &magxy*v1xL*v1yL - magxz*v1xL*v1zL - magyz*v1yL*v1zL
//                 psi2re_loc = -0.5*o49*elecxx - 0.5*o50*elecyy - 0.5*o51*eleczz&
//                      & - elecxy*v1xL*v1yL - elecxz*v1xL*v1zL - elecyz*v1yL*v1zL
//                 psi3im_loc = o31 + o32 + o33 + o34 + o35 + o36 + o37 + o38 + o39 + 5&
//                      &.d-1*magxx*v1xL*v2xL + 0.5*magxy*v2xL*v1yL + 0.5*magxy*v1xL*v2yL +&
//                      & 0.5*magyy*v1yL*v2yL + 0.5*magxz*v2xL*v1zL + 0.5*magyz*v2yL*v1zL&
//                      & + 0.5*magxz*v1xL*v2zL + 0.5*magyz*v1yL*v2zL + 0.5*magzz*v1zL*v&
//                      &2zL
//                 psi3re_loc = o40 + o41 + o42 + o43 + o44 + o45 + o46 + o47 + o48 + 5&
//                      &.d-1*elecxx*v1xL*v2xL + 0.5*elecxy*v2xL*v1yL + 0.5*elecxy*v1xL*v2&
//                      &yL + 0.5*elecyy*v1yL*v2yL + 0.5*elecxz*v2xL*v1zL + 0.5*elecyz*v&
//                      &2yL*v1zL + 0.5*elecxz*v1xL*v2zL + 0.5*elecyz*v1yL*v2zL + 0.5*ele&
//                      &czz*v1zL*v2zL
//                 psi4im_loc = o2 + o4 + o5 + o7 + o8 + o10 + o11 + o12 + o14 + o15 + &
//                      &o16 + o18 - elecxx*v2xL*v3xL - elecxy*v3xL*v2yL - elecxy*v2xL*v3yL - e&
//                      &lecyy*v2yL*v3yL - elecxz*v3xL*v2zL - elecyz*v3yL*v2zL - elecxz*v2xL*v3zL&
//                      & - elecyz*v2yL*v3zL - eleczz*v2zL*v3zL
//                 psi4re_loc = o19 + o20 + o21 + o22 + o23 + o24 + o25 + o26 + o27 + o&
//                      &28 + o29 + o30 + magxx*v2xL*v3xL + magxy*v3xL*v2yL + magxy*v2xL*v3yL +&
//                      & magyy*v2yL*v3yL + magxz*v3xL*v2zL + magyz*v3yL*v2zL + magxz*v2xL*v3zL +&
//                      & magyz*v2yL*v3zL + magzz*v2zL*v3zL
//                 o1 = psi2im_loc**2
//                 o2 = psi2re_loc**2
//                 o3 = psi3im_loc**2
//                 o4 = psi3re_loc**2
//                 o5 = psi1im_loc**2
//                 o6 = psi1re_loc**2
//                 icurvre_loc = -3.d0*o1 + 3.d0*o2 + 4.d0*psi1im_loc*psi3im_loc - psi0im_loc*psi4i&
//                      &m_loc - 4.d0*psi1re_loc*psi3re_loc + psi0re_loc*psi4re_loc
//                 icurvim_loc = psi4im_loc*psi0re_loc - 4.d0*psi3im_loc*psi1re_loc + 6.d0*psi2im_loc*psi2r&
//                      &e_loc - 4.d0*psi1im_loc*psi3re_loc + psi0im_loc*psi4re_loc
//                 jcurvre_loc = o3*psi0re_loc - o4*psi0re_loc - psi2im_loc*psi4im_loc*psi0re_loc - 2.d0*ps&
//                      &i2im_loc*psi3im_loc*psi1re_loc + 2.d0*psi1im_loc*psi4im_loc*psi1re_loc + 3.d0*o1*psi2re_loc &
//                      &- o2*psi2re_loc - 2.d0*psi1im_loc*psi3im_loc*psi2re_loc - psi0im_loc*psi4im_loc*psi2re_loc -&
//                      & 2.d0*psi1im_loc*psi2im_loc*psi3re_loc + 2.d0*psi0im_loc*psi3im_loc*psi3re_loc + 2.d0*ps&
//                      &i1re_loc*psi2re_loc*psi3re_loc + o5*psi4re_loc - o6*psi4re_loc - psi0im_loc*psi2im_loc*psi4r&
//                      &e_loc + psi0re_loc*psi2re_loc*psi4re_loc
//                 jcurvim_loc = o3*psi0im_loc - o4*psi0im_loc + o1*psi2im_loc - 3.d0*o2*psi2im_loc - 2&
//                      &.d0*psi1im_loc*psi2im_loc*psi3im_loc + o5*psi4im_loc - o6*psi4im_loc - psi0im_loc*psi2im_loc&
//                      &*psi4im_loc + psi4im_loc*psi0re_loc*psi2re_loc + 2.d0*psi3im_loc*psi1re_loc*psi2re_loc - 2.d&
//                      &0*psi3im_loc*psi0re_loc*psi3re_loc + 2.d0*psi2im_loc*psi1re_loc*psi3re_loc + 2.d0*psi1im_loc&
//                      &*psi2re_loc*psi3re_loc + psi2im_loc*psi0re_loc*psi4re_loc - 2.d0*psi1im_loc*psi1re_loc*psi4r&
//                      &e_loc + psi0im_loc*psi2re_loc*psi4re_loc
