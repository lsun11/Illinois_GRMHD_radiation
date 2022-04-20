/* electric_and_magnetic.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

   http://www.netlib.org/f2c/libf2c.zip
*/
/* System generated locals */
double tmpsqr;

/* Computing 2nd power */
tmpsqr = ggyz;
o1 = tmpsqr * tmpsqr;
o2 = -o1;
o3 = ggyy * ggzz;
o4 = o2 + o3;
/* Computing 2nd power */
tmpsqr = ggxz;
o5 = tmpsqr * tmpsqr;
o6 = -(o5 * ggyy);
o7 = ggxy * 2. * ggxz * ggyz;
o8 = -(o1 * ggxx);
/* Computing 2nd power */
tmpsqr = ggxy;
o9 = tmpsqr * tmpsqr;
o10 = -(o9 * ggzz);
o11 = ggxx * ggyy * ggzz;
o12 = o6 + o7 + o8 + o10 + o11;
o13 = 1. / o12;
/* Computing 2nd power */
tmpsqr = kkxz;
o14 = (tmpsqr * tmpsqr);
o15 = ggxz * ggyz;
o16 = -(ggxy * ggzz);
o17 = o15 + o16;
o18 = -o5;
o19 = ggxx * ggzz;
o20 = o18 + o19;
/* Computing 2nd power */
tmpsqr = kkyz;
o21 = (tmpsqr * tmpsqr);
o22 = sqrt(o12);
o23 = -(ggxz * ggyy);
o24 = ggxy * ggyz;
o25 = o23 + o24;
/* Computing 2nd power */
tmpsqr = o12;
o26 = tmpsqr * tmpsqr;
o27 = 1. / o26;
o28 = -(o20 * o25 * o27);
o29 = ggxy * ggxz;
o30 = -(ggxx * ggyz);
o31 = o29 + o30;
o32 = o17 * o27 * o31;
o33 = o28 + o32;
o34 = -(o17 * o25 * o27);
o35 = o4 * o27 * o31;
o36 = o34 + o35;
/* Computing 2nd power */
tmpsqr = o17;
o37 = tmpsqr * tmpsqr;
o38 = -(o27 * o37);
o39 = o4 * o20 * o27;
o40 = o38 + o39;
o41 = o13 * 0.5 * o25 * dygzz;
o42 = o13 * 0.5 * o17 * dzgyy;
o43 = o4 * -0.5 * o13 * dxgyz;
o44 = o4 * 0.5 * o13 * dygxz;
o45 = o4 * 0.5 * o13 * dzgxy;
o46 = o41 + o42 + o43 + o44 + o45;
o47 = o13 * 0.5 * o25 * dzgzz;
o48 = o13 * -0.5 * o17 * dygzz;
o49 = o13 * o17 * dzgyz;
o50 = o4 * -0.5 * o13 * dxgzz;
o51 = o4 * o13 * dzgxz;
o52 = o47 + o48 + o49 + o50 + o51;
o53 = o13 * 0.5 * o31 * dygzz;
o54 = o13 * 0.5 * o20 * dzgyy;
o55 = o13 * -0.5 * o17 * dxgyz;
o56 = o13 * 0.5 * o17 * dygxz;
o57 = o13 * 0.5 * o17 * dzgxy;
o58 = o53 + o54 + o55 + o56 + o57;
o59 = o13 * 0.5 * o31 * dzgzz;
o60 = o13 * -0.5 * o20 * dygzz;
o61 = o13 * o20 * dzgyz;
o62 = o13 * -0.5 * o17 * dxgzz;
o63 = o13 * o17 * dzgxz;
o64 = o59 + o60 + o61 + o62 + o63;
o65 = o13 * 0.5 * o25 * dxgzz;
o66 = o13 * 0.5 * o17 * dxgyz;
o67 = o13 * -0.5 * o17 * dygxz;
o68 = o4 * 0.5 * o13 * dzgxx;
o69 = o57 + o65 + o66 + o67 + o68;
o70 = -o9;
o71 = ggxx * ggyy;
o72 = o70 + o71;
o73 = o13 * 0.5 * o72 * dygzz;
o74 = o13 * 0.5 * o31 * dzgyy;
o75 = o13 * -0.5 * o25 * dxgyz;
o76 = o13 * 0.5 * o25 * dygxz;
o77 = o13 * 0.5 * o25 * dzgxy;
o78 = o73 + o74 + o75 + o76 + o77;
o79 = o13 * 0.5 * o72 * dzgzz;
o80 = o13 * -0.5 * o31 * dygzz;
o81 = o13 * o31 * dzgyz;
o82 = o13 * -0.5 * o25 * dxgzz;
o83 = o13 * o25 * dzgxz;
o84 = o79 + o80 + o81 + o82 + o83;
o85 = o13 * 0.5 * o31 * dxgzz;
o86 = o13 * 0.5 * o20 * dxgyz;
o87 = o13 * -0.5 * o20 * dygxz;
o88 = o13 * 0.5 * o20 * dzgxy;
o89 = o13 * 0.5 * o17 * dzgxx;
o90 = o85 + o86 + o87 + o88 + o89;
o91 = o13 * 0.5 * o72 * dxgzz;
o92 = o13 * 0.5 * o31 * dxgyz;
o93 = o13 * -0.5 * o31 * dygxz;
o94 = o13 * 0.5 * o31 * dzgxy;
o95 = o13 * 0.5 * o25 * dzgxx;
o96 = o91 + o92 + o93 + o94 + o95;
o97 = o25 * o27 * o31;
o98 = -(o17 * o27 * o72);
o99 = o97 + o98;
/* Computing 2nd power */
tmpsqr = o25;
o100 = tmpsqr * tmpsqr;
o101 = o27 * o100;
o102 = -(o4 * o27 * o72);
o103 = o101 + o102;
o104 = o17 * o25 * o27;
o105 = -(o4 * o27 * o31);
o106 = o104 + o105;
o107 = o13 * o25 * dygyz;
o108 = o13 * -0.5 * o25 * dzgyy;
o109 = o13 * 0.5 * o17 * dygyy;
o110 = o4 * -0.5 * o13 * dxgyy;
o111 = o4 * o13 * dygxy;
o112 = o107 + o108 + o109 + o110 + o111;
o113 = o13 * o31 * dygyz;
o114 = o13 * -0.5 * o31 * dzgyy;
o115 = o13 * 0.5 * o20 * dygyy;
o116 = o13 * -0.5 * o17 * dxgyy;
o117 = o13 * o17 * dygxy;
o118 = o113 + o114 + o115 + o116 + o117;
o119 = o13 * 0.5 * o25 * dxgyz;
o120 = o13 * -0.5 * o25 * dzgxy;
o121 = o13 * 0.5 * o17 * dxgyy;
o122 = o4 * 0.5 * o13 * dygxx;
o123 = o76 + o119 + o120 + o121 + o122;
o124 = o13 * o72 * dygyz;
o125 = o13 * -0.5 * o72 * dzgyy;
o126 = o13 * 0.5 * o31 * dygyy;
o127 = o13 * -0.5 * o25 * dxgyy;
o128 = o13 * o25 * dygxy;
o129 = o124 + o125 + o126 + o127 + o128;
o130 = o13 * 0.5 * o31 * dygxz;
o131 = o13 * -0.5 * o31 * dzgxy;
o132 = o13 * 0.5 * o20 * dxgyy;
o133 = o13 * 0.5 * o17 * dygxx;
o134 = o92 + o130 + o131 + o132 + o133;
o135 = o13 * 0.5 * o72 * dxgyz;
o136 = o13 * 0.5 * o72 * dygxz;
o137 = o13 * -0.5 * o72 * dzgxy;
o138 = o13 * 0.5 * o31 * dxgyy;
o139 = o13 * 0.5 * o25 * dygxx;
o140 = o135 + o136 + o137 + o138 + o139;
/* Computing 2nd power */
tmpsqr = o31;
o141 = tmpsqr * tmpsqr;
o142 = -(o27 * o141);
o143 = o20 * o27 * o72;
o144 = o142 + o143;
o145 = -(o25 * o27 * o31);
o146 = o17 * o27 * o72;
o147 = o145 + o146;
o148 = o13 * o25 * dxgxz;
o149 = o13 * -0.5 * o25 * dzgxx;
o150 = o13 * o17 * dxgxy;
o151 = o13 * -0.5 * o17 * dygxx;
o152 = o4 * 0.5 * o13 * dxgxx;
o153 = o148 + o149 + o150 + o151 + o152;
o154 = o13 * o31 * dxgxz;
o155 = o13 * -0.5 * o31 * dzgxx;
o156 = o13 * o20 * dxgxy;
o157 = o13 * -0.5 * o20 * dygxx;
o158 = o13 * 0.5 * o17 * dxgxx;
o159 = o154 + o155 + o156 + o157 + o158;
o160 = o13 * o72 * dxgxz;
o161 = o13 * -0.5 * o72 * dzgxx;
o162 = o13 * o31 * dxgxy;
o163 = o13 * -0.5 * o31 * dygxx;
o164 = o13 * 0.5 * o25 * dxgxx;
o165 = o160 + o161 + o162 + o163 + o164;
/* Computing 2nd power */
tmpsqr = kkxy;
o166 = (tmpsqr * tmpsqr);

eleczz = o4*o13*o14 + o13*o20*o21 + 2.0*o13*o17*kkxz*kkyz - o4*o13*kkxx*kkzz - 2.0*o13*o17*kkxy*kkzz - o13*o20*kkyy*kkzz - riczz_loc;

magzz = o22*o40*dxkyz + o22*o36*dxkzz - o22*o40*dykxz + o22*o33*dykzz - o22*o36*dzkxz - o22*o33*dzkyz + o22*o40*o46*kkxx + o22*o36*o52*kkxx + o22*o33*o52*kkxy + o22*o40*o58*kkxy + o22*o36*o64*kkxy - o22*o40*o69*kkxy - o22*o33*o46*kkxz - o22*o36*o69*kkxz + o22*o40*o78*kkxz + o22*o36*o84*kkxz + o22*o33*o64*kkyy - o22*o40*o90*kkyy - o22*o33*o58*kkyz + o22*o33*o84*kkyz - o22*o36*o90*kkyz - o22*o40*o96*kkyz - o22*o33*o78*kkzz - o22*o36*o96*kkzz;

elecyz = -(o13*o21*o31) + o4*o13*kkxy*kkxz + o13*o17*kkxz*kkyy - o4*o13*kkxx*kkyz - o13*o17*kkxy*kkyz - o13*o25*kkxz*kkyz + o13*o25*kkxy*kkzz + o13*o31*kkyy*kkzz - ricyz_loc;

t1 = 0.5*o22*o40*dxkyy + 0.5*o22*o36*dxkyz + 0.5*o22*o106*dxkyz + 0.5*o22*o103*dxkzz - 0.5*o22*o40*dykxy - 0.5*o22*o106*dykxz + 0.5*o22*o33*dykyz + 0.5*o22*o99*dykzz - 0.5*o22*o36*dzkxy - 0.5*o22*o103*dzkxz - 0.5*o22*o33*dzkyy - 0.5*o22*o99*dzkyz + 0.5*o22*o36*o46*kkxx + 0.5*o22*o52*o103*kkxx + 0.5*o22*o46*o106*kkxx + 0.5*o22*o40*o112*kkxx + 0.5*o22*o33*o46*kkxy + 0.5*o22*o36*o58*kkxy + 0.5*o22*o52*o99*kkxy + 0.5*o22*o64*o103*kkxy + 0.5*o22*o58*o106*kkxy - 0.5*o22*o69*o106*kkxy + 0.5*o22*o40*o118*kkxy - 0.5*o22*o40*o123*kkxy;

t2 = 0.5*o22*o36*o78*kkxz - 0.5*o22*o46*o99*kkxz - 0.5*o22*o69*o103*kkxz + 0.5*o22*o84*o103*kkxz + 0.5*o22*o78*o106*kkxz - 0.5*o22*o33*o112*kkxz - 0.5*o22*o36*o123*kkxz + 0.5*o22*o40*o129*kkxz + 0.5*o22*o33*o58*kkyy + 0.5*o22*o64*o99*kkyy - 0.5*o22*o90*o106*kkyy - 0.5*o22*o40*o134*kkyy + 0.5*o22*o33*o78*kkyz - 0.5*o22*o58*o99*kkyz + 0.5*o22*o84*o99*kkyz - 0.5*o22*o90*o103*kkyz - 0.5*o22*o96*o106*kkyz - 0.5*o22*o33*o118*kkyz - 0.5*o22*o36*o134*kkyz - 0.5*o22*o40*o140*kkyz - 0.5*o22*o78*o99*kkzz - 0.5*o22*o96*o103*kkzz - 0.5*o22*o33*o129*kkzz - 0.5*o22*o36*o140*kkzz;

magyz = t1 + t2;

elecxz = -(o13*o14*o25) - o13*o17*kkxy*kkxz - o13*o20*kkxz*kkyy + o13*o17*kkxx*kkyz + o13*o20*kkxy*kkyz - o13*o31*kkxz*kkyz + o13*o25*kkxx*kkzz + o13*o31*kkxy*kkzz - ricxz_loc;

t1 = 0.5*o22*o40*dxkxy + 0.5*o22*o36*dxkxz + 0.5*o22*o33*dxkyz + 0.5*o22*o147*dxkzz - 0.5*o22*o40*dykxx + 0.5*o22*o144*dykzz - 0.5*o22*o36*dzkxx - 0.5*o22*o33*dzkxy - 0.5*o22*o147*dzkxz - 0.5*o22*o144*dzkyz + 0.5*o22*o33*o46*kkxx + 0.5*o22*o36*o69*kkxx + 0.5*o22*o40*o123*kkxx + 0.5*o22*o52*o147*kkxx + 0.5*o22*o33*o58*kkxy + 0.5*o22*o36*o90*kkxy + 0.5*o22*o40*o134*kkxy + 0.5*o22*o52*o144*kkxy + 0.5*o22*o64*o147*kkxy - 0.5*o22*o40*o153*kkxy;

t2 = 0.5*o22*o33*o78*kkxz + 0.5*o22*o36*o96*kkxz - 0.5*o22*o33*o123*kkxz + 0.5*o22*o40*o140*kkxz - 0.5*o22*o46*o144*kkxz - 0.5*o22*o69*o147*kkxz + 0.5*o22*o84*o147*kkxz - 0.5*o22*o36*o153*kkxz + 0.5*o22*o64*o144*kkyy - 0.5*o22*o40*o159*kkyy - 0.5*o22*o33*o134*kkyz - 0.5*o22*o58*o144*kkyz + 0.5*o22*o84*o144*kkyz - 0.5*o22*o90*o147*kkyz - 0.5*o22*o36*o159*kkyz - 0.5*o22*o40*o165*kkyz - 0.5*o22*o33*o140*kkzz - 0.5*o22*o78*o144*kkzz - 0.5*o22*o96*o147*kkzz - 0.5*o22*o36*o165*kkzz;

magxz = t1 + t2;

elecyy = o13*o21*o72 + o4*o13*o166 - o4*o13*kkxx*kkyy - 2.0*o13*o25*kkxz*kkyy + 2.0*o13*o25*kkxy*kkyz - o13*o72*kkyy*kkzz - ricyy_loc;

magyy = o22*o106*dxkyy + o22*o103*dxkyz - o22*o106*dykxy + o22*o99*dykyz - o22*o103*dzkxy - o22*o99*dzkyy + o22*o46*o103*kkxx + o22*o106*o112*kkxx + o22*o46*o99*kkxy + o22*o58*o103*kkxy + o22*o106*o118*kkxy - o22*o106*o123*kkxy + o22*o78*o103*kkxz - o22*o99*o112*kkxz - o22*o103*o123*kkxz + o22*o106*o129*kkxz + o22*o58*o99*kkyy - o22*o106*o134*kkyy + o22*o78*o99*kkyz - o22*o99*o118*kkyz - o22*o103*o134*kkyz - o22*o106*o140*kkyz - o22*o99*o129*kkzz - o22*o103*o140*kkzz;

elecxy = -(o13*o17*o166) - o13*o25*kkxy*kkxz + o13*o17*kkxx*kkyy + o13*o31*kkxz*kkyy + o13*o25*kkxx*kkyz - o13*o31*kkxy*kkyz + o13*o72*kkxz*kkyz - o13*o72*kkxy*kkzz - ricxy_loc;

t1 = 0.5*o22*o106*dxkxy + 0.5*o22*o103*dxkxz + 0.5*o22*o33*dxkyy + 0.5*o22*o147*dxkyz - 0.5*o22*o106*dykxx - 0.5*o22*o33*dykxy + 0.5*o22*o99*dykxz + 0.5*o22*o144*dykyz - 0.5*o22*o103*dzkxx - 0.5*o22*o99*dzkxy - 0.5*o22*o147*dzkxy - 0.5*o22*o144*dzkyy + 0.5*o22*o69*o103*kkxx + 0.5*o22*o33*o112*kkxx + 0.5*o22*o106*o123*kkxx + 0.5*o22*o46*o147*kkxx + 0.5*o22*o69*o99*kkxy + 0.5*o22*o90*o103*kkxy + 0.5*o22*o33*o118*kkxy - 0.5*o22*o33*o123*kkxy + 0.5*o22*o106*o134*kkxy + 0.5*o22*o46*o144*kkxy + 0.5*o22*o58*o147*kkxy - 0.5*o22*o106*o153*kkxy;

t2 = 0.5*o22*o96*o103*kkxz - 0.5*o22*o99*o123*kkxz + 0.5*o22*o33*o129*kkxz + 0.5*o22*o106*o140*kkxz - 0.5*o22*o112*o144*kkxz + 0.5*o22*o78*o147*kkxz - 0.5*o22*o123*o147*kkxz - 0.5*o22*o103*o153*kkxz + 0.5*o22*o90*o99*kkyy - 0.5*o22*o33*o134*kkyy + 0.5*o22*o58*o144*kkyy - 0.5*o22*o106*o159*kkyy + 0.5*o22*o96*o99*kkyz - 0.5*o22*o99*o134*kkyz - 0.5*o22*o33*o140*kkyz + 0.5*o22*o78*o144*kkyz - 0.5*o22*o118*o144*kkyz - 0.5*o22*o134*o147*kkyz - 0.5*o22*o103*o159*kkyz - 0.5*o22*o106*o165*kkyz - 0.5*o22*o99*o140*kkzz - 0.5*o22*o129*o144*kkzz - 0.5*o22*o140*o147*kkzz - 0.5*o22*o103*o165*kkzz;

magxy = t1 + t2;

elecxx = o13*o14*o72 + o13*o20*o166 + 2.0*o13*o31*kkxy*kkxz - o13*o20*kkxx*kkyy - 2.0*o13*o31*kkxx*kkyz - o13*o72*kkxx*kkzz - ricxx_loc ;

magxx = o22*o33*dxkxy + o22*o147*dxkxz - o22*o33*dykxx + o22*o144*dykxz - o22*o147*dzkxx - o22*o144*dzkxy + o22*o33*o123*kkxx + o22*o69*o147*kkxx + o22*o33*o134*kkxy + o22*o69*o144*kkxy + o22*o90*o147*kkxy - o22*o33*o153*kkxy + o22*o33*o140*kkxz - o22*o123*o144*kkxz + o22*o96*o147*kkxz - o22*o147*o153*kkxz + o22*o90*o144*kkyy - o22*o33*o159*kkyy + o22*o96*o144*kkyz - o22*o134*o144*kkyz - o22*o147*o159*kkyz - o22*o33*o165*kkyz - o22*o140*o144*kkzz - o22*o147*o165*kkzz ;
