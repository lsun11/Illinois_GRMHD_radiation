/* ricci.f -- translated by f2c (version 20050501).
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
double sqrtmp;

/* Computing 2nd power */
sqrtmp = uyy;
t1 = sqrtmp * sqrtmp;
t2 = t1 * dygyy;
t4 = uyz * uxx;
t5 = dxgxx * dygxz;
t7 = uzz * uxy;
t8 = dxgxx * dzgyz;
t11 = uyy * uxy;
t12 = dxgxy * dygxy;
t14 = dxgxx * dygyy;
/* Computing 2nd power */
sqrtmp = uyz;
t17 = sqrtmp * sqrtmp;
t18 = t17 * dygzz;
/* Computing 2nd power */
sqrtmp = uzz;
t20 = sqrtmp * sqrtmp;
t21 = t20 * dzgzz;
t25 = uyy * uxx;
/* Computing 2nd power */
sqrtmp = dygxx;
t26 = sqrtmp * sqrtmp;
t29 = t2 * dygxx / 4. - t4 * t5 / 2. - t7 * t8 / 2. - uzz * d2zz_gxx / 
  2. - t11 * t12 - t11 * t14 / 4. + uyz * d2xy_gxz - t18 * dxgxy 
  - t21 * dxgxz / 2. - t17 * dxgxz * dzgyy + t25 * t26 / 4. - t2 * 
  dxgxy / 2.;
t30 = uyy * uzz;
/* Computing 2nd power */
sqrtmp = dygxz;
t31 = sqrtmp * sqrtmp;
/* Computing 2nd power */
sqrtmp = dxgyz;
t33 = sqrtmp * sqrtmp;
/* Computing 2nd power */
sqrtmp = dzgxy;
t35 = sqrtmp * sqrtmp;
t39 = t17 * dxgzz;
t41 = t17 * dygxz;
t45 = uzz * uxx;
/* Computing 2nd power */
sqrtmp = dzgxx;
t46 = sqrtmp * sqrtmp;
t48 = uyy * uxz;
t49 = dxgxx * dygyz;
t54 = t30 * t31 / 2. + t30 * t33 / 2. + t30 * t35 / 2. + t17 * dzgxx * 
  dzgyy / 2. + t39 * dxgyy / 2. + t41 * dzgxy + t18 * dygxx / 2. + 
  t21 * dzgxx / 4. + t45 * t46 / 4. - t48 * t49 / 2. + t17 * t33 / 
  2. + uyz * d2xz_gxy - uyy * d2yy_gxx / 2.;
/* Computing 2nd power */
sqrtmp = dxgyy;
t56 = sqrtmp * sqrtmp;
t60 = dygxx * dxgyy;
t62 = dxgxx * dxgyy;
t64 = dxgxx * dygxy;
t66 = uyz * uxz;
t67 = dzgxx * dygxz;
t69 = uyz * uzz;
t70 = dxgzz * dxgyz;
t72 = uyz * uxy;
t73 = dxgxy * dzgxy;
t75 = dxgxy * dygxz;
t77 = dxgxy * dxgyz;
t79 = uyy * uyz;
t80 = dxgyy * dxgyz;
t82 = t1 * t56 / 4. - uyy * d2xx_gyy / 2. + uyy * d2xy_gxy + t11 * 
  t60 / 4. + t25 * t62 / 4. - t25 * t64 / 2. + t66 * t67 + t69 * 
  t70 - t72 * t73 - t72 * t75 + t72 * t77 + t79 * t80;
t83 = dxgxy * dygyz;
t85 = dygzz * dzgxx;
t88 = dygzz * dxgxz;
t90 = dxgzz * dygxx;
t92 = dxgxz * dzgxy;
/* Computing 2nd power */
sqrtmp = dxgzz;
t94 = sqrtmp * sqrtmp;
t99 = dxgxx * dzgxy;
t101 = dxgxz * dxgyz;
t103 = dzgxx * dxgyy;
t105 = -t79 * t83 + t69 * t85 / 4. - uzz * d2xx_gzz / 2. - t69 * t88 / 
  2. + t66 * t90 / 2. - t66 * t92 + t20 * t94 / 4. - uyz * 
  d2yz_gxx + uzz * d2xz_gxz - t17 * t35 / 2. - t4 * t99 / 2. + 
  t66 * t101 + t72 * t103 / 2.;
t108 = dxgxx * dzgyy;
t110 = dxgxx * dygzz;
t112 = dygxx * dzgxy;
t114 = dxgxz * dygyz;
t116 = dxgxz * dzgyy;
t118 = dzgxx * dygyz;
t120 = dxgxx * dxgzz;
t122 = dzgxx * dzgxy;
t124 = dzgxx * dxgyz;
t126 = dygxx * dzgyz;
t128 = dzgxx * dzgyz;
t130 = dzgzz * dxgxy;
t132 = -t72 * t108 / 2. - t66 * t110 / 2. + t72 * t112 - t30 * t114 + t30 
  * t116 / 2. + t30 * t118 / 2. + t45 * t120 / 4. + t7 * t122 / 2. 
  + t7 * t124 / 2. + t30 * t126 / 2. + t69 * t128 / 2. - t69 * t130 
  / 2.;
t133 = dzgzz * dygxx;
t135 = dxgxy * dxgzz;
t137 = dygxx * dzgxz;
t140 = dxgxy * dzgyz;
t142 = uzz * uxz;
t143 = dxgxz * dzgxz;
t145 = dxgxz * dxgzz;
t147 = dzgxx * dzgxz;
t149 = dygzz * dygxx;
t151 = dxgxy * dzgxz;
t153 = dxgxz * dzgyz;
t156 = dzgxx * dxgzz;
t158 = t69 * t133 / 4. + t7 * t135 / 2. + t7 * t137 / 2. - t7 * t90 / 4. 
  - t30 * t140 - t142 * t143 + t142 * t145 / 2. + t142 * t147 / 2. 
  - t30 * t149 / 4. - t7 * t151 - t69 * t153 + t7 * t110 / 4. + 
  t142 * t156 / 4.;
t160 = dxgxx * dzgzz;
t163 = dzgxx * dzgyy;
t165 = dygzz * dxgxy;
t167 = dxgxx * dzgxz;
t169 = dxgxx * dxgyz;
t171 = dxgxz * dygxz;
t173 = dzgxx * dygxx;
t175 = dygxx * dygxz;
t177 = dzgxx * dygxy;
t179 = dygxx * dygyz;
t181 = dygxx * dzgyy;
t183 = dxgxy * dxgyy;
t185 = -t142 * t160 / 4. - t7 * t67 / 2. - t30 * t163 / 4. + t30 * t165 / 
  2. - t45 * t167 / 2. + t4 * t169 / 2. - t66 * t171 + t4 * t173 / 
  2. + t48 * t175 / 2. + t48 * t177 / 2. + t79 * t179 / 2. + t79 * 
  t181 / 4. + t11 * t183 / 2.;
t186 = dygxx * dygxy;
t188 = dxgxy * dzgyy;
t190 = dygxx * dxgyz;
t194 = dygyy * dxgxz;
t196 = dygyy * dzgxx;
t198 = dygxz * dzgxy;
t200 = dxgxz * dygxy;
t202 = dxgxz * dxgyy;
t207 = t11 * t186 / 2. - t79 * t188 / 2. + t48 * t190 / 2. - t48 * t112 / 
  2. - t48 * t103 / 4. - t79 * t194 / 2. + t79 * t196 / 4. - t30 * 
  t198 - t48 * t200 + t48 * t202 / 2. - uyz * d2xx_gyz + t48 * 
  t108 / 4. - t17 * t31 / 2.;
ricxx_loc = t29 + t54 + t82 + t105 + t132 + t158 + t185 + t207;
t210 = t17 * dxgyz;
/* Computing 2nd power */
sqrtmp = uxz;
t213 = sqrtmp * sqrtmp;
t214 = t213 * dzgxx;
/* Computing 2nd power */
sqrtmp = uxy;
t216 = sqrtmp * sqrtmp;
t217 = t216 * dxgxy;
t219 = t17 * dzgxy;
t222 = t216 * dxgxx;
t224 = t216 * dygxx;
t230 = t213 * dxgzz;
t234 = -t210 * dzgyy / 2. - t11 * t56 / 4. - t214 * dygxz / 2. + t217 * 
  dygxy + t219 * dygyz / 2. - t21 * dygxz / 4. + t222 * dygyy / 4. 
  - t224 * dygxy / 2. - t21 * dxgyz / 4. + t210 * dygyz / 2. - uxz *
  d2xz_gxy / 2. + uyz * d2xz_gyy / 2. - t230 * dygxx / 4. + 
  t21 * dzgxy / 4. + t39 * dygyy / 4.;
t235 = dygxz * dzgxz;
t238 = uxy * uxx;
t249 = dxgyz * dzgyz;
t251 = dzgxy * dzgyz;
t253 = dxgyy * dzgzz;
t255 = dzgxy * dzgyy;
t257 = -t142 * t235 / 2. + t66 * t116 / 2. - t238 * t26 / 4. - t224 * 
  dxgyy / 4. - t7 * t31 / 2. - uyz * d2xy_gyz / 2. - uxy * 
  d2xy_gxy + uzz * d2xz_gyz / 2. + uzz * d2yz_gxz / 2. + uxz *
  d2yz_gxx / 2. - uxz * d2xy_gxz / 2. + t45 * t90 / 4. - t69 * 
  t249 / 2. + t69 * t251 / 2. - t69 * t253 / 4. + t30 * t255 / 4.;
t259 = dygxz * dzgyz;
t261 = dxgyz * dzgxz;
t263 = dxgyy * dzgxz;
t265 = uxz * uxx;
t270 = dxgyy * dzgxy;
t272 = uxz * uxy;
t274 = dzgxy * dzgxz;
t276 = dxgyz * dzgyy;
t283 = -t69 * t259 / 2. - t142 * t261 / 2. - t7 * t263 / 2. - t265 * t173 
  / 4. + t265 * t99 / 4. + t265 * t5 / 4. - t265 * t169 / 4. - t48 *
  t270 / 4. + t272 * t73 / 2. + t142 * t274 / 2. + t30 * t276 / 4. 
  - uzz * d2zz_gxy / 2. - uzz * d2xy_gzz / 2. + t45 * t122 / 4. 
  + t45 * t67 / 4. - t45 * t124 / 4.;
t285 = dxgzz * dxgyy;
t292 = t213 * dxgxz;
t296 = t213 * dxgxx;
t301 = t17 * dxgyy;
t305 = -t45 * t137 / 2. + t7 * t285 / 4. - t7 * t126 / 2. + t7 * t149 / 
  4. + uxy * d2xx_gyy / 2. - t7 * t33 / 2. - t217 * dxgyy / 2. + 
  t292 * dzgxy / 2. + t292 * dygxz / 2. - t292 * dxgyz / 2. + t296 *
  dygzz / 4. + t20 * dxgzz * dygzz / 4. + t66 * t35 / 2. - t301 * 
  dygzz / 4. + uxy * d2yy_gxx / 2. - t7 * t116 / 2.;
t309 = dxgyz * dzgxy;
t316 = dxgyy * dygxz;
t325 = -t142 * t133 / 4. + t7 * t309 / 2. + t7 * t114 + t72 * t194 / 2. + 
  t7 * t198 / 2. + t272 * t200 - t272 * t177 / 2. - t72 * t316 / 2. 
  - t72 * t80 / 2. + t72 * t83 - t72 * t188 / 2. - t72 * t179 / 2. 
  - t7 * t118 / 2. + t7 * t163 / 2. - t272 * t202 / 2.;
t330 = dygyy * dygxx;
t332 = dxgzz * dygxy;
t340 = dygyy * dzgxy;
t342 = dygyy * dxgyz;
t346 = t238 * t64 / 2. - t272 * t175 / 2. + t272 * t49 / 2. - t238 * t62 /
  4. - t11 * t330 / 4. + t66 * t332 / 2. - t66 * t285 / 2. - t66 * 
  t309 / 2. - t4 * t112 / 4. + t4 * t190 / 4. + t66 * t118 / 2. - 
  t66 * t163 / 2. + t79 * t340 / 4. + t79 * t342 / 4. + t4 * t177 / 
  2. - t4 * t103 / 4.;
t349 = dygyy * dxgxy;
t355 = dzgxy * dygxy;
t357 = dxgyz * dygxy;
t359 = dygxz * dygxy;
t361 = dygyy * dygxz;
t363 = dxgyy * dzgyy;
t366 = dxgzz * dygyz;
t368 = dygzz * dygxz;
t370 = dygzz * dxgyz;
t372 = dygzz * dzgxy;
t374 = -uyz * d2yz_gxy / 2. + t11 * t349 / 2. - t41 * dygyz / 2. - t272 
  * t190 / 2. + uxz * d2xx_gyz / 2. - t66 * t149 / 2. + t72 * 
  t355 / 2. + t72 * t357 / 2. - t72 * t359 / 2. - t79 * t361 / 4. - 
  t79 * t363 / 4. - t4 * t175 / 4. + t69 * t366 / 2. - t69 * t368 / 
  4. + t69 * t370 / 4. + t69 * t372 / 4.;
t377 = dxgzz * dygxz;
t381 = dxgzz * dzgxy;
t390 = dygxz * dzgyy;
t392 = dxgyy * dzgyz;
t394 = dxgyy * dygzz;
t396 = uyz * d2yy_gxz / 2. - t272 * t77 / 2. + t142 * t377 / 4. - t142 *
  t70 / 4. + t272 * t75 / 2. + t142 * t381 / 4. + t142 * t88 / 2. 
  + t66 * t165 / 2. - t66 * t198 / 2. - t48 * t80 / 4. + t48 * t316 
  / 4. + t48 * t188 / 2. - t48 * t181 / 4. - t30 * t390 / 4. - t30 *
  t392 / 2. + t30 * t394 / 4.;
ricxy_loc = t234 + t257 + t283 + t305 + t325 + t346 + t374 + t396;
t414 = t217 * dzgxy / 2. - t7 * t85 / 4. + t7 * t88 / 2. + t41 * dzgyz / 
  2. + t222 * dzgyy / 4. - t18 * dxgyz / 2. + t292 * dzgxz - t2 * 
  dzgxy / 4. - t214 * dxgzz / 4. - t142 * t94 / 4. - t265 * t46 / 
  4. - t292 * dxgzz / 2. - t79 * t255 / 4. + uyz * d2xy_gzz / 2. 
  + t66 * t153;
t418 = dxgzz * dygzz;
t433 = uxz * d2zz_gxx / 2. - uyz * d2xz_gyz / 2. + uyz * d2zz_gxy / 
  2. - t69 * t418 / 4. + t301 * dzgzz / 4. - t39 * dzgyy / 4. + 
  t210 * dzgyz / 2. - t219 * dzgyz / 2. - t224 * dzgxy / 2. + t217 *
  dygxz / 2. - t217 * dxgyz / 2. + t72 * t31 / 2. - t2 * dxgyz / 
  4. + t1 * dxgyy * dzgyy / 4. - t48 * t35 / 2. - t48 * t33 / 2.;
t451 = -t214 * dzgxz / 2. - uyz * d2yz_gxz / 2. - uyy * d2yy_gxz / 2. 
  - uyy * d2xz_gyy / 2. + uyy * d2yz_gxy / 2. + uyy * 
  d2xy_gyz / 2. + uxy * d2yz_gxx / 2. - uxy * d2xy_gxz / 2. - 
  uxy * d2xz_gxy / 2. + uxy * d2xx_gyz / 2. + uxz * d2xx_gzz /
  2. - uxz * d2xz_gxz - t238 * t173 / 4. + t66 * t130 / 2. + 
  t272 * t151 - t272 * t135 / 2.;
t453 = dygxz * dxgyz;
t469 = -t272 * t137 / 2. + t48 * t453 / 2. + t48 * t140 - t48 * t165 / 2. 
  - t48 * t126 / 2. + t265 * t167 / 2. - t265 * t120 / 4. - t272 * 
  t122 / 2. - t272 * t124 / 2. + t272 * t8 / 2. - t4 * t67 / 4. + 
  t4 * t124 / 4. + t4 * t137 / 2. + t72 * t126 / 2. - t72 * t453 / 
  2. + t66 * t235 / 2.;
t476 = dzgzz * dxgyz;
t478 = dzgzz * dzgxy;
t485 = dzgzz * dxgxz;
t487 = dzgzz * dzgxx;
t491 = t66 * t261 / 2. - t66 * t274 / 2. + t79 * t390 / 4. + t79 * t392 / 
  2. + t69 * t476 / 4. - t69 * t478 / 4. - t4 * t122 / 4. + t48 * 
  t198 / 2. - t66 * t381 / 2. - t66 * t88 / 2. - t66 * t128 / 2. + 
  t142 * t485 / 2. - t142 * t487 / 4. - t66 * t70 / 2. + t72 * t263 
  / 2.;
t493 = dzgzz * dygxz;
t497 = t216 * dzgxx;
t510 = -t4 * t90 / 4. + t69 * t493 / 4. + t296 * dzgzz / 4. + t2 * dygxz /
  4. - t497 * dxgyy / 4. + t272 * t92 / 2. - t72 * t198 / 2. - t11 
  * t80 / 4. + t11 * t316 / 4. + t11 * t188 / 2. - t72 * t149 / 2. 
  + t72 * t116 / 2. - t72 * t163 / 2. + t72 * t165 / 2. - t7 * t377 
  / 4. - t72 * t285 / 2.;
t519 = dygxz * dygyz;
t529 = t272 * t171 / 2. - t272 * t101 / 2. + t11 * t270 / 4. + t238 * t99 
  / 4. + t238 * t5 / 4. - t238 * t169 / 4. + t79 * t276 / 4. + t79 *
  t519 / 2. + t25 * t112 / 4. - t25 * t190 / 4. - t11 * t196 / 4. 
  - t48 * t118 / 2. + t48 * t163 / 4. + t48 * t149 / 2. - t11 * 
  t355 / 2. - t11 * t357 / 2.;
t531 = dxgyz * dygyz;
t534 = dzgxy * dygyz;
t539 = dxgzz * dygyy;
t542 = dxgzz * dzgyy;
t550 = t11 * t359 / 2. - t79 * t531 / 2. + t25 * t175 / 4. - t79 * t534 / 
  2. - t25 * t177 / 2. + t25 * t103 / 4. - t48 * t332 / 2. - t79 * 
  t539 / 4. - t30 * t366 / 2. + t30 * t542 / 4. + t30 * t368 / 4. + 
  t30 * t370 / 4. - t30 * t372 / 4. + t48 * t285 / 4. - t7 * t70 / 
  4. + t7 * t381 / 4.;
ricxz_loc = t414 + t433 + t451 + t469 + t491 + t510 + t529 + t550;
/* Computing 2nd power */
sqrtmp = dzgyy;
t560 = sqrtmp * sqrtmp;
t563 = t213 * dygzz;
t567 = -t4 * t83 - t4 * t270 / 2. + t4 * t80 / 2. - t4 * t194 / 2. - t265 
  * t177 / 2. + t265 * t202 / 2. - t265 * t200 + t30 * t560 / 4. + 
  t21 * dzgyy / 4. + t563 * dygxx / 2. + t214 * dzgyy / 2. + t4 * 
  t196 / 4.;
t569 = t213 * dxgyz;
/* Computing 2nd power */
sqrtmp = uxx;
t577 = sqrtmp * sqrtmp;
t581 = t577 * dxgxx;
t584 = -t214 * dygyz + t569 * dzgxy + t230 * dxgyy / 2. - uxx * 
  d2xx_gyy / 2. - t230 * dygxy + t25 * t56 / 4. + t45 * t33 / 2. 
  + uxx * d2xy_gxy + t577 * t26 / 4. + t45 * t35 / 2. + t45 * t31 
  / 2. - t581 * dygxy / 2. + t581 * dxgyy / 4.;
t598 = -t4 * t181 / 4. + t4 * t179 / 2. + t4 * t316 / 2. - t45 * t163 / 
  4. + t45 * t118 / 2. + t45 * t116 / 2. - t45 * t114 + t4 * t188 / 
  2. - t45 * t309 + t238 * t186 / 2. + t238 * t183 / 2. - uxz * 
  d2yy_gxz;
t612 = t265 * t103 / 4. - uxx * d2yy_gxx / 2. + uxz * d2xy_gyz - t48 *
  t342 / 2. - t213 * t33 / 2. + t213 * t31 / 2. - t213 * t35 / 2. 
  + uxz * d2yz_gxy - uxz * d2xz_gyy - t48 * t340 / 2. + t265 * 
  t175 - t238 * t14 / 4. + t238 * t60 / 4.;
t622 = dygxy * dzgxz;
t624 = dygyy * dygzz;
t626 = dygyy * dzgyz;
t630 = -t238 * t12 + t265 * t108 / 4. - t265 * t49 / 2. + t25 * t330 / 4. 
  - t25 * t349 / 2. - t272 * t355 + t45 * t332 / 2. - t45 * t622 + 
  t30 * t624 / 4. - t30 * t626 / 2. + t7 * t539 / 4. - t66 * t534;
t631 = dygyy * dzgzz;
t634 = dygxy * dygzz;
t636 = dzgyy * dygzz;
t638 = dzgyy * dzgyz;
t640 = dygyz * dzgyz;
t642 = dygyz * dygzz;
t644 = dzgyy * dzgxz;
t650 = dygyz * dzgxz;
t652 = -t69 * t631 / 4. + t7 * t255 / 2. + t7 * t634 / 2. + t69 * t636 / 
  4. + t69 * t638 / 2. - t69 * t640 + t69 * t642 / 2. + t142 * t644 
  / 2. - t7 * t276 / 2. + t7 * t390 / 2. - t45 * t285 / 4. + t45 * 
  t263 / 2. - t142 * t650;
t656 = dygxy * dzgyz;
t660 = dygyy * dzgxz;
t663 = dzgzz * dygxy;
t670 = -t7 * t394 / 4. + t7 * t392 / 2. - t7 * t656 - t272 * t357 + t272 *
  t270 - t7 * t660 / 2. + t142 * t253 / 4. - t142 * t663 / 2. - 
  t66 * t531 + t48 * t363 / 2. + t48 * t361 / 2. + t272 * t359 + 
  uzz * d2yz_gyz;
/* Computing 2nd power */
sqrtmp = dygzz;
t673 = sqrtmp * sqrtmp;
t685 = -uzz * d2yy_gzz / 2. - uzz * d2zz_gyy / 2. + t20 * t673 / 4. - 
  t142 * t366 / 2. - t66 * t539 / 2. + t66 * t519 + t66 * t276 + 
  t66 * t394 / 2. + t142 * t368 + t142 * t542 / 4. + t272 * t181 / 
  2. - t272 * t196 / 2. - t21 * dygyz / 2.;
ricyy_loc = t567 + t584 + t598 + t612 + t630 + t652 + t670 + t685;
t702 = t216 * dygxz;
t704 = t216 * dxgyy;
t709 = -t72 * t255 / 2. - t66 * t372 / 2. - t72 * t634 / 2. - t79 * t560 /
  4. - t563 * dzgxx / 4. + t213 * dzgzz * dygxx / 4. - t230 * 
  dygxz / 2. - t213 * dzgxy * dzgxz / 2. + t569 * dzgxz / 2. + t213 
  * dygxz * dzgxz / 2. + t272 * t33 / 2. - t702 * dygxy / 2. - t704 
  * dzgxy / 2. + t216 * dzgxy * dygxy / 2. - t224 * dzgyy / 4.;
t717 = t17 * dzgyy;
t729 = t216 * dygyy * dzgxx / 4. - t4 * t31 / 2. + t581 * dxgyz / 4. - 
  t581 * dygxz / 4. + t577 * dzgxx * dygxx / 4. - t717 * dygzz / 4. 
  + t142 * t493 / 4. + t272 * t263 / 2. - t69 * t673 / 4. + t66 * 
  t259 / 2. - t48 * t255 / 4. - t72 * t392 / 2. + t72 * t656 - t72 *
  t390 / 2. - t4 * t263 / 2. + t66 * t650;
t731 = dzgzz * dzgyy;
t733 = dzgzz * dygyz;
t749 = -t69 * t731 / 4. + t69 * t733 / 2. - t66 * t644 / 2. - t66 * t366 /
  2. - t48 * t394 / 4. + t48 * t392 / 2. + t48 * t390 / 4. - t48 * 
  t276 / 4. - t142 * t418 / 4. - t66 * t251 / 2. + t66 * t249 / 2. 
  - t4 * t35 / 2. - t142 * t478 / 4. + t142 * t476 / 4. - t272 * 
  t453 / 2. - t717 * dzgyz / 2.;
t753 = t17 * dygyz;
t768 = t272 * t126 / 2. + t265 * t137 / 2. + t265 * t124 / 4. + t753 * 
  dzgyz - t265 * t122 / 4. - t753 * dygzz / 2. + t272 * t118 / 2. + 
  t238 * t190 / 4. + t17 * dygyy * dzgzz / 4. - uyz * d2yz_gyz + 
  uyz * d2yy_gzz / 2. + uxz * d2zz_gxy / 2. - t238 * t175 / 4. 
  + t72 * t531 / 2. + uyz * d2zz_gyy / 2. - t11 * t363 / 4.;
t786 = -t11 * t361 / 4. - t272 * t149 / 2. - t272 * t163 / 2. + t7 * t366 
  / 2. + t272 * t332 / 2. + t238 * t177 / 2. + t11 * t342 / 4. + 
  t72 * t534 / 2. + t11 * t340 / 4. - t72 * t519 / 2. - t272 * t309 
  / 2. - t272 * t285 / 2. + t7 * t372 / 4. - t7 * t370 / 4. - t7 * 
  t368 / 4.;
t803 = -t7 * t542 / 4. - t238 * t75 / 2. - t238 * t73 / 2. + t238 * t112 /
  4. - t265 * t110 / 4. - t238 * t108 / 4. + t25 * t270 / 4. + 
  t265 * t101 / 2. - t265 * t171 / 2. + t45 * t70 / 4. + t265 * t67 
  / 4. + t4 * t163 / 4. - t4 * t116 / 2. + t238 * t77 / 2. + t45 * 
  t85 / 4. + t45 * t377 / 4.;
t821 = -t45 * t381 / 4. - t25 * t188 / 2. - t265 * t92 / 2. - t45 * t88 / 
  2. + t4 * t285 / 2. - t79 * t624 / 4. + t79 * t626 / 2. + t72 * 
  t660 / 2. + t66 * t663 / 2. - t66 * t368 / 2. + t25 * t181 / 4. + 
  t25 * t80 / 4. + t4 * t622 + t4 * t453 / 2. + t4 * t309 / 2. - t4 
  * t332 / 2.;
t839 = -t4 * t165 / 2. - t581 * dzgxy / 4. - t25 * t316 / 4. - uxy * 
  d2xy_gyz / 2. - uxx * d2yz_gxx / 2. - uxz * d2xz_gyz / 2. - 
  uxz * d2yz_gxz / 2. + uxz * d2xy_gzz / 2. + t4 * t149 / 4. + 
  uxy * d2yy_gxz / 2. + uxx * d2xz_gxy / 2. + uxx * d2xy_gxz /
  2. - uxx * d2xx_gyz / 2. - uxy * d2yz_gxy / 2. + uxy * 
  d2xz_gyy / 2. + t216 * dxgyz * dygxy / 2.;
ricyz_loc = t709 + t729 + t749 + t768 + t786 + t803 + t821 + t839;
t854 = -t48 * t542 / 4. + t48 * t644 / 2. + t48 * t366 / 2. + t581 * 
  dxgzz / 4. - t79 * t631 / 4. - t79 * t640 + t79 * t642 / 2. + t79 
  * t638 / 2. + t72 * t542 / 2. - t272 * t133 / 2. + t11 * t394 / 
  4. + t272 * t85 / 2.;
t869 = -t11 * t392 / 2. + t79 * t636 / 4. + t272 * t377 - t11 * t656 + 
  t30 * t673 / 4. + t25 * t263 / 2. - t25 * t285 / 4. + t216 * 
  dxgzz * dxgyy / 2. + t30 * t731 / 4. - t48 * t650 + t7 * t418 / 
  2. + uyy * d2yz_gyz - t704 * dzgxz;
t883 = t72 * t251 - uyy * d2zz_gyy / 2. - t72 * t249 + t7 * t478 / 2. - 
  t7 * t476 / 2. + t1 * t560 / 4. - t72 * t259 - t7 * t493 / 2. + 
  t72 * t370 + t11 * t255 + t497 * dzgyy / 2. - t216 * t31 / 2.;
t897 = -t216 * t33 / 2. + t272 * t274 - uxy * d2zz_gxy + t702 * dxgyz - 
  t272 * t235 - t272 * t261 - t72 * t253 / 2. - uyy * d2yy_gzz / 
  2. - t224 * dzgyz + t45 * t94 / 4. + t25 * t31 / 2. - t581 * 
  dzgxz / 2. - t2 * dzgyz / 2.;
t912 = t25 * t33 / 2. + t4 * t381 / 2. - t30 * t733 / 2. + t2 * dygzz / 
  4. - t4 * t377 / 2. + t4 * t70 / 2. + t45 * t487 / 4. - t45 * 
  t485 / 2. - t4 * t85 / 4. + t4 * t128 / 2. - t4 * t153 + t4 * t88 
  / 2.;
t927 = t216 * dygzz * dygxx / 2. + t265 * t147 / 2. - t25 * t149 / 4. + 
  t25 * t35 / 2. + t265 * t145 / 2. - t265 * t143 + t4 * t133 / 4. 
  - t4 * t130 / 2. - t25 * t453 + t238 * t90 / 4. - t238 * t137 / 
  2. + t238 * t135 / 2. - t238 * t151;
t942 = -t265 * t160 / 4. + t216 * t35 / 2. + t265 * t156 / 4. + t238 * 
  t110 / 4. + uxx * d2xz_gxz - t238 * t8 / 2. + t238 * t122 + t25 
  * t126 / 2. + t25 * t165 / 2. - t25 * t140 - uxx * d2xx_gzz / 
  2. + uxy * d2yz_gxz - uxy * d2xy_gzz;
t956 = uxy * d2xz_gyz - uxx * d2zz_gxx / 2. + t577 * t46 / 4. + t48 * 
  t368 / 2. + t11 * t539 / 4. - t11 * t660 / 2. + t48 * t253 / 4. - 
  t48 * t663 / 2. - t48 * t370 / 2. + t11 * t634 / 2. + t25 * t332 /
  2. - t25 * t622 + t48 * t372 / 2.;
riczz_loc = t854 + t869 + t883 + t897 + t912 + t927 + t942 + t956;
t960 = t216 * uxy;
t961 = t960 * dygxx;
t965 = t216 * uzz;
t967 = t216 * uxx;
t972 = t960 * dxgxy;
t976 = t272 * 2. * d2xx_gyz - t961 * dygxy - t72 * 2. * d2xy_gyz - t4 
  * 2. * d2xx_gyz - t965 * 1.5 * t31 - t967 * t26 / 2. + t960 * 
  dxgxx * dygyy / 2. - t961 * dxgyy / 2. + t972 * 2. * dygxy + t72 *
  2. * d2yy_gxz - t7 * 2. * d2zz_gxy;
t977 = uzz * t577;
t980 = t17 * uyz;
t981 = t980 * dygyz;
t983 = t17 * uzz;
t985 = t980 * dzgyy;
t993 = uyy * t20;
t995 = t977 * t46 / 2. + t272 * 2. * d2yz_gxx + t981 * 2. * dzgyz - 
  t983 * t673 / 2. - t985 * dygzz / 2. + t48 * 2. * d2xy_gyz + 
  t48 * 2. * d2yz_gxy - t48 * 2. * d2yy_gxz - t66 * 2. * 
  d2yz_gxz - t30 * d2zz_gyy + t30 * 2. * d2yz_gyz + t993 * 
  t673 / 2.;
t997 = uzz * t1;
t999 = t213 * uxz;
t1000 = t999 * dzgxx;
t1005 = uxx * t1;
t1008 = uyy * t577;
t1010 = uyy * t213;
t1013 = t213 * uxx;
t1015 = t997 * t560 / 2. - t1000 * dxgzz / 2. - t45 * d2zz_gxx + t25 * 
  2. * d2xy_gxy - t25 * d2xx_gyy + t1005 * t56 / 2. - t25 * 
  d2yy_gxx + t1008 * t26 / 2. - t1010 * 1.5 * t33 - t1010 * 1.5 * 
  t35 - t1013 * t46 / 2.;
t1023 = t216 * uyy;
t1027 = uxx * t20;
t1030 = t72 * -2. * d2yz_gxy - t66 * 2. * d2xz_gyz + t72 * 2. * 
  d2xz_gyy - t272 * 2. * d2xy_gxz + t7 * 2. * d2yz_gxz + t7 * 
  2. * d2xz_gyz - t972 * dxgyy - t1023 * t56 / 2. - t965 * 1.5 * 
  t33 - t4 * 2. * d2yz_gxx + t1027 * t94 / 2. + t4 * 2. * 
  d2xz_gxy;
t1034 = t17 * uyy;
t1036 = uxx * t17;
t1038 = uyz * t577;
t1040 = uyz * dzgyy;
t1050 = uxz * dzgxx;
t1054 = -t45 * d2xx_gzz - t1034 * t560 / 2. - t1036 * 1.5 * t31 + t1038 
  * t173 + t30 * t1040 * dygzz / 2. - t981 * dygzz + t999 * dxgxx * 
  dzgzz / 2. - t30 * 2. * uxx * dygxy * dzgxz + t216 * d2yy_gxx + 
  t72 * 2. * t1050 * dygyz - t30 * d2yy_gzz;
t1055 = t999 * dxgxz;
t1058 = t213 * uzz;
t1061 = uxy * dxgxx;
t1064 = uxy * dxgxy;
t1073 = uyy * dygzz;
t1077 = -t1055 * dxgzz + t66 * 2. * d2zz_gxy - t1058 * t94 / 2. + t965 *
  t35 / 2. - t25 * t1061 * dygyy / 2. - t45 * 2. * t1064 * dzgxz + 
  t4 * 2. * d2xy_gxz + t45 * t1050 * dzgxz + t45 * 2. * 
  d2xz_gxz + t25 * 1.5 * uzz * t31 - t45 * t1073 * dygxx / 2. + 
  t1036 * t33 / 2.;
t1081 = uzz * dzgxx;
t1086 = uzz * dxgxz;
t1092 = uxz * dxgxz;
t1095 = t213 * uyz;
t1097 = uyz * dygyy;
t1103 = t1005 * t330 / 2. - t1005 * t349 - t25 * t1081 * dzgyy / 2. + t25 
  * t1081 * dygyz + t25 * t1086 * dzgyy + t1055 * 2. * dzgxz - t25 *
  2. * t1086 * dygyz + t25 * t1092 * dxgyy + t1095 * t130 + t25 * 
  t1097 * dzgxx / 2. - t25 * uzz * dygxz * dzgxy;
t1108 = uxy * dygxx;
t1113 = uyy * dygxx;
t1119 = uyy * dxgxy;
t1120 = t1119 * dzgyy;
t1122 = uyy * dxgyy;
t1123 = t1122 * dxgyz;
t1127 = uzz * dxgzz;
t1128 = t1127 * dxgyz;
t1130 = uzz * dygzz;
t1131 = t1130 * dxgxz;
t1135 = t25 * -2. * t1092 * dygxy + t25 * t1064 * dxgyy + t25 * t1108 * 
  dygxy - t25 * t1097 * dxgxz + t4 * t1113 * dzgyy / 2. + t25 * uyz 
  * dygxx * dygyz - t4 * t1120 + t4 * 2. * t1123 - t4 * 2. * t1092 *
  dzgxy + t4 * 2. * t1128 - t4 * t1131 + t4 * t1130 * dzgxx / 2.;
t1145 = uxz * dxgxx;
t1162 = t4 * -2. * t1064 * dygxz + t4 * 2. * t1064 * dxgyz - t4 * 2. * 
  t1064 * dzgxy - t4 * t1145 * dygzz + t25 * t1108 * dxgyy / 2. - 
  t25 * 2. * uyz * dxgxy * dygyz - t25 * 2. * t1064 * dygxy + t25 * 
  t1145 * dzgyy / 2. + t980 * dygyy * dzgzz / 2. - t25 * t1145 * 
  dygyz - t1027 * t485;
t1171 = uyz * dzgzz;
t1176 = uyz * dxgxz;
t1182 = uyz * t216;
t1187 = t45 * t1113 * dzgyz + t45 * t1073 * dxgxy - t45 * 2. * t1119 * 
  dzgyz + t45 * t1064 * dxgzz - t45 * t1171 * dxgxy + t45 * t1171 * 
  dygxx / 2. - t45 * 2. * t1176 * dzgyz + t45 * uyz * dzgxx * dzgyz 
  - t1182 * t270 + t17 * d2zz_gyy - t1182 * t181 / 2. + t1182 * 
  t196 / 2.;
t1197 = t17 * uxy;
t1200 = t17 * uxz;
t1202 = -t1038 * t99 - t1036 * t332 - t1038 * t5 + t1038 * t169 + t1036 * 
  2. * t622 + t1036 * t453 - t1034 * t624 / 2. + t1036 * t309 + 
  t1197 * t660 + t1034 * t626 - t1200 * t368;
t1218 = t1200 * t663 + t1036 * 1.5 * t285 - t1036 * 2. * t116 + t1036 * 
  t163 - t7 * 2. * d2xy_gzz + t4 * 2. * t1092 * dxgyz + t977 * 
  t120 / 2. - t4 * 2. * t1092 * dygxz - t4 * t1061 * dzgyy - t977 * 
  t167 - t1008 * t64 + t1008 * t62 / 2.;
t1221 = uyz * dxgyz;
t1227 = uxy * dygxz;
t1231 = uyz * dygxz;
t1234 = uxy * dygyy;
t1244 = uxz * t216;
t1246 = t48 * -2. * t1221 * dygyz - t48 * 2. * uxy * dxgyz * dygxy + t48 *
  2. * t1227 * dygxy + t66 * 2. * d2xy_gzz + t48 * 2. * t1231 * 
  dygyz - t48 * t1234 * dzgxx - t48 * 2. * uxy * dzgxy * dygxy + 
  t48 * t1221 * dzgyy - t48 * 2. * d2xz_gyy - t1200 * t372 - 
  t1244 * t103 / 2.;
t1249 = t1130 * dygxz;
t1254 = t1127 * dygyz;
t1263 = t1036 * t149 - t1036 * 2. * t165 + t48 * 2. * t1249 + t48 * t1127 
  * dzgyy / 2. - t985 * dzgyz - t48 * t1254 + t1095 * t133 / 2. - 
  t17 * 2. * d2yz_gyz - t1095 * t377 - t1095 * t85 / 2. + t72 * 
  uxz * t33 - t1244 * t112;
t1268 = uzz * dzgzz;
t1271 = uxx * dzgxx;
t1274 = uxx * dygxx;
t1275 = t1274 * dygxz;
t1277 = uyz * dzgxy;
t1280 = t1271 * dygxy;
t1282 = t1122 * dzgyz;
t1284 = uyy * dygxz;
t1292 = t66 * 2. * uxy * dxgyy * dzgxz + t66 * t1268 * dygxz + t48 * 
  t1271 * dxgyy / 2. + t48 * 2. * t1275 - t48 * 2. * t1277 * dygyz 
  - t48 * t1280 + t66 * 2. * t1282 + t66 * t1284 * dzgyy - t66 * 
  t1127 * dygzz + t213 * d2xx_gzz - t66 * 2. * t1227 * dxgyz;
t1295 = t1274 * dzgxz;
t1304 = uyy * dzgxy * dzgyy;
t1309 = uxz * t1;
t1312 = t66 * 2. * t1108 * dzgyz + t66 * 2. * t1295 - t1200 * t542 / 2. + 
  t66 * t1271 * dygxz + t66 * t1271 * dxgyz + t1200 * t253 / 2. - 
  t66 * t1304 + t1010 * t285 - t1010 * 2. * t332 + t1010 * 1.5 * 
  t149 - t1309 * t340 + t1244 * t108 / 2.;
t1331 = t1010 * -2. * t118 + t1010 * t163 + t1309 * t361 - t1309 * t342 + 
  t1309 * t363 + t272 * uyz * t31 - t1200 * t366 - t1197 * t255 - 
  t48 * uyz * dxgzz * dygyy - t1200 * 2. * t251 - t1197 * t634;
t1337 = t1271 * dzgxy;
t1340 = uxx * dxgxx;
t1349 = t17 * d2yy_gzz - t213 * 2. * d2xz_gxz + t1200 * 2. * t249 + 
  t1200 * 2. * t259 - t1200 * t370 - t66 * t1337 + t1010 * t31 / 2. 
  + t272 * t1340 * dygxz + t66 * t1268 * dxgyz - t66 * t1268 * 
  dzgxy - t272 * 2. * d2xz_gxy - t1197 * t392;
t1360 = t213 * uxy;
t1363 = t1197 * 2. * t656 - t1036 * 1.5 * t35 - t1036 * t263 - t1197 * 
  t390 + t1200 * 2. * t650 - t983 * t731 / 2. + t983 * t733 - t1200 
  * t644 + t1010 * t453 - t1360 * t137 - t1360 * t135;
t1380 = t1360 * 2. * t151 - t1360 * 2. * t101 + t1360 * t110 / 2. + t1360 
  * 2. * t171 + t1360 * 2. * t92 + t272 * 2. * t1120 + t272 * t1122 
  * dygxz - t272 * t1123 - t272 * 2. * t1231 * dzgxy + t272 * 2. * 
  t1131 + t272 * 2. * uyz * dygzz * dxgxy - t1197 * t394 / 2.;
t1391 = uyy * dygyy;
t1398 = t1197 * -2. * t519 + t72 * uxz * t35 - t1197 * t276 + t1197 * 2. *
  t534 + t1197 * t539 / 2. + t72 * 2. * t1280 + t1095 * 2. * t153 
  + t72 * t1391 * dzgxy + t72 * t1391 * dxgyz - t1095 * t381 - 
  t1095 * 2. * t274;
t1411 = t1010 * t198 + t1095 * 2. * t261 + t1095 * 2. * t235 + t1360 * t8 
  - t1360 * t124 - t1360 * t122 - t1013 * t120 / 2. + t1013 * t167 
  - t1010 * t126 - t1010 * t165 + t1010 * 2. * t140 - t1095 * t70;
t1424 = uxz * dxgyz;
t1429 = uxz * dxgzz;
t1432 = -t1058 * t487 / 2. - t1095 * t128 + t1058 * t485 - t1095 * t88 + 
  t272 * t1127 * dzgxy + t1197 * 2. * t531 - t72 * 3. * t1050 * 
  dzgyy + t72 * t1274 * dxgyz - t72 * 2. * t1424 * dzgxy + t72 * 
  t1274 * dzgxy - t72 * 3. * t1429 * dxgyy;
t1449 = uxy * t20;
t1452 = t72 * 2. * t1429 * dygxy - t1360 * t90 / 2. - t1360 * t67 - t72 * 
  t1249 + t72 * t1130 * dxgyz + t72 * 2. * t1254 - t72 * t1122 * 
  dzgyy - t72 * t1275 - t72 * t1391 * dygxz - t72 * 3. * uxz * 
  dygzz * dygxx + t1449 * t418 - t1449 * t476;
t1465 = uxx * dxgzz;
t1477 = t1449 * t478 - t7 * 2. * uxz * dygxz * dzgxz - t7 * uyz * dxgyy * 
  dzgzz + t7 * 2. * t1277 * dzgyz + t7 * t1465 * dygxx / 2. - t1000 
  * dzgxz - t7 * 2. * t1221 * dzgyz - t7 * 2. * t1231 * dzgyz + t7 *
  2. * t1304 + t7 * t1122 * dygzz / 2. - t7 * t1282;
t1491 = uxz * dzgzz;
t1495 = t7 * -2. * t1424 * dzgxz + t7 * 2. * uxz * dzgxy * dzgxz - t7 * 
  t1295 + t7 * 2. * t1337 + t1023 * t349 + t1182 * 2. * t355 - 
  t1182 * 2. * t359 + t1182 * 2. * t357 + t213 * d2zz_gxx - t1449 
  * t493 - t7 * t1491 * dygxx + t1244 * 2. * t200;
t1509 = t965 * t198 + t965 * 1.5 * t163 + t1182 * t194 + t72 * t1130 * 
  dzgxy - t1182 * t188 + t1182 * 2. * t83 - t1182 * t179 - t1182 * 
  t80 - t1182 * t316 - t1244 * t177 - t1244 * t175;
t1519 = uxy * dygxy;
t1532 = -t1244 * t202 + t967 * t64 - t967 * t62 / 2. - t1244 * t190 + 
  t1244 * t49 - t1023 * t330 / 2. + t30 * uxz * dxgyy * dzgzz / 2. 
  - t30 * 2. * t1519 * dzgyz - t30 * t1491 * dygxy - t30 * 2. * uxz 
  * dygyz * dzgxz + t30 * uxx * dxgyy * dzgxz - t30 * t1465 * dxgyy 
  / 2.;
t1537 = uyz * dygyz;
t1555 = t30 * t1040 * dzgyz + t30 * t1537 * dygzz - t30 * 2. * t1537 * 
  dzgyz + t30 * uxz * dzgyy * dzgxz + t30 * t1519 * dygzz - t30 * 
  t1097 * dzgzz / 2. - t45 * t1284 * dxgyz + t965 * t453 - t1244 * 
  2. * t77 + t1244 * 2. * t75 + t1244 * 2. * t73;
t1571 = t965 * -2. * t263 + t965 * t149 + t965 * t285 - t965 * 2. * t126 
  - t965 * t116 - t216 * 2. * d2xy_gxy + t965 * 2. * t114 + t965 *
  t309 - t272 * t1128 + t272 * t1127 * dygxz + t272 * 2. * t1176 * 
  dzgyy + t272 * t1122 * dzgxy;
t1594 = t216 * d2xx_gyy - t272 * t1340 * dxgyz + t272 * t1340 * dzgxy - 
  t272 * t1271 * dygxx - t965 * t118 - t45 * t1145 * dzgzz / 2. + 
  t45 * t1061 * dygzz / 2. + t45 * t1050 * dxgzz / 2. - t45 * 2. * 
  t1092 * dzgxz + t45 * t1092 * dxgzz + t1027 * t487 / 2. - t45 * 
  t1061 * dzgyz;
t1615 = t1036 * t198 - t993 * t733 + t25 * 1.5 * uzz * t33 + t25 * 1.5 * 
  uzz * t35 + t30 * t1465 * dygxy + t993 * t731 / 2. + t30 * uxy * 
  dxgzz * dygyy / 2. - t25 * uzz * dxgyz * dzgxy + t1010 * t309 + 
  t997 * t624 / 2. - t997 * t626 - t30 * t1234 * dzgxz;
ricscal_loc = t1030 + t1380 + t1187 + t1263 + t1246 + t976 + t995 + 
  t1015 + t1312 + t1103 + t1077 + t1162 + t1135 + t1363 + t1202 + 
  t1218 + t1292 + t1331 + t1349 + t1411 + t1571 + t1532 + t1555 + 
  t1594 + t1398 + t1509 + t1054 + t1432 + t1452 + t1477 + t1495 + 
  t1615;
