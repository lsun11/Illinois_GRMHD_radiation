fp t202;
fp t206;
fp t207;
fp t223;
fp t205;
fp t160;
fp t201;
fp t231;
fp t174;
fp t242;
fp t218;
fp t179;
fp t241;
fp t198;
fp t194;
fp t178;
fp t197;
fp t240;
fp t228;
fp t180;
fp t156;
fp t239;
fp t165;
fp t238;
fp t196;
fp t175;
fp t237;
fp t172;
fp t161;
fp t236;
fp t193;
fp t188;
fp t235;
fp t166;
fp t182;
fp t234;
fp t189;
fp t177;
fp t233;
fp t232;
fp t230;
fp t229;
fp t227;
fp t226;
fp t225;
fp t224;
fp t222;
fp t203;
fp t221;
fp t220;
fp t204;
fp t219;
fp t171;
fp t187;
fp t217;
fp t216;
fp t167;
fp t168;
fp t158;
fp t162;
fp t215;
fp t214;
fp t213;
fp t212;
fp t211;
fp t181;
fp t170;
fp t159;
fp t210;
fp t163;
fp t184;
fp t186;
fp t209;
fp t173;
fp t185;
fp t155;
fp t208;
fp t183;
fp t176;
      t202 = RATIONAL(-1.0,2.0);
      t206 = z*z;
      t207 = y*y;
      t223 = t207*t206;
      t205 = x*t223;
      t160 = t202*t205;
      t201 = RATIONAL(1.0,2.0);
      t231 = x*t201;
      t174 = t206*t231;
      t242 = t160+t174;
      t218 = t201*t207;
      t179 = x*t218;
      t241 = t179+t160;
      t198 = RATIONAL(-1.0,4.0);
      t194 = t198*t206;
      t178 = y*t194;
      t197 = RATIONAL(1.0,4.0);
      t240 = t178+t197*y;
      t228 = x*t202;
      t180 = y*t228;
      t156 = y*t174;
      t239 = t180+t156;
      t165 = t206*t180;
      t238 = y*t231+t165;
      t196 = t197*t206;
      t175 = y*t196;
      t237 = t198*y+t175;
      t172 = z*t228;
      t161 = z*t179;
      t236 = t172+t161;
      t193 = t198*t207;
      t188 = t206*t193;
      t235 = t188+t196;
      t166 = t207*t172;
      t182 = z*t231;
      t234 = t166+t182;
      t189 = t197*t207;
      t177 = t206*t189;
      t233 = t177+t194;
      t232 = y*z;
      t230 = t207*z;
      t229 = t202-x;
      t227 = y*t206;
      t226 = t201-x;
      t225 = x+t205;
      t224 = x*t232;
      t222 = z*t180+t160;
      t203 = RATIONAL(-1.0,8.0);
      t221 = t203*t207;
      t220 = t160+y*t182;
      t204 = RATIONAL(1.0,8.0);
      t219 = t204*t206;
      t171 = z*t193;
      t187 = z*t189;
      t217 = -y+t227;
      t216 = -t230+z;
      t167 = t198*t224;
      t168 = x*t175;
      t158 = x*t171;
      t162 = x*t177;
      t215 = t167+t168+t158+t162;
      t214 = t171+t197*z+t242;
      t213 = t187+t198*z+t242;
      t212 = t188+t189+t241;
      t211 = t177+t193+t241;
      t181 = y*t219;
      t170 = t204*t230;
      t159 = t197*t224;
      t210 = t181+t170+t159+t162;
      t163 = x*t178;
      t184 = t204*t232;
      t186 = z*t221;
      t209 = t163+t184+t186+t162;
      t173 = t203*t232;
      t185 = t203*t227;
      t155 = x*t187;
      t208 = t173+t185+t155+t162;
      t183 = t207*t219;
      t176 = t206*t221;
      coeffs_dx->coeff_m1_m1_m1 = t163+t176+t158+t173+t210;
      coeffs_dx->coeff_0_m1_m1 = t156+t161+t222;
      coeffs_dx->coeff_p1_m1_m1 = t158+t159+t185+t183+t209;
      coeffs_dx->coeff_m1_0_m1 = t214+t233+t236;
      coeffs_dx->coeff_0_0_m1 = t205+(-t206+t216)*x;
      coeffs_dx->coeff_p1_0_m1 = t213+t235+t236;
      coeffs_dx->coeff_m1_p1_m1 = t184+t176+t185+t170+t215;
      coeffs_dx->coeff_0_p1_m1 = t165+t161+t220;
      coeffs_dx->coeff_p1_p1_m1 = t173+t183+t181+t186+t215;
      coeffs_dx->coeff_m1_m1_0 = t211+t239+t240;
      coeffs_dx->coeff_0_m1_0 = t205+(-t207-t217)*x;
      coeffs_dx->coeff_p1_m1_0 = t212+t237+t239;
      coeffs_dx->coeff_m1_0_0 = t202+t226*t207+(t202*t207+t226)*t206+t225;
      coeffs_dx->coeff_0_0_0 = ((1.0+t223)*RATIONAL(-2.0,1.0)+(t207+t206)*
RATIONAL(2.0,1.0))*x;
      coeffs_dx->coeff_p1_0_0 = t201+t229*t207+(t218+t229)*t206+t225;
      coeffs_dx->coeff_m1_p1_0 = t211+t237+t238;
      coeffs_dx->coeff_0_p1_0 = t205+(-t207+t217)*x;
      coeffs_dx->coeff_p1_p1_0 = t212+t238+t240;
      coeffs_dx->coeff_m1_m1_p1 = t167+t181+t155+t176+t209;
      coeffs_dx->coeff_0_m1_p1 = t156+t166+t220;
      coeffs_dx->coeff_p1_m1_p1 = t163+t183+t167+t170+t208;
      coeffs_dx->coeff_m1_0_p1 = t213+t233+t234;
      coeffs_dx->coeff_0_0_p1 = t205+(-t206-t216)*x;
      coeffs_dx->coeff_p1_0_p1 = t214+t234+t235;
      coeffs_dx->coeff_m1_p1_p1 = t168+t176+t159+t186+t208;
      coeffs_dx->coeff_0_p1_p1 = t165+t166+t222;
      coeffs_dx->coeff_p1_p1_p1 = t168+t155+t183+t184+t210;
