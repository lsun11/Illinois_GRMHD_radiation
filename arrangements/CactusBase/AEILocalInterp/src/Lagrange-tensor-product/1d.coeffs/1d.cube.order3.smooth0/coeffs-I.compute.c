fp t10;
fp t9;
fp t8;
fp t7;
fp t6;
fp t5;
      t10 = x*x;
      t9 = x*t10;
      t8 = RATIONAL(-1.0,2.0);
      t7 = RATIONAL(1.0,2.0);
      t6 = RATIONAL(-1.0,6.0);
      t5 = t7*t10;
      coeffs_I->coeff_m1 = RATIONAL(-1.0,3.0)*x+t5+t6*t9;
      coeffs_I->coeff_0 = RATIONAL(1.0,1.0)+t8*x-t10+t7*t9;
      coeffs_I->coeff_p1 = x+t5+t8*t9;
      coeffs_I->coeff_p2 = t6*x+RATIONAL(1.0,6.0)*t9;