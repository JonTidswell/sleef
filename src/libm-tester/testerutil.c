//   Copyright Naoki Shibata and contributors 2010 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include <math.h>

#if defined(__MINGW32__) || defined(__MINGW64__) || defined(_MSC_VER)
#define STDIN_FILENO 0
#else
#include <unistd.h>
#include <sys/types.h>
#endif

#if defined(__MINGW32__) || defined(__MINGW64__)
#include <unistd.h>
#endif

#if defined(_MSC_VER)
#include <io.h>
#endif

#include "misc.h"

#define DENORMAL_DBL_MIN (4.9406564584124654418e-324)
#define POSITIVE_INFINITY INFINITY
#define NEGATIVE_INFINITY (-INFINITY)

#define DENORMAL_FLT_MIN (1.4012984643248170709e-45f)
#define POSITIVE_INFINITYf ((float)INFINITY)
#define NEGATIVE_INFINITYf (-(float)INFINITY)

int isnumber(double x) { return !isinf(x) && !isnan(x); }
int isPlusZero(double x) { return x == 0 && copysign(1, x) == 1; }
int isMinusZero(double x) { return x == 0 && copysign(1, x) == -1; }
double sign(double d) { return d < 0 ? -1 : 1; }
int xisnan(double x) { return x != x; }

int isnumberf(float x) { return !isinf(x) && !isnan(x); }
int isPlusZerof(float x) { return x == 0 && copysignf(1, x) == 1; }
int isMinusZerof(float x) { return x == 0 && copysignf(1, x) == -1; }
float signf(float d) { return d < 0 ? -1 : 1; }
int xisnanf(float x) { return x != x; }

int enableFlushToZero = 0;

double flushToZero(double y) {
  if (enableFlushToZero && fabs(y) < FLT_MIN) y = copysign(0.0, y);
  return y;
}

//

int readln(int fd, char *buf, int cnt) {
  int i, rcnt = 0;

  if (cnt < 1) return -1;

  while(cnt >= 2) {
    i = read(fd, buf, 1);
    if (i != 1) return i;

    if (*buf == '\n') break;

    rcnt++;
    buf++;
    cnt--;
  }

  *++buf = '\0';
  rcnt++;
  return rcnt;
}

static uint64_t xseed;

uint64_t xrand() {
  xseed = xseed * UINT64_C(6364136223846793005) + 1;
  return xseed;
}

// Fill memory with random bits
void memrand(void *p, int size) {
  uint64_t *q = (uint64_t *)p;
  int i;
  for(i=0;i<size/8;i++) *q++ = xrand();
  uint8_t *r = (uint8_t *)q;
  for(i *= 8;i<size;i++) *r++ = xrand() & 0xff;
}

void xsrand(uint64_t s) { xseed = s; }

//

#ifdef USEMPFR
#include <mpfr.h>

int cmpDenormsp(float x, mpfr_t fry) {
  float y = mpfr_get_d(fry, GMP_RNDN);
  x = flushToZero(x);
  y = flushToZero(y);
  if (xisnanf(x) && xisnanf(y)) return 1;
  if (xisnanf(x) || xisnanf(y)) return 0;
  if (isinf(x) != isinf(y)) return 0;
  if (x == POSITIVE_INFINITYf && y == POSITIVE_INFINITYf) return 1;
  if (x == NEGATIVE_INFINITYf && y == NEGATIVE_INFINITYf) return 1;
  if (y == 0) {
    if (isPlusZerof(x) && isPlusZerof(y)) return 1;
    if (isMinusZerof(x) && isMinusZerof(y)) return 1;
    return 0;
  }
  if (!xisnanf(x) && !xisnanf(y) && !isinf(x) && !isinf(y)) return signf(x) == signf(y);
  return 0;
}

int cmpDenormdp(double x, mpfr_t fry) {
  double y = mpfr_get_d(fry, GMP_RNDN);
  if (xisnan(x) && xisnan(y)) return 1;
  if (xisnan(x) || xisnan(y)) return 0;
  if (isinf(x) != isinf(y)) return 0;
  if (x == POSITIVE_INFINITY && y == POSITIVE_INFINITY) return 1;
  if (x == NEGATIVE_INFINITY && y == NEGATIVE_INFINITY) return 1;
  if (y == 0) {
    if (isPlusZero(x) && isPlusZero(y)) return 1;
    if (isMinusZero(x) && isMinusZero(y)) return 1;
    return 0;
  }
  if (!xisnan(x) && !xisnan(y) && !isinf(x) && !isinf(y)) return sign(x) == sign(y);
  return 0;
}

double countULPdp(double d, mpfr_t c) {
  mpfr_t fra, frb, frc, frd;
  mpfr_inits(fra, frb, frc, frd, NULL);

  double c2 = mpfr_get_d(c, GMP_RNDN);
  if (c2 == 0 && d != 0) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 10000;
  }
  if (isnan(c2) && isnan(d)) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 0;
  }
  if (isnan(c2) || isnan(d)) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 10001;
  }
  if (c2 == POSITIVE_INFINITY && d == POSITIVE_INFINITY) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 0;
  }
  if (c2 == NEGATIVE_INFINITY && d == NEGATIVE_INFINITY) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 0;
  }

  double v = 0;
  if (isinf(d) && !isinf(mpfr_get_d(c, GMP_RNDN))) {
    d = copysign(DBL_MAX, c2);
    v = 1;
  }

  //
  
  int e;
  frexp(mpfr_get_d(c, GMP_RNDN), &e);
  mpfr_set_ld(frb, fmaxl(ldexpl(1.0, e-53), DENORMAL_DBL_MIN), GMP_RNDN);

  mpfr_set_d(frd, d, GMP_RNDN);
  mpfr_sub(fra, frd, c, GMP_RNDN);
  mpfr_div(fra, fra, frb, GMP_RNDN);
  double u = fabs(mpfr_get_d(fra, GMP_RNDN));

  mpfr_clears(fra, frb, frc, frd, NULL);
  
  return u + v;
}

double countULP2dp(double d, mpfr_t c) {
  mpfr_t fra, frb, frc, frd;
  mpfr_inits(fra, frb, frc, frd, NULL);

  double c2 = mpfr_get_d(c, GMP_RNDN);
  if (c2 == 0 && d != 0) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 10000;
  }
  if (isnan(c2) && isnan(d)) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 0;
  }
  if (isnan(c2) || isnan(d)) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 10001;
  }
  if (c2 == POSITIVE_INFINITY && d == POSITIVE_INFINITY) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 0;
  }
  if (c2 == NEGATIVE_INFINITY && d == NEGATIVE_INFINITY) {
    mpfr_clears(fra, frb, frc, frd, NULL);
    return 0;
  }

  double v = 0;
  if (isinf(d) && !isinf(mpfr_get_d(c, GMP_RNDN))) {
    d = copysign(DBL_MAX, c2);
    v = 1;
  }

  //

  int e;
  frexp(mpfr_get_d(c, GMP_RNDN), &e);
  mpfr_set_ld(frb, fmaxl(ldexpl(1.0, e-53), DBL_MIN), GMP_RNDN);

  mpfr_set_d(frd, d, GMP_RNDN);
  mpfr_sub(fra, frd, c, GMP_RNDN);
  mpfr_div(fra, fra, frb, GMP_RNDN);
  double u = fabs(mpfr_get_d(fra, GMP_RNDN));

  mpfr_clears(fra, frb, frc, frd, NULL);
  
  return u + v;
}

double countULPsp(float d, mpfr_t c0) {
  double c = mpfr_get_d(c0, GMP_RNDN);

  d = flushToZero(d);
  float c2 = flushToZero(c);
  if (c2 == 0 && d != 0) return 10000;
  if (isnan(c2) && isnan(d)) return 0;
  if (isnan(c2) || isnan(d)) return 10001;
  if (c2 == POSITIVE_INFINITYf && d == POSITIVE_INFINITYf) return 0;
  if (c2 == NEGATIVE_INFINITYf && d == NEGATIVE_INFINITYf) return 0;

  double v = 0;
  if (isinf(d) && !isinf(c)) {
    d = copysign(FLT_MAX, c2);
    v = 1;
  }

  //

  int e;
  frexp(c, &e);

  double u = fabs(d - c) * fmin(ldexp(1.0, 24-e), 1.0 / DENORMAL_FLT_MIN);
  
  return u + v;
}

double countULP2sp(float d, mpfr_t c0) {
  double c = mpfr_get_d(c0, GMP_RNDN);

  d = flushToZero(d);
  float c2 = flushToZero(c);
  if (c2 == 0 && d != 0) return 10000;
  if (isnan(c2) && isnan(d)) return 0;
  if (isnan(c2) || isnan(d)) return 10001;
  if (c2 == POSITIVE_INFINITYf && d == POSITIVE_INFINITYf) return 0;
  if (c2 == NEGATIVE_INFINITYf && d == NEGATIVE_INFINITYf) return 0;

  double v = 0;
  if (isinf(d) && !isinf(c)) {
    d = copysign(FLT_MAX, c2);
    v = 1;
  }

  //

  int e;
  frexp(c, &e);

  double u = fabs(d - c) * fmin(ldexp(1.0, 24-e), 1.0 / FLT_MIN);
  
  return u + v;
}

//

void mpfr_sinpi(mpfr_t ret, mpfr_t arg, mpfr_rnd_t rnd) {
  mpfr_t frpi, frd;
  mpfr_inits(frpi, frd, NULL);

  mpfr_const_pi(frpi, GMP_RNDN);
  mpfr_set_d(frd, 1.0, GMP_RNDN);
  mpfr_mul(frpi, frpi, frd, GMP_RNDN);
  mpfr_mul(frd, frpi, arg, GMP_RNDN);
  mpfr_sin(ret, frd, GMP_RNDN);

  mpfr_clears(frpi, frd, NULL);
}

void mpfr_cospi(mpfr_t ret, mpfr_t arg, mpfr_rnd_t rnd) {
  mpfr_t frpi, frd;
  mpfr_inits(frpi, frd, NULL);

  mpfr_const_pi(frpi, GMP_RNDN);
  mpfr_set_d(frd, 1.0, GMP_RNDN);
  mpfr_mul(frpi, frpi, frd, GMP_RNDN);
  mpfr_mul(frd, frpi, arg, GMP_RNDN);
  mpfr_cos(ret, frd, GMP_RNDN);

  mpfr_clears(frpi, frd, NULL);
}

void mpfr_lgamma_nosign(mpfr_t ret, mpfr_t arg, mpfr_rnd_t rnd) {
  int s;
  mpfr_lgamma(ret, &s, arg, rnd);
}


static double erfinv_base(double w) {
  // polynomial from Giles erfinv gems
  double p;

  if ( w < 6.250 ) {                              
    w = w - 3.1250;                                
    p =  -3.6444120640178196996e-21;                 
    p =   -1.685059138182016589e-19 + p*w;           
    p =   1.2858480715256400167e-18 + p*w;           
    p =    1.115787767802518096e-17 + p*w;           
    p =   -1.333171662854620906e-16 + p*w;           
    p =   2.0972767875968561637e-17 + p*w;           
    p =   6.6376381343583238325e-15 + p*w;           
    p =  -4.0545662729752068639e-14 + p*w;           
    p =  -8.1519341976054721522e-14 + p*w;           
    p =   2.6335093153082322977e-12 + p*w;           
    p =  -1.2975133253453532498e-11 + p*w;           
    p =  -5.4154120542946279317e-11 + p*w;           
    p =    1.051212273321532285e-09 + p*w;           
    p =  -4.1126339803469836976e-09 + p*w;           
    p =  -2.9070369957882005086e-08 + p*w;           
    p =   4.2347877827932403518e-07 + p*w;           
    p =  -1.3654692000834678645e-06 + p*w;           
    p =  -1.3882523362786468719e-05 + p*w;           
    p =    0.0001867342080340571352 + p*w;           
    p =  -0.00074070253416626697512 + p*w;           
    p =   -0.0060336708714301490533 + p*w;           
    p =      0.24015818242558961693 + p*w;           
    p =       1.6536545626831027356 + p*w;           
  }                                                  
  else if ( w < 16.0 ) {                        
    w = sqrt(w) - 3.250;                          
    p =   2.2137376921775787049e-09;                 
    p =   9.0756561938885390979e-08 + p*w;           
    p =  -2.7517406297064545428e-07 + p*w;           
    p =   1.8239629214389227755e-08 + p*w;           
    p =   1.5027403968909827627e-06 + p*w;           
    p =   -4.013867526981545969e-06 + p*w;           
    p =   2.9234449089955446044e-06 + p*w;           
    p =   1.2475304481671778723e-05 + p*w;           
    p =  -4.7318229009055733981e-05 + p*w;           
    p =   6.8284851459573175448e-05 + p*w;           
    p =   2.4031110387097893999e-05 + p*w;           
    p =   -0.0003550375203628474796 + p*w;           
    p =   0.00095328937973738049703 + p*w;           
    p =   -0.0016882755560235047313 + p*w;           
    p =    0.0024914420961078508066 + p*w;           
    p =   -0.0037512085075692412107 + p*w;           
    p =     0.005370914553590063617 + p*w;           
    p =       1.0052589676941592334 + p*w;           
    p =       3.0838856104922207635 + p*w;           
  }                                                  
  else {                                             
    w = sqrt(w) - 5.0;                          
    p =  -2.7109920616438573243e-11;                 
    p =  -2.5556418169965252055e-10 + p*w;           
    p =   1.5076572693500548083e-09 + p*w;           
    p =  -3.7894654401267369937e-09 + p*w;           
    p =   7.6157012080783393804e-09 + p*w;           
    p =  -1.4960026627149240478e-08 + p*w;           
    p =   2.9147953450901080826e-08 + p*w;           
    p =  -6.7711997758452339498e-08 + p*w;           
    p =   2.2900482228026654717e-07 + p*w;           
    p =  -9.9298272942317002539e-07 + p*w;           
    p =   4.5260625972231537039e-06 + p*w;           
    p =  -1.9681778105531670567e-05 + p*w;           
    p =   7.5995277030017761139e-05 + p*w;           
    p =  -0.00021503011930044477347 + p*w;           
    p =  -0.00013871931833623122026 + p*w;           
    p =       1.0103004648645343977 + p*w;           
    p =       4.8499064014085844221 + p*w;           
  }                                                  
                                             
  return p;   
}


void mpfr_erfinv (mpfr_t ret, mpfr_t arg, mpfr_rnd_t rnd) {

  if ( mpfr_nan_p(arg) ) {
    mpfr_set( ret, arg, MPFR_RNDN );
    return;
  }
  if ( mpfr_zero_p(arg) ) {
    mpfr_set( ret, arg, MPFR_RNDN );
    return;
  }

  int s = mpfr_signbit(arg);
  mpfr_t parg, tmp1, tmp2, fval, fderiv, f2deriv;

  mpfr_inits( parg, tmp1, tmp2, fval, fderiv, f2deriv, NULL );
  mpfr_abs( parg, arg, MPFR_RNDN );

  int bounds = mpfr_cmp_ui( parg, 1 );
  if ( bounds == 0 ) {
    mpfr_set_inf( ret, s );
    goto done;
  }
  else if ( bounds > 0 ) {
    // abs(parg) > 1
    mpfr_set_nan( ret );
    goto done;
  }

  // our initialisation breaks down for p near zero
  // due to underflow in Float64
  // so test for it and exploit linearity ...
  mpfr_set_ui_2exp( tmp1, 1, -20, MPFR_RNDN );
  if ( mpfr_less_p( parg, tmp1 ) ) {
    mpfr_mul_d( ret, parg, 0.8862269254527583, MPFR_RNDN );
  }
  else {
    // calc these bits of giles in extended precision
    // avoids SOME numerical under/over-flow issues
    // t = (1.0 - x) * (1.0 + x);
    mpfr_sqr( tmp1, parg, MPFR_RNDN );
    mpfr_si_sub( tmp1, 1, tmp1, MPFR_RNDN );
    // w = - xlog_u05( t );    
    mpfr_log( tmp1, tmp1, MPFR_RNDN );
    mpfr_neg( tmp1, tmp1, MPFR_RNDN );

    // evaluate the polynomial in floating point
    double w = mpfr_get_d( ret, MPFR_RNDN );
    double r = erfinv_base( w );

    // result = x * poly
    mpfr_mul_d( ret, parg, r, MPFR_RNDN );
  }
 
  // now use Halley Method to refine the solution
  mpfr_erf( fval, ret, MPFR_RNDN );
  mpfr_sub( fval, fval, parg, MPFR_RNDN );

  mpfr_sqr( tmp1, ret, MPFR_RNDN );
  mpfr_neg( tmp1, tmp1, MPFR_RNDN );
  mpfr_exp( tmp1, tmp1, MPFR_RNDN );
  mpfr_mul_ui( tmp1, tmp1, 2,MPFR_RNDN );
  mpfr_const_pi( tmp2, MPFR_RNDN );
  mpfr_rec_sqrt( tmp2, tmp2, MPFR_RNDN );
  mpfr_mul( fderiv, tmp1, tmp2, MPFR_RNDN );
  
  mpfr_mul_si( f2deriv, fderiv, -2, MPFR_RNDN );
  mpfr_mul( f2deriv, f2deriv, ret, MPFR_RNDN );

  mpfr_mul( tmp1, fval, f2deriv, MPFR_RNDN );
  mpfr_sqr( tmp2, fderiv, MPFR_RNDN );
  mpfr_mul_ui( tmp2, tmp2, 2, MPFR_RNDN );
  mpfr_sub( tmp2, tmp2, tmp1, MPFR_RNDN );
  mpfr_mul( tmp1, fval, fderiv, MPFR_RNDN );
  mpfr_mul_ui( tmp1, tmp1, 2, MPFR_RNDN );
  mpfr_div( tmp1, tmp1, tmp2, MPFR_RNDN );

  mpfr_sub( ret, ret, tmp1, MPFR_RNDN );  

  mpfr_setsign( ret, ret, s, MPFR_RNDN );
done:
  mpfr_clears( parg, tmp1, tmp2, fval, fderiv, f2deriv, NULL );
  return;
}

#endif // #define USEMPFR
