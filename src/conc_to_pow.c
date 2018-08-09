#include "system_includes.h"
#include "conc_to_pow.h"
double conc_to_pow(double conc, double pow, double factorial) {
  /*
    raise conc a double representing a concentration or a count to
    a possibly non integer power, and use factorial rising if factorial > 0 or factorial falling if factorial < 0 in computing the power.
    Factorial  = 0.0 just computes conc^pow.
    Note that for factorial_falling the first term is just conc, while for factorial rising the first term is 
    (conc + factorial).
    Called by: compute_dfdke_dfdmu0,
               compute_kss,
	       rxn_likelihood,
	       rxn_likelihood_postselection,
	       lr8_approximate_ys0,
	       lr8_approximate_jacobian,
	       lr10_gradient,
	       lr11_gradient,
	       lr12_gradient,
	       lr13_gradient,
	       lr14_gradient,
	       lr4_gradient,
	       lr7_gradient,
	       lr8_gradient,
	       lr9_gradient

  */
  int64_t int_part;
  double  fract_pow;
  double  abs_pow;
  double  term;
  double  result;
  int     i;
  int     padi;
  abs_pow = pow;
  if (pow < 0.0) {
    abs_pow = 0.0 - pow;
  }
  int_part = (int64_t)abs_pow;
  result = 1.0;
  term = conc;
  fract_pow = (double)int_part - abs_pow;
  if (factorial > 0.0) {
    term += factorial;
  }
  for (i=0;i<int_part;i++) {
    result = result*term;
    term += factorial;
  }
  if (fract_pow > 0.0) {
    /*
      pow is not an integer
      we want to multiply the result by conc^fract_pos = 
      But we need conc to be strictly > 0 as we take its log.
    */
    if (conc <= 0.0) {
      fprintf(stderr,"conc_to_pow: fractional power and conc was <= 0.0, using 0.0\n");
      fflush(stderr);
      result = 0.0;
    } else {
      result = result * exp(fract_pow * log(conc));
    }
  }
  if (pow < 0) {
    if (result <= 0.0) {
      fprintf(stderr,"conc_to_pow: negative power and result was 0, using 0.0\n");
      fflush(stderr);
      result = 0.0;
    } else {
      result = 1.0/result;
    }
  }
  return(result);
}
