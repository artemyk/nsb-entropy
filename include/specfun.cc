/*
 *  File - specfun.cc
 *  
 *  Definition of special functions not present in GSL
 *
 *  Copyright (c) 2003. Ilya Nemenman. All rights reserved.
 *
 */

#include "specfun.h"

/*
  Subleading asymptotic behavior of Psi (or Digamma) function valid 
  for positive argument.
                         d
          psi_asymp(x) = --log(Gamma(x)) - log(x)
                         dx
  We aim at extremely large x, so we use the Stirling form for 
  psi, and not the Lancocs one. The stirling form has the log(x) 
  term in it, and subtraction of the logarithm can be done with no
  loss of precision.
      psi_asymp(x) = psi(x) - log(x) = -1/2x 
                - \sum_{i=1}^{\infty} B_{2i} /(2i*x^{2i}),
*/
double psi_asymp(const double& x){
  const double b[]= {		// even bernoulli coeffs; b_1=-0.5; not shown here
    0.166666666666667, -0.0333333333333333, 0.0238095238095238, -0.0333333333333333,
    0.0757575757575758, -0.253113553113553, 1.16666666666667, -7.09215686274510,
    54.9711779448622, -529.124242424242, 6192.12318840580, -86580.2531135531,	
    1425517.16666667, -27298231.0678161, 601580873.900642, -15116315767.0922,
    429614643061.167, -13711655205088.3, 4.88332318973593e+14, -19296579341940068.0};

  const int bn = 20;		// length of the above array
  const double asx = 10.0;    // value of x beyond which asymptotic is
				// believed to work
  const double xx = GSL_MAX(x, asx+ (x-floor(x))); // the value to go
				// into the asymptotic formula 
  const long recur = lround(GSL_MAX(ceil(asx-x), 0.0));	// number of
				// recursions needed to get to that value

  double f =  - series(bn, b, gsl_pow_2(1/xx), 1);
  f -= 0.5/xx; 			// adding the only odd order term

  // accounting for recursion that brought x's to asymptotic regime
  if(recur>0)
    f -= sum_ovr(x,recur-1);

  return f;
}

