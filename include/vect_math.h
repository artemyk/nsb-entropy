/*
 *  File - vect_math.h
 *  
 *  Declares necessary routines for various combinations of 
 *  fast inner products.
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef VECT_MATH_H
#define VECT_MATH_H

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>
#include <valarray>

typedef std::valarray<double> VD;

// polygamma functions
inline double psi(double x)  {return gsl_sf_psi(x);}
inline double psi1(double x) {return gsl_sf_psi_n(1,x);}
inline double psi2(double x) {return gsl_sf_psi_n(2,x);}
inline double psi3(double x) {return gsl_sf_psi_n(3,x);}
inline double psi4(double x) {return gsl_sf_psi_n(4,x);}

// 
// cumulative sums, producsts, etc 
// these all should return zero if size of arrays is 0 or less
//

/* 
   sum of elements in a double and integer arrays of length sz 
   accumulation is done in the last parameter; no overflow checking 
   is done
*/
inline double vec_sum(int sz, const double* a, double ret=0.0){
  for(int i=0; i<sz;i++)
    ret += a[i];
  return ret;
}
inline int vec_sum(int sz, const int* a, int ret=0){
  for(int i=0; i<sz;i++)
    ret += a[i];
  return ret;
}

/*
  sum of the form 1/(x+i) for i=0...N;
  summing in reverse order to have better precision
*/
inline double sum_ovr(const double& x, const long n, double sum=0.0){
  for(long i=n;i>=0;i--)
    sum+=1/(x+i);
  return sum;
}

/*
  factorial fraction; given b[0]=1/a, we have b[i]=b[i-1]/(i+a), by default c=1.0
*/
inline void fact_frac(int sz, const double& a, double* b){
  if(sz>0){
    b[0]=1/a;
    for(int i=1;i<sz;i++)
      b[i]=b[i-1]/(a+i);
  };
}


//
// various inner products
//

/*
  dot product between two regular arrays up to the given length,
  with accumulation added up to the fourth argument
*/
inline double dot(int sz, const double* a, const double* b,
		  double prod=0.0){
  for(int i=0; i < sz; i++)
    prod += a[i]*b[i];
  return prod;       
}

/*
  dot product between two vectors, the first being 
  offset by the third argument
  if vectors are of unequal size, dotting is done up to the smallest
  length 
  accumulation is added up to the lasy argument
*/
inline double dot_offset(const VD& a, const VD& b, 
			 const double offset, double prod=0.0){
  for(int i=0; i < (int) GSL_MIN(a.size(), b.size());i++)
    prod += (a[i]+offset)*b[i];
  return prod;       
}

/*
  Make f(a+a0)*b in Octave notation. For vectors a,b and a0 --
  constants this is equivalent to sum(f(a+a0).*b).  f is some function
  taking one double argument; accumulation is added up to the last
  argument, and the size of arrays is sz, the first argument.
*/
inline double dot_f_offset(int sz, const double* a, 
			   const double* b,
			   const double a0, 
			   double (*func)(double),
			   double prod=0.0){
  for(int i=0; i < sz;i++)
    prod += (*func)(a[i]+a0)*b[i];
  return prod;       
}

/*
  Make f(n,a+a0)*b in Octave notation. For vectors a,b and a0 --
  constants this is equivalent to sum(f(a+a0).*b).  f is some function
  taking one int and one double argument; accumulation is added up to
  the last argument, and the size of arrays is sz, the first argument.
*/
inline double dot_f_int_offset(int sz, const double* a, 
			       const double* b,
			       const double a0, 
			       double (*func)(int, double),
			       int n, double prod=0.0){
  for(int i=0; i < sz;i++)
    prod += (*func)(n, a[i]+a0)*b[i];
  return prod;       
}
/*
  make (a.*b)'*c in Octave notation for vectors a,b,c of the same length 
  this is equivalent to sum(a.*b.*c)
  if vectors are of unequal size, operations are done up to the smallest
  length 
  accumulation is added up to the lasy argument
*/
inline double tri_dot(const VD& a, const VD& b, const VD& c,
				 const double offset, double prod=0.0){
  for(int i=0; i < (int) GSL_MIN(GSL_MIN(a.size(), b.size()),c.size());i++)
    prod += a[i]*b[i]*c[i];
  return prod;       
}

/*
  make ((a+a0).*(b+b0))'*c in Octave notation 
  for vectors a,b,c of the same length, and a0, b0 -- constants
  this is equivalent to sum((a+a0).*(b+b0).*c)
  if vectors are of unequal size, operations are done up to the smallest
  length 
  accumulation is added up to the lasy argument
*/
inline double tri_dot_offset2(const VD& a, const VD& b, const VD& c,
			      const double a0, const double b0, 
			      double prod=0.0){
  for(int i=0; i < (int) GSL_MIN(GSL_MIN(a.size(), b.size()),c.size());i++)
    prod += (a[i]+a0)*(b[i]+b0)*c[i];
  return prod;       
}

/*
  Make (f(a+a0).*(b+b0))'*c in Octave notation. For vectors a,b,c and
  a0, b0 -- constants this is equivalent to sum(f(a+a0).*(b+b0).*c).
  f is some function taking one double argument; accumulation is added
  up to the last argument, and the size of arrays is sz, the first
  argument.
*/
inline double tri_dot_f_offset2(int sz, const double* a, 
				const double* b, const double* c,
				const double a0, const double b0, 
				double (*func)(double),
				double prod=0.0){
  for(int i=0; i < sz;i++)
    prod += (*func)(a[i]+a0)*(b[i]+b0)*c[i];
  return prod;       
}


/*
  Make f(a)*b in Octave notation. For vectors a,b. f is some function
  taking a double argument; accumulation is added up to the last
  argument, and the size of arrays is sz, the first argument.
*/
inline double dot_f(int sz, const double* a, 
		    const double* b, 
		    double (*func)(double),
		    double prod=0.0){
  for(int i=0; i < sz;i++)
    prod += (*func)(a[i])*b[i];
  return prod;       
}

/*
  Make f(n,a)*b in Octave notation. For vectors a,b and n integer. f
  is some function taking an integer and a double argument;
  accumulation is added up to the last argument, and the size of
  arrays is sz, the first argument.
*/
inline double dot_f_int(int sz, const double* a, 
			const double* b, 
			double (*func)(int, double),
			int n, double prod=0.0){
  for(int i=0; i < sz;i++)
    prod += (*func)(n, a[i])*b[i];
  return prod;       
}


//
// series summation
//
// we always some series backwards, since the last terms are usually
// smaller, and this allows not to loose precision
//

/*
  summing up a[i]*x^(i+n), where n is integer
*/
inline double series(int sz, const double* a, const double& x, int n,  
		     double result=0.0){
  
  for(int i =sz-1;i>=0;i--)
    result += a[i]*gsl_pow_int(x,i+n);

  return result;
}

/*
  summing up a[i]*x^(i+n)*(*f)(i+m,y), where n is integer, y is double
*/
inline double series_f(int sz, const double* a, const double& x, int n,
		       double (*f)(int, double), const double& y, int m,
		       double result=0.0){
  for(int i =sz-1;i>=0;i--)
    result += a[i]*gsl_pow_int(x,i+n)*((*f)(i+m,y));

  return result;
}

#endif
