/*
 *  File - BxiK_interp.cc
 *  Project - Entropy_parsing
 *
 *  Spline interpolation for b(xi) definitions.
 *
 *  Copyright (c) 2003 Ilya Nemenman. All rights reserved.
 *
 */

#include "BxiK_interp.h"


//constructor
ENTRDATA::BxiKinterp::BxiKinterp(const double& K_in, 
				 double prec_in) 
  throw(EntrData_bad_alloc): K(K_in), maxent(log(K)), ltr(1.0), prec(prec_in), 
			     slope((K-1)/K*gsl_sf_psi_n(1,1.0)),
			     ovrslope(1/slope){

  double* x;			// transformed variable
  double* y;			// xi
  double range;			// range of the transformed variable
  double d_tr;			// step of transformed varible
  double tmpB;			// temporary value for B

  // this function is a lot different from that in Octave 
  // we split the range of B=K\beta into three asymptotic regimes
  // (B<1, 1<B<K, K<B), and use the asymptotic expansion of B in each
  // of the regimes as an equispaced variable (easily analytically
  // related to B). we then build three different (overlapped for
  // better precision) splines to relate xi to this asymptotic
  // variable (and thus to B) in each of the regimes

  ltrxi = xi_KB(ltr);
  htrxi = xi_KB(K);

  if(vlevel>VERBMIN)
    std::cout<< "  Creating spline data for K=" << K << ". ";
  const long npts = (long) (1/prec);	// number of interpolation
					// points
  try{
    x= new double[npts];
    y= new double[npts];
  }catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  }

  // filling low range
  range = 1.1*low_var(ltr);
  d_tr = range/npts;
  for(long i=0; i<npts;i++){
    x[i] = i*d_tr;
    tmpB = low_B(x[i]);
    y[i] = xi_KB(tmpB);
  }

  // creating the low range spline
  try{
    low = new spline(y,x,npts);
  }catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  }

  // filling medium range
  range = (1.1*med_var(K) - 0.9*med_var(ltr));
  d_tr = range/npts;
  x[0] = 0.9*med_var(ltr);
  for(long i=1; i<npts;i++){
    x[i] = x[i-1]+d_tr;
    tmpB = med_B(x[i]);
    y[i] = xi_KB(tmpB);
  }

  // creating the medium range spline
  try{
    med = new spline(y,x,npts);
  }catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  }

  // filling high range
  range = maxent - 0.9*high_var(K);
  d_tr = range/npts;
  x[0] = 0.9*high_var(K);
  for(long i=1; i<npts;i++){
    x[i] = x[i-1]+d_tr;
    tmpB = high_B(x[i]);
    y[i] = xi_KB(tmpB);
  }

  // creating the high range spline
  try{
    high = new spline(y,x,npts);
  }catch(std::bad_alloc){
    throw EntrData_bad_alloc();
  }
  if(vlevel>VERBMIN)
    std::cout<<"Done.\n";

  delete[] x;			// cleaning up the x,y variables
  delete[] y;

}


//
// destructor
//    
ENTRDATA::BxiKinterp::~BxiKinterp(){
  delete low;
  delete med;
  delete high;
}


/*******************************************
  list of the  BxiKinterp's
*/


// destructor
ENTRDATA::BXIset::~BXIset(){
  // removing all interpolation objects
  for(ciLBXIK p=bxis.begin(); p!=bxis.end(); ++p)
    delete *p;
}

// returning a const pointer to the spline object able to generate
// B(xi) for a given value of K and precision
const ENTRDATA::BxiKinterp* ENTRDATA::BXIset::get_spline(const double& K, double prec) 
  throw (EntrData_bad_alloc){
  
  ciLBXIK p= bxis.begin();
  
  while((p!=bxis.end())&& !(((*p)->get_K()==K)&&((*p)->get_prec()==prec))){
    p++;
  }
  
  // if found the same K, prec, return the appropriate pointer
  if(p!=bxis.end()) return *p;	
  
  
  // if not found
  try{
    while(((int)size()+1)>maxlength){ // prune the set, keeping only last
				// "maxlength"-1 splines 
      delete *(bxis.begin());
      bxis.pop_front();
    }
    BxiKinterp* x = new BxiKinterp(K,prec);
    // create new spliine tables
    bxis.push_back(x);		// add them to the list

    return x;			// return the pointer to the new splines
  }
  catch (std::bad_alloc){
    throw EntrData_bad_alloc();
  }
}

