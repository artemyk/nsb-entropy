/*
 *  File - BxiK_interp.h
 *  Project - Entropy_parsing
 *
 *  Declarations of the classes used for interpolating for finding beta as 
 *  a function of xi.
 *
 *  Copyright (c) 2003 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef BXIK_INTERP_H
#define BXIK_INTERP_H

#include <iostream>
#include <list>
#include <limits>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_math.h>

#include "EntrData_except.h"
#include "EntrData.h"
#include "interp.h"
#include "vect_math.h"

extern ENTRDATA::VERBOSITY vlevel;

const double psi0_1= gsl_sf_psi(1); // necessary constant

namespace ENTRDATA{
  class BxiKinterp{
    
    spline* low;		// spline for B<=1; remember: B=K\beta
    spline* med;		// spline for 1<B<=K
    spline* high;		// spline for B>K
    double K;			// cardinality of the space
    double maxent;		// maximum entropy (log K)
    double ltr;		        // lower threshold (K is the upper threshold)
    double prec;		// interpolation precision 
				// (roughly, the number of nodes in the splines)
    double slope;		// slope for the low-B regime
    double ovrslope;		// 1/slope
    double ltrxi;		// low threshold for xi
    double htrxi;		// hi threshold for xi

    double low_var(const double& B) const { // lower range transformed variable
      return slope*B;
    }
    double low_B(const double& var) const { // reverse x-form
      return var*ovrslope;
    }

    double med_var(const double& B) const { // medium range transformed variable
      return log(B) - psi0_1;
    }
    double med_B(const double& var) const { // reverse x-form
      return exp(var+psi0_1);
    }

    double high_var(const double& B) const { // high range transfomed variable
      return maxent - K/(2*B);
    }
    double high_B(const double& var) const { // reverse x-form
      return 0.5*K/(maxent-var);
    }

  public:
    double xi_KB(const double& B) const{ // inverse function xi(B,K) 
      if(gsl_isinf(B)){ return maxent;}
      return psi(B+1.0) - psi(1.0+B/K);}
    double dxi_KB(const double& B) const // and its derivative
    {return psi1(B+1.0) - 1/K*psi1(1.0+B/K);} 
    // B as a function of xi
    double B_xiK(const double& xi) const throw (EntrData_badnum) {
      if((xi<0.0)||(xi>maxent))
	throw EntrData_badnum(NSB_ERANGE);
      if(xi<ltrxi)
	return low_B(low->eval(xi));
      else if(xi<htrxi)
	return med_B(med->eval(xi));
      else if(xi<maxent)
	return high_B(high->eval(xi));
      else
	return std::numeric_limits<double>::infinity();
    }
    double get_K(){return K;};	// figuring out the alphabet size
    double get_prec(){return prec;} // getting the precision

    BxiKinterp(const double&, double=1e-6) throw (EntrData_bad_alloc);
    ~BxiKinterp();
  };

  typedef std::list<BxiKinterp*> LBXIK;
  typedef LBXIK::const_iterator ciLBXIK;
 
  class BXIset {
  private:
    LBXIK bxis;		        // I chose list of *pointers* because I
				// don't want to implement copy
				// constructors for BxiKinterp objects
    int maxlength;		// data will not be allowed to grow beoynd maxlength
    
  public:
    BXIset(int maxl=1) throw (EntrData_bad_alloc):bxis(),maxlength(maxl){};
    ~BXIset();
    unsigned int size()const{return bxis.size();}
    int get_maxlength()const{return maxlength;}
    const BxiKinterp* get_spline(const double&, double=1e-6) throw (EntrData_bad_alloc);
    ciLBXIK begin()const{return bxis.begin();}
    ciLBXIK end()const{return bxis.end();}
  };

}

#endif
