/*
 *  File - nsb.h
 *  Project - Entropy_parsing
 *
 *  Declarations of the class used to do NSB numerics
 *
 *  Copyright (c) 2001-2010 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef NSB_H
#define NSB_H

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <limits>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>


#include "EntrData_except.h"
#include "EntrData.h"
#include "counts.h"
#include "vect_math.h"
#include "specfun.h"
#include "specfunctions/specfunctions.h"
#include "BxiK_interp.h"
#include "integration.h"

extern ENTRDATA::VERBOSITY vlevel;

namespace ENTRDATA{
  class NsbCalc;		// main NSB class
  class NsbSet;			// list of NSB objects


  const std::string BndName("bnd");
  const std::string UniName("uni");
  //  const std::string ExtName("ext");


  /* **************
     function wrapper for calculating
     the prior; this is a uniform prior on (a,b)
  */
  class Prior{			
    std::string name;
    double min;			// mini
    double max;
    double ovrZ;
  public:
    Prior(double a, double b, const std::string& n=BndName): 
      name(n), min(a<=b?a:b), max(b>a?b:a), ovrZ(1/fabs(a-b)){};
    virtual ~Prior(){};
    virtual double eval(const double& xi)const{return ovrZ;}
    //(xi>=GetMin()&&xi<=GetMax())?ovrZ:0.0;}; this function returns
    //nonzero value (ovrZ) everywhere so that later integration
    //routines can override it, if needed, and integrate beyond (a,b)
    const char* GetName() const{
      return name.c_str();
    }
    const double GetMin() const{return min;}
    const double GetMax() const{return max;}
  };

  /* ****************
     uniform prior on (0,log(K))
   */
  class UniPrior: public Prior{
  public:
    UniPrior(const double& logK): Prior(0.0,logK,UniName){}
  };


//   /* *******************
//      extrapolated prior; this class is not used anymore
//    */
//   class ExtrapPrior: public Prior{
//     double mean;		// mean of the previous estimate
//     double stddev;		// std.dev. of the previous estimate
//     double ovrZ;		// 1/normalization
//     double mpincr;		// mean plus maximum possible extropy 
// 				// increment
//   public:
//     ExtrapPrior(const double& m, const double& s, const double& ic): Prior(ExtName),
//       mean(m), stddev(s), ovrZ(1/(sqrt(2*M_PI)*stddev + ic)), mpincr(m+ic){
      
//       if((mean<0.0)||(stddev<0.0)||(ic<0.0)){
// 	std::cout<<mean<<" "<<stddev<<" "<<ic;
// 	throw EntrData_badnum(NSB_PARRANGE);
//       }
//     }
//     double eval(const double& xi)const{
//       if(xi<mean)
// 	return ovrZ*exp(-0.5*gsl_pow_2((xi-mean)/stddev));
//       else if (xi<mpincr)
// 	return ovrZ;
//       else
// 	return ovrZ*exp(-0.5*gsl_pow_2((xi-mpincr)/stddev));
//     }
//   };


  enum NSB_WARN {NSB_OK, NSB_NOCOINC, NSB_NOEVENTS, 
		 NSB_ALLCOINC, NSB_SER_B0, NSB_NR_B0_SIGN,
		 NSB_NR_B0, NSB_NR_BCL_SIGN, NSB_NR_BCL, 
		 NSB_INT_NOCONV, NSB_MANY_COINC};
  extern const char* warn_names[]; // computational error names


  typedef std::vector<NsbCalc*> LNC; // was formerly list, hence "L"
  typedef LNC::const_iterator ciLNC;
  typedef LNC::iterator iLNC;

  typedef std::vector<double> vecdoub; //vector of doubles


  


  /* *******************************************
     The main class doing (and keeping)
     all NSB calculations.
     This class is virtual and requires defining the
     function prior_xi in its children.
  */
  class NsbCalc{
  private:
    int refs;			// number of references to the object

    double* nx;			// kx bins had nx counts in them
    double* kx;			// item with nx==0 is not recorded
    double* kxng1;		// same as above but for those n
    double* nxng1;		// which are >1
    unsigned long sz;		// size of nx, k_x
    unsigned long szng1;	// size of nxng1, kxng1
				// both sizes=0 indicates that no
				// corresponding arrays are not created yet
    double  K;			// total number of bins
    double  maxent;		// maximum entropy (log(K))
    double  N;			// number of events 
    double  K1;			// number of bins with 1 or more
				// counts
    double  K2;			// -"- with two or more counts
    const BxiKinterp* interp;   // the interpolator for given K and prec
    const Prior* pr;		// pointer to the prior
    DOPART  part;		// whether the data is partitioned by
				// the sampling properties before
				// calculations

    NsbCalc* same_as;		// pointer to another obj with the
				// same data, NULL if there were none
    
    double err;			// allowed precision error
    double Snsb, dSnsb, Sml, Scl; // class, ml, and nsb values
    double Sas, dSas;		// asymptotic values
    double Bcl, xicl, dxicl;	
    double mlog;	        // value of the -log(evidence) at the saddle



    static const int maxcounter; // maximum number of newton-raphson iterations
    NSB_WARN warncode;		// error code (0 -- no errors)

    void warning(NSB_WARN warn){ // report  and record warnings
      std::cerr<< "***NSB WARNING***: " <<warn_names[warn] <<"\n";
      warncode = warn;
    }

    // writing results for K<=1; we also use this now to 
    // output the case when N==0
    void one_bin();
    // copying calculations from elsewhere or doing them
    void CopyCalcs(NsbCalc&);	
    NSB_WARN DoCalcs() throw(EntrData_badnum, EntrData_bad_alloc); 
				// in Octave -- find_nsb_entropy
    NSB_WARN DoCalcsPart(BXIset&) throw(EntrData_badnum, EntrData_bad_alloc);
				// calculation of entropies partitioning the data
    double mlog_evidence(const double&) const 
      throw(EntrData_badnum, EntrData_bad_alloc); 

    // a priori entropy and its variance (or stddev?) for given
    // of beta, and vice versa, B= B(xi)
    double xi_KB(const double& B) const {return interp->xi_KB(B);}
    double dxi_KB(const double& B) const {return interp->dxi_KB(B);} 
    double B_xiK(const double& B) const {return interp->B_xiK(B);}
    
    // mean value of entropy and entropy ^2 for give beta and counts
    double meanS(const double&) const;
    double meanS2(const double&) const throw(EntrData_bad_alloc);

    void make_ng1() throw(EntrData_bad_alloc,EntrData_badnum); // creating n>1 data
    NSB_WARN max_evidence();	// finding the saddle point

    double prior_xi(const double& xi) const{	// evaluating the prior
      return pr->eval(xi);
    }

    void CommonConstructor(BXIset&);
  public:
    NsbCalc(const counts&, const NsbSet&, BXIset&, const Prior*, DOPART=NOPART, double=1e-6) 
      throw (EntrData_bad_alloc, EntrData_badnum);
    NsbCalc(const vecdoub&, const vecdoub&, double, BXIset&,
	    const Prior*, DOPART, double)
      throw (EntrData_bad_alloc, EntrData_badnum);
    ~NsbCalc();
    
    int addref(){return ++refs;}
    int subref(){return (refs==0)?0:--refs;}

    double GetErr()const{return err;} // get allowed numerical error
    int size() const {return sz;} // returning the size of kx, nx
    NsbCalc* copied_from() const {return same_as;} // was this new data, copied?

    friend int operator==(const NsbCalc&, const NsbCalc&); // comparing
				// two structures en mass
    // integrands
    double int_1(const double&) const;
    double int_S(const double&) const;
    double int_S2(const double&) const;
    
    double get_Snsb()const{return Snsb;} // nsb entropy
    double get_dSnsb()const{return dSnsb;} // nsb posterior std dev
    double get_Sml()const{return Sml;} // maximum likelihood entropy
    double get_Scl()const{return Scl;} // value of s at the saddle
    double get_Sas()const{return Sas;} // small coincidence asymptotic
    double get_dSas()const{return Sas;} // small coincidence asymptotic for std.dev.
    double get_Bcl()const{return Bcl;} // saddle value for beta
    double get_xicl()const{return xicl;} // saddle value for a priori entropy
    double get_dxicl()const{return dxicl;} // std. dev. around the saddle

    double get_K()const{return K;} // cardinality
    double get_N()const{return N;} // # of data points7
    double get_K1()const{return K1;} // number of bins with nonzero occupancy
    double get_K2()const{return K2;} // number of bins with occupancy >1
    const Prior* get_pr()const{return pr;} // pointer to the prior
    DOPART get_part()const{return part;}// is the data partitioned for calculations?
    NSB_WARN get_warncode()const{return warncode;}; // warning code
  };

  int operator==(const NsbCalc&, const NsbCalc&); // comparing
				// two structures en mass
    
  
  // integration routines envelopes used to call the member
  // integration routines from nsb class objects; g_ in the name
  // stands for "global"
  double g_int_1(double, void*);
  double g_int_S(double, void*);
  double g_int_S2(double, void*);


  /* ******************************************
     a container class that will hold a list of 
     pointers to NsbCalc classes; if some NsbCalc's appear to be the
     same, two or more pointers will point to just one object, saving
     the storage space
  */
  class NsbSet{
  private:
    LNC data;
    BXIset* interp;
  public:
    NsbSet(int n=0): data(0), interp(NULL){data.resize(0);data.reserve(n);}
    void set_interp(BXIset* spl){interp=spl;}
    void reserve(int n){data.reserve(n);}
    unsigned int capacity(){return data.capacity();}
    void add(const counts&,  const Prior*, DOPART, double=1e-6) 
      throw (EntrData_bad_alloc, EntrData_badnum);

    ciLNC begin() const{return data.begin();}
    ciLNC end() const{return data.end();}
    const NsbCalc* operator[](int i)const{return data[i];}
    int size() const {return data.size();}

    // destructor -- removing all nsb structures
    ~NsbSet(){			// remember: delete 0;
				// is not an error
      for(iLNC p=data.begin(); p!=data.end();p++){
	if(!((*p)->subref()))	// if subtracting reference counter
				// gives zero, erase
	  delete *p;
	
	(*p)=NULL;		// in any case, the pointr iz zeroed
      }
    }
  };


  // vector of NsbSet's; actually, this is a matrix, as shown by
  // Print(). The matrix is depth * phases
  class VNsbSet:public std::vector<NsbSet>{
    int stride;			// VNsbSet is a matrix; stride is its
				// row length
    int start_depth;		// first depth recorded in the structure
    int end_depth;		// end depth recorded here
  public:
    // i -- number of elements, s -- stride, n -- rserved length of 
    // each element
    VNsbSet(BXIset& sp, int i, int s, int n, int stde, int endde):
      std::vector<NsbSet>(i),stride(s), start_depth(stde),
      end_depth(endde) 
    { 
      reserve_all(n); 
      for(int j=0; j<i;j++){
	(*this)[j].set_interp(&sp);
      }
    }

    void reserve_all(int n){
      for(int j=0; j<(int)size();j++)
	if( (int)(*this)[j].size()<n)
	  (*this)[j].reserve(n);
    }

    int Stride()const{return stride;}
    int SDepth()const{return start_depth;}
    int EDepth()const{return end_depth;}
    ~VNsbSet(){}
    void Print(const std::string&, const std::vector<std::string>&) const;
  };
}

#endif
