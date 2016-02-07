/*
 *  File - counts.h
 *  Project - Entropy_parsing
 *
 *  Declarations and inline definitions for the counter class.  The
 *  counts class contains the map (usig notation of thensb_num
 *  projects) nx -> kx, that is "event occured nx times" is mapped
 *  into "kx events like that".  
 *
 *  Copyright (c) 2001, 2002 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef COUNTS_H
#define COUNTS_H

#include <map>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "EntrData_except.h"

namespace ENTRDATA{
  typedef std::map<long, double> MCO; // map of "event occured n times" -> "m events like that"
				// this corresponds to nx->kx map in octave code
  typedef MCO::const_iterator CIMCO;
  
  
  class counts{
  private:
    double        mne;		// maximum number of possible events =
				// alphabet size ^ length; this is
				// also K in octave code 
    mutable double  ne;		// number of events in the map;
				// =sum(count*events with these
				// count) 
    mutable MCO   data;		// the data; need to declare it
				// mutable since data[0] may need to
				// be recalculated often
    mutable int   isvalid;      // is the value of ne and data[0] valid?
    mutable int   prec;		// printing precision, default 14

    void    calculate() const;	// calculating ne and data[0]

  public:
    counts(double d): mne(d), ne(0), data(), isvalid(0), prec(14) {};

    ~counts(){};
    
    int     size() const {return data.size();} // number of different counts
    CIMCO   begin() const {return data.begin();} // first data element
    CIMCO   end() const {return data.end();} // end of data
    double  operator[](const long&) const throw (EntrData_bad_alloc); 
				// subscripting operator, cannot be used 
				// to change values in data
    double& operator()(const long&) throw (EntrData_bad_alloc);
				// tricky subscripting to be used to change data
    double  NEvents() const;	// total number of events in the data
    double  Cardinality() const {return mne;};

    int precision(int newprec=0)const{ // setting and reading printing precision
      if(newprec) prec=newprec; // if 0, keep old pecision
      return prec;		// return newly set or old precision
    }
        
    friend std::ostream& operator<<(std::ostream&, const counts&); // printing the counts
				// with some precision (default 14)
  };

  std::ostream& operator<<(std::ostream&, const counts&); // printing the counts
				// with some precision (default 14)
  
}

#endif
