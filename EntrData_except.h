/*
 *  File - EntrData_except.h
 *  Project - Entropy_parsing
 *
 *  Exceptions that may arise in various data handling routines.
 *
 *  Created by nemenman on Fri Nov 02 2001.
 *  Copyright (c) 2001 Ilya Nemenman. All rights reserved.
 *
 */

#ifndef ENTRDATA_EXCEPT_H
#define ENTRDATA_EXCEPT_H

#include <exception>
#include <stdexcept>
#include <string>
#include <new>
#include <gsl/gsl_errno.h>

#include "utility.h"

namespace ENTRDATA {

  // main error class
  class EntrDataException {
  protected:
    std::string msg;
  public:
    EntrDataException(const std::string& which="Unknown error."):msg(which) {}
    virtual const char* what() const {return msg.c_str();}
    virtual ~EntrDataException() {}
  };
  
  // allocation error class; for bad memory or resource (like file)
  // allocation
  class EntrData_bad_alloc: public std::bad_alloc, public EntrDataException {
  public:
    EntrData_bad_alloc(const std::string& which="Memory allocation error."): 
      bad_alloc(), EntrDataException(which) {};
    ~EntrData_bad_alloc() throw() {};
  };
  
  // data or arguments to internal functions are out of range.    
  class EntrData_range_error: public std::range_error, 
					  public EntrDataException {
  public:
    EntrData_range_error (const std::string& what_arg): 
      range_error (what_arg), EntrDataException(what_arg) { };
    ~EntrData_range_error () throw() { };
  };
  
  // problem with translation
  class EntrData_failed_transl: public EntrDataException{
  public:
    EntrData_failed_transl(const std::string& what_arg): EntrDataException(what_arg){}
  };
  
  // numerical error codes
  enum ENTRDATA_BADNUM_EC {NSB_CALLORDER, NSB_GSL, NSB_ERANGE, NSB_PARRANGE};
  // nsb numerical error names 
  extern const char* badnum_names[];
 
  // problem with numerics
  class EntrData_badnum: public EntrDataException{
  private:
    ENTRDATA_BADNUM_EC ec;	// main error code
    int ec_ext;			// error code extension
  public:
    EntrData_badnum(ENTRDATA_BADNUM_EC ecc, int ext=0): EntrDataException("NSB numerical error. "),
					     ec(ecc), ec_ext(ext){}
    const char* what() const {
      std::string msg_ext(msg+badnum_names[ec]);
      if(ec==NSB_GSL)
	msg_ext = msg_ext + gsl_strerror (ec_ext );
      return msg_ext.c_str();
    }
    ENTRDATA_BADNUM_EC which() const {return ec;}
    int which_ext() const {return ec_ext;}
  };
}

#endif
