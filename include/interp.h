//
// a wrapper for the interpolation  structure
//

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <stdexcept>
#include <gsl/gsl_spline.h>


class spline{
private:
  gsl_spline* sp;
  int sz;
  gsl_interp_accel* acc;

public:
  // notice that we are using linear interpolation to
  // conserve memory
  spline(const double* x, const double* y, int n, 
	 const gsl_interp_type* T=gsl_interp_linear) throw(std::bad_alloc): sz(n) {
    sp = gsl_spline_alloc(T, sz);
    if (sp==NULL) throw std::bad_alloc();
    gsl_spline_init(sp,x,y,sz);
    acc = gsl_interp_accel_alloc();
    if (acc==NULL) throw std::bad_alloc();
  }

  double eval(const double& xnew) const{
    return gsl_spline_eval(sp, xnew, acc);
  }
  
  ~spline(){
    gsl_spline_free(sp);
    gsl_interp_accel_free(acc);
  };
  
  int get_size(){return sz;}
};

#endif
