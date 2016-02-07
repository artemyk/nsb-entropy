//
// a wrapper for the adaptive quadratures workspace structure
//

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <stdexcept>
#include <iostream>
#include <gsl/gsl_integration.h>


class workspace_AQ{
private:
  gsl_integration_workspace* ws;
  int sz;
  
public:
  workspace_AQ(int n=10000) throw(std::bad_alloc): ws(NULL), sz(n) {
    ws = gsl_integration_workspace_alloc(n);
    if (ws==NULL)
      throw std::bad_alloc();
  };				
  
  ~workspace_AQ(){
    gsl_integration_workspace_free (ws);
  };
  
  int get_size(){return sz;}
  gsl_integration_workspace* get_ws(){return ws;}
};

#endif
