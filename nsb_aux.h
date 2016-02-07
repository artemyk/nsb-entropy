/*
 *  File - nsb_aux.h
 *  Project - Entropy_parsing
 *
 *  Declaration of various temporary functions.
 *
 *  Copyright (c) 2003 Ilya Nemenman. All rights reserved.
 *
 */

#include "vect_math.h"

//  temporary function y = psi(x+1.0)
inline double psi_p1(const double& x){
  return psi(x+1.0);
}
