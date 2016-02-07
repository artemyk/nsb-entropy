/*
 *  File - utility.h
 *  Project - Entropy_parsing
 * 
 *  Numbers to strings conversions utilities, for GCC 3.1
 *
 *  Created by nemenman on Fri Nov 02 2001.
 *  Copyright (c) 2001 Ilya Nemenman. All rights reserved.
 *
 */



#include <sstream>


//conversion of a long to a string
inline std::string num2str (const long i){
  std::ostringstream ost;
  ost<<i;
  return ost.str();
}

//conversion of an unsigned long to a string
inline std::string num2str(const unsigned long i){
  std::ostringstream ost;
  ost<<i;
  return ost.str();
}

//conversion of an int to a string
inline std::string num2str(const int i){
  std::ostringstream ost;
  ost<<i;
  return ost.str();
}

//conversion of a double to a string
inline std::string num2str(const double d){
  std::ostringstream ost;
  ost<<d;
  return ost.str();
}
